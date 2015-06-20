/*
 * detwin.c
 * This code uses Expectation Maximization algorithm
 * We took Tom's process_hkl code and added bunch of functions to do detwin work.
 * Haiguang Liu and me wrote and modified the code..
 */

/* 
 * Assemble and process FEL Bragg intensities
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2013 Thomas White <taw@physics.org>
 *   2011      Andrew Martin <andrew.martin@desy.de>
 *   2012      Lorenzo Galli <lorenzo.galli@desy.de>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>

#include "utils.h"
#include "statistics.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "stream.h"
#include "reflist.h"
#include "image.h"
#include "crystal.h"
#include "thread-pool.h"
#include "geometry.h"

#include "detwin.h"

#include <math.h>
#include <gsl/gsl_statistics.h>

#define MAX_N_IMAGE 100000
#define MAX_REFL_PER_IMAGE 20000
#define N_TWINS 2
#define MAX_H 256
#define MAX_K 256
#define MAX_L 256

#define epsilon 1.e-15


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"read from stream file and convert the data to miller indexed FEL Bragg intensities.\n"
"\n"
"  -h, --help                Display this help message.\n"
"  -i, --input=<filename>    Specify input filename (\"-\" for stdin).\n"
"  -o, --output=<filename>   Specify output filename for merged intensities\n"
"                             Default: processed.hkl).\n"
"  -y, --symmetry=<sym>      Merge according to point group <sym>.\n"
"\n"
"      --start-after=<n>     Skip <n> crystals at the start of the stream.\n"
"      --stop-after=<n>      Stop after merging <n> crystals.\n"
"  -g, --histogram=<h,k,l>   Calculate the histogram of measurements for this\n"
"                             reflection.\n"
"  -z, --hist-parameters     Set the range for the histogram and the number of\n"
"          =<min,max,nbins>   bins. \n"
"\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model.\n"
"      --reference=<file>    Compare against intensities from <file> when\n"
"                             scaling. \n"
"      --no-polarisation     Disable polarisation correction.\n"
"      --min-measurements=<n> Require at least <n> measurements before a\n"
"                             reflection appears in the output.  Default: 2\n"
"      --min-snr=<n>         Require individual intensity measurements to\n"
"                             have I > n * sigma(I).  Default: -infinity.\n"
);
}


double compute_linear_correlation( double *x, double *y, int n_ );

void get_equivalents(int h, int k, int l, int *hs, int *ks, int *ls, int *n_twins);

static double scale_intensities(RefList *reference, RefList *new,
                              const SymOpList *sym)
{
        double s;
        double top = 0.0;
        double bot = 0.0;
        Reflection *refl;
        RefListIterator *iter;

        for ( refl = first_refl(new, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) )
        {

                double i1, i2;
                signed int hu, ku, lu;
                signed int h, k, l;
                Reflection *reference_version;

                get_indices(refl, &h, &k, &l);
                get_asymm(sym, h, k, l, &hu, &ku, &lu);
		//hu = h; ku=k;lu=l;
                reference_version = find_refl(reference, hu, ku, lu);
                if ( reference_version == NULL ) continue;

                i1 = get_intensity(reference_version);
                i2 = get_intensity(refl);

                /* Calculate LSQ estimate of scaling factor */
                top += i1 * i2;
                bot += i2 * i2;

        }

        s = top / bot;

        return s;
}

static void display_progress(int n_images, int n_crystals, int n_crystals_used)
{
        if ( !isatty(STDERR_FILENO) ) return;
        if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

        pthread_mutex_lock(&stderr_lock);
        fprintf(stderr, "\r%i images processed, %i crystals, %i crystals used.",
                n_images, n_crystals, n_crystals_used);
        pthread_mutex_unlock(&stderr_lock);

        fflush(stdout);
}

static unsigned char *flags_from_list(RefList *list)
{
        Reflection *refl;
        RefListIterator *iter;
        unsigned char *out = new_arr_flag();

        for ( refl = first_refl(list, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) ) {

                signed int h, k, l;

                get_indices(refl, &h, &k, &l);

                set_arr_flag(out, h, k, l, 1);

        }

        return out;

}

static double *intensities_from_list(RefList *list)
{
        Reflection *refl;
        RefListIterator *iter;
        double *out = new_arr_intensity();

        for ( refl = first_refl(list, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) ) {

                signed int h, k, l;
                double intensity = get_intensity(refl);

                get_indices(refl, &h, &k, &l);

                set_arr_intensity(out, h, k, l, intensity);

        }

        return out;
}



static double sym_lookup_intensity(const double *intensities,
                                   const unsigned char *flags,
                                   const SymOpList *sym,
                                   signed int h, signed int k, signed int l)
{
        int i;
        double ret = 0.0;

        for ( i=0; i<num_equivs(sym, NULL); i++ ) {

                signed int he;
                signed int ke;
                signed int le;
                double f, val;

                get_equiv(sym, NULL, i, h, k, l, &he, &ke, &le);

                f = (double)lookup_arr_flag(flags, he, ke, le);
                val = lookup_arr_intensity(intensities, he, ke, le);

                ret += f*val;

        }

        return ret;
}



static int add_crystal(RefList *model, struct image *image, Crystal *cr,
                         RefList *reference, const SymOpList *sym,
			 RefList *this_image,
                         int config_nopolar, double min_snr)
{
	Reflection *refl;
	RefListIterator *iter;

//	this_image = reflist_new();
	double scale;

	/* First, correct for polarisation */
	if ( !config_nopolar ) {
		polarisation_correction(crystal_get_reflections(cr),
		                        crystal_get_cell(cr), image);
	}

	if ( reference != NULL ) {
		scale = scale_intensities(reference,
		                          crystal_get_reflections(cr), sym);
	} else {
		scale = 1.0;
	}
	if ( isnan(scale) ) return 1;

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double refl_intensity, refl_sigma;
		signed int h, k, l;
		int model_redundancy;
		Reflection *model_version, *new_refl;
		double w;
		double temp, delta, R, mean, M2, sumweight;

		refl_intensity = scale * get_intensity(refl);
		refl_sigma = scale * get_esd_intensity(refl);
		w = 1.0;//pow(refl_sigma, -2.0);

		if ( (min_snr > -INFINITY) && isnan(refl_sigma) ) continue;
		if ( refl_intensity < min_snr * refl_sigma ) continue;

		get_indices(refl, &h, &k, &l);

		/* Put into the asymmetric unit for the target group */
		get_asymm(sym, h, k, l, &h, &k, &l);

		model_version = find_refl(model, h, k, l);
		if ( model_version == NULL ) {
			model_version = add_refl(model, h, k, l);
		}
		

		mean = get_intensity(model_version);
		sumweight = get_temp1(model_version);
		M2 = get_temp2(model_version);

		temp = w + sumweight;
		delta = refl_intensity - mean;
		R = delta * w / temp;
		set_intensity(model_version, mean + R);
		set_temp2(model_version, M2 + sumweight * delta * R);
		set_temp1(model_version, temp);

		model_redundancy = get_redundancy(model_version);
		set_redundancy(model_version, ++model_redundancy);
	//	 add this reflection to the list of this_image 
                new_refl = find_refl(this_image, h, k, l);
                if ( new_refl == NULL ) {
        //                add_refl_to_list(model_version, this_image);
			new_refl = add_refl(this_image, h, k, l);
                }
		set_intensity( new_refl, refl_intensity );
		set_redundancy( new_refl, 1);
	}
	if(num_reflections( this_image ) == 0)  { return 1;}
	else return 0;
}


static RefList** add_all(Stream *st, RefList *model, RefList *reference,
                     const SymOpList *sym,
                     double **hist_vals, signed int hist_h,
                     signed int hist_k, signed int hist_l,
                     int *hist_i, int config_nopolar, int min_measurements,
                     double min_snr, int start_after, int stop_after,
		     int *n_crystals_recorded)
{
	int rval;
	int n_images = 0;
	int n_crystals = 0;
	int n_crystals_used = 0;
	Reflection *refl;
	RefListIterator *iter;
        RefList* image_array[MAX_N_IMAGE];
//        image_array = malloc(sizeof(RefList*) * MAX_N_IMAGE);
	int i;
	for(i=0;i<MAX_N_IMAGE;i++) image_array[i] = reflist_new();

	int n_crystals_seen = 0;

	do {

		struct image image;
		int i;

		image.det = NULL;

		/* Get data from next chunk */
		rval = read_chunk(st, &image);
		if ( rval ) break;

		n_images++;

		for ( i=0; i<image.n_crystals; i++ ) {

			int r;
			Crystal *cr = image.crystals[i];

			n_crystals_seen++;
			if ( n_crystals_seen <= start_after ) continue;

			n_crystals++;
			r = add_crystal(model, &image, cr, reference, sym,
			                image_array[n_crystals_used], 
					config_nopolar, min_snr);

			if ( r == 0 ) n_crystals_used++;

			reflist_free(crystal_get_reflections(cr));
			cell_free(crystal_get_cell(cr));
			crystal_free(cr);

			if ( n_crystals_used == stop_after ) break;

		}

		free(image.filename);
		image_feature_list_free(image.features);
		free(image.crystals);

                display_progress(n_images, n_crystals_seen, n_crystals_used);

		if ( (stop_after>0) && (n_crystals_used == stop_after) ) break;

	} while ( rval == 0 );

	for ( refl = first_refl(model, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		if ( red < min_measurements ) {
			set_redundancy(refl, 0);
			continue;
		}

		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));
	}
//	STATUS("number of crystals recorded: %d \n", n_crystals_used);
	*n_crystals_recorded = n_crystals_used;
	return image_array;
}

void set_esd_for_reflist( RefList *model )
{
	RefListIterator *iter;
	Reflection *refl;

        for ( refl = first_refl(model, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) )
        {
                double var;
                int red;

                red = get_redundancy(refl);
		if( red<1 ) {set_esd_intensity( refl, 0 ); continue; }
                var = get_temp2(refl) / get_temp1(refl);
                set_esd_intensity(refl, sqrt(var)/sqrt(red));
        }
}


void get_p3_equivalents(int h, int k, int l, int *hs, int *ks, int *ls, int *n_twins)
{
	/* for PS-1, P3 */
	/*
http://www.ccp4.ac.uk/html/twinning.html
All P3i and H3:
(h,k,l) neither equivalent to (-h,-k,l) nor (k,h,-l) nor (-k,-h,-l) so we need to check all 4 possibilities. These are the only cases where tetratohedral twinning can occur:
real axes:	(a,b,c)	and	(-a,-b,c)	and	(b,a,-c)	and	(-b,-a,c)
reciprocal axes:	(a*,b*,c*)	and	(-a*,-b*,c*)	and	(b*,a*,-c*)	and	(-b*,-a*,c*)
i.e. For P3, consider reindexing (h,k,l) to (-h,-k,l) or (k,h,-l) or (-k,-h,-l).
	*/
	*n_twins = 3; 
	hs = (int*) malloc( sizeof(int)*(*n_twins) );
	ks = (int*) malloc( sizeof(int)*(*n_twins) );
	ls = (int*) malloc( sizeof(int)*(*n_twins) );
	
	hs[0] = -h; ks[0]=-k; ls[0]= l;
	hs[1] =  k; ks[1]= h; ls[1]=-l;
	hs[2] = -k; ks[2]=-h; ls[2]=-l;
	return;
}

void get_p3_ith_equivalent(int h, int k, int l, int *he, int *ke, int *le, int ith)
{
	if(ith==0) { *he= h; *ke= k; *le= l; return; }
	if(ith==1) { *he=-h; *ke=-k; *le= l; return; }
	if(ith==2) { *he= k; *ke= h; *le=-l; return; }
	if(ith==3) { *he=-k; *ke=-h; *le=-l; return; }
}

void get_p63_ith_twin( int h, int k, int l, int *he, int *ke, int *le, int ith )
{
	if(ith == 0) {*he= h; *ke= k; *le= l; return; }
	if(ith == 1) {*he= h; *ke= -k; *le=-h-l; return; } // twin operator
}

// hard coded here for P63

void stat_pearson_i_p63(RefList *list1, RefList *list2, double * val, const SymOpList *sym, SymOpList *tlaw)
{
        double *vec1, *vec2;
        double *vec3, *vec4;
        //	double *work;
        int ni = num_reflections(list1);
        int nacc = 0;
        int nacc_twin = 0;
        Reflection *refl1;
        RefListIterator *iter;

        vec1 = malloc(ni*sizeof(double));
        vec2 = malloc(ni*sizeof(double));

        vec3 = malloc(ni*sizeof(double));
        vec4 = malloc(ni*sizeof(double));
//        work = malloc(2*ni*sizeof(double));
//	printf("number of reflections %d \n", ni);

        for ( refl1 = first_refl(list1, &iter);
              refl1 != NULL;
              refl1 = next_refl(refl1, iter) )
        {
                double i1, i2;
                signed int h, k, l;
                signed int ha, ka, la;  ////Sbasu
                signed int hp, kp, lp;
                Reflection *refl2;

                get_indices(refl1, &h, &k, &l);
		if( abs(l) < 40 ) continue;
		if( h*h+k*k+l*l < 800 ) continue;
		if( h*h+k*k+l*l > 8000 ) continue;
//		get_asymm(sym, h, k, l, &hp, &kp, &lp);
		hp = h; kp=k; lp=l;
                refl2 = find_refl(list2, hp, kp, lp);
                if ( refl2 != NULL ) /* No common reflection */
		{
                	i1 = get_intensity(refl1);
                	i2 = get_intensity(refl2);
			if( i1<=0 || i2<=0 ) continue;

	                vec1[nacc] = i1;
	                vec2[nacc] = i2;
	                nacc++;
		}

	// now, the twin
		get_asymm(tlaw, ha, ka, la, &hp, &kp, &lp); // call the get_asymm function with correct twin law for merging using 'tlaw' operator...Sb
//		hp = k; kp=h; lp=-l;
                refl2 = find_refl(list2, hp, kp, lp);
                if ( refl2 != NULL ) /* No common reflection */
		{
	                i1 = get_intensity(refl1);
	                i2 = get_intensity(refl2);
			if( i1<=0 || i2<=0 ) continue;

	                vec3[nacc_twin] = i1;
	                vec4[nacc_twin] = i2;
	                nacc_twin++;
		}
        }
	//if(nacc > 10000)
	//printf("number of common reflections %d %d \n", nacc, nacc_twin);
//	printf("number of common reflections %d %d \n", nacc, nacc_twin);

	if (nacc < 2 ) val[0]=0;
	else           val[0] = gsl_stats_correlation(vec1, 1, vec2, 1, nacc);
//	else           val[0] = gsl_stats_spearman(vec1, 1, vec2, 1, nacc, work);
//	else           val[0] = compute_linear_correlation(vec1, vec2, nacc);

	if (nacc_twin < 2 ) val[1]=0;
	else        	    val[1] = gsl_stats_correlation(vec3, 1, vec4, 1, nacc_twin);
//	else           val[1] = compute_linear_correlation(vec3, vec4, nacc_twin);

        free(vec1);
        free(vec2);
        free(vec3);
        free(vec4);

        return ;
}

double compute_linear_correlation( double *x, double *y, int n_ ) {
	double mean_x_=0.0;
	double mean_y_=0.0;
	double delta_x_=0.0;
	double delta_y_=0.0;
	double numerator_ = 0.0;
	double denominator_ = 0.0;
	double sum_denominator_x_ = 0.0;
	double sum_denominator_y_ = 0.0;
	double coefficient_ = 0.0;

	int i;

        for(i=0;i<n_;i++) mean_x_ += x[i];
        for(i=0;i<n_;i++) mean_y_ += y[i];
        mean_x_ /= n_;
        mean_y_ /= n_;
        for( i=0;i<n_;i++) {
            delta_x_ = x[i] - mean_x_;
            delta_y_ = y[i] - mean_y_;
            numerator_ += delta_x_ * delta_y_;
            sum_denominator_x_ += delta_x_ * delta_x_;
            sum_denominator_y_ += delta_y_ * delta_y_;
        }
        denominator_ = sqrt(sum_denominator_x_ * sum_denominator_y_);
//	printf("np: %d, num %f, denom %f\n", n_, numerator_, denominator_);
        if (numerator_ == 0 && denominator_ == 0) {
          coefficient_ = 1;
        }
        else if (denominator_ > abs(numerator_ * epsilon)) {
          coefficient_ = numerator_ / denominator_;
	}
	return coefficient_;
}

int index_of_max_value( double* cc, int n_twins )
{
	int i;
	int winner=0;
	double max_value = cc[0];
	for(i=1;i<n_twins;i++)
	{
		if(cc[i] >= max_value) { max_value = cc[i]; winner=i; }
	}
	return winner;
}

void merge_image_at_winner_orientation( RefList* model, RefList *image, int winner, double weight, const SymOpList *sym, SymOpList *tlaw )///Sb
{
	Reflection *refl;
	Reflection *model_version;
	RefListIterator *iter;
	int h,k,l, index;
        int model_redundancy;
	double refl_intensity, refl_sigma;
        double w;
        double temp, delta, R, mean, M2, sumweight;	
	w = 1.0; // pow( refl_sigma, -2.0 );
	w = weight;

	//SymOpList *sym;
	//sym = get_pointgroup("2/m_uab");

	for ( refl = first_refl(image, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) )
	{
                get_indices(refl, &h, &k, &l);

		refl_intensity = get_intensity( refl );

		//get_p63_ith_twin(h, k, l, &h, &k, &l, winner);
        get_equiv( tlaw, NULL, 0, h, k, l, &h, &k, &l  ); ///Sb
		get_asymm(sym, h, k, l, &h, &k, &l);
                model_version = find_refl(model, h, k, l);
                if ( model_version == NULL ) {
                        model_version = add_refl(model, h, k, l);
                }

                mean = get_intensity(model_version);
                sumweight = get_temp1(model_version);
                M2 = get_temp2(model_version);
//??? weights
                temp = w + sumweight;
                delta = refl_intensity - mean;
                R = delta * w / temp;
                set_intensity(model_version, mean + R);
                set_temp2(model_version, M2 + sumweight * delta * R);
                set_temp1(model_version, temp);

                model_redundancy = get_redundancy(model_version);
                set_redundancy(model_version, ++model_redundancy);
	}
}

void compute_weights( double *cc, double *weights, int iter )
{
	double sum_w;
	// Let's do the simple case now, N_TWINS=2
	if (cc[0] == cc[1])
	{
		weights[0] = 0.5;
		weights[1] = 0.5;
	}
	else
	{
//		cc[0] += 1.0;
//		cc[1] += 1.0;
		cc[0] = exp( cc[0]*(iter+1) );
		cc[1] = exp( cc[1]*(iter+1) );
		sum_w = cc[0] + cc[1] ;
		weights[0] = cc[0] / sum_w;
		weights[1] = cc[1] / sum_w;
	}
//	printf("%f %f\n", weights[0], weights[1]);
}

void merge_image_to_model( RefList* model, RefList *image, double* cc, int iter )
{
	int twin;
	double weights[ N_TWINS ];
	//printf("%f %f\n", cc[0], cc[1]);
	compute_weights( cc, weights, iter);
	for( twin=0; twin<N_TWINS; twin++)
		merge_image_at_winner_orientation( model, image, twin, weights[twin] );	
}


RefList* make_reflections_for_uc_from_asymm( RefList* asymm, bool random_intensity, SymOpList *sym )
{
        Reflection *refl;
	Reflection *new_refl;
	Reflection *input_refl;
        RefListIterator *iter;
        int h,k,l, index;
	int he,ke,le;
	int ha,ka,la;
	int i_twin, this_N_TWINS;
        double refl_intensity; 
	RefList* this_ref_model;
	this_ref_model=reflist_new();
/* the following lines are for intensity generating for UC from Asymmetric Unit */
        double *intensities;
        unsigned char *flags;
        intensities = intensities_from_list(asymm);
        flags = flags_from_list(asymm);
        //SymOpList *sym;
        //sym = get_pointgroup("2/m_uab");
/* hard coded for P63 */
	if (random_intensity ) this_N_TWINS = N_TWINS; // generating twin pairs
	else		       this_N_TWINS = 1; // only one of the twins
	this_N_TWINS = N_TWINS;

        index = 0;
	srand(time(NULL));

        for ( refl = first_refl(asymm, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) )
        {
                get_indices(refl, &h, &k, &l);
		for(i_twin=0;i_twin<this_N_TWINS;i_twin++) {
			get_p63_ith_twin(h, k, l, &he, &ke, &le, i_twin);
			get_asymm(sym, he, ke, le, &he, &ke, &le);
			
			if( random_intensity ) 
			  refl_intensity = rand();
			else
/*			{
			  get_asymm(sym, he, ke, le, &ha, &ka, &la);
			  input_refl = find_refl(asymm, ha, ka, la);
			  if (input_refl == NULL )
			  {
				printf("could not find %d %d %d \n", ha, ka, la);
				continue;
			  }
			  refl_intensity = get_intensity( input_refl );
			}
// */			  refl_intensity = sym_lookup_intensity(intensities, flags, sym, he, ke, le);
//			printf("twins %d (%d %d %d) i_twin:= %f\n", i_twin, he, ke, le, refl_intensity );

			new_refl = find_refl( this_ref_model, he, ke, le );
			if (new_refl == NULL ) new_refl = add_refl( this_ref_model, he,ke,le );
			set_intensity( new_refl, refl_intensity );
			set_redundancy( new_refl, 1 );
		}
	}

	return this_ref_model;
	
}

void compute_partiality_for_it( RefList **input, RefList *reference, int n_image, int h, int k , int l )
{
        RefList *image_a;
        RefList *partial_list;
        int image_index;
        Reflection *refl;
        Reflection *full_refl;
        Reflection *out_refl;
        RefListIterator *iter;
        double partiality;
        double full_intensity, this_intensity;
        partial_list = reflist_new();
        
	full_refl = find_refl( reference, h, k, l );
        if (full_refl == NULL ) return;
        full_intensity = get_intensity( full_refl );
        if (full_intensity==0) return;

        for( image_index=0;image_index<n_image; image_index++)
        {
                image_a = input[ image_index ];
		refl = find_refl( image_a, h, k, l );
		if( refl == NULL ) continue;
                this_intensity = get_intensity( refl );
                partiality = this_intensity / full_intensity*100000 ;
                out_refl = add_refl( partial_list, h, k, l );
                set_intensity(out_refl, partiality);
                set_redundancy( out_refl, 1 );
        }
        write_reflist("partiality.hkl", partial_list );
}

void compute_partiality( RefList **input, RefList *reference, int n_image, const SymOpList *sym )
{
        RefList *image_a;
        RefList *partial_a;
        RefList *partial_b;
        int image_index;
	int h,k,l;
	Reflection *refl;
	Reflection *full_refl;
	Reflection *out_refl;
	RefListIterator *iter;
	double partiality;
	double full_intensity, this_intensity;
	partial_a = reflist_new();
	partial_b = reflist_new();
    signed int ha, ka, la;
	//SymOpList *sym;
	//sym = get_pointgroup("2/m_uab");

        for( image_index=0;image_index<n_image; image_index++)
        {
                image_a = input[ image_index ];

                for ( refl = first_refl(image_a, &iter);
                      refl != NULL;
                      refl = next_refl(refl, iter) )
                {
			get_indices(refl, &h, &k, &l);
			this_intensity = get_intensity( refl );
			full_refl = find_refl( reference, h, k, l );
			if (full_refl == NULL ) continue;
			full_intensity = get_intensity( full_refl );
			if (full_intensity==0) continue;
			partiality = this_intensity / full_intensity*100000 ;
			out_refl = add_refl( partial_a, h, k, l );
			set_intensity(out_refl, partiality);
			set_redundancy( out_refl, 1 );
// now the twin
			//get_asymm( sym, h, -k, -h-l, &h, &k, &l ); //for P21 pseudo-twin
            get_equiv( tlaw, NULL, 0, ha, ka, la, &h, &k, &l );
                        full_refl = find_refl( reference, h, k, l );
                        if (full_refl == NULL ) continue;
                        full_intensity = get_intensity( full_refl );
                        if (full_intensity==0) continue;
                        partiality = this_intensity / full_intensity*100000 ;
                        out_refl = add_refl( partial_b, h, k, l );
                        set_intensity(out_refl, partiality);
			set_redundancy( out_refl, 1 );
		}
	}
	write_reflist("partiality_a.hkl", partial_a );
	write_reflist("partiality_b.hkl", partial_b );

}

void compile_test_data_sets( RefList **input, RefList **output, RefList *reference, int *twin, int n_image, const SymOpList *sym, SymOpList *tlaw ) ///added sym pointer
{
	RefList *image_a;
	RefList *image_b;
	int image_index;
	int hi,ki,li;  // input h,k,l
	int he,ke,le;
	int ho,ko,lo;  // output h,k,l; twin pair of input
	int half_rand_max = (int)RAND_MAX/2;
	int i_twin;
	Reflection *refl;
	Reflection *input_refl;
	Reflection *output_refl;
	RefListIterator *iter;
	double refl_intensity;
	double threshold=1.0;
	double rand_max_double=(double) RAND_MAX;

	double *intensities;
	unsigned char *flags;

	intensities = intensities_from_list(reference);
        flags = flags_from_list(reference);

	srand( time(NULL) );
	//SymOpList *sym;
	//sym = get_pointgroup("2/m_uab");
	//printf("symmetry 2/m_uab has %d equivalence\n", num_equivs(sym, NULL) );
        int n_twin1 = 0;
	gsl_rng *rng;
	rng = gsl_rng_alloc(gsl_rng_mt19937);

	for( image_index=0;image_index<n_image; image_index++)
        {
                image_a = input[ image_index ];
		if (rand() > half_rand_max ) i_twin = 0;
		else i_twin = 1;
		twin[ image_index ] = i_twin;
		n_twin1 += i_twin;

		for ( refl = first_refl(image_a, &iter);
              	      refl != NULL;
                      refl = next_refl(refl, iter) )
        	{
                	get_indices(refl, &hi, &ki, &li);
			get_p63_ith_twin(hi, ki, li, &he, &ke, &le, i_twin);
/*
			get_asymm(sym, he, ke, le, &ho, &ko, &lo);
                        input_refl = find_refl( reference, ho, ko, lo );

                        if (input_refl == NULL ) continue;
			refl_intensity = get_intensity( input_refl );  // intensity is from twin 
*/
			refl_intensity = sym_lookup_intensity(intensities, flags, sym, he, ke, le);
			//output_refl = add_refl( image_b, hi,ki,li );   // indices are from original
			output_refl = refl;
			refl_intensity += poisson_noise( rng, refl_intensity ); //poisson noise
//			refl_intensity += sqrt( refl_intensity );
			refl_intensity *= (rand()/rand_max_double*threshold);

                        set_intensity( output_refl, refl_intensity );
                        set_redundancy( output_refl, 1 );
		}
		output[ image_index ] = image_a; 
		//output[ image_index ] = image_b; 
//		printf("number of refl: %d %d\n", num_reflections( image_a ), num_reflections( image_b ) );
	}
        printf("number of crystals in twin-1: %d \n", n_twin1);
}

void compare_sets_against_one( RefList **image_array, RefList* full_list, const SymOpList *sym,  int *twin_expected, int n_image)
{
        int image_index;
        int winner;
        RefList *image_a;
        double cc[ N_TWINS ];
        int ii;

        for( image_index=0;image_index<n_image; image_index++)
        {
                image_a = image_array[ image_index ];
                stat_pearson_i_p63(image_a, full_list, cc, sym);
		if(  index_of_max_value( cc, N_TWINS) == 0 )
			printf("%d %f %f \n", image_index, cc[0], cc[1]);
		else
			printf("%d %f %f \n", image_index, cc[1], cc[0]);
        }
}

void get_twin_info_from_reference( RefList **image_array, RefList* full_list, const SymOpList *sym,  int *twin_expected, int n_image)
{
	int image_index;
	int winner;
	RefList *image_a;
	double cc[ N_TWINS ];
	int ii;

	for( image_index=0;image_index<n_image; image_index++)
	{
		image_a = image_array[ image_index ];
		//correlation_between_two_sets( image_a, full_list, sym, hkls, x, y, cc );
		stat_pearson_i_p63(image_a, full_list, cc, sym);	
		winner = index_of_max_value( cc, N_TWINS);
		twin_expected[ image_index ] = winner;
	}	
}


RefList * emc( RefList **image_array, RefList* full_list, const SymOpList *sym,  int *twin_expected, int n_image, int max_n_iter, RefList* reference, bool WINNER_TAKES_ALL)
{
	int image_index;
	int winner, n_iter;
//	double x[ MAX_REFL_PER_IMAGE ];
//	double y[ MAX_REFL_PER_IMAGE ];
//	int hkls[3][ MAX_REFL_PER_IMAGE ];
	RefList *image_a;
	RefList *this_full_list;
	double cc[ N_TWINS ];
	double mean_cc0, mean_cc1, cc_2_ref;
	int counter[ N_TWINS ], ii;
	int n_correct;

	n_iter = 0;

	while(n_iter<max_n_iter)
	{
	mean_cc0 = 0.0;
	mean_cc1 = 0.0;
	this_full_list = reflist_new();
	n_correct = 0;
	for( ii=0; ii<N_TWINS; ii++ ) counter[ii] = 0;

	for( image_index=0;image_index<n_image; image_index++)
	{
		image_a = image_array[ image_index ];
		//correlation_between_two_sets( image_a, full_list, sym, hkls, x, y, cc );
		stat_pearson_i_p63(image_a, full_list, cc, sym);	
		/* winner takes all *
		 * find the largest cc, and assign the orientation to it
		 * add this to the diffraction volume for this round
		 */
		winner = index_of_max_value( cc, N_TWINS);
        //	printf("run %d image %d winner %d \n",n_iter, image_index, winner);

		mean_cc0 += cc[winner];
		mean_cc1 += cc[abs(winner-1)];
		counter[ winner ]++;
		if( twin_expected[ image_index ] == winner ) n_correct++;
//		STATUS("winner is %d for image %d\n", winner, image_index);
		if ( WINNER_TAKES_ALL )
			merge_image_at_winner_orientation( this_full_list, image_a, winner, 1 );
		else
			merge_image_to_model( this_full_list, image_a, cc, n_iter );
//		printf("number of reflections: %d\n", num_reflections( this_full_list ) );
//		free(cc);
		
	}	
	set_esd_for_reflist( this_full_list );
	stat_pearson_i_p63(reference, full_list, cc, sym);	
	reflist_free( full_list );
	full_list = this_full_list;
//	full_list=copy_reflist( this_full_list );
//	reflist_free( this_full_list );
	
//	cc_2_ref = stat_pearson_i( reference, full_list );

	mean_cc0 /= n_image;
	mean_cc1 /= n_image;
	//printf( "wtf: %f %d\n", mean_cc, n_image );
	printf( "mean cc, cc2ref @ iteration %d := (%f %f) (%f or %f) \n", n_iter, mean_cc0, mean_cc1, cc[0], cc[1]);
	stat_pearson_i_p63(full_list, full_list, cc, sym);	
	printf( "self cc iteration %d := %f %f or %f \n", n_iter, mean_cc0, cc[0], cc[1]);
	printf( "number of orientations" );
	for(ii=0;ii<N_TWINS;ii++) printf(" %d ", counter[ii]);
	printf( "  Correct: %d \n", n_correct );

//	if (mean_cc > 0.9 ) break;
	n_iter++;
	}	
	return full_list;
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
    char *ref_hkl = NULL;
	Stream *st;
	RefList *model;
	int config_maxonly = 0;
	int config_sum = 0;
	int config_scale = 0;
	char *sym_str = NULL;
	SymOpList *sym;
    SymOpList *tlaw;    ///Sbasu
    char *operator = NULL;  //// Sbasu
	char *histo = NULL;
	signed int hist_h, hist_k, hist_l;
	//signed int hist_nbins=50;
	//float hist_min=0.0, hist_max=0.0;
	double *hist_vals = NULL;
	int hist_i;
	//int space_for_hist = 0;
	char *histo_params = NULL;
	int config_nopolar = 0;
	char *rval;
	int min_measurements = 2;
	//int r;
	int start_after = 0;
	int stop_after = 0;
	double min_snr = -INFINITY;
        RefList** image_array=NULL;
	int max_n_iter=30;
	RefList *full_list;
	int n_crystals_recorded=0;
	bool WINNER_TAKES_ALL = false;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"max-only",           0, &config_maxonly,     1},
		{"output-every",       1, NULL,               'e'},
		{"start-after",        1, NULL,              's'},
		{"stop-after",         1, NULL,                'f'},
		{"max-niter",          1, NULL,                'm'},
		{"winner-takes-all",   0, NULL,                'w'},
		{"sum",                0, &config_sum,         1},
		{"scale",              0, &config_scale,       1},
		{"no-polarisation",    0, &config_nopolar,     1},
		{"no-polarization",    0, &config_nopolar,     1},
		{"symmetry",           1, NULL,               'y'},
		{"histogram",          1, NULL,               'g'},
		{"hist-parameters",    1, NULL,               'z'},
		{"min-measurements",   1, NULL,                2},
		{"min-snr",            1, NULL,                3},
        {"reference",          1, NULL,               'r'},
        {"operator",           1, NULL,                4}, /////Sb
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:o:y:g:s:f:m:w:z:r:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			filename = strdup(optarg);
			break;

			case 'o' :
			output = strdup(optarg);
			break;

			case 's' :
			errno = 0;
			start_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --start-after.\n");
				return 1;
			}
			break;

			case 'f' :
			errno = 0;
			stop_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --stop-after.\n");
				return 1;
			}
			break;

			case 'm' :
			max_n_iter = strtod( optarg, &rval);
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'g' :
			histo = strdup(optarg);
			break;

			case 'w' :
			WINNER_TAKES_ALL = true;
			STATUS("CONGRATULATIONS: winner takes all.\n");
			break;

			case 'z' :
			histo_params = strdup(optarg);
			break;
            
            case 'r' :
            ref_hkl = strdup(optarg);
            break;
                

			case 2 :
			errno = 0;
			min_measurements = strtol(optarg, &rval, 10);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-measurements.\n");
				return 1;
			}
			break;

			case 3 :
			errno = 0;
			min_snr = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-snr.\n");
				return 1;
			}
			ERROR("WARNING: You have used --min-snr.\n");
			ERROR("WARNING: Please read the manual carefully to "
			      "learn about possible detrimental effects of this"
			      " option.\n");
			break;
                
            case 4 :
            operator = strdup(optarg);
            break;
                
			case '?' :
			break;

			case 0 :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( filename == NULL ) {
		ERROR("Please specify filename using the -i option\n");
		return 1;
	}

	if ( output == NULL ) {
		output = strdup("processed.hkl");
	}

	if ( sym_str == NULL ) sym_str = strdup("1");
	sym = get_pointgroup(sym_str);
	free(sym_str);
    
    if ( operator ) {
        tlaw = parse_symmetry_operations(operator);
        if ( tlaw == NULL ) return 1;
        set_symmetry_name(tlaw, "Ambiguity");
    }  /// Sbasu

	/* Open the data stream */
	st = open_stream_for_read(filename);
	if ( st == NULL ) {
		ERROR("Failed to open stream.\n");
		return 1;
	}

	model = reflist_new();

	if ( histo != NULL ) {

		int r;

		r = sscanf(histo, "%i,%i,%i", &hist_h, &hist_k, &hist_l);
		if ( r != 3 ) {
			ERROR("Invalid indices for '--histogram'\n");
			return 1;
		}
		free(histo);

		/* Allocate enough space that hist_vals isn't NULL.
		 * check_hist_size will realloc it straight away */
		hist_vals = malloc(1*sizeof(double));
		STATUS("Histogramming %i %i %i -> ", hist_h, hist_k, hist_l);

		/* Put into the asymmetric cell for the target group */
		get_asymm(sym, hist_h, hist_k, hist_l,
		          &hist_h, &hist_k, &hist_l);
		STATUS("%i %i %i\n", hist_h, hist_k, hist_l);

	}

	hist_i = 0;
	image_array = add_all(st, model, NULL, sym, &hist_vals, hist_h, hist_k, hist_l,
	              &hist_i, config_nopolar, min_measurements, min_snr,
	              start_after, stop_after, &n_crystals_recorded);

	fprintf(stderr, "\n");
	if ( image_array==NULL ) {
		ERROR("Error while reading stream.\n");
		return 1;
	}

	if ( config_scale ) {
		RefList **scaled_images;
		RefList *reference;

		if ( rewind_stream(st) ) {

			ERROR("Couldn't rewind stream - scaling cannot be "
			      "performed.\n");

		} else {

			STATUS("Extra pass for scaling...\n");

			reference = model;
			model = reflist_new();

			free(hist_vals);
			hist_vals = malloc(1*sizeof(double));
			hist_i = 0;

			scaled_images = add_all(st, model, reference, sym,
				     &hist_vals, hist_h, hist_k, hist_l, &hist_i,
				     config_nopolar, min_measurements, min_snr,
				     start_after, stop_after, &n_crystals_recorded);
			fprintf(stderr, "\n");
			if ( scaled_images == NULL ) {
				ERROR("Error while reading stream.\n");
				return 1;
			}

			reflist_free(reference);

		}
		int ii=0;
		for(ii=0;ii<n_crystals_recorded;ii++) {
			reflist_free( image_array[ii] );
		}
		image_array = scaled_images;

	}

	int ii;

/*
	char out_filename[80];
	for(ii=0;ii<n_crystals_recorded;ii++)
	{
		sprintf(out_filename, "image-%d.hkl",ii );
		write_reflist(out_filename, image_array[ii]);
	}
*/

	RefList *reference;
    if ( ref_hkl == NULL ) {
        printf("%s", "i'm not gonna compare");
    }
	reference = read_reflections( ref_hkl );
	write_reflist(output, model);

	close_stream(st);

	free(output);
	free(filename);

	bool random_intensity = false;
	full_list = make_reflections_for_uc_from_asymm( model, random_intensity );
//	full_list = copy_reflist( model );
	RefList *uc_list;
	random_intensity = false;
	uc_list = make_reflections_for_uc_from_asymm( reference, random_intensity );
	double cc[N_TWINS];
	stat_pearson_i_p63(uc_list, uc_list, cc, sym);	
	printf( "Full-Full self cc refence:= %f or %f \n", cc[0], cc[1]);
	stat_pearson_i_p63(uc_list, reference, cc, sym);	
	reflist_free( uc_list );
	printf( "Reference-Full self cc refence:= %f or %f \n", cc[0], cc[1]);
	stat_pearson_i_p63(reference, reference, cc, sym);	
	printf( "Reference-refence:= %f or %f \n", cc[0], cc[1]);
	stat_pearson_i_p63(model, reference, cc, sym);	
	printf( "model-ref cc refence:= %f or %f \n", cc[0], cc[1]);

	double cc_2_ref;
	cc_2_ref = stat_pearson_i( reference, model);
	printf( "cc2ref before emc := %f \n", cc_2_ref);

	int h,k,l;
	h = 0;
	k = 0;
	l = 14;
//	compute_partiality_for_it( image_array, reference, n_crystals_recorded, h, k, l);
//	compute_partiality( image_array, reference, n_crystals_recorded);
//	return 0;

	int twin_extracted[ n_crystals_recorded ];


	int twin[ n_crystals_recorded ];
	char processed[100];
	for(ii=0;ii< n_crystals_recorded; ii++) twin[ii] = 0;
        get_twin_info_from_reference( image_array, reference, sym, twin, n_crystals_recorded );
	printf("\nnow, let's analyze data\n");
	for(ii=0;ii<10;ii++)
	{
	printf("%run %d \n",ii);
	reflist_free( full_list ); full_list=reflist_new();
	full_list = make_reflections_for_uc_from_asymm( model, true );
	full_list = emc( image_array, full_list, sym, twin, n_crystals_recorded, max_n_iter, reference, WINNER_TAKES_ALL );
//	compare_sets_against_one( image_array, image_array[0], sym, twin, n_crystals_recorded );
//	compare_sets_against_one( image_array, full_list, sym, twin, n_crystals_recorded );
	sprintf(processed, "post_process.hkl%d",ii);
	write_reflist_2(processed,full_list, sym);
	}

	return 0;
	printf("\nnow, let's compare ideally controlled data\n");
/* *** Compile my simulated crystals, easy to control twin & partiality  ***/
        RefList* my_image_array[ n_crystals_recorded ];
	compile_test_data_sets( image_array, my_image_array, reference, twin, n_crystals_recorded );
/* end of my simulations */

	reflist_free( full_list ); full_list=reflist_new();
	full_list = make_reflections_for_uc_from_asymm( model, true);
	full_list = emc( my_image_array, full_list, sym, twin, n_crystals_recorded, max_n_iter, reference, WINNER_TAKES_ALL );
	/* clean the list, release the dogs */

//	write_reflist_2("post_process.hkl", full_list, sym);
	for(ii=0;ii<n_crystals_recorded;ii++) {
		reflist_free( image_array[ii] );
//		reflist_free( my_image_array[ii] );
	}

	reflist_free(model);
	free_symoplist(sym);
	

	return 0;
}
