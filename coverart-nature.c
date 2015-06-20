/*
 * render_hkl.c
 *
 * Draw pretty renderings of reflection lists
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 */

// Modified by Shibom Basu for PSII 3D Ewald Sphere image in POV-Ray..

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <crystfel/utils.h>
#include <crystfel/symmetry.h>
#include <crystfel/render.h>
#include <crystfel/reflist.h>
#include <crystfel/reflist-utils.h>
#include <crystfel/cell.h>
#include <crystfel/cell-utils.h>


static int render_scene(UnitCell *cell, RefList *list, const SymOpList *sym,
                        double boost, double scale_top)
{
	FILE *fh;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	pid_t pid;
	int r;
	float max;
	Reflection *refl;
	RefListIterator *iter;
	SymOpMask *m;

	fh = fopen("s3.pov", "w");
	fprintf(fh, "/* POV-Ray scene */\n\n");
	fprintf(fh, "#include \"colors.inc\"\n");

	fprintf(fh, "global_settings {"
	            "  assumed_gamma 1.0"
	            "  ambient_light 3.0\n"
	            "}\n");

	fprintf(fh, "camera { location <0.0, -1.0, 0.0>"
	            " sky z direction y\n"
	            " right -x*(image_width/image_height)\n"
	            " look_at <0.0, 0.0, 0.0> }\n\n");

	//fprintf(fh, "light_source { <-60.0 -60.0 60.0> White shadowless }\n");
	//fprintf(fh, "light_source { <+60.0 -60.0 60.0> White shadowless }\n");
	fprintf(fh, "light_source { <-5.0, -10.0, 0.0> White }\n");
	fprintf(fh, "light_source { <-5.0, -10.0, 0.0> 2.0*White shadowless }\n");
	fprintf(fh, "light_source { <0.0, 10.0, 10.0> 10.0*White shadowless }\n");

	fprintf(fh, "#declare TRANS = transform { rotate <10, 0, 0> }\n");

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	max = 0.0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		float val;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		val = get_intensity(refl);
		val = (val>0.0) ? sqrt(val) : 0.0;

		if ( val > max ) max = val;

	}
	max /= boost;

	//fprintf(fh, "fog {\n");
	//fprintf(fh, "	distance 10.0\n");
	//fprintf(fh, "	color rgb<0.3, 0.3, 0.3>\n");
	//fprintf(fh, "}\n");

	/* Use manual scale top if specified */
	if ( scale_top > 0.0 ) {
		max = scale_top;
	}

	m = new_symopmask(sym);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		float radius;
		int s;
		float val, p, r, g, b;
		int j;
		int neq;

		get_indices(refl, &h, &k, &l);

		special_position(sym, m, h, k, l);
		neq = num_equivs(sym, m);

		val = get_intensity(refl);
		val = (val>0.0) ? sqrt(val) : 0.0;

		s = val / (max/6);
		p = fmod(val, max/6);
		p /= (max/6);

		r = 0;	g = 0;	b = 0;

		if ( (val < 0.0) ) {
			s = 0;
			p = 1.0;
		}
		if ( (val > max) ) {
			s = 6;
		}
		switch ( s ) {
		case 0 :   /* Black to blue */
			r = 0.0;  g = 0.0;  b = p;
			break;
		case 1 :   /* Blue to green */
			r = 0.0;  g = p;  b = 1.0-p;
			break;
		case 2 :   /* Green to red */
			r = p;  g = 1.0-p;  b = 0.0;
			break;
		case 3 :   /* Red to Orange */
			r = 1.0;  g = 0.5*p;  b = 0.0;
			break;
		case 4 :   /* Orange to Yellow */
			r = 1.0;  g = 0.5 + 0.5*p;  b = 0.0;
			break;
		case 5 :   /* Yellow to White */
			r = 1.0;  g = 1.0;  b = 1.0*p;
			break;
		case 6 :   /* Pixel has hit the maximum value */
			r = 1.0;  g = 1.0;  b = 1.0;
			break;
		}

		if ( val <= 0.0 ) continue;
		radius = 0.02 * pow(val, 1.25)/pow(max, 1.25);
		if ( radius > 0.02 ) {
			radius = 0.02;
			r = 0.0;  g = 0.0;  b = 0.4;
		}

		double res = 2.0*resolution(cell, h, k, l);
		if ( res > 0.9e9 ) continue;

		/* For each equivalent */
		for ( j=0; j<neq; j++ ) {

			signed int he, ke, le;
			float x, y, z;

			get_equiv(sym, m, j, h, k, l, &he, &ke, &le);

			x = asx*he + bsx*ke + csx*le;
			y = asy*he + bsy*ke + csy*le;
			z = asz*he + bsz*ke + csz*le;

			fprintf(fh, "/* %i %i %i */\n", h, k, l);
			fprintf(fh, "sphere { <%.5f, %.5f, %.5f>, %.5f "
			            "  texture {"
			            "    pigment {"
			            "      color rgbf <%f, %f, %f>"
			            "    } "
			            "    finish {"
			            "      ambient 0.1"
			            "      diffuse 0.3"
			            "      phong 0.7"
			            "      phong_size 10"
			            "    }"
			            "  }"
			            "  transform { TRANS }\n"
			            "}\n",
			            x/1.0e9, y/1.0e9, z/1.0e9, radius,
			            r, g, b);


		}

	}

	free_symopmask(m);

	fprintf(fh, "\n");
	fclose(fh);

	pid = fork();
	if ( !( (pid != 0) && (pid != -1) ) ) {
		if ( pid == -1 ) {
			ERROR("fork() failed.\n");
		} else {

			/* Forked successfully, child process */
			execlp("povray", "", "+W560", "+H792",
			       "+Is3.pov", "+Os3.png", "+D", "+P",
			       NULL);

		}
	}

	waitpid(pid, &r, 0);

	return 0;
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.hkl>\n\n", s);
	printf(
"Render intensity lists in raytraced 3D pictures.\n"
"\n"
"  -d, --down=<h>,<k>,<l>  Indices for the axis in the downward direction.\n"
"                           Default: 1,0,0.\n"
"  -r, --right=<h>,<k>,<l> Indices for the axis in the 'right' (roughly)\n"
"                           direction.  Default: 0,1,0.\n"
"  -o, --output=<filename> Output filename (not for POV-ray).  Default: za.pdf\n"
"      --boost=<val>       Squash colour scale by <val>.\n"
"  -p, --pdb=<file>        PDB file from which to get the unit cell.\n"
"  -y, --symmetry=<sym>    Expand reflections according to point group <sym>.\n"
"\n"
"  -c, --colscale=<scale>  Use the given colour scale.  Choose from:\n"
"                           mono    : Greyscale, black is zero.\n"
"                           invmono : Greyscale, white is zero.\n"
"                           colour  : Colour scale:\n"
"                                     black-blue-pink-red-orange-yellow-white\n"
"\n"
"  -w  --weighting=<wght>  Colour/shade the reciprocal lattice points\n"
"                           according to:\n"
"                            I      : the intensity of the reflection.\n"
"                            sqrtI  : the square root of the intensity.\n"
"                            count  : the number of measurements for the reflection.\n"
"                                     (after correcting for 'epsilon')\n"
"                            rawcts : the raw number of measurements for the\n"
"                                     reflection (no 'epsilon' correction).\n"
"\n"
"  -j <n>                  Run <n> instances of POV-ray in parallel.\n"
"  -h, --help              Display this help message.\n"
);
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	RefList *list;
	char *infile;
	char *pdb = NULL;
	int r = 0;
	char *sym_str = NULL;
	SymOpList *sym;
	char *outfile = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"pdb",                1, NULL,               'p'},
		{"symmetry",           1, NULL,               'y'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hj:p:w:c:y:d:r:o:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'p' :
			pdb = strdup(optarg);
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 0 :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( pdb == NULL ) {
		ERROR("You must specify the PDB containing the unit cell.\n");
		return 1;
	}

	if ( sym_str == NULL ) {
		sym_str = strdup("1");
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	infile = argv[optind];

	cell = load_cell_from_pdb(pdb);
	if ( cell == NULL ) {
		ERROR("Couldn't load unit cell from %s\n", pdb);
		return 1;
	}
	list = read_reflections(infile);
	if ( list == NULL ) {
		ERROR("Couldn't read file '%s'\n", infile);
		return 1;
	}
	if ( check_list_symmetry(list, sym) ) {
		ERROR("The input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(sym));
		return 1;
	}

	render_scene(cell, list, sym, 3.0, -1.0);

	free(pdb);
	free_symoplist(sym);
	reflist_free(list);
	if ( outfile != NULL ) free(outfile);

	return r;
}
