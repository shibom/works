function [ cc ] = compare_self( folder, pixel_positions, nfiles, class )
%   check the convergence of the most recent saved correlation maps at pixel positions
%   input args include:
%       folder where the correlation files are saved
%       pixel positions 
%   i.e:
%       compare_self('.',[20 50 100 150], 5, 'class0')
%   will find the most recently saved correlations for class0 (or class1 if specified in class)
%   the correlation functions at pixel pisitions [20 50 100 150] will be
%   compared at 5 latest files

if (nargin < 1)
    fprintf('default path (current folder) is used \n');
    folder='./';
end

if (nargin<=1)
   pixel_positions=[30 50 100 200];
end

if (nargin<=2)
   nfiles = 5;
end

if (nargin<=3)
    class='class1'
end

[files, nfiles] = recentfile(folder, ['*angular*' class '*h5'], nfiles);

fprintf('%s most updated file:\n\t %s\n',class,files{1} );

data = [];
max_x = [];
colors = ['b','r','g','k','c','y','m'];

for ii=1:nfiles
    data{ii} = hdf5read(files{ii},'/data/data');
    max_x(ii) = size(data{ii},1);
end


ccs=zeros(4,1);
figure;
for ii=1:4
    subplot(2,2,ii);
    hold on;
    for jj=1:nfiles
        c0=data{jj};
        this_max0 = max(c0(:,pixel_positions(ii)));
        c0(:,pixel_positions(ii) ) = c0(:,pixel_positions(ii) )/this_max0;
        if(jj<8)
            plot(c0(:,pixel_positions(ii)),'color',colors(jj));
        else
            plot(c0(:,pixel_positions(ii)),'color');            
        end
    end
    hold off;
    title(['at pixel ' int2str(pixel_positions(ii))] ,'fontsize',20);
end

text(1.2,1.2,'colors:b r g k c y m etc');

end

function [filenames, nfiles]=recentfile( path,expression, nfiles )
d = dir([path '/' expression]);
[dx dx] = sort([d.datenum],'descend');
filenames = {};
if(size(d,1) < nfiles)
    nfiles=size(d,1);
end

for ii=1:nfiles
    filenames{ii} = d(dx(ii)).name;
end
end