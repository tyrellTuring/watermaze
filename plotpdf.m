function [sh P] = plotpdf(DATA,trials,animals,varargin)

% function [H P] = plotpdf(DATA,TRIALS,ANIMALS,varargin)
%
% Plots the probability density functions for a set of trials and animals in the given data set.
%
% MANDATORY INPUTS:
% -------------------------------------------------------------------------------------------------
%
%   DATA    - DATA is a cell array that contains the raw water-maze data for a given collection of 
%             project files in the study. See 'help readwmdf' for details.
%
%   TRIALS  - A vector of indices of which trials to include in the plot. 
%
%   ANIMALS - A vector of indices of which animals to include in the plot.
%
% OUTPUT:
% -------------------------------------------------------------------------------------------------
%
%   H - A handle to the plot (a surface plot).
%
%   P - The PDF that is shown.
%
% OPTIONAL INPUTS:
% -------------------------------------------------------------------------------------------------
% Optional inputs to control plotting behaviour can be provided in parameter/value format. The optional inputs are:
%
%   'show_kldpdf' - Boolean value, determines whether to show the pdf used for the KLD calculation
%                   rather than animal's paths. Default = false.
%
<<<<<<< HEAD
%   'kldnum'    - Integer indicating which KLD struct to use when plotting the KLD pdf. Default = 1.
=======
%   'kld_num'   - To be used with 'show_kldpdf', determines which number of kld pdf to plot.
>>>>>>> a1d51b7f65ffc9ad8ce9c0194a8be3283308ed13
%
%   'show_axis' - Boolean value, determines whether the spatial axis is shown. Default = false.
%
%   'measure'   - String, defining the measure to use for combining multiple PDFs. Options are:
%                 {'mean','max','min','var'}. Default = 'mean'.
%
%   'platforms' - A K x 3 matrix of platforms to draw on the plot, where K is the number of
%                 platforms, and the columns are x-centre, y-centre, radius. Default = [];
%
%   'plat_opacity' - A scalar value indicating the opacity of the platforms. Default = 0.3.
%
%   'plat_colour' - A K x 3 matrix indicating the colors of each of the platforms. Default is all
%                   white.
%
%   'crosshair' - Boolean value, determines if cross-hair is drawn on plot. Default = false.
%
%   'cmap'      - The colour map to use for the plot. Default = 'jet'.
%
%   'climit'    - Two-element vector, [CMIN, CMAX] specifying the colour axis limits. Default =
%                 [0 max(PDF)].
%
%--------------------------------------------------------------------------------
%
% 02/2013, Frankland Lab (www.franklandlab.com)
%
% Author: Blake Richards
% Contact: blake.richards@utoronto.ca
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2013 Blake Richards (blake.richards@utoronto.ca)
%
% This file is part of the MWM Matlab Toolbox.
%
% The MWM Toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The MWM Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with the MWM Toolbox (in the file COPYING.LESSER).  If not, 
% see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% define the default optional arguments
optargs = struct('show_axis',false,...
                 'show_kldpdf',false,...
<<<<<<< HEAD
                 'kldnum',1,...
=======
                 'kld_num',1,...
>>>>>>> a1d51b7f65ffc9ad8ce9c0194a8be3283308ed13
                 'measure','mean',...
                 'platforms',[],...
                 'plat_opacity',0.3,...
                 'plat_colour',[],...
                 'crosshair',false,...
                 'cmap','jet',...
                 'climit',NaN);

% get the optional argument names
optnames = fieldnames(optargs);

% get the number of optional arguments
nargs = length(varargin)/2;

% make sure the property/value pairs are input in pairs
if round(length(varargin)/2) ~= nargs
   error('Expecting propertyName/propertyValue pairs after DATA, TRIALS and ANIMALS');
end

% step through the optional arguments, check them, and store them
for pair = reshape(varargin,2,[])

	% make it case insensitive
	inpname = lower(pair{1});

	% check whether the name matches a known option
	if any(strmatch(inpname,optnames))
		switch inpname
			case 'show_axis'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('show_axis must be a logical');
				end
			case 'show_kldpdf'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('show_kldpdf must be a logical');
				end
<<<<<<< HEAD
			case 'kldnum'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('kldnum must be an integer');
=======
			case 'kld_num'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('kld_num must be an integer');
>>>>>>> a1d51b7f65ffc9ad8ce9c0194a8be3283308ed13
				end
			case 'measure'
				if isa(pair{2},'char') && ismember(pair{2},{'mean','max','min','var'})
					optargs.(inpname) = pair{2};
				else
					error('measure must be one of {''mean'',''max'',''min'',''var''}');
				end
			case 'platforms'
				if isa(pair{2},'numeric') && size(pair{2},2) == 3
					optargs.(inpname) = pair{2};
				else
					error('platforms must be a K x 3 matrix');
				end
			case 'plat_opacity'
				if isa(pair{2},'numeric') && (pair{2} >= 0 && pair{2} <= 1)
					optargs.(inpname) = pair{2};
				else
					error('plat_opacity must be a scalar value between 0 and 1');
				end
			case 'plat_colour'
				if isa(pair{2},'numeric') && size(pair{2},2) == 3
					optargs.(inpname) = pair{2};
				else
					error('plat_colour must be a K x 3 matrix');
				end
			case 'crosshair'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('crosshair must be a Boolean value');
				end
			case 'cmap'
				if isa(pair{2},'char')
					optargs.(inpname) = pair{2};
				else
					error('cmap must be a string of the colormap name');
				end
			case 'climit'
				if isa(pair{2},numeric) && length(pair{2}) == 2
					optargs.(inpname) = pair{2};
				else
					error('climit must be a two-element vector');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THE PDF TO BE SHOWN

% get X and Y 
X = DATA{1}.PDF.x;
Y = DATA{1}.PDF.y;

% see whether we're using the PDF from the KLD
if optargs.show_kldpdf
	if isa(DATA{1}.KLD,'cell')
<<<<<<< HEAD
		P = DATA{1}.KLD{optargs.kldnum}.p;
	else
=======
		P = DATA{1}.KLD{optargs.kld_num}.p;
	else	
>>>>>>> a1d51b7f65ffc9ad8ce9c0194a8be3283308ed13
		P = DATA{1}.KLD.p;
	end
else
	% get the data from all the trials and animals
	allP = zeros(size(X,1),size(X,2),length(trials),length(animals));
	for aa = 1:length(animals)
		for tt = 1:length(trials)
			allP(:,:,tt,aa) = DATA{animals(aa)}.PDF.p(:,:,trials(tt));
		end
	end

	% take the requested measure across trials and animals
	switch optargs.measure
		case 'mean'
			P = nanmean(nanmean(allP,4),3);
		case 'max'
			P = nanmax(nanmax(allP,[],4),[],3);
		case 'min'
			P = nanmin(nanmin(allP,[],4),[],3);
		case 'var'
			varanim = nanvar(allP,[],4);
			if size(allP,3) == 1
				P = squeeze(varanim);
			else
				P = nanvar(varanim,[],3);
			end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE PDF

% calculate the color limit now if necessary
if isnan(optargs.climit)
	climit = [0 max(P(:))];
else
	climit = optargs.climit;
end

% create the figure and make it square
figure();
set(gcf,'InvertHardcopy','off');
set(gcf,'Position',[400 100 500 500],'PaperPosition',[2 2 6 6]);
set(gcf,'Color',[1 1 1]);
hold on;

% plot the PDF as a surface
sh = surf(X,Y,P);

% set the colormap
colormap(optargs.cmap);
cmap = colormap(gca);

% plot a circle of the pool edge (for a nice clean look)
circ = zeros(2001,3);
circ(:,1) = cos([0:pi/1000:2*pi])*DATA{1}.pool(3);
circ(:,2) = sin([0:pi/1000:2*pi])*DATA{1}.pool(3);
circ(:,3) = 0;
fill3(circ(:,1),circ(:,2),circ(:,3),cmap(1,:),'linewidth',2)

% do some basic formatting
set(sh,'LineStyle','none');
set(gca,'CameraPosition',[0 0 max(P(:))*1.5],'CameraTarget',[0 0 0]);
axis equal;
set(gca,'Position',[0 0 1 1]);
caxis(gca,climit);
xlim([-DATA{1}.pool(3)-5, DATA{1}.pool(3)+5])
ylim([-DATA{1}.pool(3)-5, DATA{1}.pool(3)+5])
camproj('orthographic')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMAT BASED ON THE USER'S REQUESTS

% turn off the axis if requested
if ~optargs.show_axis, axis off; end;

% draw the platforms
if ~isempty(optargs.platforms)
	
	% for each platform...
	for pp = 1:size(optargs.platforms,1)
		circ = zeros(1001,3);
		circ(:,1) = optargs.platforms(pp,1) + cos([0:pi/500:2*pi])*optargs.platforms(pp,3);
		circ(:,2) = optargs.platforms(pp,2) + sin([0:pi/500:2*pi])*optargs.platforms(pp,3);
		circ(:,3) = max(P(:))*1.1;
		ch(pp) = fill3(circ(:,1),circ(:,2),circ(:,3),[1 1 1],'linewidth',2);
		if ~isempty(optargs.plat_colour)
			set(ch(pp),'FaceColor',optargs.plat_colour(pp,:));
		else
			set(ch(pp),'FaceColor','none');
		end
		set(ch(pp),'FaceAlpha',optargs.plat_opacity,'EdgeAlpha',optargs.plat_opacity);
	end
end

% draw the crosshairs
if optargs.crosshair
	plot3([-DATA{1}.pool(3) DATA{1}.pool(3)],[0 0],[max(P(:))*1.1 max(P(:))*1.1],'w--','LineWidth',2);
	plot3([0 0],[-DATA{1}.pool(3) DATA{1}.pool(3)],[max(P(:))*1.1 max(P(:))*1.1],'w--','LineWidth',2);
end
