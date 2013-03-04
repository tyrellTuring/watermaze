function DATA = mwmpdf(DATA,varargin); 

% function DATA = mwmpdf(DATA,varargin); 
%
% Estimates spatial probability density functions from animals' search paths, using gaussian kernel 
% density estimates based on the moment-to-moment position of the animal. The obligatory 
% input structure, DATA, is a multi-level cell array that is assumed to be in the format 
% returned by readwmdf.m (see help readwmdf). The results are stored in a structure at the 
% second-level of DATA, as PDF. PDF contains the x and y coordinates where the pdf was analyzed and
% a X x Y x T matrix, p, where X and Y are the maximum number of points across the pool
% at which the pdf is evaluated (this depends on the 'resolution' argument passed, see below), T is
% the number of trials.
%
% The function utilizes Alexander Ihler's KDE toolbox for Matlab, which can be found at
% http://www.ics.uci.edu/~ihler/code/kde.html. If a bandwidth for the kernels is not specified then
% the KDE's optimal bandwidth calculation tool is used with the default spherical, uniform
% likelihood cross-validation search method. See 'help ksize' for more information. If the optimal
% bandiwdth method is used, the user can also specify to utilize the mean of the optimal bandwidths
% across all the paths, using the 'meanbw' option (see below).
%
% Note that the KDEs for each trial are also stored in the cell array PDF.kde.
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'resolution': Scalar valued spatial resolution of pdf evaluation. Default = 1 cm.
%
% - 'bandwidth' : Scalar valued spatial size of kernels for KDEs. Default = optimal estimation.
%
% - 'meanbw'    : Boolean value, if set to true and bandidth not specified, leads to mean optimal
%                 bandwidth of all data being used.
%
% 02/2013, Frankland Lab (www.franklandlab.com)
%
% Author: Blake Richards
% Contact: blake.richards@utoronto.ca
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% check the DATA argument
if ~isa(DATA,'cell')
	error('DATA should be a cell array');
end
	
% define the default optional arguments
optargs = struct('resolution',1,'bandwidth',-1,'meanbw',false);

% get the optional argument names
optnames = fieldnames(optargs);

% get the number of optional arguments
nargs = length(varargin)/2;

% make sure the property/value pairs are input in pairs
if round(length(varargin)/2) ~= nargs
   error('Expecting propertyName/propertyValue pairs after DATA');
end

% step through the optional arguments, check them, and store them
for pair = reshape(varargin,2,[])

	% make it case insensitive
	inpname = lower(pair{1});

	% check whether the name matches a known option
	if any(strmatch(inpname,optnames))
		switch inpname
			case 'resolution'
        if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
        else
					error('lambda must be a scalar value');
				end
			case 'bandwidth'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('bandwidth must be a scalar value');
				end
			case 'meanbw'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('meanbw must be a boolean value');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%  for each data set in the structure
for ff = 1:length(DATA)

	% get the points we're going to sample
	X = [-DATA{ff}.pool(3):optargs.resolution:DATA{ff}.pool(3)];
	Y = [-DATA{ff}.pool(3):optargs.resolution:DATA{ff}.pool(3)];
	[XX YY] = meshgrid(X,fliplr(Y)); 
	outofpool = sqrt(XX.^2 + YY.^2) > DATA{ff}.pool(3);
	nanXX = XX; nanXX(outofpool) = nan;
	nanYY = YY; nanYY(outofpool) = nan;

	% determine the number of points across the pool
	npoints = length(X);

	% initialize the PDF structure
	DATA{ff}.PDF.p = zeros(npoints,npoints,DATA{ff}.ntrials);
	DATA{ff}.PDF.x = nanXX;
	DATA{ff}.PDF.y = nanYY;
	DATA{ff}.PDF.res = optargs.resolution;

	% for each trial get the KDE
	for tt = 1:DATA{ff}.ntrials

		% if requested, estimate using the optimal bandwidth
		if optargs.bandwidth == -1

			% create a KDE with an arbitrary bandwidth
			tmpkde = kde(DATA{ff}.path(:,2:3,tt)', [1; 1]);

			% calculate the optimal bandwidth
			DATA{ff}.PDF.kde{tt} = ksize(tmpkde);
		else

			% just use the specified bandwidth
			DATA{ff}.PDF.kde{tt} = kde(DATA{ff}.path(:,2:3,tt)', [optargs.bandwidth;optargs.bandwidth]);
		end
	end

	% if the mean bandwidth of the KDEs was requested, then calculate it now and redo all the KDEs
	if optargs.meanbw
		sumBW = 0;
		for tt = 1:DATA{ff}.ntrials
			sumBW = sumBW + nanmean(nanmean(getBW(DATA{ff}.PDF.kde{tt})));
		end
		meanBW = sumBW/DATA{ff}.ntrials;
		for tt = 1:DATA{ff}.ntrials
			DATA{ff}.PDF.kde{tt} =  kde(DATA{ff}.path(:,2:3,tt)', [meanBW; meanBW]);
		end
	end

	% use all the KDEs to estimate the PDFs across the pool
	for tt = 1:DATA{ff}.ntrials

		% get the evaluation using the KDE
		P = reshape(evaluate(DATA{ff}.PDF.kde{tt},[XX(:)';YY(:)']),length(X),length(Y));

		% limit the data to the actual area of the pool
		P(outofpool) = nan;

		% renormalize P
		P = P./nansum(P(:)*optargs.resolution^2);

		% store the results
		DATA{ff}.PDF.p(:,:,tt) = P;
	end
end

