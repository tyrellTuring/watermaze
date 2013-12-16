function DATA = mwmkld(DATA,platforms,varargin); 

% function DATA{} = mwmkld(DATA{},PLATFORMS,varargin); 
%
% Estimates spatial Kullback-Leibler divergence values between animals' search paths and the
% distribution of a set of platforms. The "true" distribution (i.e. the platform distribution) is
% estimated using a Gaussian kernel density estimation technique (see below). PLATFORMS can be a
% cell array, in which case the KLD is calculated for each set of platforms in the cell array.
% 
% IMPORTANT: This function can only be run *after* mwmpdf.m (see help mwmpdf).
%
% The function utilizes Alexander Ihler's KDE toolbox for Matlab, which can be found at
% http://www.ics.uci.edu/~ihler/code/kde.html. If a bandwidth for the kernels is not specified then
% the KDE's optimal bandwidth calculation tool is used with the default spherical, uniform
% likelihood cross-validation search method. See 'help ksize' for more information.
%
% Note that the KDEs for the platforms as well as their estimated PDF are also stored, along with
% the actual divergence value, in the sub-structure KDE. All of the results are stored in the
% structure KLD. If multiple platform sets are given then KLD is a cell array, one structure for
% each set.
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'bandwidth' : Scalar valued spatial size of kernels for the platform KDE. Default = optimal 
%                 estimation.
%
% - 'logbase': Which logarithmic base to use for the calculation (can be 2, e, or 10). Default = e (nats).
%
% - 'weights': A vector of weights for each platforms contribution to the KDE. Default = 1 for all.
%
% - 'means'  : A cell array of vectors that can be used to calculate the KLD using the mean of
%              multiple trials' pdfs. For example, {[1 2],[4 5]} would lead to KLD being calculated
%              for pdfs averaged over the first and second, then the fourth and fifth trials. These
%              KLD values are stored in the sub-vector mkld in the order they were requested. 
%
% - 'zones' : A matrix set of zones for which the probability of each zone will be calculated for each 
%             set of platform pdfs. The zones should be an M x 3 matrix where M is the number of zones,
%             the first column is X locations, the second column is Y locations and the third column
%             is the radius of each zone. 
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

% check the DATA argument
if ~isa(DATA,'cell')
	error('DATA should be a cell array');
end

% check the platforms argument
if ~(isa(platforms,'numeric') && size(platforms,2) == 3) && ~isa(platforms,'cell')
	error('PLATFORMS should be a K x 3 matrix or a cell array');
end
	
% define the default optional arguments
optargs = struct('bandwidth',-1,...
                 'logbase',exp(1),...
                 'weights',ones(1,size(platforms,1)),...
                 'means',{{}},...
								 'zones',[0,0,0]);

% get the optional argument names
optnames = fieldnames(optargs);

% get the number of optional arguments
nargs = length(varargin)/2;

% make sure the property/value pairs are input in pairs
if round(length(varargin)/2) ~= nargs
   error('Expecting propertyName/propertyValue pairs after platforms');
end

% step through the optional arguments, check them, and store them
for pair = reshape(varargin,2,[])

	% make it case insensitive
	inpname = lower(pair{1});

	% check whether the name matches a known option
	if any(strmatch(inpname,optnames))
		switch inpname
			case 'bandwidth'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('bandwidth must be a scalar value');
				end
			case 'logbase'
				if isa(pair{2},'numeric') && ismember(pair{2},[2 exp(1) 10])
					optargs.(inpname) = pair{2};
				else
					error('logbase must be 2, e or 10');
				end
			case 'weights'
				if isa(pair{2},'numeric') && length(pair{2}) == size(platforms,1)
					optargs.(inpname) = pair{2};
				else
					error('weights must be a vector with length = size(PLATFORMS,1)');
				end
			case 'means'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('means must be a cell array');
				end
			case 'zones'
				if isa(pair{2},'numeric' && size(pair{2},2) == 3)
					optargs.(inpname) = pair{2};
				else
					error('zones must be an M x 3 matrix');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

% fix weights if necessary
if isa(platforms,'cell') && optargs.weights == ones(size(platforms,1))
	optargs.weights = [];
end

%  for each data set in the structure
for ff = 1:length(DATA)

	% get the points we're going to sample
	XX = DATA{ff}.PDF.x;
	YY = DATA{ff}.PDF.y;
	outofpool = isnan(XX);
	X = XX; X(outofpool) = 0;
	Y = YY; Y(outofpool) = 0;

	% check whether we have multiple platform sets
	if isa(platforms,'cell')

		% for each set of platforms...
		for pp = 1:length(platforms)

			% initialize the KLD{pp} structure
			DATA{ff}.KLD{pp}.kld = zeros(size(DATA{ff}.PDF.p));
			DATA{ff}.KLD{pp}.p   = zeros(size(DATA{ff}.PDF.p,1),size(DATA{ff}.PDF.p,2));

			% get the KDE
			if optargs.bandwidth == -1

				% create a KDE with an arbitrary bandwidth
				tmpkde = kde(platforms{pp}(:,1:2)', [1; 1],optargs.weights);

				% calculate the optimal bandwidth
				DATA{ff}.KLD{pp}.kde = ksize(tmpkde);
			else

				% just use the specified bandwidth
				DATA{ff}.KLD{pp}.kde = kde(platforms{pp}(:,1:2)', [optargs.bandwidth;optargs.bandwidth],optargs.weights);
			end

			% use the KDE to estimate the PDF of the platforms across the pool
			P = reshape(evaluate(DATA{ff}.KLD{pp}.kde,[X(:)';Y(:)']),size(DATA{ff}.PDF.p,1),size(DATA{ff}.PDF.p,1));

			% limit the data to the actual area of the pool
			P(outofpool) = nan;

			% renormalize P
			P = P./nansum(P(:)*DATA{ff}.PDF.res^2);

			% store the results
			DATA{ff}.KLD{pp}.p = P;

			% calculate the zone proabilities if requested
			if any(optargs.zones)
				DATA{ff}.KLD{pp}.zones = zeros(size(optargs.zones,1),1);
				for zz = 1:size(optargs.zones,1)
					DATA{ff}.KLD{pp}.zones(zz) = nansum(nansum(DATA{ff}.KLD{pp}.p.*sqrt((optargs.zones(zz,1)-XX).^2 + (optargs.zones(zz,2)-YY).^2) < optargs.zones(zz,3)));
				end
			end

			% for each trial, calcualte the KLD{pp}
			DATA{ff}.KLD{pp}.kld = zeros(DATA{ff}.ntrials,1);
			for tt = 1:DATA{ff}.ntrials
				switch optargs.logbase
					case 2
						DATA{ff}.KLD{pp}.kld(tt) = nansum(nansum(DATA{ff}.KLD{pp}.p.*log2(DATA{ff}.KLD{pp}.p./(DATA{ff}.PDF.p(:,:,tt)+exp(-100))))); 
					case exp(1)
						DATA{ff}.KLD{pp}.kld(tt) = nansum(nansum(DATA{ff}.KLD{pp}.p.*log(DATA{ff}.KLD{pp}.p./(DATA{ff}.PDF.p(:,:,tt)+exp(-100)))));
					case 10
						DATA{ff}.KLD{pp}.kld(tt) = nansum(nansum(DATA{ff}.KLD{pp}.p.*log10(DATA{ff}.KLD{pp}.p./(DATA{ff}.PDF.p(:,:,tt)+exp(-100)))));
				end
			end

			% for each set of means, calcualte the KLD{pp}
			DATA{ff}.KLD{pp}.mkld = zeros(length(optargs.means),1);
			for mm = 1:length(optargs.means)
			
				% get the mean pdf for these trials
				allP = zeros(size(DATA{ff}.PDF.p,1),size(DATA{ff}.PDF.p,2),length(optargs.means));
				for tt = 1:length(optargs.means{mm})
					allP(:,:,tt) = DATA{ff}.PDF.p(:,:,optargs.means{mm}(tt));
				end
				P = nanmean(allP,3);

				switch optargs.logbase
					case 2
						DATA{ff}.KLD{pp}.mkld(mm) = nansum(nansum(DATA{ff}.KLD{pp}.p.*log2(DATA{ff}.KLD{pp}.p./(P+exp(-100))))); 
					case exp(1)
						DATA{ff}.KLD{pp}.mkld(mm) = nansum(nansum(DATA{ff}.KLD{pp}.p.*log(DATA{ff}.KLD{pp}.p./(P+exp(-100)))));
					case 10
						DATA{ff}.KLD{pp}.mkld(mm) = nansum(nansum(DATA{ff}.KLD{pp}.p.*log10(DATA{ff}.KLD{pp}.p./(P+exp(-100)))));
				end
			end
		end
	else

		% initialize the KLD structure
		DATA{ff}.KLD.kld = zeros(size(DATA{ff}.PDF.p));
		DATA{ff}.KLD.p   = zeros(size(DATA{ff}.PDF.p,1),size(DATA{ff}.PDF.p,2));

		% get the KDE
		if optargs.bandwidth == -1

			% create a KDE with an arbitrary bandwidth
			tmpkde = kde(platforms(:,1:2)', [1; 1],optargs.weights);

			% calculate the optimal bandwidth
			DATA{ff}.KLD.kde = ksize(tmpkde);
		else

			% just use the specified bandwidth
			DATA{ff}.KLD.kde = kde(platforms(:,1:2)', [optargs.bandwidth;optargs.bandwidth],optargs.weights);
		end

		% use the KDE to estimate the PDF of the platforms across the pool
		P = reshape(evaluate(DATA{ff}.KLD.kde,[X(:)';Y(:)']),size(DATA{ff}.PDF.p,1),size(DATA{ff}.PDF.p,1));

		% limit the data to the actual area of the pool
		P(outofpool) = nan;

		% renormalize P
		P = P./nansum(P(:)*DATA{ff}.PDF.res^2);

		% store the results
		DATA{ff}.KLD.p = P;

		% calculate the zone proabilities if requested
		if any(optargs.zones)
			DATA{ff}.KLD.zones = zeros(size(optargs.zones,1),1);
			for zz = 1:size(optargs.zones,1)
				DATA{ff}.KLD.zones(zz) = nansum(nansum(DATA{ff}.KLD.p.*sqrt((optargs.zones(zz,1)-XX).^2 + (optargs.zones(zz,2)-YY).^2) < optargs.zones(zz,3)));
			end
		end

		% for each trial, calcualte the KLD
		DATA{ff}.KLD.kld = zeros(DATA{ff}.ntrials,1);
		for tt = 1:DATA{ff}.ntrials
			switch optargs.logbase
				case 2
					DATA{ff}.KLD.kld(tt) = nansum(nansum(DATA{ff}.KLD.p.*log2(DATA{ff}.KLD.p./(DATA{ff}.PDF.p(:,:,tt)+exp(-100))))); 
				case exp(1)
					DATA{ff}.KLD.kld(tt) = nansum(nansum(DATA{ff}.KLD.p.*log(DATA{ff}.KLD.p./(DATA{ff}.PDF.p(:,:,tt)+exp(-100)))));
				case 10
					DATA{ff}.KLD.kld(tt) = nansum(nansum(DATA{ff}.KLD.p.*log10(DATA{ff}.KLD.p./(DATA{ff}.PDF.p(:,:,tt)+exp(-100)))));
			end
		end

		% for each set of means, calcualte the KLD
		DATA{ff}.KLD.mkld = zeros(length(optargs.means),1);
		for mm = 1:length(optargs.means)
		
			% get the mean pdf for these trials
			allP = zeros(size(DATA{ff}.PDF.p,1),size(DATA{ff}.PDF.p,2),length(optargs.means));
			for tt = 1:length(optargs.means{mm})
				allP(:,:,tt) = DATA{ff}.PDF.p(:,:,optargs.means{mm}(tt));
			end
			P = nanmean(allP,3);

			switch optargs.logbase
				case 2
					DATA{ff}.KLD.mkld(mm) = nansum(nansum(DATA{ff}.KLD.p.*log2(DATA{ff}.KLD.p./(P+exp(-100))))); 
				case exp(1)
					DATA{ff}.KLD.mkld(mm) = nansum(nansum(DATA{ff}.KLD.p.*log(DATA{ff}.KLD.p./(P+exp(-100)))));
				case 10
					DATA{ff}.KLD.mkld(mm) = nansum(nansum(DATA{ff}.KLD.p.*log10(DATA{ff}.KLD.p./(P+exp(-100)))));
			end
		end
	end
end

