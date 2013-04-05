function DATA = mwmkld(DATA,platforms,varargin); 

% function DATA = mwmkld(DATA,PLATFORMS,varargin); 
%
% Estimates spatial Kullback-Leibler divergence values between animals' search paths and the
% distribution of a set of platforms. The "true" distribution (i.e. the platform distribution) is
% estimated using a Gaussian kernel density estimation technique (see below).
% 
% IMPORTANT: This function can only be run *after* mwmpdf.m (see help mwmpdf).
%
% The function utilizes Alexander Ihler's KDE toolbox for Matlab, which can be found at
% http://www.ics.uci.edu/~ihler/code/kde.html. If a bandwidth for the kernels is not specified then
% the KDE's optimal bandwidth calculation tool is used with the default spherical, uniform
% likelihood cross-validation search method. See 'help ksize' for more information.
%
% Note that the KDEs for the platforms as well as their estimated PDF are also stored, along with
% the actual divergence value, in the sub-structure KDE.
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
if ~isa(platforms,'numeric') && size(platforms,2) == 3
	error('PLATFORMS should be a K x 3 matrix');
end
	
% define the default optional arguments
optargs = struct('bandwidth',-1,'logbase',exp(1),'weights',ones(1,size(platforms,1)));

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
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%  for each data set in the structure
for ff = 1:length(DATA)

	% get the points we're going to sample
	XX = DATA{ff}.PDF.x;
	YY = DATA{ff}.PDF.y;
	outofpool = isnan(XX);
	X = XX; X(outofpool) = 0;
	Y = YY; Y(outofpool) = 0;

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

	% for each trial, calcualte the KLD
	DATA{ff}.KLD.kld = zeros(DATA{ff}.ntrials,1);
	for tt = 1:DATA{ff}.ntrials
		switch optargs.logbase
			case 2
				DATA{ff}.KLD.kld(tt) = nansum(nansum(DATA{ff}.KLD.p.*log2(DATA{ff}.KLD.p./DATA{ff}.PDF.p(:,:,tt)))); 
			case exp(1)
				DATA{ff}.KLD.kld(tt) = nansum(nansum(DATA{ff}.KLD.p.*log(DATA{ff}.KLD.p./DATA{ff}.PDF.p(:,:,tt))));
			case 10
				DATA{ff}.KLD.kld(tt) = nansum(nansum(DATA{ff}.KLD.p.*log10(DATA{ff}.KLD.p./DATA{ff}.PDF.p(:,:,tt))));
		end
	end
end

