function DATA = mwmentropy(DATA,varargin); 

% function DATA = mwmentropy(DATA,varargin); 
%
% Calculates the entropy measure for a Morris-Water-Maze dataset, H. The entropy
% can either be estimated using a parametric method as outlined in Maei et al. (2009) 
% or using a kernel density estimate (in which case mwmpdf must be run first). The obligatory 
% input structure, DATA, is a multi-level cell array that is assumed to be in the format 
% returned by readwmdf.m (see help readwmdf). The results are stored in the second-level of 
% DATA as either an N x 1 vector or an N x P matrix (depending on whether the 'platforms' option 
% is set, see below), where N = number of trials and P = number of platforms. 
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'parametric' : whether to use a parametric method to calculate the entropy. Default = true.
%
% - 'lambda'  : Scalar valued mixing coefficient [H = lambda*H_error + (1-lambda)*H_path]. Default = 0.5.
%
% - 'platforms': P x 3 matrix of platform x, y and radius values to use for the measurement. If this
%                argument is not set then the platform contained in the .wmdf file for each trial is
%                used.
%
%--------------------------------------------------------------------------------
%
% References
% - Maei, HR, et al. (2009). Development and Validation of a Sensitive Entropy-Based Measure for the Water
%     Maze, Frontiers in Integrative Neuroscience, 3:33.
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
	
% define the default optional arguments
optargs = struct('parametric',true,'lambda',0.5,'platforms',[]);

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
			case 'parametric'
                if isa(pair{2},'logical');
					optargs.(inpname) = pair{2};
                else
					error('parametric must be a logical');
				end
			case 'lambda'
                if isa(pair{2},'numeric') && (pair{2} >= 0 && pair{2} <= 1);
					optargs.(inpname) = pair{2};
                else
					error('lambda must be a scalar value between 0 and 1');
				end
			case 'platforms'
				if isa(pair{2},'numeric') && size(pair{2},2) == 3
					optargs.(inpname) = pair{2};
				else
					error('platforms must be a P x 3 matrix');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%  for each data set in the structure
for ff = 1:length(DATA)

	% initialize H
	if optargs.parametric && ~isempty(optargs.platforms)
		DATA{ff}.H = zeros(DATA{ff}.ntrials,size(optargs.platforms,1));
	else
		DATA{ff}.H = zeros(DATA{ff}.ntrials,1);
	end

	% for each trial
	for tt = 1:DATA{ff}.ntrials

		% get the platform data
		if isempty(optargs.platforms)
			platforms = DATA{ff}.platform(min([tt end]),1:3);
		else
			platforms = optargs.platforms;
		end

		% get the pool data
		pool = DATA{ff}.pool(min([tt end]),1:2);

		% calculate the entropy
		if optargs.parametric
	
			% for each platform...
			for pp = 1:size(platforms,1)
				H = calc_H(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platforms(pp,:), optargs.lambda);
				DATA{ff}.H(tt,pp) = H;
			end
		else
			DATA{ff}.H(tt) = -nansum(nansum(DATA{ff}.PDF.p(:,:,tt).*log(DATA{ff}.PDF.p(:,:,tt))));
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTION
%
% function H = calc_H(path, platform, lambda)
%
% Calculates the entropy measure, H, for the given platform location and the given 
% mixture coefficient, lambda.
%
function H = calc_H(path, platform, lambda)

% calculate the distance of the animal from the platform and the mean path location at each timestep
Dx = path(:,2) - platform(1);
Dy = path(:,3) - platform(2);
xd     = path(:,2) - mean(path(:,2));
yd     = path(:,3) - mean(path(:,3));

% calculate the covariance matrix of the error ellipse for the H_path measurement
mxd    = mean(xd);
myd    = mean(yd);
mxxd   = mean(xd.^2);
myyd   = mean(yd.^2);
mxyd   = mean(xd.*yd);
Sigma  = [mxxd - mxd.^2, mxyd - mxd.*myd; mxyd - mxd.*myd, myyd - myd.^2]; % note mxxd - mxd.^2 = var(x), etc.

% get the eigen values for the covariance matrix and multiply them to get sigma_a*sigma_b
eigen  = eig(Sigma);
varsig = eigen(1)*eigen(2);

% calculate the variance of the error ellipse for the H_error measurement
mD2    = mean(Dx.^2 + Dy.^2);

% calculate the entropy measure (H)
H      = lambda*2*log(mD2) + (1-lambda)*log(varsig); % H = lambda*H_error + (1-lambda)*H_path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SUB-FUNCTION
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
end
