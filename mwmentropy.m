function DATA = mwmentropy(DATA,varargin); 

% function DATA = mwmentropy(DATA,varargin); 
%
% Calculates the entropy measure for a Morris-Water-Maze dataset, H, as outlined in
% Maei et al. (2009). The obligatory input structure, DATA, is a multi-level cell array that
% is assumed to be in the format returned by readwmdf.m (see help readwmdf). The results are 
% stored in the second-level of DATA as either an N x 1 vector or an N x 4 matrix (depending
% on whether the 'allquads' option is set, see below), where N = number of trials. For example,
% if 'allquads' was not set, then DATA{i}.H(j) is the entropy for the jth trial of the ith file
% in the DATA structure.
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'lambda'  : Scalar valued mixing coefficient [H = lambda*H_error + (1-lambda)*H_path]. Default = 0.5.
%
% - 'allquads': Boolean valued flag indicating whether to simultaneously calculate the measure for all
%               four possible quadrant locations of the platform. If set to true then H is an N x 4
%               matrix where the first column corresponds to the actual platform locations. Default = false.
%
%--------------------------------------------------------------------------------
%
% References
% - Maei, HR, et al. (2009). Development and Validation of a Sensitive Entropy-Based Measure for the Water
%     Maze, Frontiers in Integrative Neuroscience, 3:33.
%
%--------------------------------------------------------------------------------
%
% 12/2011, Frankland Lab (www.franklandlab.com)
%
% Author: Hamid Maei
% Modified by: Kirill Zaslavsky & Blake Richards
% Contact: blake.richards@utoronto.ca
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% check the DATA argument
if ~isa(DATA,'cell')
	error('DATA should be a cell array');
end
	
% define the default optional arguments
optargs = struct('lambda',0.5,'allquads',false);

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
			case 'lambda'
                if isa(pair{2},'numeric') && (pair{2} >= 0 && pair{2} <= 1);
					optargs.(inpname) = pair{2};
                else
					error('lambda must be a scalar value between 0 and 1');
				end
			case 'allquads'
				if isa(pair{2},'logical');
					optargs.(inpname) = pair{2};
				else
					error('allquads must be a logical');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%  for each data set in the structure
for ff = 1:length(DATA)

	% initialize H
	if optargs.allquads
		DATA{ff}.H = zeros(DATA{ff}.ntrials,4);
	else
		DATA{ff}.H = zeros(DATA{ff}.ntrials,1);
	end

	% for each trial
	for tt = 1:DATA{ff}.ntrials

		% get the platform data
		platform = DATA{ff}.platform(min([tt end]),1:2);

		% get the pool data
		pool = DATA{ff}.pool(min([tt end]),1:2);

		% if allquads was set, calculate the positions of the equivalent platforms in each of the other quadrants
		if optargs.allquads
			plat_rel_pool = (platform - pool)';                      % location of platform relative to pool centre
			rotation_mat  = [0 -1; 1 0];                             % 90 degree rotation matrix
			platform_q2   = (rotation_mat*plat_rel_pool + pool')';   % platform rotated by 90 degrees
			platform_q3   = (rotation_mat^2*plat_rel_pool + pool')'; % platform rotated by 180 degrees
			platform_q4   = (rotation_mat^3*plat_rel_pool + pool')'; % platform rotated by 270 degrees
		end

		% calculate the entropy 
		H = calc_H(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform, optargs.lambda);
		if optargs.allquads		
			H2 = calc_H(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q2, optargs.lambda);
			H3 = calc_H(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q3, optargs.lambda);
			H4 = calc_H(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q4, optargs.lambda);
			H = [H, H2, H3, H4];
		end

		% store the result in the DATA structure
		if optargs.allquads
			DATA{ff}.H(tt,:) = H;
		else
			DATA{ff}.H(tt) = H;
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
