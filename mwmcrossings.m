function DATA = mwmcrossings(DATA,varargin); 

% function DATA = mwmcrossings(DATA,varargin); 
%
% Calculates the number of crossings measure for a Morris-Water-Maze dataset, X, as outlined in
% Maei et al. (2009). The obligatory input structure, DATA, is a multi-level cell array that
% is assumed to be in the format returned by readwmdf.m (see help readwmdf). The results are 
% stored in the second-level of DATA as either an N x 1 vector or an N x 4 matrix (depending
% on whether the 'allquads' option is set, see below), where N = number of trials. For example,
% if 'allquads' was not set, then DATA{i}.X(j) is the number of crossing for the jth trial of 
% the ith file in the DATA structure.
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'allquads': Boolean valued flag indicating whether to simultaneously calculate the measure for all
%               four possible quadrant locations of the platform. If set to true then X is a N x 4
%               matrix where the first column corresponds to the actual platform locations. Default = false.
%
%--------------------------------------------------------------------------------
%
% References
% - Maei, HR, et al. (2009). What is the Most Sensitive Measure of Water Maze Probe Test Performance?
%     Frontiers in Integrative Neuroscience, 3:4.
%
%--------------------------------------------------------------------------------
%
% 12/2011, Frankland Lab (www.franklandlab.com)
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
optargs = struct('allquads',false);

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

	% initialize X
	if optargs.allquads
		DATA{ff}.X = zeros(DATA{ff}.ntrials,4);
	else
		DATA{ff}.X = zeros(DATA{ff}.ntrials,1);
	end

	% for each trial
	for tt = 1:DATA{ff}.ntrials

		% get the platform data
		platform = DATA{ff}.platform(min([tt end]),1:3);

		% get the pool data
		pool = DATA{ff}.pool(min([tt end]),1:2);

		% if allquads was set, calculate the positions of the equivalent platforms in each of the other quadrants
		if optargs.allquads
			plat_rel_pool = (platform(1:2) - pool)';                               % location of platform relative to pool centre
			rotation_mat  = [0 -1; 1 0];                                           % 90 degree rotation matrix
			platform_q2   = [(rotation_mat*plat_rel_pool + pool')' platform(3)];   % platform rotated by 90 degrees
			platform_q3   = [(rotation_mat^2*plat_rel_pool + pool')' platform(3)]; % platform rotated by 180 degrees
			platform_q4   = [(rotation_mat^3*plat_rel_pool + pool')' platform(3)]; % platform rotated by 270 degrees
		end

		% calculate the proximity 
		X = calc_X(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform);
		if optargs.allquads		
			X2 = calc_X(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q2);
			X3 = calc_X(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q3);
			X4 = calc_X(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q4);
			X = [X, X2, X3, X4];
		end

		% store the result in the DATA structure
		if optargs.allquads
			DATA{ff}.X(tt,:) = X;
		else
			DATA{ff}.X(tt) = X;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTION
%
% function X = calc_X(path, platform)
%
% Calculates the number of crossings, X, for the given platform location.
%
function X = calc_X(path, platform)

% determine all of the path points that are within the platform area
inplat = sqrt(sum((path(:,2:3) - repmat(platform(1:2),[size(path,1) 1])).^2,2)) <= platform(3);

% count the number of non-consecutive entries into the platform
X = sum(diff(inplat) == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SUB-FUNCTION
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
end
