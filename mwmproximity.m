function DATA = mwmproximity(DATA,varargin); 

% function DATA = mwmproximity(DATA,varargin); 
%
% Calculates the proximity measure for a Morris-Water-Maze dataset, P, as outlined in
% Maei et al. (2009), which is the average distance of the animal from the platform (in cm) 
% during the trial. The obligatory input structure, DATA, is a multi-level cell array that
% is assumed to be in the format returned by readwmdf.m (see help readwmdf). The results are 
% stored in the second-level of DATA as either an N x 1 vector or an N x 4 matrix (depending
% on whether the 'allquads' option is set, see below), where N = number of trials. For example,
% if 'allquads' was not set, then DATA{i}.P(j) is the entropy for the jth trial of the ith file
% in the DATA structure.
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'allquads': Boolean valued flag indicating whether to simultaneously calculate the measure for all
%               four possible quadrant locations of the platform. If set to true then P is a N x 4
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

	% initialize P
	if optargs.allquads
		DATA{ff}.P = zeros(DATA{ff}.ntrials,4);
	else
		DATA{ff}.P = zeros(DATA{ff}.ntrials,1);
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

		% calculate the proximity 
		P = calc_P(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform);
		if optargs.allquads		
			P2 = calc_P(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q2);
			P3 = calc_P(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q3);
			P4 = calc_P(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q4);
			P = [P, P2, P3, P4];
		end

		% store the result in the DATA structure
		if optargs.allquads
			DATA{ff}.P(tt,:) = P;
		else
			DATA{ff}.P(tt) = P;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTION
%
% function P = calc_P(path, platform)
%
% Calculates the proximity measure, P, for the given platform location.
%
function P = calc_P(path, platform)

% calculate the distance of the animal from the platform at each timestep
D  = sqrt(sum((path(:,2:3) - repmat(platform,[size(path,1) 1])).^2,2));

% calculate the proximity measure (P)
P = mean(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SUB-FUNCTION
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
end
