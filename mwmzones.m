function DATA = mwmzones(DATA,varargin); 

% function DATA = mwmzones(DATA,varargin); 
%
% Calculates the zones measure for a Morris-Water-Maze dataset, Z, as outlined in
% Maei et al. (2009). This is simply the percent time spent in the "correct" zones, which
% are defined as circular zones of different radiuses around the platform.
%
% The obligatory input structure, DATA, is a multi-level cell array that is assumed to 
% be in the format returned by readwmdf.m (see help readwmdf). The results are stored 
% in a second-level variable, Z, in the cell array of DATA. Z is either an N x M matrix
% or an N x 4 x M array, depending on whether the 'allquads' option is set (see below).
% (N = number of trials and M = number of zones). For example, if 'allquads' was not set, 
% and 4 zones were defined, then DATA{i}.Z(j,k) contains the zones measures for the kth zone
% during the jth trial of the ith file in the DATA structure. 
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'zones'   : Vector indicating the radius of the zones that should be used to calculate the measure.
%               Note that the returned values, Z, will have the same number of rows as there are
%               elements in this vector. So, if the vector is [50 40 30 20 10] then Z will have 5 rows,
%		each one corresponding to the time spent in each zone of radius 50, 40, etc. centred
%               on the platform. Default = [20 15 10];
%
% - 'allquads': Boolean valued flag indicating whether to simultaneously calculate the measure for all
%               four possible quadrant locations of the platform. If set to true then Z is a 4-element
%               vector where the first element corresponds to the actual platform location. Note that
%               in this case the sum of Z = 100%. Default = false.
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
optargs = struct('zones',[20 15 10],'allquads',false);

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
			case 'zones'
				if isa(pair{2},'numeric') && all(pair{2} > 0)
					optargs.(inpname) = pair{2};
                else
					error('zones must be a positive vector');
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

	% initialize Z
	if optargs.allquads
		DATA{ff}.Z = zeros(DATA{ff}.ntrials,4,length(optargs.zones));
	else
		DATA{ff}.Z = zeros(DATA{ff}.ntrials,length(optargs.zones));
	end

	% for each trial
	for tt = 1:DATA{ff}.ntrials

		% get the platform data
		platform = DATA{ff}.platform(min([tt end]),1:2);

		% get the pool data
		pool = DATA{ff}.pool(min([tt end]),1:2);

		% if allquads was set, calculate the positions of the equivalent platforms in each of the other quadrants
		if optargs.allquads
			plat_rel_pool = (platform(1:2) - pool)';                 % location of platform relative to pool centre
			rotation_mat  = [0 -1; 1 0];                             % 90 degree rotation matrix
			platform_q2   = (rotation_mat*plat_rel_pool + pool')';   % platform rotated by 90 degrees
			platform_q3   = (rotation_mat^2*plat_rel_pool + pool')'; % platform rotated by 180 degrees
			platform_q4   = (rotation_mat^3*plat_rel_pool + pool')'; % platform rotated by 270 degrees
		end

		% calculate the zones measure 
		Z = calc_Z(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform, optargs.zones);
		if optargs.allquads		
			Z2 = calc_Z(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q2, optargs.zones);
			Z3 = calc_Z(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q3, optargs.zones);
			Z4 = calc_Z(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platform_q4, optargs.zones);
			Z = permute([Z, Z2, Z3, Z4],[2 1]);
		end

		% store the result in the DATA structure
		if optargs.allquads
			DATA{ff}.Z(tt,:,:) = Z;
		else
			DATA{ff}.Z(tt,:) = Z;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTION
%
% function Z = calc_Z(path, zones_edges)
%
% Calculates the zones measure, Z, for the given zones edges.
%
function Z = calc_Z(path, platform, zones)

% initialize Z
Z = zeros(length(zones),1);

% get the time-steps
tsteps = diff(path(:,1));

% for each zone...
for zz = 1:length(zones)

	% determine all of the path points that are within the given zone
	inzone = sqrt(sum((path(2:end,2:3) - repmat(platform(1:2),[size(path,1)-1 1])).^2,2)) <= zones(zz);

	% calculate the amount of time spent in the zones
	tzone  = sum(tsteps(inzone));

	% calculate the zones measure (Z)
	Z(zz) = (tzone./sum(tsteps)*100)'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SUB-FUNCTION
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
end
