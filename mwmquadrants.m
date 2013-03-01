function DATA = mwmquadrant(DATA,varargin); 

% function DATA = mwmquadrant(DATA,varargin); 
%
% Calculates the quadrants measure for a Morris-Water-Maze dataset, Q, as outlined in
% Maei et al. (2009). This is simply the percent time spent in the "correct" quadrant.
% The obligatory input structure, DATA, is a multi-level cell array that is assumed to 
% be in the format returned by readwmdf.m (see help readwmdf). The results are stored 
% in the second-level of the cell array DATA as either a N x 1 vector or an N x 4 matrix
% depending on whether the 'allquads' option is set (see below). So, for example, if 
% 'allquads' is not set, then DATA{i}.Q(j) is the quadrant measure for the jth trial of 
% the ith file in the DATA structure. 
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'allquads': Boolean valued flag indicating whether to simultaneously calculate the measure for all
%               four possible quadrant locations of the platform. If set to true then Q is an N x 4
%               matrix where the first column corresponds to the actual platform location. Note that
%               in this case the sum of across columns Q = 100%. Default = false.
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

	% initialize Q
	if optargs.allquads
		DATA{ff}.Q = zeros(DATA{ff}.ntrials,4);
	else
		DATA{ff}.Q = zeros(DATA{ff}.ntrials,1);
	end

	% for each trial
	for tt = 1:DATA{ff}.ntrials

		% get the platform data
		platform = DATA{ff}.platform(min([tt end]),1:2);

		% get the pool data
		pool = DATA{ff}.pool(min([tt end]),:);

		% calculate the vectors of the edges of the correct quadrant
		plat_rel_pool  = (platform - pool(1:2))';                                              % location of platform relative to pool centre
		p_rotation_mat = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];                          % 45 degree rotation matrix
		n_rotation_mat = [cos(-pi/4) -sin(-pi/4); sin(-pi/4) cos(-pi/4)];                      % -45 degree rotation matrix
		p_edge         = p_rotation_mat*plat_rel_pool;                                         % positive edge of quadrant
		n_edge         = n_rotation_mat*plat_rel_pool;                                         % negative edge of quadrant
		qedges         = [p_edge/sqrt(sum(p_edge.^2)), n_edge/sqrt(sum(n_edge.^2))]*pool(3); % quadrant edges scaled to pool size

		% if allquads was set, calculate the positions of the equivalent quadrant edges for other quadrants
		if optargs.allquads
			rotation_mat = [0 -1; 1 0];            % 90 degree rotation matrix
			qedges_q2    = rotation_mat*qedges;   % 90 degrees rotation
			qedges_q3    = rotation_mat^2*qedges; % 180 degrees rotation
			qedges_q4    = rotation_mat^3*qedges; % 270 degrees rotation
		end
	
		% make sure the path is corrected relative to the centre of the pool
		path = DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt);
		path(:,2:3) = path(:,2:3) - repmat(pool(1:2),[DATA{ff}.ntimes(tt) 1]);

		% calculate the quadrant measure 
		Q = calc_Q(path, qedges);
		if optargs.allquads		
			Q2 = calc_Q(path, qedges_q2);
			Q3 = calc_Q(path, qedges_q3);
			Q4 = calc_Q(path, qedges_q4);
			Q = [Q, Q2, Q3, Q4];
		end

		% store the result in the DATA structure
		if optargs.allquads
			DATA{ff}.Q(tt,:) = Q;
		else
			DATA{ff}.Q(tt) = Q;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTION
%
% function Q = calc_Q(path, quadrant_edges)
%
% Calculates the quadrant measure, Q, for the given quadrant edges.
%
function Q = calc_Q(path, qedges)

% get the time-steps
tsteps = diff(path(:,1));

% get the path and quadrant info in polar form
[qth qr] = cart2pol(qedges(1,:),qedges(2,:));
[pth pr] = cart2pol(path(2:end,2),path(2:end,3));

% make sure we're not in negative space
qth = mod(qth,pi*2);
pth = mod(pth,pi*2);

% determine when the animal was in the quadrant
if qth(1) > qth(2) 
	inQ = pth < qth(1) & pth > qth(2);
else
	inQ = ~(pth > qth(1) & pth < qth(2));
end

% calculate the amount of time spent in the quadrant
tQ  = sum(tsteps(inQ));

% calculate the quadrant measure (Q)
Q = tQ/sum(tsteps)*100; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SUB-FUNCTION
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
end
