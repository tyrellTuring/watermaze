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
% or an N x M x P array, depending on whether the 'platforms' option is set (see below).
% (N = number of trials, M = number of zones, P = number of platforms).
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'zones'   : Vector indicating the radius of the zones that should be used to calculate the measure.
%               Note that the returned values, Z, will have the same number of final columns as there are
%               elements in this vector. So, if the vector is [50 40 30 20 10] then Z will have a
%               size, of 5 for the third dimension, each one corresponding to the time spent in each zone 
%               of radius 50, 40, etc. centred on the platform. Default = [20 15 10];
%
% - 'platforms': P x 3 matrix of platform x, y and radius values to use for the measurement. If this
%                argument is not set then the platform contained in the .wmdf file for each trial is
%                used.
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
optargs = struct('zones',[20 15 10],'platforms',true);

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

	% initialize Z
	if isempty(optargs.platforms)
		DATA{ff}.Z = zeros(DATA{ff}.ntrials,length(optargs.zones),1);
	else
		DATA{ff}.Z = zeros(DATA{ff}.ntrials,length(optargs.zones),size(optargs.platforms,1));
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
	
		% for each platform...
		for pp = 1:size(platforms,1)

			% calculate the zones measure 
			Z = calc_Z(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platforms(pp,:), optargs.zones);
			DATA{ff}.Z(tt,:,pp) = Z;
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
