function DATA = mwmzonelat(DATA,varargin)

% function DATA = mwmzonelat(DATA,varargin); 
%
% Calculates the latency to a given set of zones in the water maze.
%
% The obligatory input structure, DATA, is a multi-level cell array that is assumed to 
% be in the format returned by readwmdf.m (see help readwmdf). The results are stored 
% in a second-level variable, Z, in the cell array of DATA. Z is either an N x M matrix
% or an N x M x P array, depending on whether the 'allquads' option is set (see below).
% (N = number of trials, M = number of zones, P = number of platforms). 
%
% Optional parameter/value inputs to the function are as follows:
%
% - 'zones'   : Vector indicating the radius of the zones that should be used to calculate the measure.
%               Note that the returned values, Z, will have the same number of final columns as there are
%               elements in this vector. So, if the vector is [50 40 30 20 10] then Z will have a
%               size, of 5 for the third dimension, each one corresponding to the latency to each zone 
%               of radius 50, 40, etc. centred on the platform. Default = [20 15 10];
%
% - 'platforms': P x 3 matrix of platform x, y and radius values to use for the measurement. If this
%                argument is not set then the platform contained in the .wmdf file for each trial is
%                used.
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
optargs = struct('zones',[20 15 10],'platforms',[]);

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
		DATA{ff}.L = zeros(DATA{ff}.ntrials,length(optargs.zones),1);
	else
		DATA{ff}.L = zeros(DATA{ff}.ntrials,length(optargs.zones),size(optargs.platforms,1));
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
			Z = calc_L(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platforms(pp,:), optargs.zones);

			% store the result in the DATA structure
			DATA{ff}.L(tt,:,pp) = Z;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-FUNCTION
%
% function L = calc_L(path, zones_edges)
%
% Calculates the latency to zone measure, L, for the given zones.
%
function L = calc_L(path, platform, zones)

% initialize L
L = nan(length(zones),1);

% get the time-steps
tsteps = path(:,1);

% for each zone...
for zz = 1:length(zones)

	% determine the first time into the given zone
	inzone = find(sqrt(sum((path(2:end,2:3) - repmat(platform(1:2),[size(path,1)-1 1])).^2,2)) <= zones(zz),1,'first');

	% calculate L
	if ~isempty(inzone), L(zz)  = tsteps(inzone); end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SUB-FUNCTION
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN FUNCTION
end
