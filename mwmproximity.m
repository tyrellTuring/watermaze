function DATA = mwmproximity(DATA,varargin); 

% function DATA{} = mwmproximity(DATA{},varargin); 
%
% Calculates the proximity measure for a Morris-Water-Maze dataset, P, as outlined in
% Maei et al. (2009), which is the average distance of the animal from the platform (in cm) 
% during the trial. The obligatory input structure, DATA, is a multi-level cell array that
% is assumed to be in the format returned by readwmdf.m (see help readwmdf). The results are 
% stored in the second-level of DATA as P: either an N x 1 vector or an N x P matrix (depending
% on whether the 'platforms' option is set, see below), where N = number of trials and P = 
% number of platforms. 
%
% Optional parameter/value inputs to the function are as follows:
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
optargs = struct('platforms',[]);

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

	% initialize P
	if isempty(optargs.platforms)
		DATA{ff}.P = zeros(DATA{ff}.ntrials,1);
	else
		DATA{ff}.P = zeros(DATA{ff}.ntrials,size(optargs.platforms,1));
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

			% calculate the proximity 
			P = calc_P(DATA{ff}.path(1:DATA{ff}.ntimes(tt),:,tt), platforms(pp,1:2));
			DATA{ff}.P(tt,pp) = P;
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
