function [GI,GVAR,GVAL] = getanovagroup(STUDY,gv,excl)

% function [GI,GVAR,GVAL] = getanovagroup(STUDY,GROUPVARS,[EXCLUDE])
%
% Returns a vector of animal indices, GI, and cell arrays of group variables, GVAR, and values, 
% GVAL, that can be used in the function anovan (which is part of the Matlab Stats package). 
% GI can be used to access data for non-excluded animals from other vectors. The variable 
% GROUPVARS should be a numeric index listing which grouping variables to include.
%
% If certain animals should be excluded they can be identified via their unique animal IDs, passed
% in as a cell array EXCLUSION. (None are excluded by default).
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

% parse the args
if ~isa(gv,'numeric')
	error('GROUPVARS should be a numeric vector');
end
if nargin < 3, excl = {}; end;

% get a logical vector of the exclusions and the animal indices
exclude = false(length(STUDY.ANIMAL.id),1);
for aa = 1:length(STUDY.ANIMAL.id)
	if ismember(STUDY.ANIMAL.id{aa},excl), exclude(aa) = true; end;
end
GI = STUDY.data_i(find(not(exclude)));

% go through each of the requested grouping variables and collect their data
for gg = 1:length(gv)
	GVAR{gg} = STUDY.ANIMAL.GROUP.vars{gv(gg)};
	GVAL{gg} = STUDY.ANIMAL.GROUP.values{gv(gg)}(find(not(exclude)));
end
