function [GI,GVAR,GVAL] = getanovagroup(STUDY,gv,dnum,excl)

% function [GI,GVAR,GVAL] = getanovagroup(STUDY,GROUPVALS,[DATANUM,EXCLUDE])
%
% Returns a vector of animal indices, GI, and cell arrays of group variables, GVAR, and values, 
% GVAL, that can be used in the function anovan (which is part of the Matlab Stats package). 
% GI can be used to access data for non-excluded animals from other vectors. GROUPVALS must be a
% cell array of cell arrays containing the group sets that are desired for the anova analyis 
% (for all grouping variables, in the order they are stored in STUDY). For example, if a study has two
% grouping variables, 'delay' and 'drug' with values '60 days', '30 days' or '1 day' and 'CNO' or 'Control'
% respectively, then GROUPVALS = {{'30 days','1 day'},{'CNO','Control'} would return variables to do
% a 2-way anova comparing the effects of '30 days' vs '1 day' and 'CNO' vs 'Control'.
%
% Note: if an array in GROUPVALS is NaN, then all values of that grouping variable are selected.
%
% DATANUM is an optional integer that can be provided to obtain indices for the animals relative to
% different data collections (see help readstudy). By default datanum = 1.
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
if ~isa(gv,'cell')
	error('GROUPVALS should be a cell array');
end
if nargin < 3, dnum = 1;  end;
if nargin < 4, excl = {}; end;

% get all the unique groups in the study
UG = getuniquegroups(STUDY);

% for each unique group get the indices of the animals and add it to our vector
GI = [];
GS = [];
for uu = 1:size(UG,1)
	G = {};
	usegroup = true;
	for gg = 1:size(UG,2)
		if ismember(UG{uu,gg},gv{gg})
			G = [G, UG{uu,gg}];
		elseif ~isa(gv{gg},'cell') && isnan(gv{gg})
			G = [G, UG{uu,gg}];
		else
			usegroup = false;
			break;
		end
	end
	if usegroup
		[gdi, gsi] = getgroup(STUDY,G,dnum,excl);
		GI = [GI, gdi];
		GS = [GS, gsi];
	end
end

% get the group variables and values for these animals
GVAR = {};
gcounter = 1;
for gg = 1:length(gv)
	if length(gv{gg}) > 1 || (~isa(gv{gg},'cell') && isnan(gv{gg}))
		GVAR = [GVAR,STUDY.ANIMAL.GROUP.vars{gg}];
		values = {};
		for vv = 1:length(GI)
			values = [values;STUDY.ANIMAL.GROUP.values{gg}{GS(vv)}];
		end
		GVAL{gcounter} = values;
		gcounter = gcounter + 1;
	end
end
