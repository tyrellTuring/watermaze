function [GD, GS] = getgroup(STUDY,gv,dnum,excl,qual,disj)

% function [GD, GS] = getgroup(STUDY,GROUPVALS,[DATANUM,EXCLUDE,MINQUALITY,DISJUNCTION])
%
% Returns two vectors, GS and GD, of indices to the data for animals in STUDY who belong to the 
% group determined by the values passed in GROUPVALS. GROUPVALS must be a cell array listing the 
% values for each of the grouping variables, in the order they are stored in STUDY. For example, 
% if a study has two grouping variables, 'delay' and 'drug' with values '30 days' or '1 day' and 
% 'CNO' or 'Control' respectively, then GROUPVALS = {'30 days','CNO'} would return the indices of all 
% animals with 'delay = 30 days' *and* 'drug = CNO'. GD are the indices for the DATA structure, GS
% are the indices for the STUDY structure.
%
% DATANUM is an optional integer that can be provided to obtain indices for the animals relative to
% different data collections (see help readstudy). By default datanum = 1.
%
% If certain animals should be excluded they can be identified via their unique animal IDs, passed
% in as a cell array EXCLUSION. (None are excluded by default).
%
% If animals should be excluded based on a minimum quality of data rating (e.g. from histology)
% pass in the minimum quality rating, QUALITY. (All quality levels are included by default.) 
%
% If the disjunction of the groups is desired, the optional argument 'DISJUNCTION' can be passed 
% in as true (it is false by default).
%
% Note: if a value in GROUPVALS is NaN, then all values of that grouping variable are selected.
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
	error('GROUPVALS should be a cell array of strings');
end
if nargin < 3, dnum = 1; end;
if nargin < 4, excl = {}; end;
if nargin < 5, qual = -1; end;
if nargin < 6, disj = false; end;

% initialize a boolean matrix for logical operations
ingroup = false(STUDY.ANIMAL.n,length(STUDY.ANIMAL.GROUP.vars));
 
% for each group variable, determine membership
for gg = 1:length(STUDY.ANIMAL.GROUP.vars)
	if isnan(gv{gg})
		ingroup(:,gg) = true;
	else
		ingroup(:,gg) = logical(ismember(STUDY.ANIMAL.GROUP.values{gg},gv{gg}));
	end
end

% get a logical vector of the exclusions
exclude = false(length(STUDY.ANIMAL.id),1);
for aa = 1:length(STUDY.ANIMAL.id)
	if ismember(STUDY.ANIMAL.id{aa},excl) || STUDY.ANIMAL.quality(aa) < qual, exclude(aa) = true; end;
end

% calculate the logical conjunction or disjunction of the groups
if disj
	GS = (find(not(exclude) & any(ingroup,2)))';
	GD = STUDY.data_i{dnum}(GS);
else
	GS = (find(not(exclude) & all(ingroup,2)))';
	GD = STUDY.data_i{dnum}(GS);
end

% remove any zeros
GS = GS(GS ~= 0);
GD = GD(GD ~= 0);
