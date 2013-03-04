function GI = getgroup(STUDY,gv,excl,disj)

% function GI = getgroup(STUDY,GROUPVALS,[EXCLUDE, DISJUNCTION])
%
% Returns a vector, GI, of indices of animals in STUDY who belong to the group determined by the
% values passed in GROUPVALS. GROUPVALS must be a cell array listing the values for each of the
% grouping variables, in the order they are stored in STUDY. For example, if a study has two
% grouping variables, 'delay' and 'drug' with values '30 days' or '1 day' and 'CNO' or 'Control'
% respectively, then GROUPVALS = {'30 days','CNO'} would return the indices of all animals with
% 'delay = 30 days' *and* 'drug = CNO'.
%
% If certain animals should be excluded they can be identified via their unique animal IDs, passed
% in as a cell array EXCLUSION. (None are excluded by default).
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

% parse the args
if ~isa(gv,'cell')
	error('GROUPVALS should be a cell array of strings');
end
if nargin < 3, excl = {}; end;
if nargin < 4, disj = false; end;

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
	if ismember(STUDY.ANIMAL.id{aa},excl), exclude(aa) = true; end;
end

% calculate the logical conjunction or disjunction of the groups
if disj
	GI = find(not(exclude) & any(ingroup,2));
else
	GI = find(not(exclude) & all(ingroup,2));
end
