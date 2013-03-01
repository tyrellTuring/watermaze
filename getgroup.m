function GI = getgroup(STUDY,gv,disj)

% function GI = getgroup(STUDY,GROUPVALS,[DISJUNCTION])
%
% Returns a vector, GI, of indices of animals in STUDY who belong to the group determined by the
% values passed in GROUPVALS. GROUPVALS must be a cell array listing the values for each of the
% grouping variables, in the order they are stored in STUDY. For example, if a study has two
% grouping variables, 'delay' and 'drug' with values '30 days' or '1 day' and 'CNO' or 'Control'
% respectively, then GROUPVALS = {'30 days','CNO'} would return the indices of all animals with
% 'delay = 30 days' *and* 'drug = CNO'. If the disjunction of the groups is desired, the optional
% argument 'DISJUNCTION' can be passed in as true (it is false by default).
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
if nargin < 3, disj = false; end;

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

% calculate the logical conjunction or disjunction of the groups
if disj
	GI = find(any(ingroup,2));
else
	GI = find(all(ingroup,2));
end
