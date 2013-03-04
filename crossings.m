function X = crossings(DATA)

% function X = crossings(DATA)
%
% Returns the mean number of platform crossings for each trial and each animal in the structure DATA (see
% below).
%
% IMPORTANT: Must be run *after* mwmcrossings.
%
% MANDATORY INPUTS:
% -------------------------------------------------------------------------------------------------
%
%   DATA - DATA is a cell array that contains the raw water-maze data for a given collection of 
%          project files in the study. See 'help readwmdf' for details.
%
% OUTPUT:
% -------------------------------------------------------------------------------------------------
%
%   X - A T x A matrix of crossings, or a T x A x 4 matrix (if 'allquads' was used, see help
%       mwmcrossings), where T is the maximum number of trials and A is the number of animals.
% 
%--------------------------------------------------------------------------------
%
% 02/2013, Frankland Lab (www.franklandlab.com)
%
% Author: Blake Richards
% Contact: blake.richards@utoronto.ca
%

% determine the maximum number of trials
ntrial = zeros(length(DATA),1);
for aa = 1:length(DATA)
	ntrial(aa) = DATA{aa}.ntrials;
end 

% initialize
X = zeros(max(ntrial),length(DATA),size(DATA{1}.X,2));

% get the latencies from each animal's trials
for aa = 1:length(DATA)
	for tt = 1:DATA{aa}.ntrials
		X(tt,aa,:) = DATA{aa}.X(tt,:);
	end
end

X = squeeze(X);
