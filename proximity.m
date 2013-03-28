function P = proximity(DATA)

% function P = proximity(DATA)
%
% Returns the mean proximity values for each trial and each animal in the structure DATA (see
% below).
%
% IMPORTANT: Must be run *after* mwmproximity.
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
%   P - A T x A matrix of proximities, or a T x A x 4 matrix (if 'allquads' was used, see help
%       mwmproximity), where T is the maximum number of trials and A is the number of animals.
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
P = zeros(max(ntrial),length(DATA),size(DATA{1}.P,2));

% get the latencies from each animal's trials
for aa = 1:length(DATA)
	for tt = 1:DATA{aa}.ntrials
		P(tt,aa,:) = DATA{aa}.P(tt,:);
	end
end

P = squeeze(P);
