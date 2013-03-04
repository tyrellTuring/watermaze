function Z = zones(DATA)

% function Z = zones(DATA)
%
% Returns the mean time spent in circular zones around the platform for each trial and each animal 
% in the structure DATA (see below).
%
% IMPORTANT: Must be run *after* mwmzones.
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
%   X - A T x A x K matrix of crossings, or a T x A x 4 x K matrix (if 'allquads' was used, see help
%       mwmzones), where T is the maximum number of trials, A is the number of animals and K is the
%       number of zones.
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
Z = zeros(max(ntrial),length(DATA),size(DATA{1}.Z,2),size(DATA{1}.Z,3));

% get the latencies from each animal's trials
for aa = 1:length(DATA)
	for tt = 1:DATA{aa}.ntrials
		Z(tt,aa,:,:) = DATA{aa}.Z(tt,:,:);
	end
end

Z = squeeze(Z);
