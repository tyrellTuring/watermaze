function L = latencies(DATA)

% function L = latencies(DATA)
%
% Returns the escape latencies in seconds for each trial and each animal in the structure DATA (see
% below).
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
%   L - A T x A matrix of latencies, where T is the maximum number of trials and A is the number of
%       animals.
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
L = zeros(max(ntrial),length(DATA));

% get the latencies from each animal's trials
for aa = 1:length(DATA)
	for tt = 1:DATA{aa}.ntrials
		L(tt,aa) = DATA{aa}.path(DATA{aa}.ntimes(tt),1,tt);
	end
end
