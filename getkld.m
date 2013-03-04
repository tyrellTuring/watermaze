function K = getkld(DATA)

% function K = getkld(DATA)
%
% Returns the Kullback-Leibler divergences from a data set. Both mwmpdf and mwmkld must have been
% run first.
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
%   K - A T x A matrix of latencies, where T is the maximum number of trials and A is the number of
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
K = zeros(max(ntrial),length(DATA));

% get the latencies from each animal's trials
for aa = 1:length(DATA)
	for tt = 1:DATA{aa}.ntrials
		K(tt,aa) = DATA{aa}.KLD.kld(tt);
	end
end
