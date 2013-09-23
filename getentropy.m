function [H, S] = getentropy(DATA)

% function [H, S] = getentropy(DATA)
%
% Returns the entropy values for each trial and each animal in the structure DATA (see
% below).
%
% IMPORTANT: Must be run *after* mwmentropy.
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
%   H - A T x A matrix of entropy estimates, or a T x A x P matrix (if P platforms were used, see help
%       mwmentropy), where T is the maximum number of trials and A is the number of animals.
%
%   S - A 2 x 2 x T x A matrix of 2D standard deviations.
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

% determine the maximum number of trials
ntrial = zeros(length(DATA),1);
for aa = 1:length(DATA)
	ntrial(aa) = DATA{aa}.ntrials;
end 

% initialize
H = zeros(max(ntrial),length(DATA),size(DATA{1}.H,2));
S = zeros(2,2,max(ntrial),length(DATA),size(DATA{1}.H,2));

% get the latencies from each animal's trials
for aa = 1:length(DATA)
	for tt = 1:DATA{aa}.ntrials
		for pp = 1:size(DATA{aa}.H,2)
			H(tt,aa,pp) = DATA{aa}.H(tt,pp);
			S(:,:,tt,aa,pp) = DATA{aa}.Sigma(:,:,tt,pp);
		end
	end
end

H = squeeze(H);
S = squeeze(S);
