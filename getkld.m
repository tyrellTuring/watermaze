function [K MK] = getkld(DATA,pindex)

% function [K MK] = getkld(DATA,[PINDEX])
%
% Returns the Kullback-Leibler divergences from a data set. Both mwmpdf and mwmkld must have been
% run first. K is the KLDs for specific trials, MK is the KLDs for the mean trial KLDs. The argument
% PINDEX can be used to select a particular platform set KLD, if multiple ones were calculated.
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
%   K - A T x A matrix of KLD values, where T is the maximum number of trials and A is the number of
%       animals.
%
%   MK - A M x A matrix of KLD values, where M is the number of trial means and A is the number of
%        animals.
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

% check arguments
if nargin < 2, pindex = 1; end;

% check whether multiple platform sets weere calculated for
if isa(DATA{1}.KLD,'cell')
	
	% determine the number of trials
	ntrial = zeros(length(DATA),1);
	for aa = 1:length(DATA)
		ntrial(aa) = DATA{aa}.ntrials;
	end 

	% determine the number of means
	nmeans = zeros(length(DATA),1);
	for aa = 1:length(DATA)
		nmeans(aa) = length(DATA{aa}.KLD{pindex}.mkld);
	end 

	% initialize 
	K  = zeros(max(ntrial),length(DATA));
	MK = zeros(max(nmeans),length(DATA));

	% get the KLD{pindex}s from each animal's trials
	for aa = 1:length(DATA)
		for tt = 1:DATA{aa}.ntrials
			K(tt,aa) = DATA{aa}.KLD{pindex}.kld(tt);
		end
	end

	% get the KLD{pindex}s from each animal's means
	for aa = 1:length(DATA)
		for mm = 1:nmeans(aa)
			MK(mm,aa) = DATA{aa}.KLD{pindex}.mkld(mm);
		end
	end
else

	% determine the number of trials
	ntrial = zeros(length(DATA),1);
	for aa = 1:length(DATA)
		ntrial(aa) = DATA{aa}.ntrials;
	end 

	% determine the number of means
	nmeans = zeros(length(DATA),1);
	for aa = 1:length(DATA)
		nmeans(aa) = length(DATA{aa}.KLD.mkld);
	end 

	% initialize 
	K  = zeros(max(ntrial),length(DATA));
	MK = zeros(max(nmeans),length(DATA));

	% get the KLDs from each animal's trials
	for aa = 1:length(DATA)
		for tt = 1:DATA{aa}.ntrials
			K(tt,aa) = DATA{aa}.KLD.kld(tt);
		end
	end

	% get the KLDs from each animal's means
	for aa = 1:length(DATA)
		for mm = 1:nmeans(aa)
			MK(mm,aa) = DATA{aa}.KLD.mkld(mm);
		end
	end
end
