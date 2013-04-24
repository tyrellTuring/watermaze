function UG = getuniquegroups(STUDY)

% Returns a 2D cell array of all unique group categories in STUDY.
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

% determine the unique values and the number of unique groups
gnum = [];
for gg = 1:length(STUDY.ANIMAL.GROUP.vars)
	unique_vals{gg} = unique(STUDY.ANIMAL.GROUP.values{gg});
	gnum = [gnum, length(unique_vals{gg})];
end

% get the group values for all unique groups
UG = cell(prod(gnum),length(unique_vals));
lastnumrem = 1;
for gg = 1:length(unique_vals)
	numrem = prod(gnum(gg+1:end));
	numpas = prod(gnum(1:gg-1));
	for pp = 1:numpas
		for vv = 1:gnum(gg)
			for rr = 1:numrem
				UG{rr + (vv-1)*numrem + (pp-1)*lastnumrem,gg} = unique_vals{gg}{vv};
			end
		end
	end
	lastnumrem = numrem;
end

