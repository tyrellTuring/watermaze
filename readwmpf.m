function PROJECT = readwmpf(fnames)

% function PROJECT = readwmpf(FILES)
%
% Reads in watermaze project information from files in the binary format produced by Actimetrics
% software (.wmpf files).
%
% MANDATORY INPUTS:
% -------------------------------------------------------------------------------------------------
%
%   FILES - A cell array with the names of all the project files to be read in.
%
% OUTPUT:
% -------------------------------------------------------------------------------------------------
% The output of the function is a cell array, PROJECT, organized as follows:
%
%   1st-level
%   -------------------
%   The first level of the cell array ranges over projects, e.g. PROJECT{3} is the project
%   information for the third project file in the list.
%
%   2nd-level
%   -------------------
%
%   PROJECT{i}.file    - The name of the file for this project (i.e. the .wmpf file).
%
%   PROJECT{i}.n       - The number of animals associated with this project.
%
%   PROJECT{i}.animals - A cell array containing the names of the animals associated with this project.
%
%   PROJECT{i}.data    - A cell array containing the names of the data file for each animal (should be
%                        animal's name according to convention).
%
%   PROJECT{i}.folder  - A string of the path to the folder containing all the data files.
%
%--------------------------------------------------------------------------------
%
% 02/2013, Frankland Lab (www.franklandlab.com)
%
% Author: Blake Richards
% Contact: blake.richards@utoronto.ca
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% check the fnames argument
if ~isa(fnames,'cell')
	error('FILES should be a cell array of file names');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP THROUGH THE FILES AND LOAD EACH ONE

% determine the number of files we're processing
if iscell(fnames), nfiles = length(fnames); else, nfiles = 1; end;

% for each file...
for ff = 1:nfiles

	% create the file identifier
	if iscell(fnames), filename = fnames{ff}; else, filename = fnames; end; 
  try
		fid = fopen(filename);
	catch
		error('Could not open file %s',filename);
	end

	% store the filename
	PROJECT{ff}.file = filename;

	% determine the size of the file in bytes
	status    = fseek(fid, 0, 'eof');
	if status < 0, error('Could not move to end of file %s',filename); end;                  
  	file_size = ftell(fid);

	% skip the first unecessary bytes of the file
	status      = fseek(fid, 0, 'bof');
	if status < 0, error('Could not move to start of file %s',filename); end;                  
	skip        = 12;
	status      = fseek(fid,skip,'bof');
	if status < 0, error('Could not move to twelfth byte in file %s',filename); end;

	% read in each of the animal's names
	aa = 0;
	alength = fread(fid,1,'uint32',0,'ieee-be'); % get the length of the first animal's name
	while alength > 0 && alength < 20

		% increment the animal counter
		aa = aa + 1;

		% read in the animal's name
		PROJECT{ff}.animal{aa} = fread(fid,alength,'*char',0,'ieee-be')';

		% store the associated filename (for easy access)
		PROJECT{ff}.data{aa}   = sprintf('%s.wmpf',PROJECT{ff}.animal{aa});

		% determine the length of the next animal's name
		alength = fread(fid,1,'uint32',0,'ieee-be');
	end

	% store the number of animals (for easy access)
	PROJECT{ff}.n = aa;

	% determine the folder name for this project
	[path, name, ext]  = fileparts(filename);
	PROJECT{ff}.folder = fullfile(path,sprintf('%s Folder',name));
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FUNCTION
end
