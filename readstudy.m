function [STUDY,DATA] = readstudy(fname,datadir,varargin)

% function [STUDY,DATA] = readstudy(STUDY_RECORDS, DATA_DIRECTORY, varargin)
%
% Reads in watermaze study information from study records in an Excel spreadsheet and loads all the
% associated data. The function also does some limited cross-referencing between the spreadsheet and
% any binary Actimetrics files to make sure the data is in broad agreement about mouse numbers, ids,
% etc.
%
% MANDATORY INPUTS:
% -------------------------------------------------------------------------------------------------
%
%   STUDY_RECORDS - A string of the name (full path) for the study records Excel spreadsheet. The
%                   spreadsheet must be a .xls file (i.e. Windows XP or earlier) and it must contain
%                   the following column headings in the following order (optional headings are
%                   indicated with '*'):
%
%                     - 'Cage #'
%                     - 'Mouse'
%                     - *'Sex'
%                     - *'DOB'
%                     - *'GROUP VAR 1'
%                     - *'GROUP VAR 2'
%                     - *'GROUP VAR ...' (there can be any number of group variables, can't be 'File')
%                     - 'File'           (.wmpf file for training or probe)
%                     - *'Notes'         (notes on this data)
%                     - ...              (there can be any number of entries in this pattern)
%
%                   Any entries that don't match these names/patterns are simply read in and stored
%                   in the EXTRAS field (see below).
%
%   DATA_DIRECTORY - The directory that contains all of the .wmpf files and their associated folders
%                    conatining .wmdf files.
%
% OUTPUT:
% -------------------------------------------------------------------------------------------------
% The output of the function is three structures, STUDY, PROJECT, and DATA with the following organization:
%
%   DATA{i} - DATA is a cell array of cell arrays. The DATA{i} contains the raw water-maze data for 
%             a given collection of project files in the study. See 'help readwmdf' for details 
%             of a single DATA entry sub-structure. The collections are determined by columns in the
%             records spreadsheet, i.e. all data from project files listed in the ith column of the
%             spreadsheet is collected into DATA{i}.
%
%   STUDY:
%
%   1st-level
%   -------------------
%
%   STUDY.record   - The file name (full path) of the study records
%
%   STUDY.data_i   - The indices for each animal within the returned DATA structure. (This is to
%                    deal with the fact that the order of animals in the spreadsheet and STUDY
%                    structure may not correspond to the order in the data structure.)
%
%   STUDY.ANIMAL   - A structure with information about each of the animals.
%
%   STUDY.DAYS     - A cell array with information about each day of water-maze for the animals.
%
%   STUDY.FILE     - A structure with information about the files associated to this study.
%
%   STUDY.PROJECT  - A cell array of cell arrays. Each cell array contains a collection of grouped
%                    projects and each cell array within those contains information about wach 
%                    water-maze project file. See 'help readwmpf' for details and sub-structure.
%
%   STUDY.EXTRAS   - A cell array with any extra entries found in the records spreadsheet. 
%
%   2nd-level
%   -------------------
%
%   STUDY.ANIMAL.n     - The total number of animals in the study.
%   STUDY.ANIMAL.id    - A cell array with each animal's unique id
%   STUDY.ANIMAL.cage  - A vector with the LAS number for each animal's cage.
%   STUDY.ANIMAL.tag   - A cell array with the tag for each animal in the cage.
%   STUDY.ANIMAL.GROUP - A structure containing info about the study groups (see below).
%   STUDY.ANIMAL.sex   - The sex of the animals. This is an optional field, see the optional input
%                         'track_sex' below.
%   STUDY.ANIMAL.dob   - The dob of the animals. This is an optional field, see the optional input
%                         'track_dob' below.
%
%   STUDY.FILE.directory   - The directory where all the data files for this study are stored (as
%                             passed in by the user).
%   STUDY.FILE.name        - A cell array of each unique .wmpf file's name.
% 
%   3rd-level
%   -------------------
%   STUDY.ANIMAL.GROUPS.vars   - A cell array with the names for each grouping variable.
%   STUDY.ANIMAL.GROUPS.values - A cell array containing n x 1 cell arrays (where n is the number
%                                 of animals in the study) with each entry in the cell arrays
%                                 indicating that animal's grouping value. Note that GROUPS.vars and
%                                 GROUPS.values can be used in the Matlab ANOVAN function. See 'help
%                                 anovan' for more details.
%
% OPTIONAL INPUTS:
% -------------------------------------------------------------------------------------------------
% Optional inputs can be provided in parameter/value format. The optional inputs are:
%
%   'cage_in_id'    - Boolean value, indicating whether the cage number is included in the animal's
%                     id in the data files. Default = true.
%
%   'centre_origin' - Boolean value, indicating whether or not to make the centre of the pool equal to the origin in the
%                     co-ordinate system (i.e. x = 0, y = 0). Default = true.
%
%   'pool_radius'   - Scalar, indicating the actual pool radius (in cm) in case scaing of the data is
%                     requested. Default = 60 cm.
%
%   'scale_data'    - Boolean value, indicating whether or not to scale the data to reflect the
%                     actual size of the pool. Default = true.
%
%   'flip_data'     - Boolean value, indicating whether or not to make the centre of the pool equal to the origin in the
%                     co-ordinate system (i.e. x = 0, y = 0). Default = true.
%
%   'track_sex'     - Boolean value, indicating whether or not to search for a 'Sex' column after
%                     the 'Mouse' column in the spreadsheet. If set to TRUE then each animal's sex is
%                     recorded and stored in the ANIMAL.sex variable. Default = false.
%
%   'track_dob'     - Boolean value, indicating whether or not to search for a 'DOB' column after
%                     the 'Mouse' column in the spreadsheet. If set to TRUE then each animal's dob is
%                     recorded and stored in the ANIMAL.dob variable. Default = false.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% check the fnames argument
if ~isa(fname,'char')
	error('STUDY_RECORDS should be a string of the study records file name');
end

% check the datadir argument and store it
if ~exist(datadir,'dir')
	error('DATA_DIRECTORY must be a valid directory');
end
STUDY.FILE.directory = datadir;	

% define the default optional arguments
optargs = struct('centre_origin',true,...
                 'cage_in_id',true,...
                 'pool_radius',60,...
                 'scale_data',true,...
                 'flip_data',true,...
                 'track_sex',false,...
                 'track_dob',false);

% get the optional argument names
optnames = fieldnames(optargs);

% get the number of optional arguments
nargs = length(varargin)/2;

% make sure the property/value pairs are input in pairs
if round(length(varargin)/2) ~= nargs
   error('Expecting propertyName/propertyValue pairs after FILE and PATH');
end

% step through the optional arguments, check them, and store them
for pair = reshape(varargin,2,[])

	% make it case insensitive
	inpname = lower(pair{1});

	% check whether the name matches a known option
	if any(strmatch(inpname,optnames))
		switch inpname
			case 'centre_origin'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('centre_origin must be a logical');
				end
			case 'cage_in_id'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('cage_in_id must be a logical');
				end
			case 'pool_radius'
				if isa(pair{2},'numeric');
					optargs.(inpname) = pair{2};
				else
					error('pool_radius must be numeric');
				end
			case 'scale_data'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('scale_data must be a logical');
				end
			case 'flip_data'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('flip_data must be a logical');
				end
			case 'track_sex'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('track_sex must be a logical');
				end
			case 'track_dob'
				if isa(pair{2},'logical')
					optargs.(inpname) = pair{2};
				else
					error('track_dob must be a logical');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE SPREADSHEET

% make sure the spreadsheet exists and store it
if ~exist(fname,'file'), error('Could not find given study records'); end;
STUDY.record = fname;

% try and read the spreadsheet
try
	[num,txt,raw] = xlsread(fname);
catch
	error('Could not read the given spreadsheet');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% step through the columns of the spreadsheet and read in the data
cc = 1;

%% BASIC INFO %%

% make sure the first column is the cages column
if strcmp(raw(1,cc),'Cage #') ~= 1, error('First column must be ''Cage #'''); end;

% store the cage numbers
STUDY.ANIMAL.cage = cell2mat(raw(2:end,cc));
cc = cc + 1;

% store the n
STUDY.ANIMAL.n = length(STUDY.ANIMAL.cage);

% make sure the second column is the mouse column
if strcmp(raw(1,cc),'Mouse') ~= 1, error('Second column must be ''Mouse'''); end;

% store the cage numbers
STUDY.ANIMAL.tag = raw(2:end,cc);
cc = cc + 1;

% create each animal's unique ID
for aa = 1:length(STUDY.ANIMAL.cage)
	if optargs.cage_in_id
		if isstr(STUDY.ANIMAL.tag{aa})
			STUDY.ANIMAL.id{aa} = sprintf('%d%s',STUDY.ANIMAL.cage(aa),STUDY.ANIMAL.tag{aa});
		else
			STUDY.ANIMAL.id{aa} = sprintf('%d%d',STUDY.ANIMAL.cage(aa),STUDY.ANIMAL.tag{aa});
		end
	else
		if isstr(STUDY.ANIMAL.tag{aa})
			STUDY.ANIMAL.id{aa} = STUDY.ANIMAL.tag{aa};
		else
			STUDY.ANIMAL.id{aa} = sprintf('%d',STUDY.ANIMAL.tag{aa});
		end
	end	
end

% if requested, get each animal's sex
if optargs.track_sex
	if strcmp(raw(1,cc),'Sex') ~= 1, error('Column after ''Mouse'' must be ''Sex'''); end;
	STUDY.ANIMAL.sex = raw(2:end,cc);
	cc = cc + 1;
end

% if requested, get each animal's date of birth
if optargs.track_dob
	if strcmp(raw(1,cc),'DOB') ~= 1, error('Column after ''Mouse'' or ''Sex'' must be ''DOB'''); end;
	for tt = 1:size(raw,1)-1
		STUDY.ANIMAL.dob{tt} = datestr(x2mdate(cell2mat(raw(tt+1,cc))),'dd mmmm yyyy');
	end
	cc = cc + 1;
end

%% GROUP INFO %%
gg = 1;
while strcmp(raw(1,min(cc,end)),'File') ~= 1
	STUDY.ANIMAL.GROUP.vars{gg}   = cell2mat(raw(1,cc));
	STUDY.ANIMAL.GROUP.values{gg} = raw(2:end,cc);
	gg = gg + 1;
	cc = cc + 1;
end

%% WATER-MAZE DATA INFO %%
ff = 1;
file_notes = {};
for ww = cc:size(raw,2)
	
	% check whether this entry is a file entry
	if strcmp(raw(1,ww),'File') == 1
	
		% get the file's info
		STUDY.FILE.name{ff} = raw(2:end,ww);
		if strcmp(raw(1,ww+1),'Notes') == 1, STUDY.FILE.notes{ff} = raw(2:end,ww+1); end;
		ff = ff + 1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE PROJECT FILES

% load each of the project file collections and make sure the animals are listed in the correct
% order in the spreadsheet
for cc = 1:length(STUDY.FILE.name)

	% get the unique project names
	[unique_projects, ind, original_order] = unique(STUDY.FILE.name{cc});

	% load the projects
	STUDY.PROJECT{cc} = readwmpf(unique_projects,STUDY.FILE.directory);

	% check that the order of the animals corresponds with the order in the spreadsheet
	aa = 1;
	for oo = 1:length(original_order)

		% reset the animal counter if necessary
		if oo > 1 && (original_order(oo) ~= original_order(oo-1)), aa = 1; end;

		% check that the animal is in the spreadsheet
		if strcmp(STUDY.ANIMAL.id{oo},STUDY.PROJECT{cc}{original_order(oo)}.animal{aa}) == 1
			aa = aa + 1;
		else
			error('The animal %s from project %s is not located in the spreadsheet where they should be',...
						STUDY.PROJECT{cc}{original_order(oo)}.animal{aa},STUDY.PROJECT{cc}{original_order(oo)}.file);
		end

		% store the animal's order in the spreadsheet
		if isfield(STUDY.PROJECT{cc}{original_order(oo)},'sheet_index')
			STUDY.PROJECT{cc}{original_order(oo)}.sheet_index = [STUDY.PROJECT{cc}{original_order(oo)}.sheet_index ,oo];
		else
			STUDY.PROJECT{cc}{original_order(oo)}.sheet_index = oo;
		end
	end
end

% store the animal's indices in the projects in an easy to access location
for cc = 1:length(STUDY.PROJECT)
	STUDY.data_i{cc} = [];
	for pp = 1:length(STUDY.PROJECT{cc})
		STUDY.data_i{cc} = [STUDY.data_i{cc},STUDY.PROJECT{cc}{pp}.sheet_index];
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE DATA FILES

% for each project collection...
for cc = 1:length(STUDY.PROJECT)

	% initialize the data collection
	DATA{cc} = {};
	for pp = 1:length(STUDY.PROJECT{cc})

		% construct the file names for the .wmdf files
		clear data_files;
		for aa = 1:length(STUDY.PROJECT{cc}{pp}.animal)
			data_files{aa} = sprintf('%s.wmdf',STUDY.PROJECT{cc}{pp}.animal{aa});
		end

		% load the data from this project
		new_data = readwmdf(data_files,STUDY.PROJECT{cc}{pp}.folder,'centre_origin',optargs.centre_origin,...
																											'pool_radius',optargs.pool_radius,...
																											'scale_data',optargs.scale_data,...
																											'flip_data',optargs.flip_data);
		% add the data to the existing collection
		DATA{cc} = [DATA{cc}, new_data];
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FUNCTION
end
