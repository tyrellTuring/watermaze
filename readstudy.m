function [STUDY,PROJECT,DATA] = readstudy(fname,datadir,varargin)

% function [STUDY,PROJECT,DATA] = readstudy(STUDY_RECORDS, DATA_DIRECTORY, varargin)
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
%                     - *'GROUP VAR ...' (there can be any number of group variables, can't be 'Date')
%                     - 'Date'              (for training or probe)
%                     - 'File'              (.wmpf file for training or probe)
%                     - 'Platform Location' (for training or probe)
%                     - 'Start Sequence'    (for training or probe)
%                     - 'Notes'             (for training or probe)
%                     - ...             (there can be any number of entries in this pattern)
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
%   PROJECT - PROJECT is a cell array with information about wach water-maze project file. See
%             'help readwmpf' for details and sub-structure.
%
%   DATA{i} - DATA is a cell array of cell arrays. The DATA{i} contains the raw water-maze data for 
%             a given collection of project files in the study. See 'help readwmdf' for details 
%             of a single DATA entry sub-structure. See the option 'collect_data' to understand
%             which project files will be included in each entry of DATA.
%
%   STUDY:
%
%   1st-level
%   -------------------
%
%   STUDY.record    - The file name (full path) of the study records
%
%   STUDY.ANIMAL   - A structure with information about each of the animals.
%
%   STUDY.DAYS      - A cell array with information about each day of water-maze for the animals.
%
%   STUDY.FILE     - A structure with information about the files associated to this study.
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
%   STUDY.DAYS{i}.date      - A cell array of the dates of the ith day in the study for each animal.
%   STUDY.DAYS{i}.file      - A cell array of the .wmpf files where the ith day's data is stored.
%   STUDY.DAYS{i}.platform  - A cell array of the platform's location on the ith day for each animal.
%   STUDY.DAYS{i}.start     - A cell array of the start sequence on the ith day for each animal.
%   STUDY.DAYS{i}.notes     - A cell array of any additional notes on the ith day's data for each animal.
%   
%   STUDY.FILE.directory   - The directory where all the data files for this study are stored (as
%                             passed in by the user).
%   STUDY.FILE.name        - A cell array of each unique .wmpf file's name.
%   STUDY.FILE.collect     - A cell array of .wmpf files that are collected together for ease of
%                             analysis (see the optional argument 'collect_data' below).
%
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
%   'collect_data'  - A cell array of cell arrays, indicating which data files should be collected
%                     together for easier joint analysis. For example, if the following is passed:
%
%                       {{'Training_Group_1.wmpf', 'Training_Group_2.wmpf'},...
%                        {'Probe_Group_1.wmpf', 'Probe_Group_2.wmpf'}}
%
%                     then the training files for Groups 1 and 2 will be collected into a single
%                     data entry, and the probe files will be collected into another.
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
                 'pool_radius',60,...
                 'scale_data',true,...
                 'flip_data',true,...
                 'track_sex',false,...
                 'track_dob',false,...
                 'collect_data',-1);

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
			case 'collect_data'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('collect_data must be a cell array');
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
	STUDY.ANIMAL.id{aa} = sprintf('%d%s',STUDY.ANIMAL.cage(aa),STUDY.ANIMAL.tag{aa});
end

% if requested, get each animal's sex
if optargs.track_sex
	if strcmp(raw(1,cc),'Sex') ~= 1, error('Column after ''Mouse'' must be ''Sex'''); end;
	STUDY.ANIMAL.sex = raw(2:end,cc);
	cc = cc + 1;
end

% if requested, get each animal's sex
if optargs.track_dob
	if strcmp(raw(1,cc),'DOB') ~= 1, error('Column after ''Mouse'' or ''Sex'' must be ''DOB'''); end;
	for tt = 1:size(raw,1)-1
		STUDY.ANIMAL.dob{tt} = datestr(x2mdate(cell2mat(raw(tt+1,cc))),'dd mmmm yyyy');
	end
	cc = cc + 1;
end

%% GROUP INFO %%
gg = 1;
while strcmp(raw(1,min(cc,end)),'Date') ~= 1
	STUDY.ANIMAL.GROUP.vars{gg}   = cell2mat(raw(1,cc));
	STUDY.ANIMAL.GROUP.values{gg} = raw(2:end,cc);
	gg = gg + 1;
	cc = cc + 1;
end

%% WATER-MAZE DATA INFO %%
dd = 1;
ff = 1;
for ww = cc:size(raw,2)
	
	% check whether this entry is a day entry
	if strcmp(raw(1,ww),'Date')   == 1 && strcmp(raw(1,ww+1),'File') == 1
	
		% get the day's info
		for tt = 1:size(raw,1)-1
			try
				STUDY.DAYS{dd}.date{tt} = datestr(x2mdate(cell2mat(raw(tt+1,ww))),'dd mmmm yyyy');
			catch
				STUDY.DAYS{dd}.date{tt} = NaN;
			end
		end
		STUDY.DAYS{dd}.file     = raw(2:end,ww+1);
		STUDY.DAYS{dd}.platform = raw(2:end,ww+2);
		STUDY.DAYS{dd}.start    = raw(2:end,ww+3);
		STUDY.DAYS{dd}.notes    = raw(2:end,ww+4);

		% if this is a new file store it in the FILE structure
		for fnew = 1:length(STUDY.DAYS{dd}.file)
			if isfield(STUDY.FILE,'name')
				if ~isnan(STUDY.DAYS{dd}.file{fnew})
					if ~ismember(STUDY.DAYS{dd}.file{fnew},STUDY.FILE.name)
						STUDY.FILE.name{ff} = STUDY.DAYS{dd}.file{fnew};
						ff = ff + 1;
					end
				end
			else
				STUDY.FILE.name{ff} = STUDY.DAYS{dd}.file{fnew};
				ff = ff + 1;
			end
		end

		% increment the counter
		dd = dd + 1;
	end
end

%% DATA COLLECTION INFO %%
if isa(optargs.collect_data,'cell')
	
	% check that each of the specified files in the collections exists in the records
	for cc = 1:length(optargs.collect_data)
		for ff = 1:length(optargs.collect_data{cc})
			infiles = false;
			for gg = 1:length(STUDY.FILE.name)
				if strcmp(optargs.collect_data{cc}{ff},STUDY.FILE.name{gg}) == 1
					infiles = true;
					break;
				end
			end
			if ~infiles
				error('''collect_data'' argument contained files not contained in study records');
			end
		end
	end
				
	% store the collection information
	STUDY.FILE.collect = optargs.collect_data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE PROJECT FILES

% load each of the project files
for pp = 1:length(STUDY.FILE.name)
	project_files{pp} = fullfile(STUDY.FILE.directory,STUDY.FILE.name{pp});
end
PROJECT = readwmpf(project_files);

% make sure each of the animals listed in the project files is listed in the study
for pp = 1:length(PROJECT)
	for aa = 1:length(PROJECT{pp}.animal)
		if ~ismember(PROJECT{pp}.animal{aa},STUDY.ANIMAL.id)
			error(sprintf('Animal %s in project %s not in records',PROJECT{pp}.animal{aa},PROJECT{pp}.file));
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE DATA FILES

% for each project
for pp = 1:length(PROJECT)

	% construct the file names for the .wmdf files
	for aa = 1:length(PROJECT{pp}.animal)
		data_files{aa} = sprintf('%s.wmdf',PROJECT{pp}.animal{aa});
	end

	% load the data from this project
	DATA{pp} = readwmdf(data_files,PROJECT{pp}.folder,'centre_origin',optargs.centre_origin,...
                                                    'pool_radius',optargs.pool_radius,...
                                                    'scale_data',optargs.scale_data,...
                                                    'flip_data',optargs.flip_data);
end

% collect the data together as requested
if isa(optargs.collect_data,'cell')
	for cc = 1:length(STUDY.FILE.collect)
		
		% determine which elements of DATA correspond to this collection
		incollection = [];
		for pp = 1:length(PROJECT)
			if ismember(PROJECT{pp}.file,STUDY.FILE.collect)
				incollection = [incollection, pp];
			end
		end
		
		% contactenate all of the members of this collection into a new DATA structure
		newDATA{cc} = {};
		for nn = 1:length(incollection)
			newDATA{cc} = {newDATA{cc}, DATA{incollection{nn}}};
		end
	end
	DATA = newDATA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FUNCTION
end
