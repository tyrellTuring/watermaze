function DATA = readwmdf(fnames, pname, varargin)

% function DATA = readwmdf(FILES, PATH, varargin)
%
% Reads in watermaze data from files in the binary format produced by Actimetrics software (.wmdf files).
%
% MANDATORY INPUTS:
% ------------------------------------------------------------------------------------------------------
%
%   FILES - A cell array with the names of all the files to be read in.
%   PATH - A string of the path to the folder containing all the files (note: it is assumed they're in the same folder).
%
% OUTPUT:
% ------------------------------------------------------------------------------------------------------
% The output of the function is a cell array, DATA, organized as follows:
%
%   1st-level
%   -------------------
%   The first level of the cell array ranges over files, e.g. DATA{3} is the data from the third file in the list.
%
%
%   2nd-level
%   -------------------
% 
%   DATA{i}.animal  - The name of the animal associated with this file (which is in fact the file name, according to
%                     Actimetrics convention).
%
%   DATA{i}.file    - The file the data is from (full path name).
%
%   DATA{i}.ntrials - The number of trials stored in this file.
%
%   DATA{i}.ntimes  - A 1 x N vector containing the number of time-steps for each trial (N = number of trials).
%
%   DATA{i}.path    - A T x 3 x N array of the path data of the animal in the maze during each trial. The first dimension
%                     ranges over time (T = number of time-steps), the second dimension is time (in seconds), x-position, 
%                     and y-position (in cm), and the third dimension ranges over trials (N = number of trials). For 
%                     example, DATA{2}.path(:,2,8) is the x-position of the path data for the eighth trial of the second 
%                     file. When there are not an equal number of time-steps in each trial, then the shorter trials are nan-
%                     padded in the first dimension.
%
%   DATA{i}.pool    - Either a 1 x 3 vector of the pool dimensions, or a N x 3 matrix of the dimensions of the pool for 
%                     each trial. An N x 3 matrix is only returend if the pool data was different on different trials. The
%                     columns contain the co-ordinates of the pool centre (x then y) and the radius of the pool (in cm).
%
%   DATA{i}.platform - Either a 1 x 3 vector of the platform dimensions, or a N x 3 matrix of the platform dimensions for 
%                      each trial. An N x 3 matrix is only returend if the platform data was different on different trials. The
%                      columns contain the co-ordinates of the platform centre (x then y) and the radius of the platform (in cm).
%
%
% OPTIONAL INPUTS:
% ------------------------------------------------------------------------------------------------------
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
%--------------------------------------------------------------------------------
%
% 12/2011, Frankland Lab (www.franklandlab.com)
%
% Author: Hamid Maei
% Modified by: Kirill Zaslavsky & Blake Richards
% Contact: blake.richards@utoronto.ca
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% check the fnames argument
if ~isa(fnames,'cell')
	error('FILES should be a cell array of file names');
end

% check the pname argument
if ~exist(pname,'dir')
	error('PATH must be a valid directory');
end
	
% define the default optional arguments
optargs = struct('centre_origin',true,...
                 'pool_radius',60,...
                 'scale_data',true,...
                 'flip_data',true);

% get the optional argument names
optnames = fieldnames(optargs);

% get the number of optional arguments
nargs = length(varargin)/2;

% make sure the property/value pairs are input in pairs
if round(length(varargin)/2) ~= nargs
   error('Expecting propertyName/propertyValue pairs after FILES and PATH');
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
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP THROUGH THE FILES AND LOAD EACH ONE

% determine the number of files we're processing
if iscell(fnames), nfiles = length(fnames); else, nfiles = 1; end;

% for each file...
for ff = 1:nfiles

	% create the file identifier
	if iscell(fnames), filename = fullfile(pname,fnames{ff}); else, filename = fullfile(pname,fnames); end; 
  	try
		fid = fopen(filename);
	catch
		error('Could not open file %s',filename);
	end

	% determine the size of the file in bytes
	try
		status    = fseek(fid, 0, 'eof');
	catch
		error('Could not move to end of file %s',filename);
	end                  
  file_size = ftell(fid);

	% go to the beginning of the first trial in the file        
	status      = fseek(fid, 0, 'bof');
	if status < 0, error('Could not move to start of file %s',filename); end;                  
	trial_start = 4;
	status      = fseek(fid,trial_start,'bof');
	if status < 0, error('Could not move to fourth byte in file %s',filename); end;                  

	% determine the number of trials in the file and the size of each trial (in bytes)
	ntrials = 0;
	status 	= 0;
	trial_size = 0; 
	trial_size(1) = fread(fid,1,'uint32',0,'ieee-be');
	while trial_size(end) > 0 && status == 0
		status                  = fseek(fid,trial_size(end),'cof');
		trial_size(ntrials + 1) = fread(fid,1,'uint32',0,'ieee-be');
		ntrials                 = ntrials + 1;
	end
	
	% store all the information about this file
	[garb1, DATA{ff}.animal, garb2] = fileparts(filename);
	DATA{ff}.file           = filename;
	DATA{ff}.ntrials        = ntrials;

	% initialize the data holders
	DATA{ff}.ntimes   = zeros(1,ntrials);
	DATA{ff}.path     = nan(10000,3,ntrials);
	DATA{ff}.pool     = zeros(ntrials,3);
	DATA{ff}.platform = zeros(ntrials,3);

	% read in each trial's data
	for tt = 1:ntrials
			
		% go to the start of the trial
		fseek(fid, trial_start, 'bof');
        	trial_size(tt) = fread(fid,1,'uint32', 0, 'ieee-be'); %number of bytes in the record for a trial

		% if we have a non-empty trial...
		if trial_size(tt) > 0

			% get the number of elements in the trial
			try
				n_elements = fread(fid,1,'uint32', 0, 'ieee-be'); 
			catch
				error('Could not read the number of elements for trial %d in file %s',tt,filename);
			end
	
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% FIRST ELEMENT OF THE TRIAL 
			
			% skip the first element (Blake says: why?)
			status = fseek(fid,12, 'cof');
			if status < 0, error('Could not move past first element for trial %d in file %s',tt,filename); end;                  

	
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% SECOND ELEMENT OF THE TRIAL (PATH)

			% get the size of the element in bytes
			try
				element_size = fread(fid, 1, 'uint32', 0, 'ieee-be');
			catch
				error('Could not read in the size of the path data for trial %d in file %s',tt,filename);
			end
	
			% get the size of the array for the path data
			try
				array_size = fread(fid, 2, 'uint32', 0, 'ieee-be');
			catch
				error('Could not read in the size of the path array for trial %d in file %s',tt,filename);
			end

			% read in the path info for this trial
			try
				thispath = fread(fid,[array_size(2) array_size(1)],'float32',0,'ieee-be');                
			catch
				error('Could not read path info for trial %d in file %s',tt,filename);
			end

			% store the number of time-steps for this trial
			DATA{ff}.ntimes(tt) = size(thispath,1);
	
			% store the path data for this trial	
			DATA{ff}.path(1:size(thispath,1),:,tt) = thispath;			

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% THIRD ELEMENT OF THE TRIAL
			
			% skip the third element (Blake says: why?)
			try
				element_size = fread(fid,1,'uint32',0,'ieee-be');
			catch
				error('Could not read in the size of the third element for trial %d in file %s',tt,filename);
			end
			status = fseek(fid,element_size,'cof');
			if status < 0, error('Could not move past third element for trial %d in file %s',tt,filename); end;                  
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% FOURTH ELEMENT OF THE TRIAL (POOL)

			% get the start position of this element in the file
			startpos = ftell(fid);

			% get the size of the element in bytes
			try
				element_size = fread(fid, 1, 'uint32', 0, 'ieee-be');
			catch
				error('Could not read in the size of the pool data for trial %d in file %s',tt,filename);
			end
		
			% get the name of the pool
			try	
				strlen    = fread(fid,1,'uint32',0,'ieee-be');
				pool_name = fread(fid,strlen,'*char',0,'ieee-be'); % we actually just ditch this right now...
			catch
				error('Could not read in the name of the pool for trial %d in file %s',tt,filename);
			end
	
			% get the coordinates and radius of the pool
			try
				x_pool = fread(fid,1,'double',0,'ieee-be');
				y_pool = fread(fid,1,'double',0,'ieee-be');
				r_pool = fread(fid,1,'double',0,'ieee-be');
			catch
				error('Could not read pool information for trial %d in file %s',tt,filename);
			end
			DATA{ff}.pool(tt,:) = [x_pool y_pool r_pool];

			% get the current position in the file
			curpos = ftell(fid);

			% move to the end of this element
			status = fseek(fid,element_size-(curpos-startpos)+4,'cof'); 
			if status < 0, error('Could not move past fourth element for trial %d in file %s',tt,filename); end;                  
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% FIFTH ELEMENT OF THE TRIAL
			
			% skip the fifth element (Blake says: why, what is this element?)
			try
				element_size = fread(fid,1,'uint32',0,'ieee-be');
			catch
				error('Could not read in the size of the fifth element for trial %d in file %s',tt,filename);
			end	
			status = fseek(fid,element_size,'cof');
			if status < 0, error('Could not move past fifth element for trial %d in file %s',tt,filename); end;                  

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% SIXTH ELEMENT OF THE TRIAL (PLATFORM)

			% get the size of the element in bytes
			try
				element_size = fread(fid, 1, 'uint32', 0, 'ieee-be');
			catch
				error('Could not read in the size of the platform data for trial %d in file %s',tt,filename);
			end
		
			% get the name of the platform
			try	
				strlen        = fread(fid,1,'uint32',0,'ieee-be');
				platform_name = fread(fid,strlen,'*char',0,'ieee-be'); % we actually just ditch this right now...
			catch
				error('Could not read in the name of the platform for trial %d in file %s',tt,filename);
			end
	
			% get the coordinates and radius of the platform
			try
				x_plat = fread(fid,1,'float64',0,'ieee-be');
				y_plat = fread(fid,1,'float64',0,'ieee-be');
				r_plat = fread(fid,1,'float64',0,'ieee-be');
			catch
				error('Could not read platform information for trial %d in file %s',tt,filename);
			end
			DATA{ff}.platform(tt,:) = [x_plat y_plat r_plat];

		end

		% move the position in the file to the start of the next trial's data
		trial_start = trial_start + trial_size(tt) + 4;
	end

	% eliminate redundant data in the pool and platform info
	if sum(diff(DATA{ff}.pool))     == 0, DATA{ff}.pool     = DATA{ff}.pool(1,:);     end;
	if sum(diff(DATA{ff}.platform)) == 0, DATA{ff}.platform = DATA{ff}.platform(1,:); end;

	% close the file identifier
	fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECT THE DATA

% for each file...
for ff = 1:length(DATA)

	% reduce the size of the path data
	DATA{ff}.path = DATA{ff}.path(1:max(DATA{ff}.ntimes),:,:);

	% put all the measurements into cm (they appear to be in half-cm)
	DATA{ff}.path(:,2:3,:) = DATA{ff}.path(:,2:3,:)*0.5;
	DATA{ff}.pool          = DATA{ff}.pool*0.5;
	DATA{ff}.platform      = DATA{ff}.platform*0.5;

	% if the centre_origin option is on, correct the spatial data for each file
	if optargs.centre_origin

		% flags to help not writing over already corrected data
		centre_set   = false;
		platform_set = false;

		% for each trial...
		for tt = 1:DATA{ff}.ntrials

			% get the centre of the pool
			if ~centre_set && size(DATA{ff}.pool,1) == 1
				centre = DATA{ff}.pool(1:2);
				centre_set = true;
			elseif size(DATA{ff}.pool,1) > 1
				centre = DATA{ff}.pool(tt,1:2);
			end

			% set the centre of the pool to zero
			DATA{ff}.pool(min([tt end]),1:2) = [0 0];

			% correct the path and platform data
			DATA{ff}.path(1:DATA{ff}.ntimes(tt),2:3,tt) = DATA{ff}.path(1:DATA{ff}.ntimes(tt),2:3,tt) - repmat(centre,[DATA{ff}.ntimes(tt) 1]);
			if ~platform_set && size(DATA{ff}.platform,1) == 1
				DATA{ff}.platform(1:2) = DATA{ff}.platform(1:2) - centre;
				platform_set = true;
			elseif size(DATA{ff}.platform,1) > 1
				DATA{ff}.platform(tt,1:2) = DATA{ff}.platform(tt,1:2) - centre;
			end
		end
	end

	% if the scale_data option is on, scale the spatial data for each file
	if optargs.scale_data

		% calculate the scale factor
		scale = optargs.pool_radius/DATA{ff}.pool(3);

		% scale the pool, platforms and paths
		DATA{ff}.pool          = DATA{ff}.pool*scale;
		DATA{ff}.platform      = DATA{ff}.platform*scale;
		DATA{ff}.path(:,2:3,:) = DATA{ff}.path(:,2:3,:)*scale;
	end

	% if the flip_data option is on, flip the spatial data for each file
	if optargs.flip_data

		% flip the platforms and paths
		DATA{ff}.platform(:,2,:) = -DATA{ff}.platform(:,2,:);
		DATA{ff}.path(:,3,:)     = -DATA{ff}.path(:,3,:);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FUNCTION
end
