MWM Matlab Toolbox
=========

Intro
-----------------------
The MWM Matlab Toolbox is a series of Matlab functions for analyzing Morris Water-Maze data (see
Morris (1984), Journal of Neuroscience Methods, 11(1):47-60 or Maei et al. (2009), Frontiers in 
Integrative Neuroscience, 3:4). To use
it, you must have Matlab version 7 or greater. Additionally, if you want to be able to use some of
the more advanced analysis methods you will need to install Alexander Ihler's Kernel Density
Estimation Toolbox (http://www.ics.uci.edu/~ihler/code/kde.html).

This toolbox is provided freely under the GNU Lesser Public License. See the file COPYING.LESSER in
the toolbox, or alternatively visit http://www.gnu.org/licenses/.

The MWM Toolbox was written by Blake Richards (blake.richards@utoronto.ca).

Installation
-----------------------
See the file INSTALL

Functions
-----------------------
Functions in the MWM Toolbox are organized into four broad groups: read functions, analysis
functions, get functions and plot functions.

### Read functions
Read functions are all named read\*.m These functions are used to read in data, either from user
designed spreadsheets or from raw water-maze data files. These functions all return multi-level
structures for storing and analyzing water-maze data. The three types of structure are:

1. STUDY   - Stores information about the study, e.g. animals, dates, groups, etc.

2. PROJECT - Stores information about water-maze "projects", as defined by the Actimetrics
	water-maze software's .wmpf files (http://www.actimetrics.com/WaterMaze/). A 
	project is essentially just a collection of water-maze trials for multiple 
	animals. The PROJECT structure is what is returned by reading a .wmpf file. PROJECT
	is a cell array, with each entry corresponding to a .wmpf project file.

3. DATA    - Stores the actual water-maze data. At its most basic this includes the animals' 
	paths, the platform locations and pool size. This is also were analyzed data is
	eventually stored, e.g. measures such as time in a quadrant, number of crossings,
	entropy, etc. DATA is a cell array of cell arrays. The first level of cell arrays
	organizes data from particular episodes. For example, DATA{1} could contain
	training data and DATA{2} could contain probe data. How this is organized depends
	on the organization of the 'records spreadsheet' (see below). The second level of 
	cell arrays iterates over each animal's data. Thus, DATA{1}{8} might be the training 
	data for animal eight for example. This ultimately corresponds to .wmdf files.

### Analysis functions
Analysis functions are all named mwm\*.m and they really provide the important features of the
toolbox. Each function uses the animal's path data to calculate a measure regarding the animal's
search path in the pool. For example, mwmquadrants provides a measurement of what percentage of
time the animal spends in each quadrant of the pool. All of these functions take a second level DATA
structure and return that structure with the measures embedded in it. For example, one could
calculate the number of crossings of a platform during the second data set with the command:

	DATA{2} = mwmcrossings(DATA{2});

The measurements are stored in each animal's cell array entry using standard names (e.g. Q for
quadrant, X for crossings or H for entropy).

### Get functions
Get functions are all named get\*.m and are there in order to extract data from the DATA structure
into something more easily analyzable (e.g. for statistics). The output of these functions are
typically matrices of the measurements that range of all animals, trials, parameters, etc.

### Plot functions
Plot functions are all named plot\*.m and provide tools for visualizing water-maze data.

Usage Requirements
-----------------------
The MWM Toolbox is made to be fairly easy to use, but it has some particular requirements. These are
stipulated below.

### Animal names
Within the actual Actimetrics water-maze project files the user must give a name to each animal in
the study. For these files to then be usable with this toolbox the names of the animals must follow
the convention (CAGE\_NUMBER)(ANIMAL). For example, if animal M1 from cage 5021 was included in a
project that animal must be named 5021M1 in the Actimetrics project file.

### Data organization
To use the MWM Toolbox your data must be organized in a particular way. Here are the requirements: 
First, the actual raw data must be from Actimetrics Water-maze program, and organized in the format 
that this program uses. This means that for each set of water-maze data there should be a .wmpf file, 
and a folder with the same name plus the suffix 'Folder' containing .wmdf files for each animal.
*This folder must be in the same directory as the .wmdf file.* For example, if you have an Actimetrics 
project file named 'My Project.wmpf' with animals 5021M1,5021M2,5021M3, etc. then there should also be a directory 
in the same folder titled 'My Project Folder' with files 5021M1.wmdf, 5021M2.wmdf, 5021M3.wmdf, etc. Any departure 
from this scheme will lead to errors in reading data.

### Study spreadsheet
One of the most important requirements that must be met to use the toolbox successfully is a
properly formatted "study spreadsheet". This is a spreadsheet that lists all the animals in the
study, their groupings (e.g. which experimental treatments they received) and the files that contain
their data. The spreadsheet must be a .xls file (otherwise a reading error may result).

*The organization of this spreadsheet must follow a specific convention.* The columns must have
specifc headings (not case-sensitive though) and they must be present in a specific order. Some columns are optional but others
aren't. The columns and the order that they must be present in are as follows (optional columns are
indicated):

- Cage #
	This is the cage number for the animal. *This must be numeric!*

- Animal
	This is the identifier for this animal in this cage.

- Sex (Optional)
	This is the animal's sex.

- DOB (Optional)
	This is the animal's date of birth.

- Grouping variables (Optional)
	These can be multiple columns of any name other than 'File'. These columns will determine the
	groups that the animal's get included in during analysis.

- File
	This is the name of each .wmpf file for the animals' data. *It must include the extension.* There
	can be multpile 'File' columns. Each column defines which data is grouped together in the output
	DATA structure (see the readstudy.m help).

- Notes (Optional)
	This contains any notes on that animal's data in the given file. Notes columns must follow File
	columns.

Usage Example
--------------
A very simple example of using the toolbox would be to use it to estimate whether there is a
difference in the amount of time spent in the correct quadrant during probe trials. So, let's say
that we had a records spreadsheet named 'wmdata.xls' in the directory '/home/timmy/data/' with the
following entries:

	Cage #, Animal, Treatment,          File,       File
	  5021,     M1, Lidocaine, Training.wmpf, Probe.wmpf
	  5021,     M2,   Vehicle, Training.wmpf, Probe.wmpf
	  5022,     M1, Lidocaine, Training.wmpf, Probe.wmpf
	  5022,     M2,   Vehicle, Training.wmpf, Probe.wmpf
	etc.

We could then determine the mean time spent in the correct quadrant during the first probe trial for each group which the
following sequence of commands:

	>> [STUDY,DATA] = readstudy('wmdata.xls','/home/timmy/data');
	>> DATA{2} = mwmquadrants(DATA{2});
	>> Q = getquadrants(DATA{2});
	>> lido = getgroup(STUDY,{'Lidocaine'});
	>> vehi = getgroup(STUDY,{'Vehicle'});
	>> mean(Q(1,1,lido)), mean(Q(1,1,mean))

Usage Help
--------
For more help on using this toolbox see the help messages contained in each of the files.

