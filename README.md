MWM Toolbox
=========

=-----------------------
Intro
=-----------------------
The MWM Matlab Toolbox is a series of Matlab functions for analyzing Morris Water-Maze data. To use
it, you must have Matlab version 7 or greater. Additionally, if you want to be able to use some of
the more advanced analysis methods you will need to install Alexander Ihler's Kernel Density
Estimation Toolbox (http://www.ics.uci.edu/~ihler/code/kde.html).

This toolbox is provided freely under the GNU Lesser Public License. See the file COPYING.LESSER in
the toolbox, or alternatively visit http://www.gnu.org/licenses/.

The MWM Toolbox was written by Blake Richards (blake.richards@utoronto.ca).

=-----------------------
Installation
=-----------------------
See the file INSTALL

=-----------------------
Functions
=-----------------------
Functions in the MWM Toolbox are organized into four broad groups: read functions, analysis
functions, get functions and plot functions.

(1) Read functions
Read functions are all named read*.m These functions are used to read in data, either from user
designed spreadsheets or from raw water-maze data files. These functions all return multi-level
structures for storing and analyzing water-maze data. The three types of structure are:

	i)   STUDY   - Stores information about the study, e.g. animals, dates, groups, etc.

	ii)  PROJECT - Stores information about water-maze "projects", as defined by the Actimetrics
	               water-maze software's .wmpf files (http://www.actimetrics.com/WaterMaze/). A 
	               project is essentially just a collection of water-maze trials for multiple 
	               animals. The PROJECT structure is what is returned by reading a .wmpf file. PROJECT
	               is a cell array, with each entry corresponding to a .wmpf project file.

	iii) DATA    - Stores the actual water-maze data. At its most basic this includes the animals' 
	               paths, the platform locations and pool size. This is also were analyzed data is
	               eventually stored, e.g. measures such as time in a quadrant, number of crossings,
	               entropy, etc. DATA is a cell array of cell arrays. The first level of cell arrays
	               organizes data from particular episodes. For example, DATA{1} could contain
	               training data and DATA{2} could contain probe data. How this is organized depends
	               on the user's input (see help readstudy and look at the 'collect_data' option to
	               learn how to control this organization). The second level of cell arrays iterates
	               over each animal's data. Thus, DATA{1}{8} might be the training data for animal
	               eight for example. This ultimately corresponds to .wmdf files.

(2) Analysis functions
Analysis functions are all named mwm*.m and they really provide the important features of the
toolbox. Each function uses the animal's path data to calculate a measure regarding the animal's
search path in the pool. For example, mwmquadrants provides a measurement of what percentage of
time the animal spends in each quadrant of the pool. All of these functions take a second level DATA
structure and return that structure with the measures embedded in it. For example, one could
calculate the number of crossings of a platform during the second data set with the command:

	DATA{2} = mwmcrossings(DATA{2});

The measurements are stored in each animal's cell array entry using standard names (e.g. Q for
quadrant, X for crossings or H for entropy).

(3) Get functions
Get functions are all named get*.m and are there in order to extract data from the DATA structure
into something more easily analyzable (e.g. for statistics). The output of these functions are
typically matrices of the measurements that range of all animals, trials, parameters, etc.

(4) Plot functions
Plot functions are all named plot*.m and provide tools for visualizing water-maze data.

=-----------------------
Usage
=-----------------------
The MWM Toolbox is made to be fairly easy to use, but it has some particular requirements. These are
stipulated below.

(1) Data organization
To use the MWM Toolbox your data must be organized in a particular way. Here are the requirements: 
First, the actual raw data must be from Actimetrics Water-maze program, and organized in the format 
that this program uses. This means that for each set of water-maze data there should be a .wmpf file, 
and a folder with the same name place the suffix 'Folder' containing .wmdf files for each animal.
For example, if you have an Actimetrics project file named 'My Project.wmpf' with animals M1,M2,M3,
etc. then there should also be a directory somewhere titled 'My Project Folder' with files M1.wmdf,
M2.wmdf, M3.wmdf, etc. Any departure from this scheme will lead to errors in reading data.

(2) Study spreadsheet
One important component of the 
