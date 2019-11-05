 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  >> make_east_objects
%
%  PURPOSE: Script to calculate mutuals, Green functions for EAST 
%	tokamak system analogous to D3D Electromagnetic Environment. 
%	Creates save file with new objects (mutuals, Greens, geometry). 
%	Note that the "minimal" set of objects which must be defined in 
%	order to save a save-set is VV, FC, and grid.
%
%  INPUTS:
%    Filenames of data files: (make any of these null: '', if want to omit)
%	vvdata_file	Filename for vvdata data (must be of form *.data)
%			(actually just has to have 4 characters after ".", and
%			 must have a name distinct from all the other *.data)
%	fcdata_file	Fcoil data file (name must be of same form as vvdata)
%	ecdata_file	Ecoil data file (name must be of same form as vvdata)
%	flzr_file	Flux loop data file (name must be "     ")
%	bpdata_file	Bprobe data file (name must be of "     ")
%	rgmin,rgmax	Min,max in major radial dimension on plasma grid [m]
%	zgmin,zgmax	Min,max in vertical dimension on plasma grid [m]
%	nr,nz  # of grid elem. in r,z dir. (must be odd; make 0 to omit grid)
%	etav  		VV resistivity, uOhm-m  (can be col vector...)
%	etaf  		Fcoil resistivity, uOhm-m (can be col vector...)
%	etae  		Ecoil resistivity, uOhm-m (can be col vector...)
%
%  OUTPUTS:
%	Save file toksys_tmp.mat with all mutuals, Greens, and geometry data.
%
%  RESTRICTIONS:
%	At present 2 E-coil segments assumed, with ecdata_file format just like
%	DIII-D 2-segment E-coil (with indices 1,2 used to identify)  9/99
%	Note that to use with DIII-D,  fldata = [flzr; zeros(4,size(flzr,2))];
%	Formats of **data files are transpose of DIII-D standard (ie nthingsx6:
%	   where the 6 cols correspond to: Z,R,dZ,dR,ang1,ang2) except for
%	   ecdata (which is nelements x 5: Z,R,dZ,dR,segment_index). The 
%	   **data objects in Matlab are the transpose of these, ie 6xn,5xn.
%	Note that the "minimal" set of objects which must be defined in 
%	order to save a save-set is VV, FC, and grid data.
%	fcnturn_file name must be such that when load it, will have fcnturn = 
%	  vector of turns (so e.g. fcnturn.mat which contains vector fcnturn,
%	  or fcnturn.dat ascii file containing vector of data).
%
%  METHOD:  
%	Loads data files, defines 
%	geometry with define_toksys_data, then uses Mike's mutind_distrib to 
%	calculate mutuals. bgreens_distrib to calculate fields at B-probes.

%  WRITTEN BY:  Dave Humphreys  ON	3/31/04
%
%  MODIFICATION HISTORY:
%	3/31/04  DAH  Derived from make_iter_objects.m
%	3/04/05  MLW  "functionify" the code
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicit Inputs:
make_tok_inputs = struct( ...
'tokamak', 'east', ...
'datadir', [gatools_root '/tokamaks/east/make/'], ...
'limdata_file', 'limdata', ...
'vvdata_file', 'vvdata3', ...     
'fcdata_file', 'fcdata.mat', ...
'fcnames_file', 'fcnames', ...
'ecdata_file', '', ...
'flzr_file',   'fldata', ...	  
'flnames_file', 'flnames', ...
'bpdata_file', 'bpdata', ...
'bpnames_file', 'bpnames', ...
'fcnturn_file', 'fcnturn.mat', ... %ascii file with EAST PF turns 
'fcturn_file', 'fcturn.mat', ...
'turnfc_file', 'turnfc.mat', ...
'fcid_file', 'fcid.mat', ...
'rgmin', 1.20, ...     		%rgrid1 EFIT variable
'rgmax', 1.20+1.4, ... 		%rgrid1 + EFIT variable xdim
'zgmin', 0-2.4/2, ...		%zmid - zdim/2 EFIT variables
'zgmax', 0+2.4/2, ...		%zmid + zdim/2 EFIT variables
'nr', 65, ...			%EFIT var nw
'nz', 65, ...			%EFIT var nh
'etav', 7.4e-1, ... 		%VV res uOhm-m (316SS)
'etaf', 1.7e-2, ... 
'etae', 1.7e-2, ... %Ecoil resistivity, uOhm-m
'nminvv', 2, ...     %# of conductors in *min* dimension of VV
'nminfc', 3, ...     %# of conductors in *min* dimension of FC
'nmingg', 2, ...     %# of conductors in *min* dimension of GG
'drfl', 0.001, ...   %width of FL for mutuals to points (like grid ggdata)
'dzfl', 0.001)   %height of FL for mutuals to points (like grid ggdata)

make_tok_objects(make_tok_inputs,'plot_east_geo')

