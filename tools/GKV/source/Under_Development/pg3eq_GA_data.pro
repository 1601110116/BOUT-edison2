;
; Reads data from the the netcdf files (*.x.nc, *.h.nc and *.c.nc) 
; created by pg3eq using Dan Shumaker's NETCDF_pg3eq package,and  
; uses this data to create GKV data objects
;
;  Written by W.M. Nevins 1/17/99
;



Function pg3eq_data, path=file_path, file=file_name, var=var_name_in
;
; Use netcdf_pg3eq to read data in to structure array
;
data_in=netcdf_pg3eq(path=file_path, file=file_name, var=var_name_in)

file_name=data_in.file_name					; Get name of data file

separator="/"								; Path separator for unix devices
IF (!D.Name EQ "MAC") then separator=":"			; or for MAC's
substrings=STRSPLIT(file_name, separator, /Extract)	; Break file_name into substrings
n_strings=N_ELEMENTS(substrings)
filename=substrings[n_strings-1]				; Strip leading directory names of filename
separator="."
subsubstrings=STRSPLIT(filename, separator, /Extract)	; Break filename at dots
Run_name=subsubstrings[0]						; Run_name is leading piece of "filename"

num_vars=data_in.num_vars
objects=OBJARR(num_vars)						; Create object array
FOR ivar=0,num_vars-1 DO BEGIN
	var_structure=data_in.vars[ivar]
	num_dims=var_structure.num_dims
	IF (num_dims eq 1) THEN BEGIN				; Get signal data from data_in structure
		signal_data={GKVs1D}					; Create named structure
		signal_data.mnemonic	= var_structure.var_name
		signal_data.Title	= var_structure.plot_title
;		signal_data.units	= ???				; NEED to get signal units from pg3eq!!!
		signal_data.values	= var_structure.values
		signal_data.CodeName	= "pg3eq"			; CodeName is ALWAYS pg3eq for this function
		signal_data.CodePI 	= "A. Dimits"		; pg3eq PI is Dimits
		signal_data.RunID 	= Run_name			; RunID stripped out of "filename"
		signal_data.FileID	= filename			; FileID is filename (with directories stripped out)		

		xmin=MIN(*var_structure.grids[0].grid)
		xmax=MAX(*var_structure.grids[0].grid)
		imax=N_ELEMENTS(*var_structure.grids[0].grid) - 1
		signal_data.xtitle	= var_structure.grids[0].grid_label
		signal_data.xunits	= var_structure.grids[0].grid_units
		signal_data.xgrid	= var_structure.grids[0].grid
		signal_data.xbndry	= "open"			; time-dimension is "open"		
		signal_data.xmin	= xmin
		signal_data.xmax	= xmax
		signal_data.imin	= 0
		signal_data.imax	= imax
;
;  Register object, and store in object array
;
		objects[ivar]=Obj_New("GKVs1D", signal_data)
	ENDIF 
	IF (num_dims eq 2) THEN BEGIN
		signal_data={GKVs2D}					; Get signal data from data_in structure
		signal_data.mnemonic	= var_structure.var_name
		signal_data.Title	= var_structure.plot_title
;		signal_data.units	= ???				; NEED to get signal units from pg3eq!!!
		signal_data.values	= var_structure.values
		signal_data.CodeName	= "pg3eq"			; CodeName is ALWAYS pg3eq for this function
		signal_data.CodePI 	= "A. Dimits"		; pg3eq PI is Dimits
		signal_data.RunID 	= Run_name			; RunID stripped out of "filename"
		signal_data.FileID	= filename			; FileID is filename (with directories stripped out)		

		xmin=MIN(*var_structure.grids[0].grid)
		xmax=MAX(*var_structure.grids[0].grid)
		imax=N_ELEMENTS(*var_structure.grids[0].grid) - 1
		signal_data.xtitle	= var_structure.grids[0].grid_label
		signal_data.xunits	= var_structure.grids[0].grid_units
		signal_data.xgrid	= var_structure.grids[0].grid
		signal_data.xbndry	= "periodic"		; x-dimension is periodic in pg3eq		
		signal_data.xmin	= xmin
		signal_data.xmax	= xmax
		signal_data.imin	= 0
		signal_data.imax	= imax
		
		ymin=MIN(*var_structure.grids[1].grid)
		ymax=MAX(*var_structure.grids[1].grid)
		jmax=N_ELEMENTS(*var_structure.grids[1].grid) - 1		
		signal_data.ytitle	= var_structure.grids[1].grid_label
		signal_data.yunits	= var_structure.grids[1].grid_units
		signal_data.ygrid	= var_structure.grids[1].grid
		signal_data.ybndry	= "open"			; time-dimension is "open"
		signal_data.ymin	= ymin
		signal_data.ymax	= ymax
		signal_data.jmin	= 0
		signal_data.jmax	= jmax
;
;  Register object, and store in object array
;
		objects[ivar]=Obj_New("GKVs2D", signal_data)
	ENDIF 
ENDFOR
;
; Pack data objects into structure using mnemonics as tag names
;
command_string="structure={ "						; set up command_string
FOR ivar=0,num_vars-1 DO BEGIN						; ivar_string is a string containing
	ivar_string=STRING(FORMAT='(I2)',ivar)			; value of ivar (in ascii)
	command_string=	command_string + 			$	
				data_in.vars[ivar].var_name + 	$	; tag_name:objects[ivar]
				":" + "objects[" + ivar_string + "]"	
	IF(ivar LT num_vars-1) THEN 				$	; If not last object, then
				command_string= command_string + ", "	; add comma (as the delimiter between tag definitions) 
ENDFOR
command_string=command_string + "}"					; close bracket at end of structure definition statement
OK=EXECUTE(command_string)							; Execute command_string (thereby creating structure)
	

RETURN, structure
END 
	