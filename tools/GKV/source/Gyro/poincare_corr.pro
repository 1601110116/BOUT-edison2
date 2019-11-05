FUNCTION Poincare_Full_Corr, filename

CD, CURRENT=current_working_directory

IF(N_ELEMENTS(fileName) EQ 0) THEN fileIn=DIALOG_PICKFILE(Path=path, Filter='*.log', Get_Path=NetCDF_DIR) ELSE fileIn=fileName
ok = FINDFILE(fileIn, Count=nfiles) ; Check if filename points to a valid file
IF (nfiles NE 1) THEN BEGIN     ; 	we didn't find a file... 
    MESSAGE, "Couldn't find .log file", /INFORMATIONAL
    RETURN, 0
ENDIF

fileName=fileIN

found:
n_iterations= 144000
output = FLTARR(4,3,n_iterations)
i=0L
j=-1L

GET_LUN, ioUnit
OPENR, ioUnit, fileName, ERROR=err

IF(err NE 0) THEN BEGIN
    PRINT, "Error opening run.out file:"
    PRINT, "Error Message = ", !ERR_STRING
    RETURN, output
ENDIF

thisLine=""
ReadLine:	IF( EOF(ioUnit) ) THEN GOTO, Done
READF, ioUnit, thisLine
IF (strpos(thisLine,'orbit') GT 0) THEN BEGIN
    i = 0L
    j = j+1
    GOTO, ReadLine
ENDIF

thisBr = strmid(thisLine, 1,12) + 0.
thisBPhi = strmid(thisLine, 14,12) + 0.
thisTheta = strmid(thisLine, 26,12) + 0.


output[j,0,i] = thisBR
output[j,1,i] = thisBPhi
output[j,2,i] = thistheta
i = i+1



GOTO, ReadLine

DONE:

FREE_LUN, ioUnit                

CD, current_working_directory

n_theta=1
arr = output[0,0,*]
;n_iterations = 48000
output2 = fltarr(n_iterations)
FOR j = 0L, n_iterations-1 DO BEGIN
    sum = complex(0,0)
    FOR i = 0L, (n_iterations-1)/48 DO BEGIN
        arg = j+i*n_theta
        IF arg GT n_iterations-1 THEN BEGIN
            sum = sum
        ENDIF ELSE BEGIN
            sum = sum + arr[i*n_theta]*arr[i*n_theta+j]        
        ENDELSE
    ENDFOR
    output2[j] = sum/n_iterations    
    IF (j mod 4800 EQ 0) THEN BEGIN
        PRINT, "j = ", j
    ENDIF
    ;PRINT, "j = ", j, " sum = ", sum   
ENDFOR

print, "j = ", j

gridStructure = {Grid}
gridStructure.Mnemonic = "Steps"
gridStructure.Title = "Steps"
gridValues= output[0,2,*]
gridValues = REFORM(gridvalues)
gridStructure.Values= PTR_NEW(INDGEN(n_iterations, /long)/48.) ;gridValues
nPoints = N_ELEMENTS(*gridStructure.Values)
gridStructure.Range=[0,npoints/48]
gridStructure.irange=[0,npoints-1]

objStructure = {GKVs1D}
objStructure.Grid1 = gridStructure
objStructure.mnemonic = 'Corr'
objStructure.title = 'Corr'
objStructure.indices = PTR_NEW('*')
objStructure.units = ''
objValues = output2[*]/output2[0]
objValues = REFORM(objValues)
objStructure.values = PTR_NEW(objValues)
objStructure.vrange = [MIN(objValues), MAX(objValues)]
obj = obj_new("gkvs1d", objstructure)


RETURN, obj
END

FUNCTION GKVs1D::Poincare_Corr, j0, j_max

arr = *self.values
n_iterations = N_ELEMENTS(arr)
sum = complex(0,0)
i_c = complex(0,1)
;normarr = FFT(arr,1)
n_theta=1
output = complexarr(j_max+1)
FOR j = 0, j_max DO BEGIN
    sum = complex(0,0)
    FOR i = 0, n_iterations-1 DO BEGIN
        arg = j0+j+i*n_theta
        IF arg GT n_iterations-1 THEN BEGIN
            sum = sum
        ENDIF ELSE BEGIN
        sum = sum + arr[j0+i*n_theta]*arr[j0+i*n_theta+j]        
        ENDELSE
    ENDFOR
    output[j] = sum/n_iterations    
    ;PRINT, "j = ", j, " sum = ", sum   
ENDFOR

objOutput = self -> MakeCopy()
PTR_FREE, objOutput.values
objOutput.values = ptr_new(output)
objOutput.title = "Corr"
objOutput.mnemonic = "Corr"
objOutput.units=''

objOutput.vrange = [MIN(output), MAX(output)]
RETURN, objOutput
END


FUNCTION Poincare_Corr, j0, j_max, arr

n_theta=8
n_iterations=100
sum = complex(0,0)
i_c = complex(0,1)
j = 0
output = complexarr(j_max+1)

FOR j = 0, j_max DO BEGIN
    FOR i = 0, n_iterations-1 DO BEGIN
        l_bar = FLOOR((j0 + j+0.) / n_theta)
        j_bar = (j0+j) MOD n_theta
        sum = sum + arr[j0,i]*CONJ(arr[j_bar,i])*exp(-2*!pi*i_c*l_bar*i/n_iterations)
    ENDFOR
    output[j] = sum    
     PRINT, "j = ", j, " sum = ", sum
    sum = complex(0,0)
   
ENDFOR

;output = 1/n_iterations^2 * sum

RETURN, output
END


FUNCTION Poincare_Corr2, j0, j_max, arr

n_theta=8
n_iterations=1000
i_c = complex(0,1)
j = 0
output = complexarr(j_max+1)

normarr = FFT(arr,1)

arr1d = REFORM(arr,800)

FOR j = 0, j_max DO BEGIN
    sum = complex(0,0)
    FOR i = 0, n_iterations-1 DO BEGIN
        ;l_bar = FLOOR((j0 + j+0.) / n_theta)
        ;j_bar = (j0+j) MOD n_theta
        ;a = (j0) MOD n_theta
        ;b = ((j0+i)/n_theta)
        ;c = (a+j) MOD n_theta
        ;d = FLOOR((j0+j)/n_theta) + i
        ;PRINT, a, b, c, d
        arg = j0+i*n_theta+j
        IF arg GT 7999 THEN BEGIN
            sum = sum
        ENDIF ELSE BEGIN
        sum = sum + arr1d[j0+i*n_theta]*arr1d[j0+i*n_theta+j]        
        ENDELSE
    ENDFOR
    output[j] = sum/n_iterations    
    ;PRINT, "j = ", j, " sum = ", sum   
ENDFOR

;output = 1/n_iterations^2 * sum

RETURN, output
END
