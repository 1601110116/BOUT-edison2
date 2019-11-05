FUNCTION GKVs1D::Make_Spline, _extra=extra
;
; Returns GKVsnD second derivative of self.  Should be called only
; once when attempting to create a cubic spline interpolation of
; object.  
;
; SAMPLE Calls:  sp = gkv_object -> make_spline(axis=1)
;                print, gkv_object -> splint(x=value, axis=1)
;                print, gkv_object -> splint(x=value2, axis=1)
;                print, gkv_object -> splint(x=value3, axis=1)
; NOTES:         Be sure axis in splint matches axis in make_spline.  
;

nDims = self -> NumDims()
iaxis=0
result = GetKeyWord('axis', extra)
IF(Query_Integer(result)) THEN iaxis = result ; command line of form axis = axisnumber
IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result) ; command line of the form axis = 'mnemnic'
IF(iaxis EQ 0)  THEN BEGIN					; 'iaxis' is not yet set, so try to get axis ID from 'extra'
	IF(typeOf(extra) NE 8) THEN BEGIN			; No more unparsed keywords in Extra, so print error message and return.
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
	ENDIF
	axisInfo = self -> GetAxis(extra)
	FOR i=0, nDims-1 DO BEGIN
		iaxis = i+1
		IF(typeOf(axisInfo.(i)) NE 7) THEN BEGIN
			gridValues = axisInfo.(i)
			GOTO, DONE1
		ENDIF
	ENDFOR
;
; an axis mnemonic did not appear on the command line, so print error message and return
;
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
ENDIF
DONE1:

CASE iaxis OF
    1 : xaxis=self.grid1
    2 : xaxis=self.grid2
    3 : xaxis=self.grid3
    4 : xaxis=self.grid4
ENDCASE

;
; Note:xlength, xaxis and xunits refer to the user defined axis.  This will
; be separate from the variable xrange- used for grid1 only.
;

xLength = xaxis.range[1] - xaxis.range[0]
xUnits = xaxis.irange[1] - xaxis.irange[0]

result = self -> makecopy()
values = *(self.values)
x= *(xaxis.values)
dx = x[1]-x[0]
result_values = *(result.values)

CASE ndims OF
    1 :  BEGIN
        u = *(xaxis.values)
        f0=values[0]
        f1=values[1]
        f2=values[2]
        fN=values[xUnits]
        fN1=values[xUnits-1]
        fN2=values[xUnits-2]
        
        yp1=(-f2+8*f1-8*fN+fN1)/(12*dx)
        ypn=-(-f1+8*f0-8*fN1+fN2)/(12*dx)
        
        (result_values)[0]=-0.5
        u[0]=(3./(x[1]-x[0]))*((values[1]-values[0])/(x[1]-x[0])-yp1)
                
        FOR i = 1,xUnits-1 DO BEGIN
            sig = 1;(x[i]-x[i-1])/(x[i+1]-x[i])
            p = sig*result_values[i-1] + 2
            result_values[i] = 0; (sig - 1.)/p
            u[i] = (6 *((values[i+1]-values[i])/(x[i+1]-x[i])-(values[i]-values[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p
        ENDFOR
        
        Qn=.5
;        UN=(3/(x[xUnits]-x[xUnits-1]))*(YPN-(values[xUnits]-values[xUnits-1])/(x[xUnits]-x[xUnits-1]))        
;        result_values[xUnits]=(UN-qn*u[xUnits-1])/(qn*result_values[xUnits-1]+1)
        sig = 1
        result_values[xUnits]=0
        p=sig*result_values[xUnits-1]+2
        u[xUnits]= (6*((values[0]-values[xUnits])/dx-(values[xUnits]-values[xUnits-1])/dx)/(2*dx)-sig*u[xUnits-1])/p        
       
        FOR I= 0, xUnits-1 DO BEGIN
            k= xUnits-I-1
            result_values[k] =result_values[k]*result_values[k+1]+u[k]
        ENDFOR
    END
    2:  BEGIN
        xrange = self.grid1.irange[1]- self.grid1.irange[0]+1
        yrange = self.grid2.irange[1] - self.grid2.irange[0]+1
        CASE iaxis OF            
            1: BEGIN                
                u = fltarr(xrange,yrange); need to be a vector [xgrid,ygrid,etc]
                f0= values[0:*]
                f1=values[1,*]
                f2=values[2,*]
                fN=values[xUnits,*]
                fN1=values[xUnits-1,*]
                fN2=values[xUnits-2,*]

                yp1=(-f2+8*f1-8*fN+fN1)/(12*dx)
                ypn=(-f1+8*f0-8*fN1+fN2)/(12*dx)
                
                (result_values)[0,*]=-.5
                u[0,*]=(3./(x[1]-x[0]))*((values[1,*]-values[0,*])/(x[1]-x[0])-yp1)
                FOR i = 1,xUnits-1 DO BEGIN
                    sig = (x[i]-x[i-1])/(x[i+1]-x[i])
                    p = sig*result_values[i-1,*] + 2  
                    result_values[i,*] = (sig - 1.)/p
                    u[i,*] = (6 *((values[i+1,*]-values[i,*])/(x[i+1]-x[i])-(values[i,*]-values[i-1,*])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1,*])/p 
                ENDFOR
                
                Qn=.5
                ;UN=(3/(x[xUnits]-x[xUnits-1]))*(YPN-(values[xUnits,*]-values[xUnits-1,*])/(x[xUnits]-x[xUnits-1]))
                ;result_values[xUnits,*]=(UN-qn*u[xUnits-1,*])/(qn*result_values[xUnits-1,*]+1) 
                sig = 1
                result_values[xUnits,*]=0
                p=sig*result_values[xUnits-1,*]+2
                u[xUnits,*]= (6*((values[0,*]-values[xUnits,*])/dx-(values[xUnits,*]-values[xUnits-1,*])/dx)/(2*dx)-sig*u[xUnits-1,*])/p        
                

                FOR I= 0, xUnits-1 DO BEGIN
                    k= xUnits-I-1
                    result_values[k,*] =result_values[k,*]*result_values[k+1,*]+u[k,*]
                ENDFOR
            END
            2: BEGIN
                u = fltarr(xrange,yrange)
                f0= values[*,0]
                f1=values[*,1]
                f2=values[*,2]
                fN=values[*,xUnits]
                fN1=values[*,xUnits-1]
                fN2=values[*,xUnits-2]

                yp1=(-f2+8*f1-8*fN+fN1)/(12*dx)
                ypn=-(-f1+8*f0-8*fN1+fN2)/(12*dx)
                
                (result_values)[*,0]=-.5
                u[*,0]=(3./(x[1]-x[0]))*((values[*,1]-values[*,0])/(x[1]-x[0])-yp1) 
                FOR i = 1,xUnits-1 DO BEGIN
                    sig = (x[i]-x[i-1])/(x[i+1]-x[i])
                    p = sig*result_values[*,i-1] + 2  
                    result_values[*,i] = (sig - 1.)/p
                    u[i,*] = (6 *((values[*,i+1]-values[*,i])/(x[i+1]-x[i])-(values[*,i]-values[*,i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[*,i-1])/p 
                ENDFOR
                
                Qn=.5
                UN=(3/(x[xUnits]-x[xUnits-1]))*(YPN-(values[*,xUnits]-values[*,xUnits-1])/(x[xUnits]-x[xUnits-1]))
                result_values[*,xUnits]=(UN-qn*u[*,xUnits-1])/(qn*result_values[*,xUnits-1]+1)   ; 
                
                FOR I= 0, xUnits-1 DO BEGIN
                    k= xUnits-I-1
                    result_values[*,k] =result_values[*,k]*result_values[*,k+1]+u[*,k]
                ENDFOR
            END
        ENDCASE
    END
ENDCASE
PTR_FREE, result.values
result.values = PTR_NEW(result_values)
vmin = GKVsd_Min(result_values, Max=vmax)
result.vrange=[vmin,vmax]

RETURN, result
END;   ***********  GKVs1d::Make_Spline ************

FUNCTION GKVs1D::Splint,x=x,spline=spline,_extra=extra
;
;
;

nDims = self -> NumDims()
;
; get axis to be interpolated
;
iaxis=0
result = GetKeyWord('axis', extra)
IF(Query_Integer(result)) THEN iaxis = result ; command line of form axis = axisnumber
IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result) ; command line of the form axis = 'mnemnic'
IF(iaxis EQ 0)  THEN BEGIN					; 'iaxis' is not yet set, so try to get axis ID from 'extra'
	IF(typeOf(extra) NE 8) THEN BEGIN			; No more unparsed keywords in Extra, so print error message and return.
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
	ENDIF
	axisInfo = self -> GetAxis(extra)
	FOR i=0, nDims-1 DO BEGIN
		iaxis = i+1
		IF(typeOf(axisInfo.(i)) NE 7) THEN BEGIN
			gridValues = axisInfo.(i)
			GOTO, DONE1
		ENDIF
	ENDFOR
;
; an axis mnemonic did not appear on the command line, so print error message and return
;
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
ENDIF
DONE1:


CASE iaxis OF
    1 : xaxis=self.grid1
    2 : xaxis=self.grid2
    3 : xaxis=self.grid3
    4 : xaxis=self.grid4
ENDCASE

xLength = xaxis.range[1] - xaxis.range[0]
xUnits = xaxis.irange[1] - xaxis.irange[0]
xgrid = *(xaxis.values)
values=*(self.values)

KLo=0
KHi=xUnits

result= GetKeyWord('kHi',extra)
IF(Query_Integer(result)) THEN kHi=result
result= GetKeyWord('kLo',extra)
IF(Query_Integer(result)) THEN kLo=result
IF (kHi LT kLo) THEN BEGIN
    MESSAGE, 'input kHi must be larger than input kLo'/INFORMATIONAL
    RETURN, 0
ENDIF 

IF(kHi GT xUnits) THEN kHi=xUnits
IF(kLo LT 0)THEN kLo =0

IF (xgrid[kHi] LT x) THEN BEGIN
    kHi = xUnits    
ENDIF  
IF (xgrid[kLo] GT x) THEN BEGIN
    kLo = 0
ENDIF

RETRY1: 
IF(kHi- kLo GT 1) THEN BEGIN       
    k= (khi+klo)/2
    IF (xgrid[k] GT x) THEN BEGIN
        kHi = k
    ENDIF ELSE BEGIN
        kLo = k
    ENDELSE
    GOTO, RETRY1
ENDIF



nDims = self -> NumDims()
spvalues = *(spline.values)

H=xgrid[khi]-xgrid[klo]
A= (xgrid[khi]-x)/H
B=(x-xgrid[klo])/h
output= a*values[klo]+b*values[khi]+((a^3-a)*spvalues[klo]+(B^3-b)*spvalues[khi])*(h^2)/6

RETURN, output
END;************ GKVs1D::Splint****************

FUNCTION GKVs1D::PlotSpline, _extra=extra

nDims = self -> NumDims()

;
; get axis to be interpolated
;
iaxis=0
result = GetKeyWord('axis', extra)
IF(Query_Integer(result)) THEN iaxis = result ; command line of form axis = axisnumber
IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result) ; command line of the form axis = 'mnemnic'
IF(iaxis EQ 0)  THEN BEGIN					; 'iaxis' is not yet set, so try to get axis ID from 'extra'
	IF(typeOf(extra) NE 8) THEN BEGIN			; No more unparsed keywords in Extra, so print error message and return.
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
	ENDIF
	axisInfo = self -> GetAxis(extra)
	FOR i=0, nDims-1 DO BEGIN
		iaxis = i+1
		IF(typeOf(axisInfo.(i)) NE 7) THEN BEGIN
			gridValues = axisInfo.(i)
			GOTO, DONE1
		ENDIF
	ENDFOR
;
; an axis mnemonic did not appear on the command line, so print error message and return
;
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
ENDIF
DONE1:

n=2
result = GetKeyWord('n', extra)
IF(Query_Integer(result)) THEN n = result

CASE iaxis OF
    1 : xaxis=self.grid1
    2 : xaxis=self.grid2
    3 : xaxis=self.grid3
    4 : xaxis=self.grid4
ENDCASE

spline = self ->make_spline(axis=iaxis)

xLength = xaxis.range[1] - xaxis.range[0]
xUnits = xaxis.irange[1] - xaxis.irange[0]
xgrid = *(xaxis.values)
selfvalues=*(self.values)
dx = xgrid[2]-xgrid[1]

finalRange= n*xunits-1

finalGrid = gkvsd_gridcopy(xaxis)
PTR_FREE, finalGrid.values
finalGrid.irange[1] = 2*xUnits -1
finalgridvalues = FLTARR(2*xUnits+1)
values = FLTARR(2*xUnits+1)


FOR i = 0, 2*xUnits DO BEGIN
IF (i MOD 2) EQ 0 THEN BEGIN
finalGridvalues[i] = xgrid[i/2]
values[i] = selfvalues[i/2]
ENDIF ELSE BEGIN
xPosition = xgrid[(i-1)/2]+dx/2
finalGridvalues[i]= xPosition
values[i] = self -> splint(axis=iaxis, spline=spline,x=xPosition)
ENDELSE
ENDFOR

output = self -> makecopy()
PTR_FREE, finalGrid.values
finalGrid.values =PTR_NEW(finalgridvalues)
finalGrid.irange[1] = xaxis.irange[1]+ n * (xUnits-1)
output.grid1 = finalGrid
PTR_FREE, output.values
output.values = PTR_NEW(values)

RETURN, output
END;


FUNCTION GKVs1D::DummySpline, x=x
ndims = self -> numdims()
IF (ndims NE 1) THEN BEGIN
    MESSAGE, 'dummyspline only works on 1D objects', /INFORMATIONAL
    RETURN, 0
ENDIF

yvalues = *(self.values)
xvalues = *(self.grid1.values)
output = SPLINE(xvalues,yvalues,x)

RETURN, output
END;  ***  GKVs1D::DummySpline  ***


FUNCTION GKVs1D::Splint_Vector,delx=delx,spline=spline,_extra=extra
;
;
;

IF (delx EQ 0) THEN BEGIN
output = *(self.values)
GOTO, DONE2
ENDIF

nDims = self -> NumDims()
;
; get axis to be interpolated
;
iaxis=0
result = GetKeyWord('axis', extra)
IF(Query_Integer(result)) THEN iaxis = result ; command line of form axis = axisnumber
IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result) ; command line of the form axis = 'mnemnic'
IF(iaxis EQ 0)  THEN BEGIN					; 'iaxis' is not yet set, so try to get axis ID from 'extra'
	IF(typeOf(extra) NE 8) THEN BEGIN			; No more unparsed keywords in Extra, so print error message and return.
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
	ENDIF
	axisInfo = self -> GetAxis(extra)
	FOR i=0, nDims-1 DO BEGIN
		iaxis = i+1
		IF(typeOf(axisInfo.(i)) NE 7) THEN BEGIN
			gridValues = axisInfo.(i)
			GOTO, DONE1
		ENDIF
	ENDFOR
;
; an axis mnemonic did not appear on the command line, so print error message and return
;
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN,0
ENDIF
DONE1:


CASE iaxis OF
    1 : xaxis=self.grid1
    2 : xaxis=self.grid2
    3 : xaxis=self.grid3
    4 : xaxis=self.grid4
ENDCASE

xLength = xaxis.range[1] - xaxis.range[0]
xUnits = xaxis.irange[1] - xaxis.irange[0]
xgrid = *(xaxis.values)
values=*(self.values)

KLo=0
KHi=1

output = fltarr(xUnits+1)

spvalues = *(spline.values)

dx = xgrid[1]-xgrid[0]
b=dx-delx

CASE ndims OF
    1: BEGIN 
        output = delx*values[*]+(dx-delx)*SHIFT(values,1)+((delx^3-delx)*spvalues+(b^3-b)*SHIFT(spvalues,1))*dx^2/6    
    END
    2:  BEGIN
        CASE iaxis OF
            1: BEGIN
                output = delx*values[*,*]+(dx-delx)*SHIFT(values,1,0)+((delx^3-delx)*spvalues+(b^3-b)*SHIFT(spvalues,1,0))*dx^2/6
                PRINT, self.values - output
            END
            2: BEGIN
                output = delx*values[*,*]+(dx-delx)*SHIFT(values,0,1)+((delx^3-delx)*spvalues+(b^3-b)*SHIFT(spvalues,0,1))*dx^2/6
            END
        ENDCASE
    END
    3: BEGIN
        CASE iaxis OF
            1: BEGIN
                output = delx*values[*,*,*]+(dx-delx)*SHIFT(values,1,0,0)+((delx^3-delx)*spvalues+(b^3-b)*SHIFT(spvalues,1,0,0))*dx^2/6
            END
            2: BEGIN
                output = delx*values[*,*,*]+(dx-delx)*SHIFT(values,0,1,0)+((delx^3-delx)*spvalues+(b^3-b)*SHIFT(spvalues,0,1,0))*dx^2/6
            END
            3: BEGIN
                output = delx*values[*,*,*]+(dx-delx)*SHIFT(values,0,0,1)+((delx^3-delx)*spvalues+(b^3-b)*SHIFT(spvalues,0,0,1))*dx^2/6
            END
        ENDCASE
    END
ENDCASE
;PTR_FREE, xgrid
;PTR_FREE, values

;H=xgrid[khi]-xgrid[klo]
;A= (xgrid[khi]-x)/H
;B=(x-xgrid[klo])/h
;output= a*values[klo]+b*values[khi]+((a^3-a)*spvalues[klo]+(B^3-b)*spvalues[khi])*(h^2)/6


DONE2:
RETURN, output
END;************ GKVs1D::Splint_vector****************
