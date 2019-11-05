
FUNCTION GKVs1D::Velocity_Transform, v=v, _extra=extra 
;
; FUNCTION form of radial transform.
; 
; Sample call:  transformed_data = data -> Velocity_Transform(vel=3,axis=2)
;
; Written by E. Wang 
;     7/21/09
;
copy = self -> makeCopy()
copy -> Velocity_Transform, v=v, _extra=extra
RETURN, copy
END;   ******** Velocity_Transform *******

PRO GKVs1D::Velocity_Transform, v=v, _extra=extra
;
; Adjusts values of self by specified velocity along specified axis
;
; NOTE:  does not work if time grid is not last grid.  Assumes
; supplied axis is uniform.  Also assumses periodic(open) boundary conditions.
;
; Sample calls: data -> Velocity_Transform,vel=3,axis=1
;               data -> Velocity_Transform,vel=-2,axis='r'
;
; Written by E. Wang
;     7/21/09

nDims = self -> NumDims()
IF (ndims LT 2) THEN BEGIN
    MESSAGE, 'Radial_Transform called with less than 2 dimensions'
    RETURN
ENDIF

;
; get axis to be transformed
;
iaxis=0
result = GetKeyWord('axis', extra)
IF(Query_Integer(result)) THEN iaxis = result ; command line of form axis = axisnumber
IF(typeOf(result) EQ 7)  THEN iaxis = self -> AxisNumber(result) ; command line of the form axis = 'mnemnic'
IF(iaxis EQ 0)  THEN BEGIN					; 'iaxis' is not yet set, so try to get axis ID from 'extra'
	IF(typeOf(extra) NE 8) THEN BEGIN			; No more unparsed keywords in Extra, so print error message and return.
		MESSAGE, 'No Valid axis ID', /INFORMATIONAL
		RETURN
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
		RETURN
ENDIF
DONE1:

IF iaxis EQ nDims THEN BEGIN  ; Make sure input axis is not time axis (time assumed to be last grid)
    MESSAGE, 'input axis cannot be time axis', /INFORMATIONAL
RETURN
ENDIF

CASE iaxis OF
    1 : rgrid=self.grid1 
    2 : rgrid=self.grid2
    3 : rgrid=self.grid3 
    4 : rgrid=self.grid4 
ENDCASE

;
; Assume last grid is time grid.
;

CASE ndims OF
    1 : tgrid=self.grid1
    2 : tgrid=self.grid2
    3 : tgrid=self.grid3
    4 : tgrid=self.grid4
ENDCASE

;IF(rgrid.uniform EQ 0) THEN BEGIN
;    MESSAGE, 'nonuniform radial grid, will not adjust', /INFORMATIONAL
;    RETURN
;ENDIF

;
; Basics about radial and time grid, get values
;

rLength = rgrid.range[1] - rgrid.range[0]
rUnits = rgrid.irange[1] - rgrid.irange[0]
dr = rLength / rUnits
rLength = rLength + dr

tLength = tgrid.range[1] - tgrid.range[0]
tUnits = tgrid.irange[1] - tgrid.irange[0]
dt = tLength/tUnits

values = *(self.values)


;
; Begin modifying values.
;
FOR i = 0,tUnits DO BEGIN
tStep = i*dt
distance = v * tStep MOD rLength + (v LT 0)*rLength
steps = FIX(distance/dr)
overstep = (distance- steps*dr)/dr
CASE ndims OF
    2: BEGIN
        offset = values[*,i]
        temp1 = SHIFT(offset, steps)
        temp2 = SHIFT(offset, steps+1)
    END
    3: BEGIN 
        offset = values[*,*,i]
        CASE iaxis OF
            1: BEGIN
                temp1 = SHIFT(offset,[steps,0])
                temp2 = SHIFT(offset,[steps+1,0])
            END
            2: BEGIN
                temp1 = SHIFT(offset,[0,steps])
                temp2 = SHIFT(offset,[0,steps+1])
            END
        ENDCASE
    END
    4: BEGIN 
        offset = values[*,*,*,i]
        CASE iaxis OF
            1: BEGIN
                temp1 = SHIFT(offset,steps,0,0)
                temp2 = SHIFT(offset,[steps+1,0,0])
            END
            2: BEGIN
                temp1 = SHIFT(offset,[0,steps,0])
                temp2 = SHIFT(offset,[0,steps+1,0])
            END
            3:  BEGIN
                temp1 = SHIFT(offset,[0,0,steps])
                temp2 = SHIFT(offset,[0,0,steps+1])
            END
        ENDCASE
    END
ENDCASE

final = (1-overstep)*temp1 + overstep*temp2
CASE ndims OF
    2: BEGIN        
        values[*,i] = final
    END
    3: BEGIN
        values[*,*,i] = final
    END
    4: BEGIN
        values[*,*,*,i] = final
    END
ENDCASE

ENDFOR

;
; Save new values
;

PTR_FREE, self.values
ptr = PTR_NEW(values)
self.values= ptr

RETURN
END;  *******  GKVs1D::Velocity_Transform *****
