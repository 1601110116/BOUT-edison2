FUNCTION GKV_Create_Object2, data, tValues=tValues, title=title, mnemonic=mnemonic
;
;  Create a GKV object out of an array (data) whose dimensions
;  are [nx,ny=1,nz, nt]
;
;	tValues		Should be set to a 1-D vector of time values
;			If there are more than 'nt' elements in tValues,
;			then the last nt values are used.
;
;	title		Should be set to a string variable containing
;			(in Hershey vector font, if desired) the 'title'
;			for the GKV object
;
;	mnemonic	Should be set to a string variable containing
;			the 'mnemonic' for the GKV object
;
;
;------------------------------------------



structure={gkvs3d}
info=SIZE(data)
nx=info[1]
nz=info[3]
nt=info[4]

xGrid={grid}
xGrid.mnemonic='r'
xGrid.title='r'
xgrid.units=''
xGrid.values=PTR_NEW(FINDGEN(nx))
xGrid.irange=[0,nx-1]
xGrid.range=[0,nx-1]

zGrid={grid}
zGrid.mnemonic='zeta'
zGrid.title='!4f!X'
zgrid.units=''
zGrid.values=PTR_NEW(FINDGEN(nz))
zGrid.irange=[0,nz-1]
zGrid.range=[0,nz-1]

tGrid={grid}
tGrid.mnemonic='t'
tGrid.title='t'
tgrid.units=''
ntValues=N_ELEMENTS(tValues)
IF(ntValues EQ 0) THEN ntValues = FINDGEN(nt)
tGrid.values=PTR_NEW(tValues[ntValues-nt:ntValues-1])
tGrid.irange=[0,nt-1]
tGrid.range=[tValues[ntValues-nt],tValues[ntValues-1]]

structure.values=PTR_NEW(REFORM(data))
structure.grid1=xGrid
structure.grid2=zGrid
structure.grid3=tGrid
structure.indices=PTR_NEW(['*', '*', '*'])

IF(N_ELEMENTS(title) NE 0) THEN structure.title=title
IF(N_ELEMENTS(mnemonic) NE 0) THEN structure.mnemonic=mnemonic

obj=OBJ_NEW('GKVs3D', structure)
return, obj
end


FUNCTION GKV_Create_Object1, data, tValues=tValues, title=title, mnemonic=mnemonic
;
;  Create a GKV object out of an array (data) whose dimensions
;  are [nx=1,ny,nz, nt]
;
;	tValues		Should be set to a 1-D vector of time values
;			If there are more than 'nt' elements in tValues,
;			then the last nt values are used.
;
;	title		Should be set to a string variable containing
;			(in Hershey vector font, if desired) the 'title'
;			for the GKV object
;
;	mnemonic	Should be set to a string variable containing
;			the 'mnemonic' for the GKV object
;
;
;------------------------------------------



structure={gkvs3d}
info=SIZE(data)
ny=info[2]
nz=info[3]
nt=info[4]

yGrid={grid}
yGrid.mnemonic='theta'
yGrid.title='!4h!x'
ygrid.units=''
yGrid.values=PTR_NEW(FINDGEN(ny))
yGrid.irange=[0,ny-1]
yGrid.range=[0,ny-1]

zGrid={grid}
zGrid.mnemonic='zeta'
zGrid.title='!4f!X'
zgrid.units=''
zGrid.values=PTR_NEW(FINDGEN(nz))
zGrid.irange=[0,nz-1]
zGrid.range=[0,nz-1]

tGrid={grid}
tGrid.mnemonic='t'
tGrid.title='t'
tgrid.units=''
ntValues=N_ELEMENTS(tValues)
IF(ntValues EQ 0) THEN ntValues = FINDGEN(nt)
tGrid.values=PTR_NEW(tValues[ntValues-nt:ntValues-1])
tGrid.irange=[0,nt-1]
tGrid.range=[tValues[ntValues-nt],tValues[ntValues-1]]

structure.values=PTR_NEW(REFORM(data))
structure.grid1=yGrid
structure.grid2=zGrid
structure.grid3=tGrid
structure.indices=PTR_NEW(['*', '*', '*'])

IF(N_ELEMENTS(title) NE 0) THEN structure.title=title
IF(N_ELEMENTS(mnemonic) NE 0) THEN structure.mnemonic=mnemonic

obj=OBJ_NEW('GKVs3D', structure)
return, obj
end


FUNCTION GKV_Create_Object, data, tValues=tValues, title=title, mnemonic=mnemonic
;
;  Create a GKV object out of an array (data) whose dimensions
;  are [nx=1,ny,nz, nt] or [nx, ny=1, nz, nt]
;
;	tValues		Should be set to a 1-D vector of time values
;			If there are more than 'nt' elements in tValues,
;			then the last nt values are used.
;
;	title		Should be set to a string variable containing
;			(in Hershey vector font, if desired) the 'title'
;			for the GKV object
;
;	mnemonic	Should be set to a string variable containing
;			the 'mnemonic' for the GKV object
;
;
;------------------------------------------

info=SIZE(data)

CASE 1 OF 

    info[1] EQ 1: RETURN, GKV_Create_Object1(data, tValues=tValues, title=title, mnemonic=mnemonic)

    info[2] EQ 1: RETURN, GKV_Create_Object2(data, tValues=tValues, title=title, mnemonic=mnemonic)

    ELSE : MESSAGE, "GKV_Create_Object called with wrong argument", /INFORMATIONAL

ENDCASE

RETURN, 0
END
