PRO XGKC_Check_Names_CleanUp, tlb
;
; Nothing to be done.  In particular, DON'T free up
;	info.checkPtr, as then there would be no way
;	to retrieve the 'check' data for return to 
;	the calling routine.

END ; ****** XGKC_Check_Names_CleanUp ****** ;

PRO XGKC_Check_Names_Quit, event 
   Widget_control, event.top, /Destroy    ; Destroy widget structure 
END ; ****** XGKC_Check_Names_Quit ****** ; 

PRO XGKC_Check_Names_All, event 

WIDGET_CONTROL, event.top, Get_UValue=info, /No_Copy        	; Get user info 
nNames = N_ELEMENTS(info.names)
(*info.checksPtr) = REPLICATE(1, nNames)
WIDGET_CONTROL, event.top, Set_UValue=info, /No_Copy        	; Set user info 
WIDGET_CONTROL, event.top, /Destroy    ; Destroy widget structure 
END ; ****** XGKC_Check_Names_All ****** ; 


PRO XGKV_Check_Names_Resize, event
;
; Event processor to resize the XGKV_Check_Names widget
;
; Begin by getting user info
;
WIDGET_CONTROL,   event.top, Get_UValue=info, /No_Copy
;
; Get 'geometry' information for TLG, and the two buttons
;
tlbGeometry  = WIDGET_INFO(event.top,   /GEOMETRY)
allGeometry  = WIDGET_INFO(info.allID,  /GEOMETRY)
doneGeometry = WIDGET_INFO(info.doneID, /GEOMETRY)
;
; Dimensions of resized tlb widget window are in the structure 'event'
;
newXsize = event.x > 100
newYsize = event.y - allGeometry.SCR_YSIZE - doneGeometry.SCR_YSIZE
;
; Compute new ysize
;
newysize=newysize - 2*tlbgeometry.ypad - 2*tlbgeometry.space
;
; Reset width of 'ALL' button
;
WIDGET_CONTROL, info.allID,  SCR_XSize=newXsize
;
; Reset width and height of list of names
;
WIDGET_CONTROL, info.listID, SCR_XSize=newXsize, SCR_YSize=newYsize
;
; Reset width of 'DONE' button
;
WIDGET_CONTROL, info.doneID, SCR_XSize=newXsize
;
; Return user info
;
WIDGET_CONTROL, event.top, Set_UValue=info, /No_Copy 

END ; ****** XGKV_Check_Names_Resize ****** ;


PRO XGKV_Check_Names_Event, event
;
; Event processor for list widget in XGKV_Check_Names
;
; Set character for 'checkmark'
;
checkmark = '+'
;
; Get user info
;
Widget_Control, event.top, Get_UValue=info, /No_Copy 
;
; Toggle element of 'checks' array touched by user
;      	
(*info.checksPtr)[event.INDEX] = (*info.checksPtr)[event.INDEX] EQ 0
;
; Make array of check marks, chectTxt
;
nNames = N_ELEMENTS(info.names)
checkTxt = REPLICATE(" ", nNames)
FOR i=0, nNames-1 DO BEGIN
	IF( (*info.checksPtr)[i] EQ 1) THEN checkTxt[i]=checkmark
ENDFOR
;
; Make list with appropriate names 'checked'
;
checkedList = checkTxt + info.names
;
; Update elemenets of list widget
;
WIDGET_CONTROL, info.listID, SET_VALUE=checkedList
;
; Return user info
WIDGET_CONTROL, event.top, SET_UValue=info, /No_Copy

RETURN
END ; ****** XGKV_Check_Names_Event ****** ;

FUNCTION XGKV_Check_Names, parent, Names=names,  Window_Title=windowTitle 
;
;
;	Purpose:
;			
;		This widget program queries user regarding list of names supplied as first argument
;
;	Input Arguments:
;
;		parent		The widget ID of the parent widget.
;					(Optional)
;
;	Input Keywords:
;
;		Names			a string array containing names to be queried
;					(Required)
;
;		Window_Title	a string variable used as the title of the window.
;					(Optinoal)
;
;	Output:
;
;		Checks		A byte array whose value is '1' if the corresponding element of 
;					'names' has been selected by the user, and '0' otherwise.
; 	 
;  Written by W.M. Nevins
;	8/16/00
;
IF(TypeOF(names) NE 7) THEN BEGIN
	MESSAGE, "No 'names' provided'", /INFORMATIONAL
	RETURN, 0
ENDIF
nNames = N_ELEMENTS(names)
CASE N_PARAMS() OF
	0:	tlb = WIDGET_BASE(        COLUMN=1, TITle=windowTitle, TLB_FRAME_ATTR=8, /TLB_SIZE_EVENTS)
	1:	tlb = WIDGET_BASE(parent, COLUMN=1, TITle=windowTitle, TLB_FRAME_ATTR=8, /TLB_SIZE_EVENTS)
ENDCASE
allID  = WIDGET_BUTTON(tlb, VALUE='ALL',  EVENT_PRO='XGKC_Check_Names_All' )
checkTxt = REPLICATE(' ', nNames)
checkedList = checkTxt + names
checks = BYTARR(nNames)
checksPtr = PTR_NEW(checks)
listID   = WIDGET_LIST(tlb, /MULTIPLE, VALUE=checkedList, EVENT_PRO='XGKV_Check_Names_Event', SCR_XSIZE=100, SCR_YSIZE=300)
doneID = WIDGET_BUTTON(tlb, VALUE='DONE', EVENT_PRO='XGKC_Check_Names_Quit')
info={listID:listID, allID:allID, doneID:doneID, names:names, checksPtr:checksPtr}
Widget_Control, tlb, Set_UValue=info, /No_Copy ; Save info structure 
Widget_Control, tlb, /Realize			; Put widget on screen
XMANAGER, 'XGKV_Check_Names', tlb, EVENT_HANDLER='XGKV_Check_Names_Resize', CLEANUP="GKC_Check_Names_CleanUp"
checks = *checksPtr
PTR_FREE, checksPTR
RETURN, checks
END ; ****** XGKV_Check_Names ****** ;

