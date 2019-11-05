Pro GKV_View_Resize, event  Widget_Control,   event.top, Get_UValue=info, /No_Copy Widget_Control,   info.drawID, Draw_XSize=event.x > 200, $                   Draw_YSize=event.y > 200 WSet, info.wID TvLct, info.r, info.g, info.b    ; Reload previous color table info.object -> Draw, _Extra=*info.extra Widget_Control, event.top, Set_UValue=info, /No_Copy end  Pro GKV_View_CleanUp, tlb ; ; final duties...this is called as widget structure melts away... ; (see call to XManager in Pro XDisplay) ;    Widget_Control, tlb, Get_UValue=info, /No_Copy    if N_Elements(info) ne 0 then $ 		PTR_FREE, info.extra;end Pro GKV_View_Quit, event    Widget_control, event.top, /Destroy    ; Destroy widget structure end  Pro GKV_View_Colors, event ; ; Change color table ; Widget_Control, event.top, Get_UValue=info, /No_Copy      ; Get UValues  thisEvent=Tag_Names(event, /Structure_Name) CASE thisevent of    'WIDGET_BUTTON'   : Begin       XColors, NColors=info.ncolors,  $  ; Change Color Table                Group_Leader=event.top, NotifyID=[event.id, event.top], $                Title='ChangeColors' + '(' + StrTrim(info.wID,2) + ')'                         end ; XColors will notify 'color change button' ; (i.e., THIS proceedure!) when colors are changed. ; This will result in an 'XCOLORS_LOAD' event    'XCOLORS_LOAD' : begin          info.r=event.r    ; load new color table          info.b=event.b    ; into 'info' for          info.g=event.g    ; color protection          WSet, info.wID      info.object -> Draw, _Extra=*info.extra                        end ENDCASE    Widget_Control, event.top, Set_UValue=info, /No_Copy ; return info to UValue  end Pro GKV_View_PostScript, event; ; output image as postscript file ;Widget_Control, event.top, Get_UValue=info, /No_Copy		; Get user info keys=ps_form(Cancel=cancelled, filename='GKV_View.ps')	; use PS_Form to define PostScript file if cancelled then return thisDevice=!D.name set_plot, 'ps'                                              ; Make current device Device, _extra=keys                                         ; a PostScript file TvLct, info.r, info.g, info.b    ; Reload previous color table info.object -> draw, _Extra=*info.extraDevice, /close_file              ; Close PostScript file set_plot, thisDevice             ; Replace original 'device' WSet, info.wID                   ; Return to original window (is this necessary???) Widget_Control, event.top, Set_Uvalue=info, /No_Copy        ; Return user info end  Pro GKV_View_GIF, event; ; output image as GIF file; filename=Dialog_PickFile(/write, File='GKV_View.gif')      ; Get filename if filename eq "" then return                               	; Return if no filename supplied Widget_Control, event.top, Get_UValue=info, /No_Copy        	; Get user info Device, Get_Visual_Depth=thisDepth                          	; Check visual depth if thisDepth gt 8 then begin    image24=TvRD(True=1)    image2D=Color_Quan(image24,1,r,g,b) ENDIF ELSE BEGIN    image2d=TvRd()   image2d=TvRd()    r=info.r    g=info.g    b=info.b ENDELSE Write_Gif, filename, image2D, r, g, b Widget_Control, event.top, Set_Uvalue=info, /No_Copy        ; Return user info end   Pro GKV_View_JPEG, event; ; output image as JPEG file; filename=Dialog_PickFile(/write, File='GKV_View.jpg')      ; Get filename if filename eq "" then return                               	; Return if no filename supplied Widget_Control, event.top, Get_UValue=info, /No_Copy        	; Get user info Device, Get_Visual_Depth=thisDepth                          	; Check visual depth if thisDepth gt 8 then begin    image24=TvRD(True=1)                                     	; 24 (or 32) bit system endif ELSE Begin   Snap=TvRD()                                              ; bit system    s=Size(snap, /Dimensions)    image24=BytArr( 3, s[0], s[1])    image24[0,*,*]=info.r[snap]    image24[1,*,*]=info.g[snap]    image24[2,*,*]=info.b[snap] ENDELSE Write_JPEG, filename, image24, True=1, Quality=75           ; write figure to JPEG file Widget_Control, event.top, Set_Uvalue=info, /No_Copy        ; Return user info end Pro GKV_View_Print, event ; ; Send processed image to printer...but doesn't really work ; because printer doesn't understand !P.MULTI ; ok=Dialog_PrinterSetup() if ok eq 0 then return Widget_Control, event.top, Get_UValue=info, /No_Copy        ; Get user info Widget_Control, event.top, /HourGlass thisDevice=!D.nameWset, info.Wid key=PsWindow() Set_Plot, 'printer' Device, _Extra=key TvLct, info.r, info.g, info.b    ; Reload previous color table info.object -> draw, _Extra=*info.extraDevice, /Close_document Set_Plot, thisDevice Widget_Control, event.top, Set_Uvalue=info, /No_Copy        ; Return user info end    Pro GKV_View, object, NColors=ncolors, _Extra=extra, Group_Leader=group; ; Resizable graphics window for any object which ; has a "Draw" method ; if N_Elements(object) eq 0 then begin    ok=dialog_message('GKV_View:  must pass an Object as the parameter')    Return endif ;  Column sets number of columns ;  /TLB_Size_Events for resizable graphis window ;  MBar for menu bar tlb=Widget_Base(column=1, /TLB_SIZE_Events, MBar=menuID, Title='GKV_View') ; ; Create pull-down menus ;    fileID=Widget_Button(menuID, Value='File')    ; add items to file menu here...     ;  put a "save as" button on FILE menu          saveID=Widget_Button(fileID, Value='Save as ...', /Menu)          ; Add "children" to 'save as' button..          postscriptID=Widget_Button(saveID, Value='PostScript', Event_Pro='GKV_View_Postscript')          gifID=Widget_Button(saveID, Value='GIF', Event_Pro='GKV_View_GIF')          jpegID=Widget_Button(saveID, Value='JPEG', Event_Pro='GKV_View_JPEG')      ;  put a "print" button in the FILE menu          printID=Widget_Button(fileID, Value='Print...', Event_Pro='GKV_View_Print')    ;  put 'quit' button on FIlE menu          quitID=Widget_Button(FileID, Value='Quit', /Separator, Event_Pro='GKV_View_Quit')    ;    colorID=Widget_Button(menuID, Value='Color')       image_ID=Widget_Button(colorID, Value='Image Colors', Event_Pro='GKV_View_Colors') ; drawID=Widget_Draw(tlb, XSize=400, YSize=400) ; ; Draw object into draw widget ; Widget_Control, tlb, /Realize			; Put widget on screen Widget_Control, drawID, Get_Value=wID		; Get window ID of Draw Widget WSet, wID                             	; Make this current graphics window loadct, 5, NColors=ncolors             	; Load color table #1 object -> Draw,  _Extra=extra			; Draw object into Draw Widget if N_elements(ncolors) eq 0 then ncolors=!D.TABLE_SIZE ; ; Make info structure to save information needed to redisplay object ; TvLct, r,g,b, /Get      ; Get current color table info={object:object, drawID:drawID, wID:wID, $       r:r, g:g, b:b, ncolors:ncolors, Extra:PTR_NEW(extra)}  Widget_Control, tlb, Set_UValue=info, /No_Copy ; Save info structure  Xmanager,   'GKV_View', tlb, Event_Handler='GKV_View_Resize', $             /No_Block, Group_Leader=group, Cleanup="GKV_View_CleanUp"  end