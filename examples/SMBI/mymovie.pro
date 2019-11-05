PRO mymovie, data, uedge, mpeg=mpeg, delay=delay

  on_error,2  ; If an error occurs, return to caller

  
  mg=uedge
  val = reform(data[*,*,0,*])

  s = SIZE(val)
 
  wcur = !D.WINDOW ; Get current window ID
  xsize = !D.X_size ; Size of the current window
  ysize = !D.Y_size 

  IF NOT KEYWORD_SET(DELAY) THEN delay=1.0

  IF KEYWORD_SET(mpeg) THEN BEGIN
	 mpegfile = mpeg+".mpeg"
	 mpegid = MPEG_OPEN([xsize,ysize],file=mpegfile)
         PRINT, "Writing to: " +mpegfile
  ENDIF

     
      nx = s[1]
      nz = s[2]
      nt = s[3]
      g=uedge
     x=(g.psixy[*,32]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
   ;   loadct, 39
    ;  device, decomposed=0
                                ;safe_colors, /first
   ;   TVLCT, r, g, b, /get        
  
      FOR i=0, nt-1 DO BEGIN
          ;  contour, val[*,*,j],mg.rxy,mg.zxy, /fill, nlevels=50
           ;plot,x,val[*,32,i],thick=5,charsize=2,charthick=2,psym=-4,title="Time "+strtrim(i,2)+" Ta"   
            contour, val[*,*,i],mg.rxy,mg.zxy, /fill, nlevels=50

               
      ;;        IF delay LT 0.0 THEN BEGIN
       ;;           cursor, x, y, /down
        ;;      ENDIF ELSE WAIT, delay
           WAIT, delay
          IF KEYWORD_SET(mpeg) THEN BEGIN
              PRINT, "  frame "+strtrim(string(i+1),1)+"/"+strtrim(string(nt),1)
	      frame = tvrd(true=1,/order)
              MPEG_PUT, mpegid, image=frame, frame=i
              WAIT, delay
          ENDIF

       ENDFOR

      IF KEYWORD_SET(mpeg) THEN BEGIN
	PRINT, "Saving MPEG file: "+mpegfile
	MPEG_SAVE,mpegid
	MPEG_CLOSE,mpegid
      ENDIF
END
    
