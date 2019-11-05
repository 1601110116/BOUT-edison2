PRO GKVS1D::WriteASCII, FileName=name, Path=path
;
; This proceedure writes contents of 'self' into an ascii file
;
; Keywords
;
;	FileName	Set to an text string containing
;			the desired 'root' for the
;			name of the file to be written. 
;			FileName will then be root.txt
;			Defaults to self.mnemonic
;
;	Path		Set to the desired path to the 
;			directory in which ASCII file
;			is to be written.
;			Defaults to the current working
;			directory.
;
; Written by W.M. Nevins
;  12/1/03
;
CD, CURRENT=CurrentWorkingDirectory
IF(TypeOF(Path) EQ 7) THEN CD, path

FileName=self.mnemonic
IF(TypeOF(name) EQ 7) THEN fileName=name
FileName=FileName + ".txt"
;
; Get I/O unit, and open file for writting
;
GET_LUN, IOunit
OPENW, IOunit, FileName
;
; Write headers to file
;
Independent = self.Grid1.mnemonic
dependent   = self.mnemonic
PRINTF, IOunit, FORMAT='(2A15)', independent, dependent
;
; Write values to file
;
irange      =self.Grid1.irange
independent = *self.Grid1.values
dependent   = *self.values
for i=irange[0], irange[1] DO PRINTF, IOunit, FORMAT='("  ", G13.6, "  ", G13.6)', 	$
			independent[i], dependent[i]
;
; Final duties
;
CLOSE, IOunit
FREE_LUN, IOunit
CD, CurrentWorkingDirectory
RETURN
END   ; ****** GKVS1D::WriteASCII ****** ;


PRO GKVS2D::WriteASCII, FileName=name, Path=path
;
; This proceedure writes contents of 'self' into an ascii file
;
; Keywords
;
;	FileName	Set to an text string containing
;			the desired 'root' for the
;			name of the file to be written. 
;			FileName will then be root.txt
;			Defaults to self.mnemonic
;
;	Path		Set to the desired path to the 
;			directory in which ASCII file
;			is to be written.
;			Defaults to the current working
;			directory.
;
; Written by W.M. Nevins
;  12/1/03
;
; This routine is just a 'trap' to prevent
; calling GKVS1D::WriteASCII on higher dimensional
; objects.
;
MESSAGE, "WriteASCII only implimented in 1-D", /INFORMATIONAL
RETURN
END   ; ****** GKVS2D::WriteASCII ****** ;


