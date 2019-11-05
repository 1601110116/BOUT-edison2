path="data"
print, 'data collected from path=./data'
;;g=file_import("data/uedge.grd_Up_Ni_Tei_2d.nc")
g=file_import("data/circle.grd.hl2a.nc")

print, 'Collecting Ni'
ni=collect(path=path,var="Ni")
print, 'Collecting Ti'
Ti=collect(path=path,var="Ti")
print, 'Collecting Vi'
vi=collect(path=path,var="Vi")
print, 'Collecting Te'
Te=collect(path=path,var="Te")
print, 'Collecting Nn'
nn=collect(path=path,var="Nn")
print, 'Collecting Tn'
tn=collect(path=path,var="Tn")
print, 'Collecting Vn'
Vn=collect(path=path,var="Vn")
print, 'Collecting Nm'
nm=collect(path=path,var="Nm")
print, 'Collecting Vmx'
vmx=collect(path=path,var="Vmx")

;;IF KEYWORD_SET(output) THEN BEGIN
  ;;  SET_PLOT, 'PS'
   ;; DEVICE, file="testplot.ps", /color, /landscape
;;ENDIF

safe_colors, /first

ntime=30

;!P.MULTI=[0,2,2,0,0]
;;Ni

plot, g.rxy[*,32],ni[*,12,0,ntime], color=1, xtitle="X", title="Ni"
oplot,g.rxy[*,32], ni[*,22,0,ntime], color=2
oplot,g.rxy[*,32], ni[*,32,0,ntime], color=4


;;Ti

;plot, Ti[*,12,0,ntime], color=1, xtitle="X", title="Ti"
;oplot, Ti[*,22,0,ntime], color=2
;oplot, Ti[*,32,0,ntime], color=4
;;Vi

;plot, vi[*,12,0,ntime], color=1, xtitle="X", title="vi"
;oplot, vi[*,22,0,ntime], color=2
;oplot, vi[*,32,0,ntime], color=4

;;Te

;plot, te[*,12,0,ntime], color=1, xtitle="X", title="Te"
;oplot, te[*,22,0,ntime], color=2
;oplot, te[*,32,0,ntime], color=4

;;IF KEYWORD_SET(output) THEN BEGIN
 ;;   DEVICE, /close
 ;;   SET_PLOT, 'X'
;;ENDIF

;!P.MULTI=0

