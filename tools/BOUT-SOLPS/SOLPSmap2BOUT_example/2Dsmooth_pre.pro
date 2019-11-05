.r smooth_bout2D.pro
filename='cmod-1100212023_260x64y_0.9psi_8PFR_v1_bout.grd.nc'
;filename='d3d_144981.03175_572_v2_260x64y_bout.grd.nc'
g=file_import(filename)
print,max(g.rxy[g.nx-1,*], ypeak), ypeak
iter=100

profne=g.neexp
profne=smooth_bout2D(profne, filename, iter=iter)

profni=g.niexp
profni=smooth_bout2D(profni, filename, iter=iter)

profTe=g.Teexp
profTe=smooth_bout2D(profTe, filename, iter=iter)

profTi=g.Tiexp
profTi=smooth_bout2D(profTi, filename, iter=iter)

;nimpbar=max(g.n_imp[*,*])
;profn_imp=10*g.n_imp/nimpbar
;profn_imp=smooth_bout2D(profn_imp, filename, iter=iter)
;profn_imp=profn_imp*nimpbar/10

;profpressure_s=g.pressure_s
;profpressure_s=smooth_bout2D(profpressure_s, filename, iter=iter)

;profE_r=g.E_r
;profE_r=smooth_bout2D(profE_r, filename, iter=iter)

;profbxcvz=g.bxcvz
;profbxcvz=smooth_bout2D(profbxcvz, filename, iter=iter)

;profsinty=g.sinty
;profsinty=smooth_bout2D(profsinty, filename, iter=iter)

;set_plot,'ps'
;device,/color,bits=8

;safe_colors
su=surface(profne[*,*])
su=surface(profni[*,*])
su=surface(profTe[*,*])
su=surface(profTi[*,*])
;su=surface(profn_imp[*,*])
;su=surface(profpressure_s[*,*])
;su=surface(profE_r[*,*])
;su=surface(profbxcvz[*,*])



;window, 2
;plot, profne[*,38]
;window, 3
;plot, profne[350,*]

;device,/close
;set_plot,'x'

handle = file_open(filename, /write)
  s = file_write(handle, 'Niexp', profni)
  s = file_write(handle, 'Neexp', profne)
  s = file_write(handle, 'Tiexp', profTi)
  s = file_write(handle, 'Teexp', profTe)
 ; s = file_write(handle, 'pressure_s', profpressure_s)
 ; s = file_write(handle, 'N_imp', profn_imp)
 ; s = file_write(handle, 'E_r', profE_r)
 ; s = file_write(handle, 'bxcvz', profbxcvz)
 ; s = file_write(handle, 'sinty', profsinty)
  file_close, handle



;handle = file_open(filename, /write)
;s = file_write(handle, 'Neexp', profne)
;file_close, handle
