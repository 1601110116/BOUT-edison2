pro fitfunc_ga, x, c, f, pder
  z=2.*(c[0]-x)/c[1]
  pz1=1.+c[4]*z+c[5]*z*z+c[6]*z*z*z
  pz2=1.+c[7]*z
  f=0.5*(c[2]-c[3])*(pz1*exp(z)-pz2*exp(-z))/(exp(z)+exp(-z))+0.5*(c[2]+c[3])
  
  pz1_z = c[4]+2.*c[5]*z+3.*c[6]*z*z
  pz2_z = c[7]
  f_pz1 = 0.5*(c[2]-c[3])/(exp(z)+exp(-z))*exp(z)
  f_pz2 = -0.5*(c[2]-c[3])/(exp(z)+exp(-z))*exp(-z)
  pder_z=-0.5*(c[2]-c[3])*( pz1_z*exp(2.*z)-pz2_z*exp(-2.*z)+pz1_z-pz2_z+2.*pz1+2.*pz2 )/(exp(z)+exp(-z))^2.

  pder0 = pder_z*2./c[1]
  pder1 = pder_z*(-z/c[1])
  pder2 = 0.5*(pz1*exp(z)-pz2*exp(-z))/(exp(z)+exp(-z))+0.5
  pder3 = -0.5*(pz1*exp(z)-pz2*exp(-z))/(exp(z)+exp(-z))+0.5
  pder4 = f_pz1*z
  pder5 = f_pz1*z*z
  pder6 = f_pz1*z*z*z
  pder7 = f_pz2*z

  pder = [[pder0],[pder1],[pder2],[pder3],[pder4],[pder5],[pder6],[pder7]]
;  help,z,pder
end

function fit_prof, psn, prof
  aa=[0.98244772,0.034941133,0.51108721,0.029467707,0.28668591,-0.0069066899,9.0527052e-05,-0.015914946]
;  aa=[1,1,1,1,1,1,1,1]
;  weights = replicate(1.0,n_elements(psn))
;  weights = psn
  weights = replicate(0.0,n_elements(psn))
;  help,aa,weights
  result = mpcurvefit(psn,prof,weights,aa,sigma,function_name='fitfunc_ga')
;  result = curvefit(psn,prof,weights,aa,sigma,function_name='fitfunc_ga')

  result = result*(prof[0]/result[0])

  plot,psn,prof
  oplot,psn,result,col=2
  oplot,psn,deriv(deriv(deriv(result))),col=3

  return,result
end
