Function CrossCoherence, Arg1, Arg2, lagWindow=lw
;
; Computes cross coherence and cross phase between the two objects Arg1 and Arg2
;
arg1 -> scaleAxis, 't', /uniform
narg2 = arg2 -> interpolate(arg1)
Spect1 = arg1 -> xspect(lw=lw)
spect2 = narg2 -> xspect(lw=lw)
xspect = arg1 -> xspect(lw=lw, ref=narg2)
narg2 -> trash
;

Spect1 -> get, values=ptr1
Spect2 -> get, values=ptr2
xspect -> get, values=ptr3
s1 = *ptr1
s2 = *ptr2
xs = *ptr3
temp = xs/SQRT(s1*s2)
xCoh = ABS(temp)
xPhaseValues = IMAGINARY(ALOG(temp))

xCoherence = xspect -> MakeCopy(/noValues)
arg1 -> get, title=title1, mnemonic=mnemonic1
arg2 -> get, title=title2, mnemonic=mnemonic2
mnemonic='xCoh{' + mnemonic1 + ',' + mnemonic2 + '}'
title = 'xCoh{' + title1 + ',' + title2 + '}'
xCoherence -> set, values=PTR_NEW(xCoh), units='', mnemonic=mnemonic, title=title
xCoherence -> set, vrange=[0,1]

xPhase = xspect -> MakeCopy(/noValues)
mnemonic='xPhase{' + mnemonic1 + ',' + mnemonic2 + '}'
title = 'xPhase{' + title1 + ',' + title2 + '}'
xPhase -> set, values=PTR_NEW(xPhaseValues), units='', mnemonic=mnemonic, title=title
xPhase -> set, vrange=[-!PI, !PI]

spect1 -> trash
spect2 -> trash
xspect -> trash

result= { Name:"CrossCoherence", xCoherence:xCoherence, xPhase:xPhase }
RETURN, result
END ; ****** CrossCoherence ****** ;