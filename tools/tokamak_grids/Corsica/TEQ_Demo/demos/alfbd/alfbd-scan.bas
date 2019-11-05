# Vary fuzzy-marker weight (alfbd) over wide range and look at sensitivity
# on coil cost function and plasma shape deviation.

# Save base-case plasma boundary...
real rls0 = rls
real zls0 = zls

# Save base-case coil cost-function...
real cost0 = sum(rc*abs(cc))

function shape_dev

  # Reurn average deviation [cm] of plasma boundary versus base-case.
  real dsq = (rls - rls0)**2 + (zls - zls0)**2

  return sqrt(sum(dsq))/mls

endf  # shape_dev

package eq

integer i
real cost(0), dev(0), zalfbd(0)
do i = -3, 3
  alfbd = 16*10.0**i
  zalfbd := alfbd(1)
  run
  cost := sum(rc*abs(cc))
  dev := shape_dev
  pbg
enddo

oframe
titles ":F2:alfbd scan for ITER" ":F2:alfbd" ":F2:Coil cost & shape deviation"
attr labels yes
frame ,,0,2
plot cost/cost0 zalfbd scale loglin labels="cost"
plot dev zalfbd labels="dev"
plot [0, max(dev)],[16,16] style=dotted color=gray70 labels=" "
cframe


