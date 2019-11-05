# Start-up with the reference SOF equilibrium
chameleon basename = trim(probid)
integer i
if (plcm < 15) then
  i = index(basename, "]")
  basename = substr(basename, 1, i)//" @ "//format(plcm, 0, 1, 1)//" MA"
endif

# ceq settings...
package ceq
nctot = 2
vo = ["li(3)", "betap(1)"]
vo0 = [ li(3), 0.1]
vi = ["alfa(0)", "betaj"]
x0 = [alfa(0), betaj]
betp(0) = 1
ihy = 99; factor = 0.01; lop0 = 1

# Storage for results...
real zli = 0.01*fromone(iota(70, 120, 5))
integer n = length(zli)
real zq0(n), zufc6(n), zufc9(n)

# Survey loop...
logical addq0 = false, dropq0 = false
win
do i = 1, n
  vo0(1) = zli(i)
  probid = basename//" li="//format(vo0(1), 0, 2, 1)
  probid
  run
  if (qsrf(1) < 0.9 & ~addq0) then
    <<return<<"Adding q(0) constraint for "<<trim(probid)
    addq0 = true
    nctot = 3
    vo(nctot) = "qsrf(1)"
    vo0(nctot) = 0.9
    vi(nctot) = "betp(0)"
    x0(nctot) = betp(0)
    run
  elseif (betp(0) > 4 & ~dropq0) then
    <<return<<"Dropping q(0) constraint for "<<trim(probid)
    dropq0 = true
    nctot = 2
    run
  endif
  layout; profiles; pufc
  chameleon sname = "sof-li="//format(vo0(1), 0, 2, 1)//".sav"
#   saveq(sname)
  zq0(i) = qsrf(1)
  zufc6(i) = ufc(6)
  zufc9(i) = ufc(9)
enddo

oframe
attr labels=yes
chameleon title = ":F2:li-survey @ SOF"
if (plcm < 15) then
  title = ":F2:li-survey @ "//format(plcm, 0, 1, 1)//" MA"
endif
titles title ":F2:li(3)" ":F2:q(0) and UFC"
plot zq0, zli labels="q(0)"
plot zufc6, zli labels="PF6"
plot zufc9, zli labels="CS1"
cframe

<<return<<"Summary..."<<return
<<"   li(3)    q(0)  ufc(PF6)  ufc(CS1)"
do i = 1, n
  <<format(zli(i), 8, 2, 1)<<format(zq0(i), 8, 2, 1) \
   <<format(zufc6(i), 10, 2, 1)<<format(zufc9(i), 10, 2, 1)
enddo

