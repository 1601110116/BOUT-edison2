# Mal-formed ceq problem: trying to use pressure parameters to control the
# current profile instead of J0 parameters.

package ceq
nctot = 3
vo = ["betap(1)", "li(3)", "qsrf(1)"]
vo0 = [0.5, 1, 0.9]
vi = ["betaj", "alfa(1)", "betp(1)"]
x0 = [betaj, alfa(1), betp(1)]
ihy = 20
run
