# Demonstrate ceq with more complex constaints for an ITER EOB state.

# Here, we want to...
# (1) push the CS1 coils to their B-I limit-line: max(ufc(9:10)) == 1, and
# (2) push the CS2 coils to the vertical repulsion limit: csfz_rep == 120, while
# (3) holding betap, li and q(0) at nominal values.

package ceq

# Turn on coil diagnostics...
lop0 = 1

# Enable expression parsing in vo list, e.g. "max(ufc(9:10))"...
iequa = -2

# Force coil diagnostics be updated at each object function evaluation...
kbfc = 1

# Initialize vltf to the present value of flux linking plasma from coils and
# turn off flux constraint per Ejima scaling...
vltf = vltsnd(1); cejima = 0

# The CS2 currents will now be controled by HYBRD, not the G-S solver, so...
ic(7:12) = [7, 0, 8, 8, 0, 9]

nctot = 6
vo = ["betap(1)", "li(3)", "qsrf(1)", "max(ufc(9:10))", "csfz_rep", "cc(11)/cc(8)"]
vo0 = [0.65, 0.9, 0.9, 1, 120, 1]
vi = ["betaj", "alfa(0)", "betp(0)", "vltf", "cc(8)", "cc(11)"]
x0 = [betaj, alfa(0), betp(0), vltf, cc(8), cc(11)]
ihy = 20
run
