mdef Wait()=
  pause($1)
  call basisexe("/usr/bin/clear")
mend

restore iter_tfc.pfb
chameleon probid0 = probid
chameleon note0 = note

<<"# Demonstrate Corsica's graphics.bas routine"

<<"# Open an EZN Xgks window..."
<<"corsica> win"
win

Wait("scripts")
<<"# List all script files which have been loaded..."
<<"corsica> scripts"
scripts

Wait("graphics")
<<"# Display the graphics.bas intro..."
<<"corsica> graphics"
graphics

Wait("colors")
<<"# colors command..."
<<"corsica> colors"
colors

Wait("fonts")
<<"# fonts command..."
<<"corsica> fonts"
fonts

Wait("list-fonts")
<<"# list fonts..."
<<"corsica> fonts(""list"")"
fonts("list")

Wait("greek-fonts")
<<"# Display greek (F8) fonts..."
<<"corsica> fonts(""f8"")"
fonts("f8")

Wait("use-greek")
<<"# Use greek fonts..."
<<"corsica> probid = "":F2:This looks :F8:greek:F2: to me"""
<<"corsica> note = ""note:: there's a colon in this string"""
<<"corsica> layout"
probid = ":F2:This looks :F8:greek:F2: to me"
note = "note:: there's a colon in this string"""
layout
probid = probid0
note = note0

Wait("layout-0")
<<"# Layout style 0"
<<"corsica> layout(0,0)"
note = "layout style 0"
layout(0,0)

Wait("layout-1")
<<"# Layout style 1: show coils in proportion to their current"
<<"corsica> layout(1,0)"
note = "layout style 1"
layout(1,0)

Wait("layout-2")
<<"# Layout style 2: show filaments"
<<"corsica> layout(2,0)"
note = "layout style 2"
layout(2,0)

Wait("layout-3")
<<"# Layout style 3: show coil indices"
<<"corsica> layout(3,0)"
note = "layout style 3"
layout(3,0)

Wait("layout-4")
<<"# Layout style 4: show filaments and R-Z grid"
<<"corsica> layout(4,0)"
note = "layout style 4"
layout(4,0)

Wait("zoom-1")
<<"# Change plot scale to match R-Z grid"
<<"corsica> zoom; layout(0,0)"
note = "Layout on R-Z grid domain"
zoom; layout(0,0)

Wait("zoom-2")
<<"# Layout near divertor"
<<"corsica> zoom(4, 6.5, -4.5, -2.5); layout(0,0)"
note = "Near divertor"
zoom(4, 6.5, -4.5, -2.5); layout(0,0)

zoom("reset")
note = " "

Wait("pbx")
<<"# Plot all x-point surfaces"
<<"corsica> pbx"
pbx

Wait("pfbd")
<<"# Plot fuzzy-boundary points with index values"
<<"corsica> pfbd"
pfbd

Wait("pls")
<<"# Plot plasma boundary with index values"
<<"corsica> pls"
pls

Wait("pplates")
<<"# Plot first-wall plates with plate numbers"
<<"corsica> zoom; pplates"
zoom; pplates
zoom("reset")

Wait("profiles")
<<"# Plot plasma profiles"
<<"corsica> profiles"
profiles

Wait("pwires")
<<"# Plot passive structure ""wires"""
<<"corsica> pwires"
pwires

Wait("quit")
quit
