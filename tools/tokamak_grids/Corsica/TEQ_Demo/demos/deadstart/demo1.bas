echo = yes

# Create an equilibrium (actually, two) using the tokamak deadstart procedure
# instead of using a save-file.
#
# Start-up the session with...
#
#     caltrans -probname fdf2 demo1.bas
#
# (1) Load script defining deadstart procedure...

      read tokamak.bas

# (2) Execute the desdstart procedure...

      ds("tokamak.inp", "coils.inp")

      # NOTE: save-file gets named with name in tokamak.inp file as in
      #       tolower(<name>)//".sav"

# (3) Impose up/down symmetry...

      set_symmetric

      # NOTE: save-file gets named using caltrans probname as in
      #       <probname>//"-sym.sav"

# (4) Rename the save-files...

      call basisexe("mv fdf2.sav fdf2-asym.sav")
      call basisexe("mv "//trim(probname)//"-sym.sav fdf2.sav")
