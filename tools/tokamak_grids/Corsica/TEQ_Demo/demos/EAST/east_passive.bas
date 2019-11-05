nwires = 80
integer i, io = basopen("east_vv.data", "r")
do i = 1,nwires
  io>>rwires(i)>>zwires(i)>>drwires(i)>>dzwires(i)>>awires(i)>>awires2(i)
enddo
close(io)
awires=awires*pi/180
awires2=awires2*pi/180
rhwires = 0.74e-06 # per Humphrys/make_east_objects.m (etav entry)
idwires(1:40) = "VVi"
idwires(41:) = "VVo"
