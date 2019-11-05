!+
!  Interface routine to allow IDL to call procedure is a sharealble library
!  via CALL_EXTERNAL. Intended for linux OS.
!  Created: 02/3/2010 by rjg    
!-

      ! Interface routine to SLATEC DFC routine (from http://www.netlib.org/slatec/src/)
      SUBROUTINE IDL_DFC(ARGC,ARGV)
      IMPLICIT NONE
      INTEGER*4 ARGC
      INTEGER*4 ARGV(*)
      


      ! Here is the normal fortran call 
      ! All TYPE REAL variables are DOUBLE PRECISION

!      call SUBROUTINE DFC (           
!                           NDATA,              ! in: # data points             
!                           XDATA,              ! in: array of independent variables
!                           YDATA,              ! in: Array of dependent variables
!                           SDDATA,             ! in: Standard deviations of ydata
!                           NORD,               ! in: order of spline (The order of the spline
!                                               !     is one more than the  degree of the piecewise
!                                               !     polynomial defined on each interval.  This is 
!                                               !     consistent with the  B-spline package
!                                               !     convention.  For example, NORD=4 when we are 
!                                               !     using piecewise cubics.)
!                           NBKPT,              ! in: number of knots in b-spline
!                           BKPT,               ! in: array of knots
!                           NCONST,             ! in: number of constraints
!                           XCONST,             ! in: x-location of constraint
!                           YCONST,             ! in: value of constraint 
!                           NDERIV,             ! in: type of constraint and derivative value
!                           MODE,               ! in/out: on input, mode for type of lsq solution
!                                                         on output, status of the fit
!                           COEFF,              ! out: array of coefficients from fit if mode = 0 or 1;
!                                                      size is nbkt-nord
!                           W,                  ! in/out: storage for floats
!                           IW)                 ! in/out: storage for intergers; in addition,
!                                               !         IW(1),IW(2) must contain the amounts of
!                                                         working storage actually  allocated for
!                                                         the working arrays W(*) and IW(*).
      print *,' Hello from IDL_DFC '     

      ! Here is the CALL_EXTERNAL call 
      call dfc (
                 %val(argv(1)),       &  ! ndata
                 %val(argv(2)),       &  ! xdata
                 %val(argv(3)),       &  ! ydata
                 %val(argv(4)),       &  ! sddata
                 %val(argv(5)),       &  ! nord
                 %val(argv(6)),       &  ! nbkpt
                 %val(argv(7)),       &  ! nconst
                 %val(argv(8)),       &  ! xconst
                 %val(argv(9)),       &  ! yconst
                 %val(argv(10)),      &  ! nderiv
                 %val(argv(11)),      &  ! mode
                 %val(argv(12)),      &  ! coeff
                 %val(argv(13)),      &  ! w
                 %val(argv(14))   )      ! iw



      print *, ' IDL_DFC: return from call to dfc'

      end
!
!
