!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           14 Aug 25 15:39:00 *
!**************************************************************
! User     : Full professional version
! Notebook : consistentproduct
! Evaluation time                 : 6 s     Mode  : Optimal
! Number of formulae              : 6       Method: Automatic
! Subroutine                      : consistentproduct size: 148
! Total size of Mathematica  code : 148 subexpressions
! Total size of Fortran code      : 535 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE consistentproduct(v,a,b,prod,dprodda,dproddb)
USE SMSUtility
IMPLICIT NONE
DOUBLE PRECISION v(36),a(6),b(6),prod,dprodda(6),dproddb(6)
v(31)=2d0*b(6)
v(30)=2d0*b(5)
v(29)=2d0*b(4)
prod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*v(29)+a(5)*v(30)+a(6)*v(31)
dprodda(1)=b(1)
dprodda(2)=b(2)
dprodda(3)=b(3)
dprodda(4)=v(29)
dprodda(5)=v(30)
dprodda(6)=v(31)
dproddb(1)=a(1)
dproddb(2)=a(2)
dproddb(3)=a(3)
dproddb(4)=2d0*a(4)
dproddb(5)=2d0*a(5)
dproddb(6)=2d0*a(6)
END
