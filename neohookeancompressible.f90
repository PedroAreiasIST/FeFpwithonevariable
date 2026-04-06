!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           13 Oct 25 16:00:11 *
!**************************************************************
! User     : Full professional version
! Notebook : neohookeancompressible
! Evaluation time                 : 1 s     Mode  : Optimal
! Number of formulae              : 59      Method: Automatic
! Subroutine                      : neohookeancompressible size: 999
! Total size of Mathematica  code : 999 subexpressions
! Total size of Fortran code      : 2761 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE neohookeancompressible(v,young,poiss,ce,psi,se,dsedce)
USE SMSUtility
IMPLICIT NONE
DOUBLE PRECISION v(167),young,poiss,ce(6),psi,se(6),dsedce(6,6)
v(154)=ce(5)**2
v(153)=ce(4)**2
v(152)=2d0*poiss
v(151)=young/(2d0+v(152))
v(150)=ce(2)*ce(3)-ce(6)**2
v(149)=(-2d0)*ce(6)
v(148)=2d0*ce(5)
v(147)=2d0*ce(4)
v(155)=ce(5)*v(147)
v(146)=ce(6)*v(147)-ce(2)*v(148)
v(145)=ce(1)*ce(3)-v(154)
v(144)=ce(1)*ce(2)-v(153)
v(74)=-(ce(3)*v(147))+ce(6)*v(148)
v(76)=ce(1)*v(149)+v(155)
v(63)=v(151)/2d0
v(54)=((-2d0)*v(151)+young/(1d0-v(152)))/3d0
v(162)=v(54)/2d0
v(56)=ce(1)*v(150)-ce(3)*v(153)-ce(2)*v(154)+ce(6)*v(155)
v(157)=2d0/v(56)
v(84)=sqrt(v(56))
v(156)=v(162)/v(84)**2
v(78)=1d0/v(56)**2
v(83)=-(v(76)*v(78))
v(82)=-(v(146)*v(78))
v(81)=-(v(74)*v(78))
v(103)=v(156)*v(76)
v(101)=v(146)*v(156)
v(99)=v(156)*v(74)
v(58)=dlog(v(84))
v(60)=-v(151)+v(54)*v(58)
v(129)=v(157)*v(60)
v(158)=(v(156)*v(157)-v(129)*v(56)*v(78))/4d0
v(109)=(v(103)/v(56)+v(60)*v(83))/2d0
v(161)=2d0*v(109)
v(108)=(v(101)/v(56)+v(60)*v(82))/2d0
v(160)=2d0*v(108)
v(107)=(v(60)*v(81)+v(99)/v(56))/2d0
v(159)=2d0*v(107)
v(105)=v(145)*v(158)
v(104)=v(150)*v(158)
v(64)=v(60)/(2d0*v(56))
v(137)=v(109)*v(150)+v(149)*v(64)
v(133)=v(108)*v(145)-v(148)*v(64)
v(127)=v(107)*v(144)-v(147)*v(64)
v(111)=2d0*(v(104)*v(145)+ce(3)*v(64))
v(112)=2d0*(v(104)*v(144)+ce(2)*v(64))
v(114)=2d0*(v(105)*v(144)+ce(1)*v(64))
v(116)=v(150)*v(159)
v(117)=v(145)*v(159)
v(119)=v(150)*v(160)
v(121)=v(144)*v(160)
v(123)=v(145)*v(161)
v(124)=v(144)*v(161)
v(130)=(ce(6)*v(129)+v(160)*v(74))/2d0
v(131)=(ce(5)*v(129)+v(161)*v(74))/2d0
v(136)=(ce(4)*v(129)+v(146)*v(161))/2d0
psi=v(151)*(((-3d0)+ce(1)+ce(2)+ce(3))/2d0-v(58))+v(162)*(v(58)*v(58))
se(1)=2d0*(v(63)+v(150)*v(64))
se(2)=2d0*(v(63)+v(145)*v(64))
se(3)=2d0*(v(63)+v(144)*v(64))
se(4)=v(64)*v(74)
se(5)=v(146)*v(64)
se(6)=v(64)*v(76)
dsedce(1,1)=2d0*v(104)*v(150)
dsedce(1,2)=v(111)
dsedce(1,3)=v(112)
dsedce(1,4)=v(116)
dsedce(1,5)=v(119)
dsedce(1,6)=2d0*v(137)
dsedce(2,1)=v(111)
dsedce(2,2)=2d0*v(105)*v(145)
dsedce(2,3)=v(114)
dsedce(2,4)=v(117)
dsedce(2,5)=2d0*v(133)
dsedce(2,6)=v(123)
dsedce(3,1)=v(112)
dsedce(3,2)=v(114)
dsedce(3,3)=2d0*(v(144)*v(144))*v(158)
dsedce(3,4)=2d0*v(127)
dsedce(3,5)=v(121)
dsedce(3,6)=v(124)
dsedce(4,1)=v(116)/2d0
dsedce(4,2)=v(117)/2d0
dsedce(4,3)=v(127)
dsedce(4,4)=(-(ce(3)*v(129))+v(159)*v(74))/2d0
dsedce(4,5)=v(130)
dsedce(4,6)=v(131)
dsedce(5,1)=v(119)/2d0
dsedce(5,2)=v(133)
dsedce(5,3)=v(121)/2d0
dsedce(5,4)=v(130)
dsedce(5,5)=(-(ce(2)*v(129))+v(146)*v(160))/2d0
dsedce(5,6)=v(136)
dsedce(6,1)=v(137)
dsedce(6,2)=v(123)/2d0
dsedce(6,3)=v(124)/2d0
dsedce(6,4)=v(131)
dsedce(6,5)=v(136)
dsedce(6,6)=(-(ce(1)*v(129))+v(161)*v(76))/2d0
END
