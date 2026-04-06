!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           30 Mar 26 22:17:55 *
!**************************************************************
! User     : Full professional version
! Notebook : yld_ghorbel
! Evaluation time                 : 34 s    Mode  : Optimal
! Number of formulae              : 72      Method: Automatic
! Subroutine                      : yld_ghorbel size: 1051
! Total size of Mathematica  code : 1051 subexpressions
! Total size of Fortran code      : 3457 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE yld_ghorbel(v,rsigma,st0,mandel,seq,dseqdmandel,norm,dnormdmandel)
USE SMSUtility
IMPLICIT NONE
DOUBLE PRECISION v(181),rsigma,st0,mandel(6),seq,dseqdmandel(6),norm(6),dnormdmandel(6,6)
v(167)=2d0*mandel(6)**2
v(166)=2d0*mandel(5)**2
v(165)=2d0*mandel(4)**2
v(164)=0.5625d-2*st0**2
v(163)=mandel(1)+mandel(2)+mandel(3)
v(162)=1d0/rsigma
v(161)=(-1d0)+rsigma
v(173)=0.5d0*v(161)
v(97)=v(162)*v(173)
v(73)=2d0*v(163)
v(57)=(v(163)*v(163))
v(56)=v(164)+v(57)
v(74)=1d0/sqrt(v(56))
v(171)=v(73)*v(74)
v(75)=(-0.5d0)*v(171)/v(56)
v(86)=v(163)*v(75)
v(87)=v(162)*(0.15d1+v(161)*(0.75d0+0.75d0*v(163)*v(74)))
v(55)=1d0/sqrt(v(56))
v(100)=v(55)+v(86)
v(96)=v(100)*v(97)
v(49)=(-1d0/3d0)*v(163)
v(65)=mandel(3)+v(49)
v(63)=mandel(2)+v(49)
v(60)=mandel(1)+v(49)
v(50)=(v(165)+v(166)+v(167)+(v(60)*v(60))+(v(63)*v(63))+(v(65)*v(65)))/2d0
v(170)=sqrt(3d0*v(50))
v(76)=1d0/v(170)
v(169)=0.15d1*v(76)
v(168)=3d0*v(76)
v(85)=mandel(6)*v(168)
v(105)=v(85)*v(96)
v(83)=mandel(5)*v(168)
v(104)=v(83)*v(96)
v(81)=mandel(4)*v(168)
v(103)=v(81)*v(96)
v(79)=v(169)*v(65)
v(78)=v(169)*v(63)
v(77)=v(169)*v(60)
v(99)=v(170)*((v(171)*v(57))/v(56)**2+(3d0-v(57)/v(56))*v(75))*v(97)
v(102)=v(79)*v(96)+v(99)
v(101)=v(78)*v(96)+v(99)
v(98)=v(77)*v(96)+v(99)
v(90)=1d0/v(170)**2
v(172)=-(v(87)*v(90))
v(95)=v(172)*v(85)
v(94)=v(172)*v(83)
v(93)=v(172)*v(81)
v(89)=v(169)*(v(74)+v(86))*v(97)
v(92)=v(172)*v(79)+v(89)
v(91)=v(172)*v(78)+v(89)
v(88)=v(172)*v(77)+v(89)
v(62)=v(76)*v(87)
v(130)=(2d0/3d0)*v(62)
v(124)=(-1d0/3d0)*v(62)
v(176)=v(124)+v(98)
v(175)=v(102)+v(124)
v(174)=v(101)+v(124)
v(115)=2d0*v(62)
v(58)=v(170)*v(96)
v(61)=v(58)+v(60)*v(62)
v(64)=v(58)+v(62)*v(63)
v(66)=v(58)+v(62)*v(65)
v(67)=mandel(4)*v(115)
v(68)=mandel(5)*v(115)
v(69)=mandel(6)*v(115)
v(146)=mandel(4)*v(94)
v(147)=mandel(4)*v(95)
v(152)=mandel(5)*v(95)
seq=v(162)*v(170)*(1d0+v(173)*(1d0+v(163)*v(55)))
dseqdmandel(1)=v(61)
dseqdmandel(2)=v(64)
dseqdmandel(3)=v(66)
dseqdmandel(4)=v(67)
dseqdmandel(5)=v(68)
dseqdmandel(6)=v(69)
norm(1)=v(61)
norm(2)=v(64)
norm(3)=v(66)
norm(4)=v(67)/2d0
norm(5)=v(68)/2d0
norm(6)=v(69)/2d0
dnormdmandel(1,1)=v(130)+v(60)*v(88)+v(98)
dnormdmandel(1,2)=v(174)+v(60)*v(91)
dnormdmandel(1,3)=v(175)+v(60)*v(92)
dnormdmandel(1,4)=v(103)+v(60)*v(93)
dnormdmandel(1,5)=v(104)+v(60)*v(94)
dnormdmandel(1,6)=v(105)+v(60)*v(95)
dnormdmandel(2,1)=v(176)+v(63)*v(88)
dnormdmandel(2,2)=v(101)+v(130)+v(63)*v(91)
dnormdmandel(2,3)=v(175)+v(63)*v(92)
dnormdmandel(2,4)=v(103)+v(63)*v(93)
dnormdmandel(2,5)=v(104)+v(63)*v(94)
dnormdmandel(2,6)=v(105)+v(63)*v(95)
dnormdmandel(3,1)=v(176)+v(65)*v(88)
dnormdmandel(3,2)=v(174)+v(65)*v(91)
dnormdmandel(3,3)=v(102)+v(130)+v(65)*v(92)
dnormdmandel(3,4)=v(103)+v(65)*v(93)
dnormdmandel(3,5)=v(104)+v(65)*v(94)
dnormdmandel(3,6)=v(105)+v(65)*v(95)
dnormdmandel(4,1)=mandel(4)*v(88)
dnormdmandel(4,2)=mandel(4)*v(91)
dnormdmandel(4,3)=mandel(4)*v(92)
dnormdmandel(4,4)=(v(115)+2d0*mandel(4)*v(93))/2d0
dnormdmandel(4,5)=v(146)
dnormdmandel(4,6)=v(147)
dnormdmandel(5,1)=mandel(5)*v(88)
dnormdmandel(5,2)=mandel(5)*v(91)
dnormdmandel(5,3)=mandel(5)*v(92)
dnormdmandel(5,4)=v(146)
dnormdmandel(5,5)=(v(115)+2d0*mandel(5)*v(94))/2d0
dnormdmandel(5,6)=v(152)
dnormdmandel(6,1)=mandel(6)*v(88)
dnormdmandel(6,2)=mandel(6)*v(91)
dnormdmandel(6,3)=mandel(6)*v(92)
dnormdmandel(6,4)=v(147)
dnormdmandel(6,5)=v(152)
dnormdmandel(6,6)=(v(115)+2d0*mandel(6)*v(95))/2d0
END
