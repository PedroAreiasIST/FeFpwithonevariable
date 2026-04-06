!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           6 Apr 26 20:47:21  *
!**************************************************************
! User     : Full professional version
! Notebook : hosfordcriterionsymmpd
! Evaluation time                 : 66 s    Mode  : Optimal
! Number of formulae              : 65      Method: Automatic
! Subroutine                      : hosfordcriterionsymmpd size: 902
! Total size of Mathematica  code : 902 subexpressions
! Total size of Fortran code      : 3045 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE hosfordcriterionsymmpd(v,n,mandel,phi,dphidmandel,norm,dnormdmandel)
USE SMSUtility
IMPLICIT NONE
DOUBLE PRECISION v(153),n,mandel(6),phi,dphidmandel(6),norm(6),dnormdmandel(6,6)
v(139)=6d0*mandel(6)**n
v(138)=6d0*mandel(5)**n
v(137)=6d0*mandel(4)**n
v(133)=n/2d0
v(132)=mandel(1)-mandel(3)
v(131)=(-1d0)+n
v(136)=v(132)**v(131)
v(130)=mandel(2)-mandel(3)
v(134)=v(130)**v(131)
v(129)=mandel(1)-mandel(2)
v(135)=v(129)**v(131)
v(128)=1d0/n
v(141)=1d0*v(128)
v(66)=(-1d0)+v(128)
v(68)=v(133)*(v(134)-v(135))
v(71)=v(133)*(-v(134)-v(136))
v(65)=v(133)*(v(135)+v(136))
v(73)=mandel(4)**v(131)
v(144)=0.15d1*v(73)
v(76)=mandel(5)**v(131)
v(148)=0.15d1*v(76)
v(79)=mandel(6)**v(131)
v(47)=(v(129)**n+v(130)**n+v(132)**n+v(137)+v(138)+v(139))/2d0
v(69)=v(47)**((-2d0)+v(128))
v(142)=v(66)*v(69)
v(140)=3d0*n*v(142)
v(81)=v(140)*v(79)
v(87)=v(141)*v(81)
v(78)=v(140)*v(76)
v(86)=v(141)*v(78)
v(75)=v(140)*v(73)
v(85)=v(141)*v(75)
v(72)=v(142)*v(71)
v(147)=0.15d1*v(72)
v(70)=v(142)*v(68)
v(146)=0.15d1*v(70)
v(83)=v(141)*v(70)
v(67)=v(142)*v(65)
v(145)=0.15d1*v(67)
v(82)=v(141)*v(67)
v(61)=v(47)**v(66)
v(143)=3d0*v(61)
v(54)=v(141)*v(61)
v(117)=v(131)*v(143)
v(90)=(-1d0)+v(131)
v(89)=v(131)*v(133)
v(92)=v(129)**v(90)*v(89)
v(91)=v(132)**v(90)*v(89)
v(88)=v(130)**v(90)*v(89)
v(53)=v(54)*v(65)
v(57)=v(54)*v(68)
v(60)=v(54)*v(71)
v(62)=v(143)*v(73)
v(63)=v(143)*v(76)
v(64)=v(143)*v(79)
v(94)=v(68)*v(82)-v(54)*v(92)
v(95)=v(71)*v(82)-v(54)*v(91)
v(97)=v(71)*v(83)-v(54)*v(88)
v(112)=v(144)*v(78)
v(113)=v(144)*v(81)
v(119)=v(148)*v(81)
phi=v(47)**v(128)
dphidmandel(1)=v(53)
dphidmandel(2)=v(57)
dphidmandel(3)=v(60)
dphidmandel(4)=v(62)
dphidmandel(5)=v(63)
dphidmandel(6)=v(64)
norm(1)=v(53)
norm(2)=v(57)
norm(3)=v(60)
norm(4)=v(62)/2d0
norm(5)=v(63)/2d0
norm(6)=v(64)/2d0
dnormdmandel(1,1)=v(65)*v(82)+v(54)*(v(91)+v(92))
dnormdmandel(1,2)=v(94)
dnormdmandel(1,3)=v(95)
dnormdmandel(1,4)=v(65)*v(85)
dnormdmandel(1,5)=v(65)*v(86)
dnormdmandel(1,6)=v(65)*v(87)
dnormdmandel(2,1)=v(94)
dnormdmandel(2,2)=v(68)*v(83)+v(54)*(v(88)+v(92))
dnormdmandel(2,3)=v(97)
dnormdmandel(2,4)=v(68)*v(85)
dnormdmandel(2,5)=v(68)*v(86)
dnormdmandel(2,6)=v(68)*v(87)
dnormdmandel(3,1)=v(95)
dnormdmandel(3,2)=v(97)
dnormdmandel(3,3)=v(141)*v(71)*v(72)+v(54)*(v(88)+v(91))
dnormdmandel(3,4)=v(71)*v(85)
dnormdmandel(3,5)=v(71)*v(86)
dnormdmandel(3,6)=v(71)*v(87)
dnormdmandel(4,1)=v(145)*v(73)
dnormdmandel(4,2)=v(146)*v(73)
dnormdmandel(4,3)=v(147)*v(73)
dnormdmandel(4,4)=(mandel(4)**v(90)*v(117)+3d0*v(73)*v(75))/2d0
dnormdmandel(4,5)=v(112)
dnormdmandel(4,6)=v(113)
dnormdmandel(5,1)=v(145)*v(76)
dnormdmandel(5,2)=v(146)*v(76)
dnormdmandel(5,3)=v(148)*v(72)
dnormdmandel(5,4)=v(112)
dnormdmandel(5,5)=(mandel(5)**v(90)*v(117)+3d0*v(76)*v(78))/2d0
dnormdmandel(5,6)=v(119)
dnormdmandel(6,1)=v(145)*v(79)
dnormdmandel(6,2)=v(146)*v(79)
dnormdmandel(6,3)=v(147)*v(79)
dnormdmandel(6,4)=v(113)
dnormdmandel(6,5)=v(119)
dnormdmandel(6,6)=(mandel(6)**v(90)*v(117)+3d0*v(79)*v(81))/2d0
END
