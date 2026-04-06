!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
subroutine mainelast(ce,psi,se,dsedce)
  use basfun
  implicit real(8)(a-h,o-z)
  real(8)::young,poiss,num,mum,km
  real(8),dimension(6)::ce,se,dsdvv
  real(8),dimension(3)::a0
  real(8),dimension(6,3)::dseda0
  real(8),dimension(6,6)::dsedce
!**************
  young=2.34044d09
  poiss=0.38d00  
  call NEOHOOKEANINCOMPRESSIBLE(YOUNG,POISS,CE,PSI,SE,DSEDCE)  
!!  call neohookeancompressible(young,poiss,ce,psi,se,dsedce)  
  return
  A0(1)=1.0D00
  A0(2)=0.0d00
  A0(3)=0.0d00  
  CALL VECTNORMALIZE(3,A0)
  CALL GUO2008DRIVINGFORCE(A0,CE,FORCE)
!*** shear modulus of matrix:
  MUM=0.93d9!1.1508d9
!*** poisson coefficient of matrix
  NUM=0.38d00!4864D00
!*** bulk modulus of matrix
  KM=MUM*2.0D00*(1.0D00+NUM)/(3.0D00-6.0D00*NUM)
  young=mum*(2.0d00*(1.0d00+num))
!*** initial void fraction
  g=16.971d6
  smax=47.0d6
  e0=0.01916d00
  ebar=0.00175d00  
  VV0=0.06781D00
  VVO=VV0
  A0(1)=1.0D00
  A0(2)=0.0d00
  A0(3)=0.0d00
  CALL VECTNORMALIZE(3,A0)
  CALL GUO2008DRIVINGFORCE(A0,CE,FORCE)
  CALL GUO2008VOIDFROMFORCE(VV0,VVO,VVN,E0,EBAR,FORCE)
  vvn=min(0.99d00,vvn)  
  VVMAX=VV0
  IF(FORCE.LE.1.02D00)THEN
     CALL GUO2008MODIFIED(VVMAX,MUM,KM,A0,CE,PSI,SE,DSEDCE,DSEDA0,DSDVV)
  ELSEIF(FORCE.LE.1.07d00)THEN
     CALL GUO2008MODIFIED(VVMAX,0.4D00*MUM,0.4D00*KM,A0,CE,PSI,SE,DSEDCE,DSEDA0,DSDVV)
  ELSE
     CALL GUO2008MODIFIED(VVMAX,0.2D00*MUM,0.2D00*KM,A0,CE,PSI,SE,DSEDCE,DSEDA0,DSDVV)
  END IF
end subroutine mainelast

!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
subroutine mainplast(mandel,seq,dseqdmandel,norm,dnormdmandel)
  implicit real(8)(a-h,o-z)
  real(8)::n=8.0d00
  real(8),dimension(6)::mandel,dseqdmandel,norm
  real(8),dimension(6,6)::dnormdmandel
  rsigma=1.387d00
  st0=62.0d06
!  call druckerprager(rsigma,mandel,seq,dseqdmandel,norm,dnormdmandel)
!  call unsymmetric(rsigma,mandel,seq,dseqdmandel,norm,dnormdmandel)
  call yld_ghorbel(rsigma,st0,mandel,seq,dseqdmandel,norm,dnormdmandel)
!  call hosfordcriterionsymm(8.0d00,mandel,seq,dseqdmandel,norm,dnormdmandel)  
end subroutine mainplast
!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
subroutine mainhard(kn,kn1,y,dydk)
  implicit real(8)(a-h,o-z)
  real(8)::n,kn,kn1,y,dydk
  st0=62.0d06
  ssat=78.0d06
  b=14.0d00
  n=1.30d00
  y=ssat-(ssat-st0)*exp(-b*kn1**n)
  dydk=1.0d00*n*(ssat-st0)*(b*kn1**(n-1.0d00))*exp(-b*kn1**n)
  return
  if (kn <= 0.03d00) then
! Linear hardening region
     Y = 37.0d0 + 20.0d0 * kn1
     dYdk = 20.0d0
  else 
     Y = 37.0d00+20.0d00*0.03d00
     dYdk =0.0d00!2156.25d0 * kn1**2 - 1318.75d0 * kn1 + 120.3125d0
  end if
  y=1.0d6*y
  dydk=1.0d6*dydk
end subroutine mainhard
real(8) function vectnorm2(n, v)
  implicit real(8)(a-h,o-z)
  integer, intent(in) :: n
  real(8), intent(in) :: v(*)
  real(8) :: nrm
  integer :: i
  nrm = 0.0d0
  do i=1,n
     nrm = nrm + v(i)*v(i)
  end do
  nrm = sqrt(nrm)
  vectnorm2=nrm
end function vectnorm2

module intmod
  type plas1
     real(8),dimension(6)::norm,normfixed
     real(8)::gamma,seq,gammamax
     real(8),dimension(6,6)::depdcetr

!*** c
     real(8),dimension(6)::c,dfdep
!*** end c

!*** ce
     real(8),dimension(6)::ce,cetr,dcedgamma,ep,dgammadcetr
     real(8),dimension(6,6)::dcedcetr,dcednorm,dcedep,dcetrdc
!*** end ce

     real(8),dimension(6,6)::dnormdmandel,dnormdce,dnormfixedcetr

!*** se
     real(8),dimension(6)::se
     real(8),dimension(6,6)::dsedce,dsedcetr
!*** end se

!*** residual
!*** end residual

!*** y
     real(8)::y
     real(8)::dydk
     real(8),dimension(6)::dydep,dydcetr
!*** end y

!*** f
     real(8)::f
     real(8)::dfdgamma
     real(8),dimension(6)::dfdnorm,dfdcetr
!*** end f

!*** mandel
     real(8),dimension(6)::mandel
     real(8),dimension(6)::dmandelgamma
     real(8),dimension(6,6)::dmandeldnorm,dmandeldcetr,dmandeldce,dmandeldep
!*** end mandel
  end type plas1
end module intmod

!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
SUBROUTINE substepping1(NSUBI,SUBELAST,SUBPLAST,SUBHARD,GN,GN1,KN,KN1,STRAIN,PSI,S,DSTRESSDSTRAIN)
  use intmod
  implicit real(8)(a-h,o-z)
  real(8),parameter::small=1.0d-12
  external::subelast,subplast,subhard
  real(8)::kn,kn1,psi,kntemp
  real(8),dimension(3,3)::gn,gn1,gntemp
  real(8),dimension(6)::strain,s,straintemp
  real(8),dimension(6,6)::dsdc,DSTRESSDSTRAIN
  type(plas1)::p
!-----------------------------------
!*** STORE OLD HISTORICAL VARIABLES
!-----------------------------------  
  NSUB=ABS(NSUBI)
!-----------------
!*** TIME SCALING
!-----------------
  KNTEMP=KN
  GNTEMP=GN
  IF(NSUBI.LT.0)THEN
     DO ISUB=1,NSUB
        COMB=1.0D00*ISUB/(1.0D00*NSUB)
        STRAINTEMP=STRAIN
        CALL integrate1(SUBELAST,SUBPLAST,SUBHARD,GN,GN1,KN,KN1,STRAINTEMP,PSI,S,DSDC,P)
        KN=KN1
        GN=GN1
     END DO
  ELSE
     DO ISUB=1,NSUB
        COMB=1.0D00*ISUB/(1.0D00*NSUB)
        STRAINTEMP=STRAIN*COMB
        CALL integrate1(SUBELAST,SUBPLAST,SUBHARD,GN,GN1,KN,KN1,STRAINTEMP,PSI,S,DSDC,P)
        KN=KN1
        GN=GN1        
        DO JSUB=1,10
           CALL integrate1(SUBELAST,SUBPLAST,SUBHARD,GN,GN1,KN,KN1,STRAINTEMP,PSI,S,DSDC,P)
           KN=KN1
           GN=GN1
        END DO
     END DO
  END IF
!-------------------
!*** RESTORE VALUES
!-------------------
  KN=KNTEMP
  GN=GNTEMP
!-------------------  
!*** DSTRESSDSTRAIN
!-------------------
  DO J=1,6
     DO I=1,6
        IF(J.LE.3)THEN
           DSTRESSDSTRAIN(I,J)=2.0D00*DSDC(I,J)
        ELSE
           DSTRESSDSTRAIN(I,J)=DSDC(I,J)
        END IF
     END DO
  END DO
!  write(*,*)"Dstressdstrain(4,4)",dstressdstrain(4,4)
END SUBROUTINE substepping1

!--------------------
!*** integrate in 1d
!--------------------
!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
subroutine integrate1(subelast,subplast,subhard,gn,gn1,kn,kn1,strain,psi,s,dsdc,p)
  use intmod
  use basfun
  implicit real(8)(a-h,o-z)
  type(plas1)::p
  real(8),parameter::small=1.0d-10
  external::subelast,subplast,subhard
  real(8)::kn,kn1,psi
  real(8),dimension(3,3)::gn,gn1
  real(8),dimension(6)::strain,s
  real(8),dimension(6,6)::dsdc,dsdctrial
!------------------------
!*** total strain
!*** beware of tolerances
!------------------------
  if(vectnorm2(6,strain).le.small)then
     strain(3)=-1.0d-9
     strain(6)=1.0d-9
  end if
!----------------------------------------------
!*** c right cauchy-green strain
!*** attention: strain is in engineering format
!----------------------------------------------
  do iv=1,3
     p%c(iv)=2.0d00*strain(iv)+1.0d00
  end do
  do iv=4,6
     p%c(iv)=strain(iv)
  end do
!------------
!*** check gn
!------------
  if(vectnorm2(9,gn).le.small)then
     do id=1,3
        gn(id,id)=1.0d00
     end do
  end if
!----------------------------------
!*** ensure that history is correct
!*** for early exits
!----------------------------------
  kn1=kn
  gn1=gn
  p%gamma=0.0d00
!------------
!*** ce trial
!------------
  call cefromc2(p%c,gn,p%cetr,p%dcetrdc)!attention dcetrdc is needed in the end!
  p%ep=0.0d00
  call ffromgamma(subelast,subplast,subhard,kn,kn1,psi,p)
  call dettotal(gn,gn1,p%se,p%dsedce,p%dcetrdc,s,dsdc)
  dsdctrial=dsdc
  p%normfixed=p%norm
  p%dnormfixedcetr=p%dnormdce
!*** ep now should be null
  p%ep=p%gamma*p%normfixed
  if(p%f.gt.0.0d00)then
     if(p%f.le.0.0d00)then
        write(*,*)"f violation=",p%f,"gamma=",p%gamma
     end if
     emax=2.0d00
     ns=2000
     p%gammamax=emax/vectnorm2(6,p%norm)
     dgamma=p%gammamax/(1.0d00*ns)
     do i=1,ns
        p%gamma=dgamma*i
        p%ep=p%gamma*p%normfixed
        kn1=kn
        call ffromgamma(subelast,subplast,subhard,kn,kn1,psi,p)
        if(p%f.le.0.0d00)then
           p%gammamax=p%gamma
           exit
        end if
     enddo
!*** ensure that f is negative or zero
     if(p%f.gt.0.0d00)then
        write(*,*)"f>0 equals",p%f
        p%gamma=1.0d10
        rn=vectnorm2(6,p%normfixed)
        p%normfixed=p%mandel*rn/vectnorm2(6,p%mandel)
        p%ep=p%gamma*p%normfixed
        call ffromgamma(subelast,subplast,subhard,kn,kn1,psi,p)
        if(p%f.gt.0)stop "ERROR !!!!"
     else
        p%gamma=p%gammamax
     endif
     call iterate(subelast,subplast,subhard,kn,kn1,psi,p)
     p%dgammadcetr=-p%dfdcetr/p%dfdgamma
     do iv=1,6
        do jv=1,6
           p%depdcetr(iv,jv)=p%normfixed(iv)*p%dgammadcetr(jv)+p%gamma*p%dnormfixedcetr(iv,jv)
        end do
     end do
     do iv=1,6
        do jv=1,6
           p%dsedcetr(iv,jv)=0.0d00
           do kv=1,6
              p%dsedcetr(iv,jv)=p%dsedcetr(iv,jv)+p%dsedce(iv,kv)*p%dcedcetr(kv,jv)
           end do
        end do
     end do
     CALL tangentmodulus(GN,GN1,p%EP,p%SE,p%DEPDCETR,p%DCETRDC,p%DSEDCETR,S,DSDC)
  else
     CALL dettotal(GN,GN1,p%SE,p%DSEDCE,p%DCETRDC,S,DSDC)
  end if
  dsdc=0.4d00*dsdc+0.6d00*dsdctrial
end subroutine integrate1

!----------------------------------------------------
!*** determines kn1
!*** and y and derivatives of y with respect to:
!*** gamma, norm and cetr (dydep, dydcetr)
!----------------------------------------------------
subroutine hardeningsubproblem(subhard,kn,kn1,p)
  use intmod
  implicit real(8)(a-h,o-z)
  type(plas1)::p
  real(8)::kn,kn1
  external::subhard
  real(8)::res
  real(8),dimension(6)::dresdep
  real(8),dimension(6)::dresdcetr   
  real(8),dimension(6)::dproddmandel,dproddep
  kn1=kn
  call consistentproduct(p%mandel,p%ep,prod,dproddmandel,dproddep)
  if(prod.gt.0.0d00)then
     nitr=10000
     do itr=1,nitr
        call subhard(kn,kn1,p%y,p%dydk)
        res=(kn1-kn)*p%y-prod
        dres=p%y+(kn1-kn)*p%dydk
        if(abs(dres).gt.1.0d-20)then
           kn1=kn1-0.25d00*res/dres
        else
           exit
        end if
        kn1=max(kn1,1.0d-20)
        if(prod.lt.0.0d00)stop "prod negative"
        if(abs(res).le.1.0d-5*abs(prod)+1.0d-5)exit
     end do
     if((itr.gt.nitr).or.(kn1.lt.kn))then
        kn1=kn
        call subhard(kn,kn1,p%y,p%dydk)
     end if
     do iv=1,6
        dresdep(iv)=-dproddep(iv)
        do jv=1,6
           dresdep(iv)=dresdep(iv)-dproddmandel(jv)*p%dmandeldep(jv,iv)
        end do
     end do
     do iv=1,6
        dresdcetr(iv)=0.0d00
        do jv=1,6
           dresdcetr(iv)=dresdcetr(iv)-dproddmandel(jv)*p%dmandeldcetr(jv,iv)         
        end do
     end do
     if(abs(dres).gt.1.0d-20)then
        p%dydep=-(dydk/dres)*dresdep
        p%dydcetr=-(dydk/dres)*dresdcetr
     else
        p%dydep=0.0d00
        p%dydcetr=0.0d00
     end if
  else
     call subhard(kn,kn1,p%y,p%dydk)
     p%dydep=0.0d00
     p%dydcetr=0.0d00
  end if
end subroutine hardeningsubproblem

!--------------------------------------------------
!*** given gamma, calculate mandel and derivatives
!*** with respect to gamma, norm and cetr
!--------------------------------------------------
!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
subroutine mandelfromep(subelast,psi,p)
  use intmod
  implicit real(8)(a-h,o-z)
  type(plas1)::p
  real(8)::length,dtime,psi
  external::subelast
  call doublesideed(p%cetr,p%ce,p%ep,p%dcedep,p%dcedcetr)
  call subelast(p%ce,psi,p%se,p%dsedce)
  call detmandelsymmetric(p%ce,p%se,p%dsedce,p%mandel,p%dmandeldce)
  call mandelderivatives(p%dmandeldce,p%dcedep,p%dcedcetr,p%dmandeldep,p%dmandeldcetr)
end subroutine mandelfromep

!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
subroutine fanddf(subelast,subplast,subhard,kn,kn1,psi,p)
  use intmod
  implicit real(8)(a-h,o-z)
  type(plas1)::p
  real(8)::kn,kn1
  external::subelast,subplast,subhard  
  call ffromgamma(subelast,subplast,subhard,kn,kn1,psi,p)     
  p%dfdgamma=dot_product(p%dfdep,p%normfixed)
end subroutine fanddf

subroutine iterate(subelast,subplast,subhard,kn,kn1,psi,p)
  use intmod
  implicit real(8)(a-h,o-z)
  type(plas1)::p
  real(8)::psi,kn,kn1
  external::subelast,subplast,subhard
  p%gammamax=p%gamma
  p%ep=p%gamma*p%normfixed
  if(p%f.gt.0.0d00)then
     write(*,*)"wait a minute...."
     stop
  end if
  call rootg(0.0d00,p%gammamax,funcinternal,gammas,ifl)
  p%gamma=gammas
contains
  subroutine funcinternal(x,fx,dfx)
    implicit real(8)(a-h,o-z)
    p%gamma=x
    p%ep=p%gamma*p%normfixed
    call fanddf(subelast,subplast,subhard,kn,kn1,psi,p)
    fx=p%f
    dfx=p%dfdgamma
  end subroutine funcinternal
end subroutine iterate

!---------------------------------------------
!*** determines f from gamma and derivatives:
!*** with respect to gamma
!*** with respect to norm
!*** with respect to cetr
!*** AND updates kn1. Beware of that!
!---------------------------------------------
!DIR$ OPTIMIZE:0
!DIR$ NOINLINE
!seems ok
subroutine ffromgamma(subelast,subplast,subhard,kn,kn1,psi,p)
  use intmod
  use basfun
  implicit real(8)(a-h,o-z)
  type(plas1)::p
  external::subelast,subplast,subhard
  real(8)::kn,kn1,psi
  real(8),dimension(6)::dseqdmandel
!-------------------------------
!*** mandel stress
!*** equivalent stress
!*** hardening function and kn1
!-------------------------------
  call mandelfromep(subelast,psi,p)
  call subplast(p%mandel,p%seq,dseqdmandel,p%norm,p%dnormdmandel)
!*** dnormdce
  do i=1,6
     do j=1,6
        p%dnormdce(i,j)=0.0d00
        do k=1,6
           p%dnormdce(i,j)=p%dnormdce(i,j)+p%dnormdmandel(i,k)*p%dmandeldce(k,j)
        end do
     end do
  end do
  call hardeningsubproblem(subhard,kn,kn1,p)
  do id=1,6
     p%dfdep(id)=0.0d00
     do jd=1,6
        p%dfdep(id)=p%dfdep(id)+dseqdmandel(jd)*p%dmandeldep(jd,id)
     end do
  end do
  do id=1,6
     p%dfdcetr(id)=0.0d00
     do jd=1,6
        p%dfdcetr(id)=p%dfdcetr(id)+dseqdmandel(jd)*p%dmandeldcetr(jd,id)
     end do
  end do
!----------------------- 
!*** hardening function
!-----------------------
  p%f=p%seq-p%y
  p%dfdep=p%dfdep-p%dydep
  p%dfdcetr=p%dfdcetr-p%dydcetr
end subroutine ffromgamma

!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           24 Oct 25 11:35:16 *
!**************************************************************
! User     : Full professional version
! Notebook : doublesideed
! Evaluation time                 : 261 s   Mode  : Optimal
! Number of formulae              : 337     Method: Automatic
! Subroutine                      : doublesideed size: 6862
! Total size of Mathematica  code : 6862 subexpressions
! Total size of Fortran code      : 16692 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE doublesideed(cetr,ce,deltaep,dcedep,dcedcetr)
  IMPLICIT NONE
  DOUBLE PRECISION v(525),cetr(6),ce(6),deltaep(6),dcedep(6,6),dcedcetr(6,6)
  v(501)=deltaep(3)**2
  v(500)=v(501)/16d0
  v(495)=(-12d0)+deltaep(3)
  v(492)=(-0.25d0)*deltaep(6)
  v(490)=(deltaep(4)*deltaep(5)+(deltaep(2)+deltaep(3))*deltaep(6))/16d0
  v(489)=deltaep(2)**2
  v(488)=v(489)/16d0
  v(483)=deltaep(1)**2
  v(481)=v(483)/16d0
  v(480)=(-0.25d0)*deltaep(4)
  v(479)=(-0.25d0)*deltaep(3)
  v(478)=(deltaep(1)*deltaep(4)+deltaep(2)*deltaep(4)+deltaep(5)*deltaep(6))/16d0
  v(477)=(-0.25d0)*deltaep(5)
  v(475)=(-0.25d0)*deltaep(2)
  v(474)=(deltaep(1)*deltaep(5)+deltaep(3)*deltaep(5)+deltaep(4)*deltaep(6))/16d0
  v(473)=deltaep(6)**2/16d0
  v(472)=deltaep(6)/8d0
  v(470)=deltaep(6)/16d0
  v(471)=(-0.25d0)*v(470)
  v(468)=deltaep(5)**2/16d0
  v(467)=(-0.25d0)*deltaep(1)
  v(466)=deltaep(5)/8d0
  v(464)=deltaep(5)/16d0
  v(465)=(-0.25d0)*v(464)
  v(463)=deltaep(4)**2/16d0
  v(462)=deltaep(4)/8d0
  v(460)=deltaep(4)/16d0
  v(469)=(-1d0/24d0)*v(460)
  v(461)=(-0.25d0)*v(460)
  v(459)=deltaep(3)/16d0
  v(458)=deltaep(2)/16d0
  v(457)=deltaep(1)/16d0
  v(154)=v(457)+v(458)
  v(207)=v(458)+v(459)
  v(147)=v(457)+v(459)
  v(169)=deltaep(2)*v(461)
  v(167)=v(460)/2d0
  v(155)=deltaep(4)*v(461)
  v(160)=v(155)/6d0
  v(165)=v(462)*v(467)
  v(214)=deltaep(5)*v(469)
  v(152)=6d0*v(214)
  v(164)=v(152)/3d0
  v(187)=deltaep(3)*v(465)
  v(183)=v(464)/2d0
  v(149)=deltaep(5)*v(465)
  v(161)=v(149)/6d0
  v(184)=v(466)*v(467)
  v(186)=deltaep(6)*v(469)
  v(171)=(-1d0/24d0)*(deltaep(6)*v(464))
  v(219)=deltaep(3)*v(471)
  v(215)=v(470)/2d0
  v(189)=deltaep(6)*v(471)
  v(210)=v(189)/6d0
  v(174)=(-(deltaep(6)*v(147))+24d0*v(215)-deltaep(4)*v(466)-deltaep(2)*v(470))/24d0
  v(157)=deltaep(4)*v(471)
  v(212)=v(157)/3d0
  v(185)=v(183)+(v(157)+v(184)+v(187))/6d0
  v(150)=deltaep(5)*v(471)
  v(166)=v(167)+(v(150)+v(165)+v(169))/6d0
  v(216)=v(472)*v(475)
  v(491)=v(152)+v(216)
  v(217)=v(215)+(v(219)+v(491))/6d0
  v(175)=(-0.25d0)*v(474)
  v(188)=v(183)+(v(175)+v(187))/6d0
  v(176)=v(183)+(v(157)+v(175)+v(464)*v(475))/6d0
  v(151)=v(175)+v(147)*v(477)
  v(476)=v(151)+v(157)
  v(250)=v(464)+(v(476)+v(466)*v(479))/6d0
  v(163)=v(464)+(v(184)+v(476))/6d0
  v(107)=v(474)*v(477)
  v(168)=(-0.25d0)*v(478)
  v(191)=v(167)+(v(150)+v(168)+v(460)*v(479))/6d0
  v(170)=v(167)+(v(168)+v(169))/6d0
  v(156)=v(168)+v(154)*v(480)
  v(514)=v(150)+v(156)
  v(100)=v(478)*v(480)
  v(493)=6d0+v(100)
  v(91)=v(463)+v(468)+v(481)
  v(484)=6d0+v(91)
  v(482)=(-6d0)-v(473)-v(91)
  v(172)=(-0.25d0)*v(91)
  v(190)=(12d0*v(147)-deltaep(3)*v(147)-deltaep(5)*v(466)+v(482))/24d0
  v(173)=(12d0*v(154)-deltaep(2)*v(154)-deltaep(4)*v(462)+v(482))/24d0
  v(159)=(-0.25d0)+v(457)+(v(149)+v(155)+v(172)-v(483)/32d0)/6d0
  v(92)=(-(deltaep(1)*v(484))+4d0*(v(107)+v(493)+3d0*v(91)))/24d0
  v(487)=2d0*v(92)
  v(95)=(-(deltaep(6)*v(474))+12d0*v(478)-deltaep(2)*v(478)-deltaep(4)*v(484))/24d0
  v(485)=2d0*v(95)
  v(274)=v(191)*v(95)
  v(235)=v(160)*v(95)
  v(182)=v(176)*v(485)
  v(181)=v(174)*v(485)
  v(180)=v(173)*v(485)
  v(179)=v(171)*v(485)
  v(178)=v(170)*v(485)
  v(177)=v(166)*v(485)
  v(113)=(v(95)*v(95))
  v(96)=(-(deltaep(6)*v(478))-deltaep(5)*v(484)-v(474)*v(495))/24d0
  v(486)=2d0*v(96)
  v(271)=v(161)*v(96)
  v(261)=v(176)*v(96)
  v(197)=v(191)*v(486)
  v(203)=v(182)+v(197)+v(164)*v(487)
  v(196)=v(190)*v(486)
  v(202)=v(181)+v(196)+v(163)*v(487)
  v(195)=v(174)*v(486)
  v(194)=v(188)*v(486)
  v(200)=v(179)+v(194)+v(161)*v(487)
  v(199)=v(178)+v(186)*v(486)+v(160)*v(487)
  v(192)=v(185)*v(486)
  v(198)=v(177)+v(192)+v(159)*v(487)
  v(115)=(v(96)*v(96))
  v(110)=v(113)+v(115)+(v(92)*v(92))
  v(512)=2d0*v(110)
  v(98)=v(463)+v(473)+v(488)
  v(494)=6d0+v(98)
  v(221)=(-0.25d0)*v(98)
  v(222)=((-6d0)+12d0*v(207)-deltaep(3)*v(207)-v(468)-deltaep(6)*v(472)-v(98))/24d0
  v(209)=(-0.25d0)+v(458)+(v(155)+v(189)+v(221)-v(489)/32d0)/6d0
  v(218)=(-0.25d0)*v(490)
  v(220)=v(215)+(v(218)+v(219))/6d0
  v(208)=v(218)+v(207)*v(492)
  v(213)=v(470)+(v(208)+v(491))/6d0
  v(106)=v(490)*v(492)
  v(111)=(-(deltaep(2)*v(494))+4d0*(v(106)+v(493)+3d0*v(98)))/24d0
  v(498)=2d0*v(111)
  v(496)=v(111)+v(92)
  v(102)=(-(deltaep(5)*v(478))-deltaep(6)*v(494)-v(490)*v(495))/24d0
  v(497)=2d0*v(102)
  v(273)=v(102)*v(174)
  v(258)=v(102)*v(210)
  v(240)=v(102)*v(191)+v(176)*v(496)+(v(164)+v(213))*v(95)+v(222)*v(96)
  v(239)=v(102)*v(190)+v(174)*v(496)+(v(163)+v(212))*v(95)+v(191)*v(96)
  v(238)=v(177)+v(178)+v(261)+v(273)+v(173)*v(496)
  v(237)=v(102)*v(188)+v(171)*v(496)+(v(161)+v(210))*v(95)+v(220)*v(96)
  v(236)=v(102)*v(186)+v(235)+v(170)*v(496)+v(209)*v(95)+v(217)*v(96)
  v(234)=v(102)*v(185)+v(235)+v(166)*v(496)+v(159)*v(95)+v(214)*v(96)
  v(228)=v(222)*v(497)
  v(233)=v(182)+v(228)+v(213)*v(498)
  v(227)=v(191)*v(497)
  v(232)=v(181)+v(227)+v(212)*v(498)
  v(226)=v(176)*v(497)
  v(225)=v(220)*v(497)
  v(230)=v(179)+v(225)+v(210)*v(498)
  v(224)=v(217)*v(497)
  v(229)=v(178)+v(224)+v(209)*v(498)
  v(116)=(v(102)*v(102))
  v(120)=(v(111)*v(111))+v(113)+v(116)
  v(509)=v(110)+v(120)
  v(506)=2d0*v(120)
  v(101)=v(496)*v(95)+v(102)*v(96)
  v(499)=2d0*v(101)
  v(289)=v(101)*v(199)
  v(246)=v(240)*v(499)
  v(245)=v(239)*v(499)
  v(244)=v(238)*v(499)
  v(243)=v(237)*v(499)
  v(242)=v(236)*v(499)
  v(241)=v(234)*v(499)
  v(119)=(v(101)*v(101))
  v(105)=v(468)+v(473)+v(500)
  v(248)=(-0.25d0)+v(459)+((-0.25d0)*v(105)+v(149)+v(189)-v(501)/32d0)/6d0
  v(112)=(-(deltaep(3)*(6d0+v(105)))+4d0*(6d0+3d0*v(105)+v(106)+v(107)))/24d0
  v(504)=2d0*v(112)
  v(503)=v(111)+v(112)
  v(502)=v(112)+v(92)
  v(275)=v(192)+v(194)+v(273)+v(274)+v(190)*v(502)
  v(272)=v(102)*v(171)+v(271)+v(188)*v(502)+v(220)*v(95)+v(248)*v(96)
  v(270)=v(102)*v(170)+v(186)*v(502)+v(217)*v(95)+(v(160)+v(210))*v(96)
  v(269)=v(102)*v(166)+v(271)+v(185)*v(502)+v(214)*v(95)+v(159)*v(96)
  v(262)=v(224)+v(225)+v(261)+v(274)+v(222)*v(503)
  v(260)=v(102)*(v(212)+v(250))+v(191)*v(503)+v(190)*v(95)+v(174)*v(96)
  v(343)=v(101)*v(260)
  v(259)=v(102)*v(248)+v(258)+v(220)*v(503)+v(188)*v(95)+v(171)*v(96)
  v(257)=v(102)*v(209)+v(258)+v(217)*v(503)+v(186)*v(95)+v(170)*v(96)
  v(256)=v(102)*(v(160)+v(161))+v(214)*v(503)+v(185)*v(95)+v(166)*v(96)
  v(255)=v(197)+v(228)+(v(470)+(v(152)+v(208)+v(472)*v(479))/6d0)*v(504)
  v(252)=v(194)+v(225)+v(248)*v(504)
  v(122)=(v(112)*v(112))+v(115)+v(116)
  v(511)=2d0*v(122)
  v(508)=v(110)+v(122)
  v(507)=v(120)+v(122)
  v(117)=v(102)*v(503)+v(95)*v(96)
  v(505)=2d0*v(117)
  v(340)=v(117)*v(230)
  v(305)=v(117)*v(239)
  v(268)=v(262)*v(505)
  v(337)=v(246)+v(268)+v(233)*v(506)
  v(267)=v(260)*v(505)
  v(336)=v(245)+v(267)+v(232)*v(506)
  v(266)=v(240)*v(505)
  v(335)=v(244)+v(266)+v(506)*(v(180)+v(226)+v(498)*(v(460)+(v(462)*v(475)+v(514))/6d0))
  v(265)=v(259)*v(505)
  v(334)=v(243)+v(265)+v(230)*v(506)
  v(388)=cetr(6)*v(334)
  v(264)=v(257)*v(505)
  v(333)=v(242)+v(264)+v(229)*v(506)
  v(125)=(v(117)*v(117))
  v(108)=v(102)*v(95)+v(502)*v(96)
  v(513)=2d0*v(108)
  v(510)=v(108)*v(260)
  v(342)=v(108)*v(240)
  v(344)=v(264)+v(265)+v(342)+v(343)+v(262)*v(507)
  v(350)=cetr(6)*v(344)
  v(341)=v(108)*v(237)+v(117)*v(252)+v(101)*v(272)+v(340)+v(259)*v(507)
  v(347)=cetr(6)*v(341)
  v(339)=v(117)*v(229)+v(108)*v(236)+v(101)*v(270)+v(340)+v(257)*v(507)
  v(346)=cetr(6)*v(339)
  v(338)=v(117)*(v(199)+v(200))+v(108)*v(234)+v(101)*v(269)+v(256)*v(507)
  v(345)=cetr(6)*v(338)
  v(307)=v(117)*v(240)+v(108)*(v(203)+v(255))+v(101)*v(262)+v(260)*v(508)
  v(419)=cetr(4)*v(307)
  v(355)=cetr(3)*v(307)
  v(349)=cetr(6)*v(307)
  v(313)=cetr(5)*v(307)
  v(303)=v(108)*v(200)
  v(304)=v(117)*v(237)+v(108)*v(252)+v(101)*v(259)+v(303)+v(272)*v(508)
  v(310)=cetr(5)*v(304)
  v(302)=v(108)*(v(199)+v(230))+v(117)*v(236)+v(101)*v(257)+v(270)*v(508)
  v(309)=cetr(5)*v(302)
  v(301)=v(108)*v(198)+v(117)*v(234)+v(101)*v(256)+v(303)+v(269)*v(508)
  v(308)=cetr(5)*v(301)
  v(294)=v(101)*(v(203)+v(233))+v(117)*v(260)+v(108)*v(262)+v(240)*v(509)
  v(402)=cetr(2)*v(294)
  v(371)=cetr(1)*v(294)+cetr(4)*v(337)+cetr(5)*v(344)
  v(368)=cetr(5)*v(294)
  v(357)=cetr(6)*v(337)+cetr(3)*v(344)+v(368)
  v(348)=cetr(6)*v(294)
  v(300)=cetr(4)*v(294)
  v(363)=v(300)+cetr(2)*v(337)+v(350)
  v(293)=v(101)*(v(202)+v(232))+v(117)*v(275)+v(239)*v(509)+v(510)
  v(401)=cetr(1)*v(293)
  v(370)=v(313)+cetr(4)*v(336)+v(401)
  v(323)=cetr(6)*v(293)
  v(311)=cetr(5)*v(293)
  v(356)=v(311)+cetr(6)*v(336)+v(355)
  v(299)=cetr(4)*v(293)
  v(362)=v(299)+cetr(2)*v(336)+v(349)
  v(292)=v(241)+v(242)+v(305)+v(342)+v(238)*v(509)
  v(369)=cetr(1)*v(292)+cetr(4)*v(335)+v(368)
  v(354)=cetr(5)*v(292)+cetr(3)*v(294)+cetr(6)*v(335)
  v(298)=cetr(4)*v(292)
  v(361)=v(298)+cetr(2)*v(335)+v(348)
  v(291)=v(101)*(v(200)+v(230))+v(108)*v(259)+v(117)*v(272)+v(237)*v(509)
  v(367)=cetr(1)*v(291)+cetr(4)*v(334)+cetr(5)*v(341)
  v(353)=cetr(5)*v(291)+cetr(3)*v(341)+v(388)
  v(297)=cetr(4)*v(291)
  v(360)=v(297)+cetr(2)*v(334)+v(347)
  v(290)=v(101)*v(229)+v(108)*v(257)+v(117)*v(270)+v(289)+v(236)*v(509)
  v(366)=cetr(1)*v(290)+cetr(4)*v(333)+cetr(5)*v(339)
  v(352)=cetr(5)*v(290)+cetr(6)*v(333)+cetr(3)*v(339)
  v(296)=cetr(4)*v(290)
  v(359)=v(296)+cetr(2)*v(333)+v(346)
  v(288)=v(101)*v(198)+v(108)*v(256)+v(117)*v(269)+v(289)+v(234)*v(509)
  v(295)=cetr(4)*v(288)
  v(281)=2d0*v(510)
  v(375)=v(268)+v(281)+v(255)*v(511)
  v(287)=v(246)+v(281)+v(203)*v(512)
  v(332)=cetr(1)*v(287)+v(300)+v(313)
  v(326)=cetr(4)*v(287)+v(349)+v(402)
  v(319)=cetr(5)*v(287)+v(348)+v(355)
  v(280)=v(275)*v(513)
  v(374)=v(267)+v(280)+(v(196)+v(227)+v(250)*v(504))*v(511)
  v(286)=v(245)+v(280)+v(202)*v(512)
  v(279)=v(239)*v(513)
  v(373)=v(266)+v(279)+((2d0/3d0)*v(112)*v(150)+v(195)+v(226))*v(511)
  v(285)=v(244)+v(279)+v(512)*(v(180)+v(195)+((v(165)+6d0*v(460)+v(514))*v(92))/3d0)
  v(330)=cetr(1)*v(285)+v(298)+v(311)
  v(324)=cetr(4)*v(285)+cetr(2)*v(292)+v(323)
  v(317)=cetr(5)*v(285)+cetr(6)*v(292)+cetr(3)*v(293)
  v(278)=v(272)*v(513)
  v(372)=v(265)+v(278)+v(252)*v(511)
  v(284)=v(243)+v(278)+v(200)*v(512)
  v(381)=cetr(5)*v(284)
  v(329)=cetr(1)*v(284)+v(297)+v(310)
  v(322)=cetr(4)*v(284)+cetr(2)*v(291)+cetr(6)*v(304)
  v(316)=cetr(6)*v(291)+cetr(3)*v(304)+v(381)
  v(283)=v(242)+v(199)*v(512)+v(270)*v(513)
  v(364)=cetr(4)*v(283)
  v(365)=cetr(1)*v(288)+cetr(5)*v(338)+v(364)
  v(358)=cetr(2)*v(283)+v(295)+v(345)
  v(351)=cetr(6)*v(283)+cetr(5)*v(288)+cetr(3)*v(338)
  v(328)=cetr(1)*v(283)+v(296)+v(309)
  v(321)=cetr(2)*v(290)+cetr(6)*v(302)+v(364)
  v(315)=cetr(5)*v(283)+cetr(6)*v(290)+cetr(3)*v(302)
  v(276)=v(269)*v(513)
  v(306)=v(276)+v(278)+v(305)+v(343)+v(275)*v(508)
  v(325)=cetr(4)*v(286)+cetr(2)*v(293)+cetr(6)*v(306)
  v(318)=cetr(5)*v(286)+cetr(3)*v(306)+v(323)
  v(312)=cetr(5)*v(306)
  v(331)=cetr(1)*v(286)+v(299)+v(312)
  v(282)=v(241)+v(276)+v(198)*v(512)
  v(327)=cetr(1)*v(282)+v(295)+v(308)
  v(320)=cetr(4)*v(282)+cetr(2)*v(288)+cetr(6)*v(301)
  v(314)=cetr(5)*v(282)+cetr(6)*v(288)+cetr(3)*v(301)
  v(124)=(v(108)*v(108))
  v(109)=(v(110)*v(110))+v(119)+v(124)
  v(114)=v(108)*v(117)+v(101)*v(509)
  v(515)=v(109)*v(114)
  v(133)=cetr(4)*v(114)
  v(118)=v(101)*v(117)+v(108)*v(508)
  v(517)=2d0*v(118)
  v(516)=v(109)*v(118)
  v(139)=cetr(5)*v(118)
  v(130)=cetr(5)*v(109)+cetr(6)*v(114)+cetr(3)*v(118)
  v(416)=v(130)*v(307)
  v(129)=cetr(4)*v(109)+cetr(2)*v(114)+cetr(6)*v(118)
  v(415)=v(129)*v(294)
  v(128)=cetr(1)*v(109)+v(133)+v(139)
  v(411)=v(128)*v(293)
  v(121)=v(119)+(v(120)*v(120))+v(125)
  v(518)=v(114)*v(121)
  v(123)=v(101)*v(108)+v(117)*v(507)
  v(520)=v(121)*v(123)
  v(519)=2d0*v(123)
  v(140)=cetr(6)*v(123)
  v(137)=cetr(5)*v(114)+cetr(6)*v(121)+cetr(3)*v(123)
  v(136)=cetr(2)*v(121)+v(133)+v(140)
  v(135)=cetr(1)*v(114)+cetr(4)*v(121)+cetr(5)*v(123)
  v(408)=v(135)*v(293)
  v(126)=(v(122)*v(122))+v(124)+v(125)
  v(380)=cetr(4)*v(118)+cetr(2)*v(123)+cetr(6)*v(126)
  v(379)=cetr(1)*v(118)+cetr(4)*v(123)+cetr(5)*v(126)
  v(378)=cetr(3)*v(126)+v(139)+v(140)
  v(425)=(v(114)*v(114))
  v(426)=(v(118)*v(118))
  v(429)=v(114)*v(118)
  v(431)=(v(123)*v(123))
  v(433)=v(114)*v(123)
  v(436)=v(118)*v(123)
  v(443)=v(109)*v(123)+v(429)
  v(444)=v(118)*v(121)+v(433)
  v(449)=v(114)*v(126)+v(436)
  ce(1)=v(109)*v(128)+v(114)*v(129)+v(118)*v(130)
  ce(2)=v(114)*v(135)+v(121)*v(136)+v(123)*v(137)
  ce(3)=v(126)*v(378)+v(118)*v(379)+v(123)*v(380)
  ce(4)=v(114)*v(128)+v(121)*v(129)+v(123)*v(130)
  ce(5)=v(118)*v(128)+v(123)*v(129)+v(126)*v(130)
  ce(6)=v(118)*v(135)+v(123)*v(136)+v(126)*v(137)
  dcedep(1,1)=v(128)*v(282)+v(129)*v(288)+v(130)*v(301)+v(118)*v(314)+v(114)*v(320)+v(109)*v(327)
  dcedep(1,2)=v(128)*v(283)+v(129)*v(290)+v(130)*v(302)+v(118)*v(315)+v(114)*v(321)+v(109)*v(328)
  dcedep(1,3)=v(128)*v(284)+v(129)*v(291)+v(130)*v(304)+v(118)*v(316)+v(114)*v(322)+v(109)*v(329)
  dcedep(1,4)=v(128)*v(285)+v(129)*v(292)+v(130)*v(293)+v(118)*v(317)+v(114)*v(324)+v(109)*v(330)
  dcedep(1,5)=v(128)*v(286)+v(129)*v(293)+v(130)*v(306)+v(118)*v(318)+v(114)*v(325)+v(109)*v(331)
  dcedep(1,6)=v(128)*v(287)+v(118)*v(319)+v(114)*v(326)+v(109)*v(332)+v(415)+v(416)
  dcedep(2,1)=v(136)*v(283)+v(135)*v(288)+v(137)*v(338)+v(123)*v(351)+v(121)*v(358)+v(114)*v(365)
  dcedep(2,2)=v(135)*v(290)+v(136)*v(333)+v(137)*v(339)+v(123)*v(352)+v(121)*v(359)+v(114)*v(366)
  dcedep(2,3)=v(135)*v(291)+v(136)*v(334)+v(137)*v(341)+v(123)*v(353)+v(121)*v(360)+v(114)*v(367)
  dcedep(2,4)=v(135)*v(292)+v(137)*v(294)+v(136)*v(335)+v(123)*v(354)+v(121)*v(361)+v(114)*v(369)
  dcedep(2,5)=v(137)*v(307)+v(136)*v(336)+v(123)*v(356)+v(121)*v(362)+v(114)*v(370)+v(408)
  dcedep(2,6)=v(135)*v(294)+v(136)*v(337)+v(137)*v(344)+v(123)*v(357)+v(121)*v(363)+v(114)*v(371)
  dcedep(3,1)=v(123)*(cetr(6)*v(284)+cetr(4)*v(301)+cetr(2)*v(338))+v(126)*(cetr(3)*v(284)+v(308)+v(345))+v(284)*v(378)+v&
       &(301)*v(379)+v(338)*v(380)+v(118)*(cetr(1)*v(301)+cetr(4)*v(338)+v(381))
  dcedep(3,2)=v(118)*(cetr(1)*v(302)+cetr(5)*v(334)+cetr(4)*v(339))+v(126)*(v(309)+cetr(3)*v(334)+v(346))+v(334)*v(378)+v&
       &(302)*v(379)+v(339)*v(380)+v(123)*(cetr(4)*v(302)+cetr(2)*v(339)+v(388))
  dcedep(3,3)=v(126)*(v(310)+v(347)+cetr(3)*v(372))+v(118)*(cetr(1)*v(304)+cetr(4)*v(341)+cetr(5)*v(372))+v(123)*(cetr(4&
       &)*v(304)+cetr(2)*v(341)+cetr(6)*v(372))+v(372)*v(378)+v(304)*v(379)+v(341)*v(380)
  dcedep(3,4)=v(126)*(v(311)+v(348)+cetr(3)*v(373))+v(373)*v(378)+v(293)*v(379)+v(294)*v(380)+v(118)*(v(300)+cetr(5)*v&
       &(373)+v(401))+v(123)*(v(299)+cetr(6)*v(373)+v(402))
  dcedep(3,5)=v(126)*(v(312)+v(349)+cetr(3)*v(374))+v(123)*(cetr(4)*v(306)+cetr(2)*v(307)+cetr(6)*v(374))+v(374)*v(378)+v&
       &(306)*v(379)+v(307)*v(380)+v(118)*(cetr(1)*v(306)+cetr(5)*v(374)+v(419))
  dcedep(3,6)=v(126)*(v(313)+v(350)+cetr(3)*v(375))+v(118)*(cetr(1)*v(307)+cetr(4)*v(344)+cetr(5)*v(375))+v(375)*v(378)+v&
       &(307)*v(379)+v(344)*v(380)+v(123)*(cetr(2)*v(344)+cetr(6)*v(375)+v(419))
  dcedep(4,1)=v(129)*v(283)+v(128)*v(288)+v(123)*v(314)+v(121)*v(320)+v(114)*v(327)+v(130)*v(338)
  dcedep(4,2)=v(128)*v(290)+v(123)*v(315)+v(121)*v(321)+v(114)*v(328)+v(129)*v(333)+v(130)*v(339)
  dcedep(4,3)=v(128)*v(291)+v(123)*v(316)+v(121)*v(322)+v(114)*v(329)+v(129)*v(334)+v(130)*v(341)
  dcedep(4,4)=v(128)*v(292)+v(130)*v(294)+v(123)*v(317)+v(121)*v(324)+v(114)*v(330)+v(129)*v(335)
  dcedep(4,5)=v(123)*v(318)+v(121)*v(325)+v(114)*v(331)+v(129)*v(336)+v(411)+v(416)
  dcedep(4,6)=v(128)*v(294)+v(123)*v(319)+v(121)*v(326)+v(114)*v(332)+v(129)*v(337)+v(130)*v(344)
  dcedep(5,1)=v(130)*v(284)+v(128)*v(301)+v(126)*v(314)+v(123)*v(320)+v(118)*v(327)+v(129)*v(338)
  dcedep(5,2)=v(128)*v(302)+v(126)*v(315)+v(123)*v(321)+v(118)*v(328)+v(130)*v(334)+v(129)*v(339)
  dcedep(5,3)=v(128)*v(304)+v(126)*v(316)+v(123)*v(322)+v(118)*v(329)+v(129)*v(341)+v(130)*v(372)
  dcedep(5,4)=v(126)*v(317)+v(123)*v(324)+v(118)*v(330)+v(130)*v(373)+v(411)+v(415)
  dcedep(5,5)=v(128)*v(306)+v(129)*v(307)+v(126)*v(318)+v(123)*v(325)+v(118)*v(331)+v(130)*v(374)
  dcedep(5,6)=v(128)*v(307)+v(126)*v(319)+v(123)*v(326)+v(118)*v(332)+v(129)*v(344)+v(130)*v(375)
  dcedep(6,1)=v(137)*v(284)+v(135)*v(301)+v(136)*v(338)+v(126)*v(351)+v(123)*v(358)+v(118)*v(365)
  dcedep(6,2)=v(135)*v(302)+v(137)*v(334)+v(136)*v(339)+v(126)*v(352)+v(123)*v(359)+v(118)*v(366)
  dcedep(6,3)=v(135)*v(304)+v(136)*v(341)+v(126)*v(353)+v(123)*v(360)+v(118)*v(367)+v(137)*v(372)
  dcedep(6,4)=v(136)*v(294)+v(126)*v(354)+v(123)*v(361)+v(118)*v(369)+v(137)*v(373)+v(408)
  dcedep(6,5)=v(135)*v(306)+v(136)*v(307)+v(126)*v(356)+v(123)*v(362)+v(118)*v(370)+v(137)*v(374)
  dcedep(6,6)=v(135)*v(307)+v(136)*v(344)+v(126)*v(357)+v(123)*v(363)+v(118)*v(371)+v(137)*v(375)
  dcedcetr(1,1)=(v(109)*v(109))
  dcedcetr(1,2)=v(425)
  dcedcetr(1,3)=v(426)
  dcedcetr(1,4)=2d0*v(515)
  dcedcetr(1,5)=2d0*v(516)
  dcedcetr(1,6)=v(114)*v(517)
  dcedcetr(2,1)=v(425)
  dcedcetr(2,2)=(v(121)*v(121))
  dcedcetr(2,3)=v(431)
  dcedcetr(2,4)=2d0*v(518)
  dcedcetr(2,5)=v(114)*v(519)
  dcedcetr(2,6)=2d0*v(520)
  dcedcetr(3,1)=v(426)
  dcedcetr(3,2)=v(431)
  dcedcetr(3,3)=(v(126)*v(126))
  dcedcetr(3,4)=v(123)*v(517)
  dcedcetr(3,5)=v(126)*v(517)
  dcedcetr(3,6)=v(126)*v(519)
  dcedcetr(4,1)=v(515)
  dcedcetr(4,2)=v(518)
  dcedcetr(4,3)=v(436)
  dcedcetr(4,4)=v(109)*v(121)+v(425)
  dcedcetr(4,5)=v(443)
  dcedcetr(4,6)=v(444)
  dcedcetr(5,1)=v(516)
  dcedcetr(5,2)=v(433)
  dcedcetr(5,3)=v(118)*v(126)
  dcedcetr(5,4)=v(443)
  dcedcetr(5,5)=v(109)*v(126)+v(426)
  dcedcetr(5,6)=v(449)
  dcedcetr(6,1)=v(429)
  dcedcetr(6,2)=v(520)
  dcedcetr(6,3)=v(123)*v(126)
  dcedcetr(6,4)=v(444)
  dcedcetr(6,5)=v(449)
  dcedcetr(6,6)=v(121)*v(126)+v(431)
END SUBROUTINE doublesideed


subroutine consistentproduct(a,b,prod,dprodda,dproddb)
  implicit none
  double precision v(36),a(6),b(6),prod,dprodda(6),dproddb(6)
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
end subroutine consistentproduct

!------------------
!*** Voigt product
!------------------
subroutine consistentproductplain(a,b,prod)
  implicit none
  double precision v(36),a(6),b(6),prod,dprodda(6),dproddb(6)
  v(31)=2d0*b(6)
  v(30)=2d0*b(5)
  v(29)=2d0*b(4)
  prod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*v(29)+a(5)*v(30)+a(6)*v(31)
end subroutine consistentproductplain

!--------------
!*** ce from c
!--------------
subroutine cefromc2(c,g,ce,dcedc)
  implicit none
  double precision v(121),c(6),g(3,3),ce(6),dcedc(6,6)
  v(116)=2d0*g(2,3)
  v(115)=2d0*g(2,2)
  v(114)=2d0*g(2,1)
  v(113)=c(1)*g(1,2)+c(4)*g(2,2)+c(5)*g(3,2)
  v(112)=c(4)*g(1,2)+c(2)*g(2,2)+c(6)*g(3,2)
  v(111)=c(5)*g(1,2)+c(6)*g(2,2)+c(3)*g(3,2)
  v(110)=c(1)*g(1,1)+c(4)*g(2,1)+c(5)*g(3,1)
  v(109)=c(4)*g(1,1)+c(2)*g(2,1)+c(6)*g(3,1)
  v(108)=c(5)*g(1,1)+c(6)*g(2,1)+c(3)*g(3,1)
  ce(1)=g(3,1)*v(108)+g(2,1)*v(109)+g(1,1)*v(110)
  ce(2)=g(3,2)*v(111)+g(2,2)*v(112)+g(1,2)*v(113)
  ce(3)=g(3,3)*(c(5)*g(1,3)+c(6)*g(2,3)+c(3)*g(3,3))+g(1,3)*(c(1)*g(1,3)+c(4)*g(2,3)+c(5)*g(3,3))+g(2,3)*(c(4)*g(1,3)+c(2&
       &)*g(2,3)+c(6)*g(3,3))
  ce(4)=g(3,2)*v(108)+g(2,2)*v(109)+g(1,2)*v(110)
  ce(5)=g(3,3)*v(108)+g(2,3)*v(109)+g(1,3)*v(110)
  ce(6)=g(3,3)*v(111)+g(2,3)*v(112)+g(1,3)*v(113)
  dcedc(1,1)=g(1,1)**2
  dcedc(1,2)=g(2,1)**2
  dcedc(1,3)=g(3,1)**2
  dcedc(1,4)=g(1,1)*v(114)
  dcedc(1,5)=2d0*g(1,1)*g(3,1)
  dcedc(1,6)=g(3,1)*v(114)
  dcedc(2,1)=g(1,2)**2
  dcedc(2,2)=g(2,2)**2
  dcedc(2,3)=g(3,2)**2
  dcedc(2,4)=g(1,2)*v(115)
  dcedc(2,5)=2d0*g(1,2)*g(3,2)
  dcedc(2,6)=g(3,2)*v(115)
  dcedc(3,1)=g(1,3)**2
  dcedc(3,2)=g(2,3)**2
  dcedc(3,3)=g(3,3)**2
  dcedc(3,4)=g(1,3)*v(116)
  dcedc(3,5)=2d0*g(1,3)*g(3,3)
  dcedc(3,6)=g(3,3)*v(116)
  dcedc(4,1)=g(1,1)*g(1,2)
  dcedc(4,2)=g(2,1)*g(2,2)
  dcedc(4,3)=g(3,1)*g(3,2)
  dcedc(4,4)=g(1,2)*g(2,1)+g(1,1)*g(2,2)
  dcedc(4,5)=g(1,2)*g(3,1)+g(1,1)*g(3,2)
  dcedc(4,6)=g(2,2)*g(3,1)+g(2,1)*g(3,2)
  dcedc(5,1)=g(1,1)*g(1,3)
  dcedc(5,2)=g(2,1)*g(2,3)
  dcedc(5,3)=g(3,1)*g(3,3)
  dcedc(5,4)=g(1,3)*g(2,1)+g(1,1)*g(2,3)
  dcedc(5,5)=g(1,3)*g(3,1)+g(1,1)*g(3,3)
  dcedc(5,6)=g(2,3)*g(3,1)+g(2,1)*g(3,3)
  dcedc(6,1)=g(1,2)*g(1,3)
  dcedc(6,2)=g(2,2)*g(2,3)
  dcedc(6,3)=g(3,2)*g(3,3)
  dcedc(6,4)=g(1,3)*g(2,2)+g(1,2)*g(2,3)
  dcedc(6,5)=g(1,3)*g(3,2)+g(1,2)*g(3,3)
  dcedc(6,6)=g(2,3)*g(3,2)+g(2,2)*g(3,3)
end subroutine cefromc2

!----------------------
!*** derivatives of ce
!----------------------
subroutine cederivatives(norm,gamma,dcedep,ce,dcednorm,dcedgamma)
  implicit real(8)(a-h,o-z)
  real(8)::norm(6),ce(6),dcedep(6,6),dcednorm(6,6),dcedgamma(6)
  do id=1,6
     do jd=1,6
        dcednorm(id,jd)=dcedep(id,jd)*gamma
     end do
     dcedgamma(id)=0.0d00
     do jd=1,6
        dcedgamma(id)=dcedgamma(id)+dcedep(id,jd)*norm(jd)
     end do
  end do
end subroutine cederivatives

!-----------------------------------------
!*** calculates symmetrized mandel stress
!-----------------------------------------
subroutine detmandelsymmetric(ce,se,dse,mandel,dmandeldce)
  implicit none
  double precision v(190),ce(6),se(6),dse(6,6),mandel(6),dmandeldce(6,6)
  v(185)=se(2)+se(3)
  v(184)=ce(2)+ce(3)
  v(183)=se(1)+se(3)
  v(182)=ce(1)+ce(3)
  v(181)=se(1)+se(2)
  v(180)=ce(1)+ce(2)
  v(179)=ce(6)*dse(6,6)+se(6)
  v(178)=ce(6)*dse(6,5)
  v(177)=ce(6)*dse(6,4)
  v(176)=ce(6)*dse(6,3)
  v(175)=ce(6)*dse(6,2)
  v(174)=ce(6)*dse(6,1)
  v(173)=ce(5)*dse(5,6)
  v(172)=ce(5)*dse(5,5)+se(5)
  v(171)=ce(5)*dse(5,4)
  v(170)=ce(5)*dse(5,3)
  v(169)=ce(5)*dse(5,2)
  v(168)=ce(5)*dse(5,1)
  v(167)=ce(4)*dse(4,6)
  v(166)=ce(4)*dse(4,5)
  v(165)=ce(4)*dse(4,4)+se(4)
  v(164)=ce(4)*dse(4,3)
  v(163)=ce(4)*dse(4,2)
  v(162)=ce(4)*dse(4,1)
  v(161)=ce(6)*se(6)
  v(160)=ce(5)*se(5)
  v(159)=ce(4)*se(4)
  mandel(1)=ce(1)*se(1)+v(159)+v(160)
  mandel(2)=ce(2)*se(2)+v(159)+v(161)
  mandel(3)=ce(3)*se(3)+v(160)+v(161)
  mandel(4)=(ce(6)*se(5)+ce(5)*se(6)+se(4)*v(180)+ce(4)*v(181))/2d0
  mandel(5)=(ce(6)*se(4)+ce(4)*se(6)+se(5)*v(182)+ce(5)*v(183))/2d0
  mandel(6)=(ce(5)*se(4)+ce(4)*se(5)+se(6)*v(184)+ce(6)*v(185))/2d0
  dmandeldce(1,1)=ce(1)*dse(1,1)+se(1)+v(162)+v(168)
  dmandeldce(1,2)=ce(1)*dse(1,2)+v(163)+v(169)
  dmandeldce(1,3)=ce(1)*dse(1,3)+v(164)+v(170)
  dmandeldce(1,4)=ce(1)*dse(1,4)+v(165)+v(171)
  dmandeldce(1,5)=ce(1)*dse(1,5)+v(166)+v(172)
  dmandeldce(1,6)=ce(1)*dse(1,6)+v(167)+v(173)
  dmandeldce(2,1)=ce(2)*dse(2,1)+v(162)+v(174)
  dmandeldce(2,2)=ce(2)*dse(2,2)+se(2)+v(163)+v(175)
  dmandeldce(2,3)=ce(2)*dse(2,3)+v(164)+v(176)
  dmandeldce(2,4)=ce(2)*dse(2,4)+v(165)+v(177)
  dmandeldce(2,5)=ce(2)*dse(2,5)+v(166)+v(178)
  dmandeldce(2,6)=ce(2)*dse(2,6)+v(167)+v(179)
  dmandeldce(3,1)=ce(3)*dse(3,1)+v(168)+v(174)
  dmandeldce(3,2)=ce(3)*dse(3,2)+v(169)+v(175)
  dmandeldce(3,3)=ce(3)*dse(3,3)+se(3)+v(170)+v(176)
  dmandeldce(3,4)=ce(3)*dse(3,4)+v(171)+v(177)
  dmandeldce(3,5)=ce(3)*dse(3,5)+v(172)+v(178)
  dmandeldce(3,6)=ce(3)*dse(3,6)+v(173)+v(179)
  dmandeldce(4,1)=(ce(4)*(dse(1,1)+dse(2,1))+ce(6)*dse(5,1)+ce(5)*dse(6,1)+se(4)+dse(4,1)*v(180))/2d0
  dmandeldce(4,2)=(ce(4)*(dse(1,2)+dse(2,2))+ce(6)*dse(5,2)+ce(5)*dse(6,2)+se(4)+dse(4,2)*v(180))/2d0
  dmandeldce(4,3)=(ce(4)*(dse(1,3)+dse(2,3))+ce(6)*dse(5,3)+ce(5)*dse(6,3)+dse(4,3)*v(180))/2d0
  dmandeldce(4,4)=(ce(4)*(dse(1,4)+dse(2,4))+ce(6)*dse(5,4)+ce(5)*dse(6,4)+dse(4,4)*v(180)+v(181))/2d0
  dmandeldce(4,5)=(ce(4)*(dse(1,5)+dse(2,5))+ce(6)*dse(5,5)+ce(5)*dse(6,5)+se(6)+dse(4,5)*v(180))/2d0
  dmandeldce(4,6)=(ce(4)*(dse(1,6)+dse(2,6))+ce(6)*dse(5,6)+ce(5)*dse(6,6)+se(5)+dse(4,6)*v(180))/2d0
  dmandeldce(5,1)=(ce(5)*(dse(1,1)+dse(3,1))+ce(6)*dse(4,1)+ce(4)*dse(6,1)+se(5)+dse(5,1)*v(182))/2d0
  dmandeldce(5,2)=(ce(5)*(dse(1,2)+dse(3,2))+ce(6)*dse(4,2)+ce(4)*dse(6,2)+dse(5,2)*v(182))/2d0
  dmandeldce(5,3)=(ce(5)*(dse(1,3)+dse(3,3))+ce(6)*dse(4,3)+ce(4)*dse(6,3)+se(5)+dse(5,3)*v(182))/2d0
  dmandeldce(5,4)=(ce(5)*(dse(1,4)+dse(3,4))+ce(6)*dse(4,4)+ce(4)*dse(6,4)+se(6)+dse(5,4)*v(182))/2d0
  dmandeldce(5,5)=(ce(5)*(dse(1,5)+dse(3,5))+ce(6)*dse(4,5)+ce(4)*dse(6,5)+dse(5,5)*v(182)+v(183))/2d0
  dmandeldce(5,6)=(ce(5)*(dse(1,6)+dse(3,6))+ce(6)*dse(4,6)+ce(4)*dse(6,6)+se(4)+dse(5,6)*v(182))/2d0
  dmandeldce(6,1)=(ce(6)*(dse(2,1)+dse(3,1))+ce(5)*dse(4,1)+ce(4)*dse(5,1)+dse(6,1)*v(184))/2d0
  dmandeldce(6,2)=(ce(6)*(dse(2,2)+dse(3,2))+ce(5)*dse(4,2)+ce(4)*dse(5,2)+se(6)+dse(6,2)*v(184))/2d0
  dmandeldce(6,3)=(ce(6)*(dse(2,3)+dse(3,3))+ce(5)*dse(4,3)+ce(4)*dse(5,3)+se(6)+dse(6,3)*v(184))/2d0
  dmandeldce(6,4)=(ce(6)*(dse(2,4)+dse(3,4))+ce(5)*dse(4,4)+ce(4)*dse(5,4)+se(5)+dse(6,4)*v(184))/2d0
  dmandeldce(6,5)=(ce(6)*(dse(2,5)+dse(3,5))+ce(5)*dse(4,5)+ce(4)*dse(5,5)+se(4)+dse(6,5)*v(184))/2d0
  dmandeldce(6,6)=(ce(6)*(dse(2,6)+dse(3,6))+ce(5)*dse(4,6)+ce(4)*dse(5,6)+dse(6,6)*v(184)+v(185))/2d0
end subroutine detmandelsymmetric

!-----------------------------------------
!*** mandel derivatives
!*** from dmandeldce, dcedep and dcedcetr
!*** dmandeldep, dmandeldcetr
!-----------------------------------------
subroutine mandelderivatives(dmandeldce,dcedep,dcedcetr,dmandeldep,dmandeldcetr)
  implicit real(8)(a-h,o-z)
  real(8),dimension(6,6)::dmandeldce,dcedep,dcedcetr,dmandeldep,dmandeldcetr
  dmandeldep=0.0d00
  dmandeldcetr=0.0d000
  do id=1,6
     do jd=1,6
        do kd=1,6
           dmandeldep(id,jd)=dmandeldep(id,jd)+dmandeldce(id,kd)*dcedep(kd,jd)
           dmandeldcetr(id,jd)=dmandeldcetr(id,jd)+dmandeldce(id,kd)*dcedcetr(kd,jd)
        end do
     end do
  end do
end subroutine mandelderivatives



!**************************************************************
!* AceGen    6.808 Linux (6 Sep 16)                           *
!*           Co. J. Korelc  2013           25 Oct 25 14:28:52 *
!**************************************************************
! User     : Full professional version
! Notebook : tangentmodulus
! Evaluation time                 : 141 s   Mode  : Optimal
! Number of formulae              : 473     Method: Automatic
! Subroutine                      : tangentmodulus size: 18513
! Total size of Mathematica  code : 18513 subexpressions
! Total size of Fortran code      : 38144 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE tangentmodulus(gn,gn1,ep,se,depdcetr,dcetrdc,dsedcetr,s,dsdc)
  IMPLICIT NONE
  DOUBLE PRECISION v(784),gn(3,3),gn1(3,3),ep(6),se(6),depdcetr(6,6),dcetrdc(6,6),dsedcetr(6,6),s(6),dsdc(6,6)
  v(767)=ep(3)/8d0
  v(761)=(-0.25d0)*ep(6)
  v(760)=ep(2)/8d0
  v(756)=ep(1)/8d0
  v(754)=(-0.25d0)*ep(4)
  v(752)=(-0.25d0)*ep(5)
  v(751)=ep(6)/8d0
  v(750)=ep(5)/8d0
  v(749)=ep(4)/8d0
  v(748)=dcetrdc(1,6)*dsedcetr(6,1)+dcetrdc(2,6)*dsedcetr(6,2)+dcetrdc(3,6)*dsedcetr(6,3)+dcetrdc(4,6)*dsedcetr(6,4)&
       &+dcetrdc(5,6)*dsedcetr(6,5)+dcetrdc(6,6)*dsedcetr(6,6)
  v(747)=dcetrdc(1,5)*dsedcetr(6,1)+dcetrdc(2,5)*dsedcetr(6,2)+dcetrdc(3,5)*dsedcetr(6,3)+dcetrdc(4,5)*dsedcetr(6,4)&
       &+dcetrdc(5,5)*dsedcetr(6,5)+dcetrdc(6,5)*dsedcetr(6,6)
  v(746)=dcetrdc(1,4)*dsedcetr(6,1)+dcetrdc(2,4)*dsedcetr(6,2)+dcetrdc(3,4)*dsedcetr(6,3)+dcetrdc(4,4)*dsedcetr(6,4)&
       &+dcetrdc(5,4)*dsedcetr(6,5)+dcetrdc(6,4)*dsedcetr(6,6)
  v(745)=dcetrdc(1,3)*dsedcetr(6,1)+dcetrdc(2,3)*dsedcetr(6,2)+dcetrdc(3,3)*dsedcetr(6,3)+dcetrdc(4,3)*dsedcetr(6,4)&
       &+dcetrdc(5,3)*dsedcetr(6,5)+dcetrdc(6,3)*dsedcetr(6,6)
  v(744)=dcetrdc(1,2)*dsedcetr(6,1)+dcetrdc(2,2)*dsedcetr(6,2)+dcetrdc(3,2)*dsedcetr(6,3)+dcetrdc(4,2)*dsedcetr(6,4)&
       &+dcetrdc(5,2)*dsedcetr(6,5)+dcetrdc(6,2)*dsedcetr(6,6)
  v(743)=dcetrdc(1,1)*dsedcetr(6,1)+dcetrdc(2,1)*dsedcetr(6,2)+dcetrdc(3,1)*dsedcetr(6,3)+dcetrdc(4,1)*dsedcetr(6,4)&
       &+dcetrdc(5,1)*dsedcetr(6,5)+dcetrdc(6,1)*dsedcetr(6,6)
  v(742)=dcetrdc(1,6)*dsedcetr(5,1)+dcetrdc(2,6)*dsedcetr(5,2)+dcetrdc(3,6)*dsedcetr(5,3)+dcetrdc(4,6)*dsedcetr(5,4)&
       &+dcetrdc(5,6)*dsedcetr(5,5)+dcetrdc(6,6)*dsedcetr(5,6)
  v(741)=dcetrdc(1,5)*dsedcetr(5,1)+dcetrdc(2,5)*dsedcetr(5,2)+dcetrdc(3,5)*dsedcetr(5,3)+dcetrdc(4,5)*dsedcetr(5,4)&
       &+dcetrdc(5,5)*dsedcetr(5,5)+dcetrdc(6,5)*dsedcetr(5,6)
  v(740)=dcetrdc(1,4)*dsedcetr(5,1)+dcetrdc(2,4)*dsedcetr(5,2)+dcetrdc(3,4)*dsedcetr(5,3)+dcetrdc(4,4)*dsedcetr(5,4)&
       &+dcetrdc(5,4)*dsedcetr(5,5)+dcetrdc(6,4)*dsedcetr(5,6)
  v(739)=dcetrdc(1,3)*dsedcetr(5,1)+dcetrdc(2,3)*dsedcetr(5,2)+dcetrdc(3,3)*dsedcetr(5,3)+dcetrdc(4,3)*dsedcetr(5,4)&
       &+dcetrdc(5,3)*dsedcetr(5,5)+dcetrdc(6,3)*dsedcetr(5,6)
  v(738)=dcetrdc(1,2)*dsedcetr(5,1)+dcetrdc(2,2)*dsedcetr(5,2)+dcetrdc(3,2)*dsedcetr(5,3)+dcetrdc(4,2)*dsedcetr(5,4)&
       &+dcetrdc(5,2)*dsedcetr(5,5)+dcetrdc(6,2)*dsedcetr(5,6)
  v(737)=dcetrdc(1,1)*dsedcetr(5,1)+dcetrdc(2,1)*dsedcetr(5,2)+dcetrdc(3,1)*dsedcetr(5,3)+dcetrdc(4,1)*dsedcetr(5,4)&
       &+dcetrdc(5,1)*dsedcetr(5,5)+dcetrdc(6,1)*dsedcetr(5,6)
  v(736)=dcetrdc(1,6)*dsedcetr(4,1)+dcetrdc(2,6)*dsedcetr(4,2)+dcetrdc(3,6)*dsedcetr(4,3)+dcetrdc(4,6)*dsedcetr(4,4)&
       &+dcetrdc(5,6)*dsedcetr(4,5)+dcetrdc(6,6)*dsedcetr(4,6)
  v(735)=dcetrdc(1,5)*dsedcetr(4,1)+dcetrdc(2,5)*dsedcetr(4,2)+dcetrdc(3,5)*dsedcetr(4,3)+dcetrdc(4,5)*dsedcetr(4,4)&
       &+dcetrdc(5,5)*dsedcetr(4,5)+dcetrdc(6,5)*dsedcetr(4,6)
  v(734)=dcetrdc(1,4)*dsedcetr(4,1)+dcetrdc(2,4)*dsedcetr(4,2)+dcetrdc(3,4)*dsedcetr(4,3)+dcetrdc(4,4)*dsedcetr(4,4)&
       &+dcetrdc(5,4)*dsedcetr(4,5)+dcetrdc(6,4)*dsedcetr(4,6)
  v(733)=dcetrdc(1,3)*dsedcetr(4,1)+dcetrdc(2,3)*dsedcetr(4,2)+dcetrdc(3,3)*dsedcetr(4,3)+dcetrdc(4,3)*dsedcetr(4,4)&
       &+dcetrdc(5,3)*dsedcetr(4,5)+dcetrdc(6,3)*dsedcetr(4,6)
  v(732)=dcetrdc(1,2)*dsedcetr(4,1)+dcetrdc(2,2)*dsedcetr(4,2)+dcetrdc(3,2)*dsedcetr(4,3)+dcetrdc(4,2)*dsedcetr(4,4)&
       &+dcetrdc(5,2)*dsedcetr(4,5)+dcetrdc(6,2)*dsedcetr(4,6)
  v(731)=dcetrdc(1,1)*dsedcetr(4,1)+dcetrdc(2,1)*dsedcetr(4,2)+dcetrdc(3,1)*dsedcetr(4,3)+dcetrdc(4,1)*dsedcetr(4,4)&
       &+dcetrdc(5,1)*dsedcetr(4,5)+dcetrdc(6,1)*dsedcetr(4,6)
  v(730)=dcetrdc(1,6)*dsedcetr(3,1)+dcetrdc(2,6)*dsedcetr(3,2)+dcetrdc(3,6)*dsedcetr(3,3)+dcetrdc(4,6)*dsedcetr(3,4)&
       &+dcetrdc(5,6)*dsedcetr(3,5)+dcetrdc(6,6)*dsedcetr(3,6)
  v(729)=dcetrdc(1,5)*dsedcetr(3,1)+dcetrdc(2,5)*dsedcetr(3,2)+dcetrdc(3,5)*dsedcetr(3,3)+dcetrdc(4,5)*dsedcetr(3,4)&
       &+dcetrdc(5,5)*dsedcetr(3,5)+dcetrdc(6,5)*dsedcetr(3,6)
  v(728)=dcetrdc(1,4)*dsedcetr(3,1)+dcetrdc(2,4)*dsedcetr(3,2)+dcetrdc(3,4)*dsedcetr(3,3)+dcetrdc(4,4)*dsedcetr(3,4)&
       &+dcetrdc(5,4)*dsedcetr(3,5)+dcetrdc(6,4)*dsedcetr(3,6)
  v(727)=dcetrdc(1,3)*dsedcetr(3,1)+dcetrdc(2,3)*dsedcetr(3,2)+dcetrdc(3,3)*dsedcetr(3,3)+dcetrdc(4,3)*dsedcetr(3,4)&
       &+dcetrdc(5,3)*dsedcetr(3,5)+dcetrdc(6,3)*dsedcetr(3,6)
  v(726)=dcetrdc(1,2)*dsedcetr(3,1)+dcetrdc(2,2)*dsedcetr(3,2)+dcetrdc(3,2)*dsedcetr(3,3)+dcetrdc(4,2)*dsedcetr(3,4)&
       &+dcetrdc(5,2)*dsedcetr(3,5)+dcetrdc(6,2)*dsedcetr(3,6)
  v(725)=dcetrdc(1,1)*dsedcetr(3,1)+dcetrdc(2,1)*dsedcetr(3,2)+dcetrdc(3,1)*dsedcetr(3,3)+dcetrdc(4,1)*dsedcetr(3,4)&
       &+dcetrdc(5,1)*dsedcetr(3,5)+dcetrdc(6,1)*dsedcetr(3,6)
  v(724)=dcetrdc(1,6)*dsedcetr(2,1)+dcetrdc(2,6)*dsedcetr(2,2)+dcetrdc(3,6)*dsedcetr(2,3)+dcetrdc(4,6)*dsedcetr(2,4)&
       &+dcetrdc(5,6)*dsedcetr(2,5)+dcetrdc(6,6)*dsedcetr(2,6)
  v(723)=dcetrdc(1,5)*dsedcetr(2,1)+dcetrdc(2,5)*dsedcetr(2,2)+dcetrdc(3,5)*dsedcetr(2,3)+dcetrdc(4,5)*dsedcetr(2,4)&
       &+dcetrdc(5,5)*dsedcetr(2,5)+dcetrdc(6,5)*dsedcetr(2,6)
  v(722)=dcetrdc(1,4)*dsedcetr(2,1)+dcetrdc(2,4)*dsedcetr(2,2)+dcetrdc(3,4)*dsedcetr(2,3)+dcetrdc(4,4)*dsedcetr(2,4)&
       &+dcetrdc(5,4)*dsedcetr(2,5)+dcetrdc(6,4)*dsedcetr(2,6)
  v(721)=dcetrdc(1,3)*dsedcetr(2,1)+dcetrdc(2,3)*dsedcetr(2,2)+dcetrdc(3,3)*dsedcetr(2,3)+dcetrdc(4,3)*dsedcetr(2,4)&
       &+dcetrdc(5,3)*dsedcetr(2,5)+dcetrdc(6,3)*dsedcetr(2,6)
  v(720)=dcetrdc(1,2)*dsedcetr(2,1)+dcetrdc(2,2)*dsedcetr(2,2)+dcetrdc(3,2)*dsedcetr(2,3)+dcetrdc(4,2)*dsedcetr(2,4)&
       &+dcetrdc(5,2)*dsedcetr(2,5)+dcetrdc(6,2)*dsedcetr(2,6)
  v(719)=dcetrdc(1,1)*dsedcetr(2,1)+dcetrdc(2,1)*dsedcetr(2,2)+dcetrdc(3,1)*dsedcetr(2,3)+dcetrdc(4,1)*dsedcetr(2,4)&
       &+dcetrdc(5,1)*dsedcetr(2,5)+dcetrdc(6,1)*dsedcetr(2,6)
  v(718)=dcetrdc(1,6)*dsedcetr(1,1)+dcetrdc(2,6)*dsedcetr(1,2)+dcetrdc(3,6)*dsedcetr(1,3)+dcetrdc(4,6)*dsedcetr(1,4)&
       &+dcetrdc(5,6)*dsedcetr(1,5)+dcetrdc(6,6)*dsedcetr(1,6)
  v(717)=dcetrdc(1,5)*dsedcetr(1,1)+dcetrdc(2,5)*dsedcetr(1,2)+dcetrdc(3,5)*dsedcetr(1,3)+dcetrdc(4,5)*dsedcetr(1,4)&
       &+dcetrdc(5,5)*dsedcetr(1,5)+dcetrdc(6,5)*dsedcetr(1,6)
  v(716)=dcetrdc(1,4)*dsedcetr(1,1)+dcetrdc(2,4)*dsedcetr(1,2)+dcetrdc(3,4)*dsedcetr(1,3)+dcetrdc(4,4)*dsedcetr(1,4)&
       &+dcetrdc(5,4)*dsedcetr(1,5)+dcetrdc(6,4)*dsedcetr(1,6)
  v(715)=dcetrdc(1,3)*dsedcetr(1,1)+dcetrdc(2,3)*dsedcetr(1,2)+dcetrdc(3,3)*dsedcetr(1,3)+dcetrdc(4,3)*dsedcetr(1,4)&
       &+dcetrdc(5,3)*dsedcetr(1,5)+dcetrdc(6,3)*dsedcetr(1,6)
  v(714)=dcetrdc(1,2)*dsedcetr(1,1)+dcetrdc(2,2)*dsedcetr(1,2)+dcetrdc(3,2)*dsedcetr(1,3)+dcetrdc(4,2)*dsedcetr(1,4)&
       &+dcetrdc(5,2)*dsedcetr(1,5)+dcetrdc(6,2)*dsedcetr(1,6)
  v(713)=dcetrdc(1,1)*dsedcetr(1,1)+dcetrdc(2,1)*dsedcetr(1,2)+dcetrdc(3,1)*dsedcetr(1,3)+dcetrdc(4,1)*dsedcetr(1,4)&
       &+dcetrdc(5,1)*dsedcetr(1,5)+dcetrdc(6,1)*dsedcetr(1,6)
  v(712)=dcetrdc(1,6)*depdcetr(6,1)+dcetrdc(2,6)*depdcetr(6,2)+dcetrdc(3,6)*depdcetr(6,3)+dcetrdc(4,6)*depdcetr(6,4)&
       &+dcetrdc(5,6)*depdcetr(6,5)+dcetrdc(6,6)*depdcetr(6,6)
  v(711)=dcetrdc(1,5)*depdcetr(6,1)+dcetrdc(2,5)*depdcetr(6,2)+dcetrdc(3,5)*depdcetr(6,3)+dcetrdc(4,5)*depdcetr(6,4)&
       &+dcetrdc(5,5)*depdcetr(6,5)+dcetrdc(6,5)*depdcetr(6,6)
  v(710)=dcetrdc(1,4)*depdcetr(6,1)+dcetrdc(2,4)*depdcetr(6,2)+dcetrdc(3,4)*depdcetr(6,3)+dcetrdc(4,4)*depdcetr(6,4)&
       &+dcetrdc(5,4)*depdcetr(6,5)+dcetrdc(6,4)*depdcetr(6,6)
  v(709)=dcetrdc(1,3)*depdcetr(6,1)+dcetrdc(2,3)*depdcetr(6,2)+dcetrdc(3,3)*depdcetr(6,3)+dcetrdc(4,3)*depdcetr(6,4)&
       &+dcetrdc(5,3)*depdcetr(6,5)+dcetrdc(6,3)*depdcetr(6,6)
  v(708)=dcetrdc(1,2)*depdcetr(6,1)+dcetrdc(2,2)*depdcetr(6,2)+dcetrdc(3,2)*depdcetr(6,3)+dcetrdc(4,2)*depdcetr(6,4)&
       &+dcetrdc(5,2)*depdcetr(6,5)+dcetrdc(6,2)*depdcetr(6,6)
  v(707)=dcetrdc(1,1)*depdcetr(6,1)+dcetrdc(2,1)*depdcetr(6,2)+dcetrdc(3,1)*depdcetr(6,3)+dcetrdc(4,1)*depdcetr(6,4)&
       &+dcetrdc(5,1)*depdcetr(6,5)+dcetrdc(6,1)*depdcetr(6,6)
  v(706)=dcetrdc(1,6)*depdcetr(5,1)+dcetrdc(2,6)*depdcetr(5,2)+dcetrdc(3,6)*depdcetr(5,3)+dcetrdc(4,6)*depdcetr(5,4)&
       &+dcetrdc(5,6)*depdcetr(5,5)+dcetrdc(6,6)*depdcetr(5,6)
  v(705)=dcetrdc(1,5)*depdcetr(5,1)+dcetrdc(2,5)*depdcetr(5,2)+dcetrdc(3,5)*depdcetr(5,3)+dcetrdc(4,5)*depdcetr(5,4)&
       &+dcetrdc(5,5)*depdcetr(5,5)+dcetrdc(6,5)*depdcetr(5,6)
  v(704)=dcetrdc(1,4)*depdcetr(5,1)+dcetrdc(2,4)*depdcetr(5,2)+dcetrdc(3,4)*depdcetr(5,3)+dcetrdc(4,4)*depdcetr(5,4)&
       &+dcetrdc(5,4)*depdcetr(5,5)+dcetrdc(6,4)*depdcetr(5,6)
  v(703)=dcetrdc(1,3)*depdcetr(5,1)+dcetrdc(2,3)*depdcetr(5,2)+dcetrdc(3,3)*depdcetr(5,3)+dcetrdc(4,3)*depdcetr(5,4)&
       &+dcetrdc(5,3)*depdcetr(5,5)+dcetrdc(6,3)*depdcetr(5,6)
  v(702)=dcetrdc(1,2)*depdcetr(5,1)+dcetrdc(2,2)*depdcetr(5,2)+dcetrdc(3,2)*depdcetr(5,3)+dcetrdc(4,2)*depdcetr(5,4)&
       &+dcetrdc(5,2)*depdcetr(5,5)+dcetrdc(6,2)*depdcetr(5,6)
  v(701)=dcetrdc(1,1)*depdcetr(5,1)+dcetrdc(2,1)*depdcetr(5,2)+dcetrdc(3,1)*depdcetr(5,3)+dcetrdc(4,1)*depdcetr(5,4)&
       &+dcetrdc(5,1)*depdcetr(5,5)+dcetrdc(6,1)*depdcetr(5,6)
  v(700)=dcetrdc(1,6)*depdcetr(4,1)+dcetrdc(2,6)*depdcetr(4,2)+dcetrdc(3,6)*depdcetr(4,3)+dcetrdc(4,6)*depdcetr(4,4)&
       &+dcetrdc(5,6)*depdcetr(4,5)+dcetrdc(6,6)*depdcetr(4,6)
  v(699)=dcetrdc(1,5)*depdcetr(4,1)+dcetrdc(2,5)*depdcetr(4,2)+dcetrdc(3,5)*depdcetr(4,3)+dcetrdc(4,5)*depdcetr(4,4)&
       &+dcetrdc(5,5)*depdcetr(4,5)+dcetrdc(6,5)*depdcetr(4,6)
  v(698)=dcetrdc(1,4)*depdcetr(4,1)+dcetrdc(2,4)*depdcetr(4,2)+dcetrdc(3,4)*depdcetr(4,3)+dcetrdc(4,4)*depdcetr(4,4)&
       &+dcetrdc(5,4)*depdcetr(4,5)+dcetrdc(6,4)*depdcetr(4,6)
  v(697)=dcetrdc(1,3)*depdcetr(4,1)+dcetrdc(2,3)*depdcetr(4,2)+dcetrdc(3,3)*depdcetr(4,3)+dcetrdc(4,3)*depdcetr(4,4)&
       &+dcetrdc(5,3)*depdcetr(4,5)+dcetrdc(6,3)*depdcetr(4,6)
  v(696)=dcetrdc(1,2)*depdcetr(4,1)+dcetrdc(2,2)*depdcetr(4,2)+dcetrdc(3,2)*depdcetr(4,3)+dcetrdc(4,2)*depdcetr(4,4)&
       &+dcetrdc(5,2)*depdcetr(4,5)+dcetrdc(6,2)*depdcetr(4,6)
  v(695)=dcetrdc(1,1)*depdcetr(4,1)+dcetrdc(2,1)*depdcetr(4,2)+dcetrdc(3,1)*depdcetr(4,3)+dcetrdc(4,1)*depdcetr(4,4)&
       &+dcetrdc(5,1)*depdcetr(4,5)+dcetrdc(6,1)*depdcetr(4,6)
  v(694)=dcetrdc(1,6)*depdcetr(3,1)+dcetrdc(2,6)*depdcetr(3,2)+dcetrdc(3,6)*depdcetr(3,3)+dcetrdc(4,6)*depdcetr(3,4)&
       &+dcetrdc(5,6)*depdcetr(3,5)+dcetrdc(6,6)*depdcetr(3,6)
  v(693)=dcetrdc(1,5)*depdcetr(3,1)+dcetrdc(2,5)*depdcetr(3,2)+dcetrdc(3,5)*depdcetr(3,3)+dcetrdc(4,5)*depdcetr(3,4)&
       &+dcetrdc(5,5)*depdcetr(3,5)+dcetrdc(6,5)*depdcetr(3,6)
  v(692)=dcetrdc(1,4)*depdcetr(3,1)+dcetrdc(2,4)*depdcetr(3,2)+dcetrdc(3,4)*depdcetr(3,3)+dcetrdc(4,4)*depdcetr(3,4)&
       &+dcetrdc(5,4)*depdcetr(3,5)+dcetrdc(6,4)*depdcetr(3,6)
  v(691)=dcetrdc(1,3)*depdcetr(3,1)+dcetrdc(2,3)*depdcetr(3,2)+dcetrdc(3,3)*depdcetr(3,3)+dcetrdc(4,3)*depdcetr(3,4)&
       &+dcetrdc(5,3)*depdcetr(3,5)+dcetrdc(6,3)*depdcetr(3,6)
  v(690)=dcetrdc(1,2)*depdcetr(3,1)+dcetrdc(2,2)*depdcetr(3,2)+dcetrdc(3,2)*depdcetr(3,3)+dcetrdc(4,2)*depdcetr(3,4)&
       &+dcetrdc(5,2)*depdcetr(3,5)+dcetrdc(6,2)*depdcetr(3,6)
  v(689)=dcetrdc(1,1)*depdcetr(3,1)+dcetrdc(2,1)*depdcetr(3,2)+dcetrdc(3,1)*depdcetr(3,3)+dcetrdc(4,1)*depdcetr(3,4)&
       &+dcetrdc(5,1)*depdcetr(3,5)+dcetrdc(6,1)*depdcetr(3,6)
  v(688)=dcetrdc(1,6)*depdcetr(2,1)+dcetrdc(2,6)*depdcetr(2,2)+dcetrdc(3,6)*depdcetr(2,3)+dcetrdc(4,6)*depdcetr(2,4)&
       &+dcetrdc(5,6)*depdcetr(2,5)+dcetrdc(6,6)*depdcetr(2,6)
  v(687)=dcetrdc(1,5)*depdcetr(2,1)+dcetrdc(2,5)*depdcetr(2,2)+dcetrdc(3,5)*depdcetr(2,3)+dcetrdc(4,5)*depdcetr(2,4)&
       &+dcetrdc(5,5)*depdcetr(2,5)+dcetrdc(6,5)*depdcetr(2,6)
  v(686)=dcetrdc(1,4)*depdcetr(2,1)+dcetrdc(2,4)*depdcetr(2,2)+dcetrdc(3,4)*depdcetr(2,3)+dcetrdc(4,4)*depdcetr(2,4)&
       &+dcetrdc(5,4)*depdcetr(2,5)+dcetrdc(6,4)*depdcetr(2,6)
  v(685)=dcetrdc(1,3)*depdcetr(2,1)+dcetrdc(2,3)*depdcetr(2,2)+dcetrdc(3,3)*depdcetr(2,3)+dcetrdc(4,3)*depdcetr(2,4)&
       &+dcetrdc(5,3)*depdcetr(2,5)+dcetrdc(6,3)*depdcetr(2,6)
  v(684)=dcetrdc(1,2)*depdcetr(2,1)+dcetrdc(2,2)*depdcetr(2,2)+dcetrdc(3,2)*depdcetr(2,3)+dcetrdc(4,2)*depdcetr(2,4)&
       &+dcetrdc(5,2)*depdcetr(2,5)+dcetrdc(6,2)*depdcetr(2,6)
  v(683)=dcetrdc(1,1)*depdcetr(2,1)+dcetrdc(2,1)*depdcetr(2,2)+dcetrdc(3,1)*depdcetr(2,3)+dcetrdc(4,1)*depdcetr(2,4)&
       &+dcetrdc(5,1)*depdcetr(2,5)+dcetrdc(6,1)*depdcetr(2,6)
  v(682)=dcetrdc(1,6)*depdcetr(1,1)+dcetrdc(2,6)*depdcetr(1,2)+dcetrdc(3,6)*depdcetr(1,3)+dcetrdc(4,6)*depdcetr(1,4)&
       &+dcetrdc(5,6)*depdcetr(1,5)+dcetrdc(6,6)*depdcetr(1,6)
  v(681)=dcetrdc(1,5)*depdcetr(1,1)+dcetrdc(2,5)*depdcetr(1,2)+dcetrdc(3,5)*depdcetr(1,3)+dcetrdc(4,5)*depdcetr(1,4)&
       &+dcetrdc(5,5)*depdcetr(1,5)+dcetrdc(6,5)*depdcetr(1,6)
  v(680)=dcetrdc(1,4)*depdcetr(1,1)+dcetrdc(2,4)*depdcetr(1,2)+dcetrdc(3,4)*depdcetr(1,3)+dcetrdc(4,4)*depdcetr(1,4)&
       &+dcetrdc(5,4)*depdcetr(1,5)+dcetrdc(6,4)*depdcetr(1,6)
  v(679)=dcetrdc(1,3)*depdcetr(1,1)+dcetrdc(2,3)*depdcetr(1,2)+dcetrdc(3,3)*depdcetr(1,3)+dcetrdc(4,3)*depdcetr(1,4)&
       &+dcetrdc(5,3)*depdcetr(1,5)+dcetrdc(6,3)*depdcetr(1,6)
  v(678)=dcetrdc(1,2)*depdcetr(1,1)+dcetrdc(2,2)*depdcetr(1,2)+dcetrdc(3,2)*depdcetr(1,3)+dcetrdc(4,2)*depdcetr(1,4)&
       &+dcetrdc(5,2)*depdcetr(1,5)+dcetrdc(6,2)*depdcetr(1,6)
  v(677)=dcetrdc(1,1)*depdcetr(1,1)+dcetrdc(2,1)*depdcetr(1,2)+dcetrdc(3,1)*depdcetr(1,3)+dcetrdc(4,1)*depdcetr(1,4)&
       &+dcetrdc(5,1)*depdcetr(1,5)+dcetrdc(6,1)*depdcetr(1,6)
  v(676)=ep(3)**2/16d0
  v(675)=(-12d0)+ep(3)
  v(672)=(ep(4)*ep(5)+(ep(2)+ep(3))*ep(6))/16d0
  v(762)=(-0.25d0)*v(672)
  v(671)=ep(2)**2/16d0
  v(669)=ep(1)**2/16d0
  v(668)=(ep(1)*ep(4)+ep(2)*ep(4)+ep(5)*ep(6))/16d0
  v(755)=(-0.25d0)*v(668)
  v(667)=(ep(1)*ep(5)+ep(3)*ep(5)+ep(4)*ep(6))/16d0
  v(753)=(-0.25d0)*v(667)
  v(666)=ep(6)**2/16d0
  v(665)=ep(5)**2/16d0
  v(664)=ep(4)**2/16d0
  v(203)=v(667)*v(752)
  v(196)=v(668)*v(754)
  v(673)=6d0+v(196)
  v(187)=v(664)+v(665)+v(669)
  v(670)=6d0+v(187)
  v(188)=(-(ep(1)*v(670))+4d0*(3d0*v(187)+v(203)+v(673)))/24d0
  v(759)=2d0*v(188)
  v(191)=(-(ep(6)*v(667))+12d0*v(668)-ep(2)*v(668)-ep(4)*v(670))/24d0
  v(757)=2d0*v(191)
  v(209)=(v(191)*v(191))
  v(192)=(-(ep(6)*v(668))-ep(5)*v(670)-v(667)*v(675))/24d0
  v(758)=2d0*v(192)
  v(211)=(v(192)*v(192))
  v(206)=(v(188)*v(188))+v(209)+v(211)
  v(774)=2d0*v(206)
  v(194)=v(664)+v(666)+v(671)
  v(674)=6d0+v(194)
  v(202)=v(672)*v(761)
  v(207)=(4d0*(3d0*v(194)+v(202)+v(673))-ep(2)*v(674))/24d0
  v(765)=v(188)+v(207)
  v(764)=2d0*v(207)
  v(198)=(-(ep(5)*v(668))-ep(6)*v(674)-v(672)*v(675))/24d0
  v(763)=2d0*v(198)
  v(212)=(v(198)*v(198))
  v(216)=(v(207)*v(207))+v(209)+v(212)
  v(777)=2d0*v(216)
  v(775)=v(206)+v(216)
  v(197)=v(192)*v(198)+v(191)*v(765)
  v(766)=2d0*v(197)
  v(215)=(v(197)*v(197))
  v(201)=v(665)+v(666)+v(676)
  v(768)=6d0+v(201)
  v(208)=(4d0*(6d0+3d0*v(201)+v(202)+v(203))-ep(3)*v(768))/24d0
  v(772)=v(188)+v(208)
  v(770)=v(207)+v(208)
  v(769)=2d0*v(208)
  v(218)=(v(208)*v(208))+v(211)+v(212)
  v(779)=2d0*v(218)
  v(778)=v(216)+v(218)
  v(776)=v(206)+v(218)
  v(213)=v(191)*v(192)+v(198)*v(770)
  v(771)=2d0*v(213)
  v(221)=(v(213)*v(213))
  v(204)=v(191)*v(198)+v(192)*v(772)
  v(773)=2d0*v(204)
  v(220)=(v(204)*v(204))
  v(205)=(v(206)*v(206))+v(215)+v(220)
  v(210)=v(204)*v(213)+v(197)*v(775)
  v(214)=v(197)*v(213)+v(204)*v(776)
  v(238)=gn(3,1)*v(205)+gn(3,2)*v(210)+gn(3,3)*v(214)
  v(233)=gn(2,1)*v(205)+gn(2,2)*v(210)+gn(2,3)*v(214)
  v(217)=v(215)+(v(216)*v(216))+v(221)
  v(219)=v(197)*v(204)+v(213)*v(778)
  v(240)=gn(3,1)*v(210)+gn(3,2)*v(217)+gn(3,3)*v(219)
  v(234)=gn(2,1)*v(210)+gn(2,2)*v(217)+gn(2,3)*v(219)
  v(222)=(v(218)*v(218))+v(220)+v(221)
  v(236)=gn(3,1)*v(214)+gn(3,2)*v(219)+gn(3,3)*v(222)
  v(626)=se(3)*v(236)+se(5)*v(238)+se(6)*v(240)
  v(625)=se(6)*v(236)+se(4)*v(238)+se(2)*v(240)
  v(624)=se(5)*v(236)+se(1)*v(238)+se(4)*v(240)
  v(232)=gn(2,1)*v(214)+gn(2,2)*v(219)+gn(2,3)*v(222)
  v(241)=se(6)*v(232)+se(4)*v(233)+se(2)*v(234)
  v(239)=se(5)*v(232)+se(1)*v(233)+se(4)*v(234)
  v(237)=se(3)*v(232)+se(5)*v(233)+se(6)*v(234)
  v(223)=gn(1,1)*v(214)+gn(1,2)*v(219)+gn(1,3)*v(222)
  v(224)=gn(1,1)*v(205)+gn(1,2)*v(210)+gn(1,3)*v(214)
  v(225)=gn(1,1)*v(210)+gn(1,2)*v(217)+gn(1,3)*v(219)
  v(229)=se(6)*v(223)+se(4)*v(224)+se(2)*v(225)
  v(228)=se(5)*v(223)+se(1)*v(224)+se(4)*v(225)
  v(227)=se(3)*v(223)+se(5)*v(224)+se(6)*v(225)
  v(316)=v(695)*v(749)
  v(317)=v(696)*v(749)
  v(318)=v(697)*v(749)
  v(319)=v(698)*v(749)
  v(320)=v(699)*v(749)
  v(321)=v(700)*v(749)
  v(322)=v(701)*v(750)
  v(323)=v(702)*v(750)
  v(324)=v(703)*v(750)
  v(325)=v(704)*v(750)
  v(326)=v(705)*v(750)
  v(327)=v(706)*v(750)
  v(328)=v(707)*v(751)
  v(329)=v(708)*v(751)
  v(330)=v(709)*v(751)
  v(331)=v(710)*v(751)
  v(332)=v(711)*v(751)
  v(333)=v(712)*v(751)
  v(334)=(ep(5)*(v(677)+v(689))+ep(6)*v(695)+ep(1)*v(701)+ep(3)*v(701)+ep(4)*v(707))/16d0
  v(335)=(ep(5)*(v(678)+v(690))+ep(6)*v(696)+ep(1)*v(702)+ep(3)*v(702)+ep(4)*v(708))/16d0
  v(336)=(ep(5)*(v(679)+v(691))+ep(6)*v(697)+ep(1)*v(703)+ep(3)*v(703)+ep(4)*v(709))/16d0
  v(337)=(ep(5)*(v(680)+v(692))+ep(6)*v(698)+ep(1)*v(704)+ep(3)*v(704)+ep(4)*v(710))/16d0
  v(338)=(ep(5)*(v(681)+v(693))+ep(6)*v(699)+ep(1)*v(705)+ep(3)*v(705)+ep(4)*v(711))/16d0
  v(339)=(ep(5)*(v(682)+v(694))+ep(6)*v(700)+ep(1)*v(706)+ep(3)*v(706)+ep(4)*v(712))/16d0
  v(340)=v(334)*v(752)+v(701)*v(753)
  v(341)=v(335)*v(752)+v(702)*v(753)
  v(342)=v(336)*v(752)+v(703)*v(753)
  v(343)=v(337)*v(752)+v(704)*v(753)
  v(344)=v(338)*v(752)+v(705)*v(753)
  v(345)=v(339)*v(752)+v(706)*v(753)
  v(346)=(ep(4)*(v(677)+v(683))+ep(1)*v(695)+ep(2)*v(695)+ep(6)*v(701)+ep(5)*v(707))/16d0
  v(347)=(ep(4)*(v(678)+v(684))+ep(1)*v(696)+ep(2)*v(696)+ep(6)*v(702)+ep(5)*v(708))/16d0
  v(348)=(ep(4)*(v(679)+v(685))+ep(1)*v(697)+ep(2)*v(697)+ep(6)*v(703)+ep(5)*v(709))/16d0
  v(349)=(ep(4)*(v(680)+v(686))+ep(1)*v(698)+ep(2)*v(698)+ep(6)*v(704)+ep(5)*v(710))/16d0
  v(350)=(ep(4)*(v(681)+v(687))+ep(1)*v(699)+ep(2)*v(699)+ep(6)*v(705)+ep(5)*v(711))/16d0
  v(351)=(ep(4)*(v(682)+v(688))+ep(1)*v(700)+ep(2)*v(700)+ep(6)*v(706)+ep(5)*v(712))/16d0
  v(352)=v(346)*v(754)+v(695)*v(755)
  v(353)=v(347)*v(754)+v(696)*v(755)
  v(354)=v(348)*v(754)+v(697)*v(755)
  v(355)=v(349)*v(754)+v(698)*v(755)
  v(356)=v(350)*v(754)+v(699)*v(755)
  v(357)=v(351)*v(754)+v(700)*v(755)
  v(358)=v(316)+v(322)+v(677)*v(756)
  v(359)=v(317)+v(323)+v(678)*v(756)
  v(360)=v(318)+v(324)+v(679)*v(756)
  v(361)=v(319)+v(325)+v(680)*v(756)
  v(362)=v(320)+v(326)+v(681)*v(756)
  v(363)=v(321)+v(327)+v(682)*v(756)
  v(364)=(-(ep(1)*v(358))+4d0*(v(340)+v(352)+3d0*v(358))-v(670)*v(677))/24d0
  v(365)=(-(ep(1)*v(359))+4d0*(v(341)+v(353)+3d0*v(359))-v(670)*v(678))/24d0
  v(366)=(-(ep(1)*v(360))+4d0*(v(342)+v(354)+3d0*v(360))-v(670)*v(679))/24d0
  v(367)=(-(ep(1)*v(361))+4d0*(v(343)+v(355)+3d0*v(361))-v(670)*v(680))/24d0
  v(368)=(-(ep(1)*v(362))+4d0*(v(344)+v(356)+3d0*v(362))-v(670)*v(681))/24d0
  v(369)=(-(ep(1)*v(363))+4d0*(v(345)+v(357)+3d0*v(363))-v(670)*v(682))/24d0
  v(370)=(-(ep(6)*v(334))+12d0*v(346)-ep(2)*v(346)-ep(4)*v(358)-v(668)*v(683)-6d0*v(695)-v(187)*v(695)-v(667)*v(707))&
       &/24d0
  v(371)=(-(ep(6)*v(335))+12d0*v(347)-ep(2)*v(347)-ep(4)*v(359)-v(668)*v(684)-6d0*v(696)-v(187)*v(696)-v(667)*v(708))&
       &/24d0
  v(372)=(-(ep(6)*v(336))+12d0*v(348)-ep(2)*v(348)-ep(4)*v(360)-v(668)*v(685)-6d0*v(697)-v(187)*v(697)-v(667)*v(709))&
       &/24d0
  v(373)=(-(ep(6)*v(337))+12d0*v(349)-ep(2)*v(349)-ep(4)*v(361)-v(668)*v(686)-6d0*v(698)-v(187)*v(698)-v(667)*v(710))&
       &/24d0
  v(374)=(-(ep(6)*v(338))+12d0*v(350)-ep(2)*v(350)-ep(4)*v(362)-v(668)*v(687)-6d0*v(699)-v(187)*v(699)-v(667)*v(711))&
       &/24d0
  v(375)=(-(ep(6)*v(339))+12d0*v(351)-ep(2)*v(351)-ep(4)*v(363)-v(668)*v(688)-6d0*v(700)-v(187)*v(700)-v(667)*v(712))&
       &/24d0
  v(376)=v(370)*v(757)
  v(377)=v(371)*v(757)
  v(378)=v(372)*v(757)
  v(379)=v(373)*v(757)
  v(380)=v(374)*v(757)
  v(381)=v(375)*v(757)
  v(382)=(12d0*v(334)-ep(3)*v(334)-ep(6)*v(346)-ep(5)*v(358)-v(667)*v(689)-6d0*v(701)-v(187)*v(701)-v(668)*v(707))/24d0
  v(383)=(12d0*v(335)-ep(3)*v(335)-ep(6)*v(347)-ep(5)*v(359)-v(667)*v(690)-6d0*v(702)-v(187)*v(702)-v(668)*v(708))/24d0
  v(384)=(12d0*v(336)-ep(3)*v(336)-ep(6)*v(348)-ep(5)*v(360)-v(667)*v(691)-6d0*v(703)-v(187)*v(703)-v(668)*v(709))/24d0
  v(385)=(12d0*v(337)-ep(3)*v(337)-ep(6)*v(349)-ep(5)*v(361)-v(667)*v(692)-6d0*v(704)-v(187)*v(704)-v(668)*v(710))/24d0
  v(386)=(12d0*v(338)-ep(3)*v(338)-ep(6)*v(350)-ep(5)*v(362)-v(667)*v(693)-6d0*v(705)-v(187)*v(705)-v(668)*v(711))/24d0
  v(387)=(12d0*v(339)-ep(3)*v(339)-ep(6)*v(351)-ep(5)*v(363)-v(667)*v(694)-6d0*v(706)-v(187)*v(706)-v(668)*v(712))/24d0
  v(388)=v(382)*v(758)
  v(389)=v(383)*v(758)
  v(390)=v(384)*v(758)
  v(391)=v(385)*v(758)
  v(392)=v(386)*v(758)
  v(393)=v(387)*v(758)
  v(394)=v(376)+v(388)+v(364)*v(759)
  v(395)=v(377)+v(389)+v(365)*v(759)
  v(396)=v(378)+v(390)+v(366)*v(759)
  v(397)=v(379)+v(391)+v(367)*v(759)
  v(398)=v(380)+v(392)+v(368)*v(759)
  v(399)=v(381)+v(393)+v(369)*v(759)
  v(400)=v(316)+v(328)+v(683)*v(760)
  v(401)=v(317)+v(329)+v(684)*v(760)
  v(402)=v(318)+v(330)+v(685)*v(760)
  v(403)=v(319)+v(331)+v(686)*v(760)
  v(404)=v(320)+v(332)+v(687)*v(760)
  v(405)=v(321)+v(333)+v(688)*v(760)
  v(406)=(ep(6)*(v(683)+v(689))+ep(5)*v(695)+ep(4)*v(701)+ep(2)*v(707)+ep(3)*v(707))/16d0
  v(407)=(ep(6)*(v(684)+v(690))+ep(5)*v(696)+ep(4)*v(702)+ep(2)*v(708)+ep(3)*v(708))/16d0
  v(408)=(ep(6)*(v(685)+v(691))+ep(5)*v(697)+ep(4)*v(703)+ep(2)*v(709)+ep(3)*v(709))/16d0
  v(409)=(ep(6)*(v(686)+v(692))+ep(5)*v(698)+ep(4)*v(704)+ep(2)*v(710)+ep(3)*v(710))/16d0
  v(410)=(ep(6)*(v(687)+v(693))+ep(5)*v(699)+ep(4)*v(705)+ep(2)*v(711)+ep(3)*v(711))/16d0
  v(411)=(ep(6)*(v(688)+v(694))+ep(5)*v(700)+ep(4)*v(706)+ep(2)*v(712)+ep(3)*v(712))/16d0
  v(412)=v(406)*v(761)+v(707)*v(762)
  v(413)=v(407)*v(761)+v(708)*v(762)
  v(414)=v(408)*v(761)+v(709)*v(762)
  v(415)=v(409)*v(761)+v(710)*v(762)
  v(416)=v(410)*v(761)+v(711)*v(762)
  v(417)=v(411)*v(761)+v(712)*v(762)
  v(418)=(-(ep(2)*v(400))+4d0*(v(352)+3d0*v(400)+v(412))-v(674)*v(683))/24d0
  v(419)=(-(ep(2)*v(401))+4d0*(v(353)+3d0*v(401)+v(413))-v(674)*v(684))/24d0
  v(420)=(-(ep(2)*v(402))+4d0*(v(354)+3d0*v(402)+v(414))-v(674)*v(685))/24d0
  v(421)=(-(ep(2)*v(403))+4d0*(v(355)+3d0*v(403)+v(415))-v(674)*v(686))/24d0
  v(422)=(-(ep(2)*v(404))+4d0*(v(356)+3d0*v(404)+v(416))-v(674)*v(687))/24d0
  v(423)=(-(ep(2)*v(405))+4d0*(v(357)+3d0*v(405)+v(417))-v(674)*v(688))/24d0
  v(424)=(-(ep(5)*v(346))-ep(6)*v(400)+12d0*v(406)-ep(3)*v(406)-v(672)*v(689)-v(668)*v(701)-6d0*v(707)-v(194)*v(707))&
       &/24d0
  v(425)=(-(ep(5)*v(347))-ep(6)*v(401)+12d0*v(407)-ep(3)*v(407)-v(672)*v(690)-v(668)*v(702)-6d0*v(708)-v(194)*v(708))&
       &/24d0
  v(426)=(-(ep(5)*v(348))-ep(6)*v(402)+12d0*v(408)-ep(3)*v(408)-v(672)*v(691)-v(668)*v(703)-6d0*v(709)-v(194)*v(709))&
       &/24d0
  v(427)=(-(ep(5)*v(349))-ep(6)*v(403)+12d0*v(409)-ep(3)*v(409)-v(672)*v(692)-v(668)*v(704)-6d0*v(710)-v(194)*v(710))&
       &/24d0
  v(428)=(-(ep(5)*v(350))-ep(6)*v(404)+12d0*v(410)-ep(3)*v(410)-v(672)*v(693)-v(668)*v(705)-6d0*v(711)-v(194)*v(711))&
       &/24d0
  v(429)=(-(ep(5)*v(351))-ep(6)*v(405)+12d0*v(411)-ep(3)*v(411)-v(672)*v(694)-v(668)*v(706)-6d0*v(712)-v(194)*v(712))&
       &/24d0
  v(430)=v(424)*v(763)
  v(431)=v(425)*v(763)
  v(432)=v(426)*v(763)
  v(433)=v(427)*v(763)
  v(434)=v(428)*v(763)
  v(435)=v(429)*v(763)
  v(436)=v(376)+v(430)+v(418)*v(764)
  v(437)=v(377)+v(431)+v(419)*v(764)
  v(438)=v(378)+v(432)+v(420)*v(764)
  v(439)=v(379)+v(433)+v(421)*v(764)
  v(440)=v(380)+v(434)+v(422)*v(764)
  v(441)=v(381)+v(435)+v(423)*v(764)
  v(442)=v(198)*v(382)+v(191)*(v(364)+v(418))+v(192)*v(424)+v(370)*v(765)
  v(443)=v(198)*v(383)+v(191)*(v(365)+v(419))+v(192)*v(425)+v(371)*v(765)
  v(444)=v(198)*v(384)+v(191)*(v(366)+v(420))+v(192)*v(426)+v(372)*v(765)
  v(445)=v(198)*v(385)+v(191)*(v(367)+v(421))+v(192)*v(427)+v(373)*v(765)
  v(446)=v(198)*v(386)+v(191)*(v(368)+v(422))+v(192)*v(428)+v(374)*v(765)
  v(447)=v(198)*v(387)+v(191)*(v(369)+v(423))+v(192)*v(429)+v(375)*v(765)
  v(448)=v(442)*v(766)
  v(449)=v(443)*v(766)
  v(450)=v(444)*v(766)
  v(451)=v(445)*v(766)
  v(452)=v(446)*v(766)
  v(453)=v(447)*v(766)
  v(454)=v(322)+v(328)+v(689)*v(767)
  v(455)=v(323)+v(329)+v(690)*v(767)
  v(456)=v(324)+v(330)+v(691)*v(767)
  v(457)=v(325)+v(331)+v(692)*v(767)
  v(458)=v(326)+v(332)+v(693)*v(767)
  v(459)=v(327)+v(333)+v(694)*v(767)
  v(460)=(-(ep(3)*v(454))+4d0*(v(340)+v(412)+3d0*v(454))-v(689)*v(768))/24d0
  v(461)=(-(ep(3)*v(455))+4d0*(v(341)+v(413)+3d0*v(455))-v(690)*v(768))/24d0
  v(462)=(-(ep(3)*v(456))+4d0*(v(342)+v(414)+3d0*v(456))-v(691)*v(768))/24d0
  v(463)=(-(ep(3)*v(457))+4d0*(v(343)+v(415)+3d0*v(457))-v(692)*v(768))/24d0
  v(464)=(-(ep(3)*v(458))+4d0*(v(344)+v(416)+3d0*v(458))-v(693)*v(768))/24d0
  v(465)=(-(ep(3)*v(459))+4d0*(v(345)+v(417)+3d0*v(459))-v(694)*v(768))/24d0
  v(466)=v(388)+v(430)+v(460)*v(769)
  v(467)=v(389)+v(431)+v(461)*v(769)
  v(468)=v(390)+v(432)+v(462)*v(769)
  v(469)=v(391)+v(433)+v(463)*v(769)
  v(470)=v(392)+v(434)+v(464)*v(769)
  v(471)=v(393)+v(435)+v(465)*v(769)
  v(472)=v(192)*v(370)+v(191)*v(382)+v(198)*(v(418)+v(460))+v(424)*v(770)
  v(473)=v(192)*v(371)+v(191)*v(383)+v(198)*(v(419)+v(461))+v(425)*v(770)
  v(474)=v(192)*v(372)+v(191)*v(384)+v(198)*(v(420)+v(462))+v(426)*v(770)
  v(475)=v(192)*v(373)+v(191)*v(385)+v(198)*(v(421)+v(463))+v(427)*v(770)
  v(476)=v(192)*v(374)+v(191)*v(386)+v(198)*(v(422)+v(464))+v(428)*v(770)
  v(477)=v(192)*v(375)+v(191)*v(387)+v(198)*(v(423)+v(465))+v(429)*v(770)
  v(478)=v(472)*v(771)
  v(479)=v(473)*v(771)
  v(480)=v(474)*v(771)
  v(481)=v(475)*v(771)
  v(482)=v(476)*v(771)
  v(483)=v(477)*v(771)
  v(484)=v(198)*v(370)+v(191)*v(424)+v(192)*(v(364)+v(460))+v(382)*v(772)
  v(485)=v(198)*v(371)+v(191)*v(425)+v(192)*(v(365)+v(461))+v(383)*v(772)
  v(486)=v(198)*v(372)+v(191)*v(426)+v(192)*(v(366)+v(462))+v(384)*v(772)
  v(487)=v(198)*v(373)+v(191)*v(427)+v(192)*(v(367)+v(463))+v(385)*v(772)
  v(488)=v(198)*v(374)+v(191)*v(428)+v(192)*(v(368)+v(464))+v(386)*v(772)
  v(489)=v(198)*v(375)+v(191)*v(429)+v(192)*(v(369)+v(465))+v(387)*v(772)
  v(490)=v(484)*v(773)
  v(491)=v(485)*v(773)
  v(492)=v(486)*v(773)
  v(493)=v(487)*v(773)
  v(494)=v(488)*v(773)
  v(495)=v(489)*v(773)
  v(496)=v(448)+v(490)+v(394)*v(774)
  v(497)=v(449)+v(491)+v(395)*v(774)
  v(498)=v(450)+v(492)+v(396)*v(774)
  v(499)=v(451)+v(493)+v(397)*v(774)
  v(500)=v(452)+v(494)+v(398)*v(774)
  v(501)=v(453)+v(495)+v(399)*v(774)
  v(502)=v(197)*(v(394)+v(436))+v(204)*v(472)+v(213)*v(484)+v(442)*v(775)
  v(503)=v(197)*(v(395)+v(437))+v(204)*v(473)+v(213)*v(485)+v(443)*v(775)
  v(504)=v(197)*(v(396)+v(438))+v(204)*v(474)+v(213)*v(486)+v(444)*v(775)
  v(505)=v(197)*(v(397)+v(439))+v(204)*v(475)+v(213)*v(487)+v(445)*v(775)
  v(506)=v(197)*(v(398)+v(440))+v(204)*v(476)+v(213)*v(488)+v(446)*v(775)
  v(507)=v(197)*(v(399)+v(441))+v(204)*v(477)+v(213)*v(489)+v(447)*v(775)
  v(508)=v(213)*v(442)+v(204)*(v(394)+v(466))+v(197)*v(472)+v(484)*v(776)
  v(509)=v(213)*v(443)+v(204)*(v(395)+v(467))+v(197)*v(473)+v(485)*v(776)
  v(510)=v(213)*v(444)+v(204)*(v(396)+v(468))+v(197)*v(474)+v(486)*v(776)
  v(511)=v(213)*v(445)+v(204)*(v(397)+v(469))+v(197)*v(475)+v(487)*v(776)
  v(512)=v(213)*v(446)+v(204)*(v(398)+v(470))+v(197)*v(476)+v(488)*v(776)
  v(513)=v(213)*v(447)+v(204)*(v(399)+v(471))+v(197)*v(477)+v(489)*v(776)
  v(514)=gn(3,1)*v(496)+gn(3,2)*v(502)+gn(3,3)*v(508)
  v(515)=gn(3,1)*v(497)+gn(3,2)*v(503)+gn(3,3)*v(509)
  v(516)=gn(3,1)*v(498)+gn(3,2)*v(504)+gn(3,3)*v(510)
  v(517)=gn(3,1)*v(499)+gn(3,2)*v(505)+gn(3,3)*v(511)
  v(518)=gn(3,1)*v(500)+gn(3,2)*v(506)+gn(3,3)*v(512)
  v(519)=gn(3,1)*v(501)+gn(3,2)*v(507)+gn(3,3)*v(513)
  v(520)=gn(2,1)*v(496)+gn(2,2)*v(502)+gn(2,3)*v(508)
  v(521)=gn(2,1)*v(497)+gn(2,2)*v(503)+gn(2,3)*v(509)
  v(522)=gn(2,1)*v(498)+gn(2,2)*v(504)+gn(2,3)*v(510)
  v(523)=gn(2,1)*v(499)+gn(2,2)*v(505)+gn(2,3)*v(511)
  v(524)=gn(2,1)*v(500)+gn(2,2)*v(506)+gn(2,3)*v(512)
  v(525)=gn(2,1)*v(501)+gn(2,2)*v(507)+gn(2,3)*v(513)
  v(526)=v(448)+v(478)+v(436)*v(777)
  v(527)=v(449)+v(479)+v(437)*v(777)
  v(528)=v(450)+v(480)+v(438)*v(777)
  v(529)=v(451)+v(481)+v(439)*v(777)
  v(530)=v(452)+v(482)+v(440)*v(777)
  v(531)=v(453)+v(483)+v(441)*v(777)
  v(532)=v(204)*v(442)+v(213)*(v(436)+v(466))+v(197)*v(484)+v(472)*v(778)
  v(533)=v(204)*v(443)+v(213)*(v(437)+v(467))+v(197)*v(485)+v(473)*v(778)
  v(534)=v(204)*v(444)+v(213)*(v(438)+v(468))+v(197)*v(486)+v(474)*v(778)
  v(535)=v(204)*v(445)+v(213)*(v(439)+v(469))+v(197)*v(487)+v(475)*v(778)
  v(536)=v(204)*v(446)+v(213)*(v(440)+v(470))+v(197)*v(488)+v(476)*v(778)
  v(537)=v(204)*v(447)+v(213)*(v(441)+v(471))+v(197)*v(489)+v(477)*v(778)
  v(538)=gn(3,1)*v(502)+gn(3,2)*v(526)+gn(3,3)*v(532)
  v(539)=gn(3,1)*v(503)+gn(3,2)*v(527)+gn(3,3)*v(533)
  v(540)=gn(3,1)*v(504)+gn(3,2)*v(528)+gn(3,3)*v(534)
  v(541)=gn(3,1)*v(505)+gn(3,2)*v(529)+gn(3,3)*v(535)
  v(542)=gn(3,1)*v(506)+gn(3,2)*v(530)+gn(3,3)*v(536)
  v(543)=gn(3,1)*v(507)+gn(3,2)*v(531)+gn(3,3)*v(537)
  v(544)=gn(2,1)*v(502)+gn(2,2)*v(526)+gn(2,3)*v(532)
  v(545)=gn(2,1)*v(503)+gn(2,2)*v(527)+gn(2,3)*v(533)
  v(546)=gn(2,1)*v(504)+gn(2,2)*v(528)+gn(2,3)*v(534)
  v(547)=gn(2,1)*v(505)+gn(2,2)*v(529)+gn(2,3)*v(535)
  v(548)=gn(2,1)*v(506)+gn(2,2)*v(530)+gn(2,3)*v(536)
  v(549)=gn(2,1)*v(507)+gn(2,2)*v(531)+gn(2,3)*v(537)
  v(550)=v(478)+v(490)+v(466)*v(779)
  v(551)=v(479)+v(491)+v(467)*v(779)
  v(552)=v(480)+v(492)+v(468)*v(779)
  v(553)=v(481)+v(493)+v(469)*v(779)
  v(554)=v(482)+v(494)+v(470)*v(779)
  v(555)=v(483)+v(495)+v(471)*v(779)
  v(556)=gn(3,1)*v(508)+gn(3,2)*v(532)+gn(3,3)*v(550)
  v(557)=gn(3,1)*v(509)+gn(3,2)*v(533)+gn(3,3)*v(551)
  v(558)=gn(3,1)*v(510)+gn(3,2)*v(534)+gn(3,3)*v(552)
  v(559)=gn(3,1)*v(511)+gn(3,2)*v(535)+gn(3,3)*v(553)
  v(560)=gn(3,1)*v(512)+gn(3,2)*v(536)+gn(3,3)*v(554)
  v(561)=gn(3,1)*v(513)+gn(3,2)*v(537)+gn(3,3)*v(555)
  v(562)=gn(2,1)*v(508)+gn(2,2)*v(532)+gn(2,3)*v(550)
  v(563)=gn(2,1)*v(509)+gn(2,2)*v(533)+gn(2,3)*v(551)
  v(564)=gn(2,1)*v(510)+gn(2,2)*v(534)+gn(2,3)*v(552)
  v(565)=gn(2,1)*v(511)+gn(2,2)*v(535)+gn(2,3)*v(553)
  v(566)=gn(2,1)*v(512)+gn(2,2)*v(536)+gn(2,3)*v(554)
  v(567)=gn(2,1)*v(513)+gn(2,2)*v(537)+gn(2,3)*v(555)
  v(568)=se(4)*v(520)+se(2)*v(544)+se(6)*v(562)+v(234)*v(719)+v(233)*v(731)+v(232)*v(743)
  v(569)=se(4)*v(521)+se(2)*v(545)+se(6)*v(563)+v(234)*v(720)+v(233)*v(732)+v(232)*v(744)
  v(570)=se(4)*v(522)+se(2)*v(546)+se(6)*v(564)+v(234)*v(721)+v(233)*v(733)+v(232)*v(745)
  v(571)=se(4)*v(523)+se(2)*v(547)+se(6)*v(565)+v(234)*v(722)+v(233)*v(734)+v(232)*v(746)
  v(572)=se(4)*v(524)+se(2)*v(548)+se(6)*v(566)+v(234)*v(723)+v(233)*v(735)+v(232)*v(747)
  v(573)=se(4)*v(525)+se(2)*v(549)+se(6)*v(567)+v(234)*v(724)+v(233)*v(736)+v(232)*v(748)
  v(574)=se(1)*v(520)+se(4)*v(544)+se(5)*v(562)+v(233)*v(713)+v(234)*v(731)+v(232)*v(737)
  v(575)=se(1)*v(521)+se(4)*v(545)+se(5)*v(563)+v(233)*v(714)+v(234)*v(732)+v(232)*v(738)
  v(576)=se(1)*v(522)+se(4)*v(546)+se(5)*v(564)+v(233)*v(715)+v(234)*v(733)+v(232)*v(739)
  v(577)=se(1)*v(523)+se(4)*v(547)+se(5)*v(565)+v(233)*v(716)+v(234)*v(734)+v(232)*v(740)
  v(578)=se(1)*v(524)+se(4)*v(548)+se(5)*v(566)+v(233)*v(717)+v(234)*v(735)+v(232)*v(741)
  v(579)=se(1)*v(525)+se(4)*v(549)+se(5)*v(567)+v(233)*v(718)+v(234)*v(736)+v(232)*v(742)
  v(580)=se(5)*v(520)+se(6)*v(544)+se(3)*v(562)+v(232)*v(725)+v(233)*v(737)+v(234)*v(743)
  v(581)=se(5)*v(521)+se(6)*v(545)+se(3)*v(563)+v(232)*v(726)+v(233)*v(738)+v(234)*v(744)
  v(582)=se(5)*v(522)+se(6)*v(546)+se(3)*v(564)+v(232)*v(727)+v(233)*v(739)+v(234)*v(745)
  v(583)=se(5)*v(523)+se(6)*v(547)+se(3)*v(565)+v(232)*v(728)+v(233)*v(740)+v(234)*v(746)
  v(584)=se(5)*v(524)+se(6)*v(548)+se(3)*v(566)+v(232)*v(729)+v(233)*v(741)+v(234)*v(747)
  v(585)=se(5)*v(525)+se(6)*v(549)+se(3)*v(567)+v(232)*v(730)+v(233)*v(742)+v(234)*v(748)
  v(586)=gn(1,1)*v(508)+gn(1,2)*v(532)+gn(1,3)*v(550)
  v(587)=gn(1,1)*v(509)+gn(1,2)*v(533)+gn(1,3)*v(551)
  v(588)=gn(1,1)*v(510)+gn(1,2)*v(534)+gn(1,3)*v(552)
  v(589)=gn(1,1)*v(511)+gn(1,2)*v(535)+gn(1,3)*v(553)
  v(590)=gn(1,1)*v(512)+gn(1,2)*v(536)+gn(1,3)*v(554)
  v(591)=gn(1,1)*v(513)+gn(1,2)*v(537)+gn(1,3)*v(555)
  v(592)=gn(1,1)*v(496)+gn(1,2)*v(502)+gn(1,3)*v(508)
  v(593)=gn(1,1)*v(497)+gn(1,2)*v(503)+gn(1,3)*v(509)
  v(594)=gn(1,1)*v(498)+gn(1,2)*v(504)+gn(1,3)*v(510)
  v(595)=gn(1,1)*v(499)+gn(1,2)*v(505)+gn(1,3)*v(511)
  v(596)=gn(1,1)*v(500)+gn(1,2)*v(506)+gn(1,3)*v(512)
  v(597)=gn(1,1)*v(501)+gn(1,2)*v(507)+gn(1,3)*v(513)
  v(598)=gn(1,1)*v(502)+gn(1,2)*v(526)+gn(1,3)*v(532)
  v(599)=gn(1,1)*v(503)+gn(1,2)*v(527)+gn(1,3)*v(533)
  v(600)=gn(1,1)*v(504)+gn(1,2)*v(528)+gn(1,3)*v(534)
  v(601)=gn(1,1)*v(505)+gn(1,2)*v(529)+gn(1,3)*v(535)
  v(602)=gn(1,1)*v(506)+gn(1,2)*v(530)+gn(1,3)*v(536)
  v(603)=gn(1,1)*v(507)+gn(1,2)*v(531)+gn(1,3)*v(537)
  v(604)=se(6)*v(586)+se(4)*v(592)+se(2)*v(598)+v(225)*v(719)+v(224)*v(731)+v(223)*v(743)
  v(605)=se(6)*v(587)+se(4)*v(593)+se(2)*v(599)+v(225)*v(720)+v(224)*v(732)+v(223)*v(744)
  v(606)=se(6)*v(588)+se(4)*v(594)+se(2)*v(600)+v(225)*v(721)+v(224)*v(733)+v(223)*v(745)
  v(607)=se(6)*v(589)+se(4)*v(595)+se(2)*v(601)+v(225)*v(722)+v(224)*v(734)+v(223)*v(746)
  v(608)=se(6)*v(590)+se(4)*v(596)+se(2)*v(602)+v(225)*v(723)+v(224)*v(735)+v(223)*v(747)
  v(609)=se(6)*v(591)+se(4)*v(597)+se(2)*v(603)+v(225)*v(724)+v(224)*v(736)+v(223)*v(748)
  v(610)=se(5)*v(586)+se(1)*v(592)+se(4)*v(598)+v(224)*v(713)+v(225)*v(731)+v(223)*v(737)
  v(611)=se(5)*v(587)+se(1)*v(593)+se(4)*v(599)+v(224)*v(714)+v(225)*v(732)+v(223)*v(738)
  v(612)=se(5)*v(588)+se(1)*v(594)+se(4)*v(600)+v(224)*v(715)+v(225)*v(733)+v(223)*v(739)
  v(613)=se(5)*v(589)+se(1)*v(595)+se(4)*v(601)+v(224)*v(716)+v(225)*v(734)+v(223)*v(740)
  v(614)=se(5)*v(590)+se(1)*v(596)+se(4)*v(602)+v(224)*v(717)+v(225)*v(735)+v(223)*v(741)
  v(615)=se(5)*v(591)+se(1)*v(597)+se(4)*v(603)+v(224)*v(718)+v(225)*v(736)+v(223)*v(742)
  v(616)=se(3)*v(586)+se(5)*v(592)+se(6)*v(598)+v(223)*v(725)+v(224)*v(737)+v(225)*v(743)
  v(617)=se(3)*v(587)+se(5)*v(593)+se(6)*v(599)+v(223)*v(726)+v(224)*v(738)+v(225)*v(744)
  v(618)=se(3)*v(588)+se(5)*v(594)+se(6)*v(600)+v(223)*v(727)+v(224)*v(739)+v(225)*v(745)
  v(619)=se(3)*v(589)+se(5)*v(595)+se(6)*v(601)+v(223)*v(728)+v(224)*v(740)+v(225)*v(746)
  v(620)=se(3)*v(590)+se(5)*v(596)+se(6)*v(602)+v(223)*v(729)+v(224)*v(741)+v(225)*v(747)
  v(621)=se(3)*v(591)+se(5)*v(597)+se(6)*v(603)+v(223)*v(730)+v(224)*v(742)+v(225)*v(748)
  s(1)=v(223)*v(227)+v(224)*v(228)+v(225)*v(229)
  s(2)=v(232)*v(237)+v(233)*v(239)+v(234)*v(241)
  s(3)=v(238)*v(624)+v(240)*v(625)+v(236)*v(626)
  s(4)=v(227)*v(232)+v(228)*v(233)+v(229)*v(234)
  s(5)=v(227)*v(236)+v(228)*v(238)+v(229)*v(240)
  s(6)=v(236)*v(237)+v(238)*v(239)+v(240)*v(241)
  gn1(1,1)=v(224)
  gn1(1,2)=v(225)
  gn1(1,3)=v(223)
  gn1(2,1)=v(233)
  gn1(2,2)=v(234)
  gn1(2,3)=v(232)
  gn1(3,1)=v(238)
  gn1(3,2)=v(240)
  gn1(3,3)=v(236)
  dsdc(1,1)=v(227)*v(586)+v(228)*v(592)+v(229)*v(598)+v(225)*v(604)+v(224)*v(610)+v(223)*v(616)
  dsdc(1,2)=v(227)*v(587)+v(228)*v(593)+v(229)*v(599)+v(225)*v(605)+v(224)*v(611)+v(223)*v(617)
  dsdc(1,3)=v(227)*v(588)+v(228)*v(594)+v(229)*v(600)+v(225)*v(606)+v(224)*v(612)+v(223)*v(618)
  dsdc(1,4)=v(227)*v(589)+v(228)*v(595)+v(229)*v(601)+v(225)*v(607)+v(224)*v(613)+v(223)*v(619)
  dsdc(1,5)=v(227)*v(590)+v(228)*v(596)+v(229)*v(602)+v(225)*v(608)+v(224)*v(614)+v(223)*v(620)
  dsdc(1,6)=v(227)*v(591)+v(228)*v(597)+v(229)*v(603)+v(225)*v(609)+v(224)*v(615)+v(223)*v(621)
  dsdc(2,1)=v(239)*v(520)+v(241)*v(544)+v(237)*v(562)+v(234)*v(568)+v(233)*v(574)+v(232)*v(580)
  dsdc(2,2)=v(239)*v(521)+v(241)*v(545)+v(237)*v(563)+v(234)*v(569)+v(233)*v(575)+v(232)*v(581)
  dsdc(2,3)=v(239)*v(522)+v(241)*v(546)+v(237)*v(564)+v(234)*v(570)+v(233)*v(576)+v(232)*v(582)
  dsdc(2,4)=v(239)*v(523)+v(241)*v(547)+v(237)*v(565)+v(234)*v(571)+v(233)*v(577)+v(232)*v(583)
  dsdc(2,5)=v(239)*v(524)+v(241)*v(548)+v(237)*v(566)+v(234)*v(572)+v(233)*v(578)+v(232)*v(584)
  dsdc(2,6)=v(239)*v(525)+v(241)*v(549)+v(237)*v(567)+v(234)*v(573)+v(233)*v(579)+v(232)*v(585)
  dsdc(3,1)=v(514)*v(624)+v(538)*v(625)+v(556)*v(626)+v(238)*(se(1)*v(514)+se(4)*v(538)+se(5)*v(556)+v(238)*v(713)+v(240&
       &)*v(731)+v(236)*v(737))+v(240)*(se(4)*v(514)+se(2)*v(538)+se(6)*v(556)+v(240)*v(719)+v(238)*v(731)+v(236)*v(743))+v(236&
       &)*(se(5)*v(514)+se(6)*v(538)+se(3)*v(556)+v(236)*v(725)+v(238)*v(737)+v(240)*v(743))
  dsdc(3,2)=v(515)*v(624)+v(539)*v(625)+v(557)*v(626)+v(238)*(se(1)*v(515)+se(4)*v(539)+se(5)*v(557)+v(238)*v(714)+v(240&
       &)*v(732)+v(236)*v(738))+v(240)*(se(4)*v(515)+se(2)*v(539)+se(6)*v(557)+v(240)*v(720)+v(238)*v(732)+v(236)*v(744))+v(236&
       &)*(se(5)*v(515)+se(6)*v(539)+se(3)*v(557)+v(236)*v(726)+v(238)*v(738)+v(240)*v(744))
  dsdc(3,3)=v(516)*v(624)+v(540)*v(625)+v(558)*v(626)+v(238)*(se(1)*v(516)+se(4)*v(540)+se(5)*v(558)+v(238)*v(715)+v(240&
       &)*v(733)+v(236)*v(739))+v(240)*(se(4)*v(516)+se(2)*v(540)+se(6)*v(558)+v(240)*v(721)+v(238)*v(733)+v(236)*v(745))+v(236&
       &)*(se(5)*v(516)+se(6)*v(540)+se(3)*v(558)+v(236)*v(727)+v(238)*v(739)+v(240)*v(745))
  dsdc(3,4)=v(517)*v(624)+v(541)*v(625)+v(559)*v(626)+v(238)*(se(1)*v(517)+se(4)*v(541)+se(5)*v(559)+v(238)*v(716)+v(240&
       &)*v(734)+v(236)*v(740))+v(240)*(se(4)*v(517)+se(2)*v(541)+se(6)*v(559)+v(240)*v(722)+v(238)*v(734)+v(236)*v(746))+v(236&
       &)*(se(5)*v(517)+se(6)*v(541)+se(3)*v(559)+v(236)*v(728)+v(238)*v(740)+v(240)*v(746))
  dsdc(3,5)=v(518)*v(624)+v(542)*v(625)+v(560)*v(626)+v(238)*(se(1)*v(518)+se(4)*v(542)+se(5)*v(560)+v(238)*v(717)+v(240&
       &)*v(735)+v(236)*v(741))+v(240)*(se(4)*v(518)+se(2)*v(542)+se(6)*v(560)+v(240)*v(723)+v(238)*v(735)+v(236)*v(747))+v(236&
       &)*(se(5)*v(518)+se(6)*v(542)+se(3)*v(560)+v(236)*v(729)+v(238)*v(741)+v(240)*v(747))
  dsdc(3,6)=v(519)*v(624)+v(543)*v(625)+v(561)*v(626)+v(238)*(se(1)*v(519)+se(4)*v(543)+se(5)*v(561)+v(238)*v(718)+v(240&
       &)*v(736)+v(236)*v(742))+v(240)*(se(4)*v(519)+se(2)*v(543)+se(6)*v(561)+v(240)*v(724)+v(238)*v(736)+v(236)*v(748))+v(236&
       &)*(se(5)*v(519)+se(6)*v(543)+se(3)*v(561)+v(236)*v(730)+v(238)*v(742)+v(240)*v(748))
  dsdc(4,1)=v(228)*v(520)+v(229)*v(544)+v(227)*v(562)+v(234)*v(604)+v(233)*v(610)+v(232)*v(616)
  dsdc(4,2)=v(228)*v(521)+v(229)*v(545)+v(227)*v(563)+v(234)*v(605)+v(233)*v(611)+v(232)*v(617)
  dsdc(4,3)=v(228)*v(522)+v(229)*v(546)+v(227)*v(564)+v(234)*v(606)+v(233)*v(612)+v(232)*v(618)
  dsdc(4,4)=v(228)*v(523)+v(229)*v(547)+v(227)*v(565)+v(234)*v(607)+v(233)*v(613)+v(232)*v(619)
  dsdc(4,5)=v(228)*v(524)+v(229)*v(548)+v(227)*v(566)+v(234)*v(608)+v(233)*v(614)+v(232)*v(620)
  dsdc(4,6)=v(228)*v(525)+v(229)*v(549)+v(227)*v(567)+v(234)*v(609)+v(233)*v(615)+v(232)*v(621)
  dsdc(5,1)=v(228)*v(514)+v(229)*v(538)+v(227)*v(556)+v(240)*v(604)+v(238)*v(610)+v(236)*v(616)
  dsdc(5,2)=v(228)*v(515)+v(229)*v(539)+v(227)*v(557)+v(240)*v(605)+v(238)*v(611)+v(236)*v(617)
  dsdc(5,3)=v(228)*v(516)+v(229)*v(540)+v(227)*v(558)+v(240)*v(606)+v(238)*v(612)+v(236)*v(618)
  dsdc(5,4)=v(228)*v(517)+v(229)*v(541)+v(227)*v(559)+v(240)*v(607)+v(238)*v(613)+v(236)*v(619)
  dsdc(5,5)=v(228)*v(518)+v(229)*v(542)+v(227)*v(560)+v(240)*v(608)+v(238)*v(614)+v(236)*v(620)
  dsdc(5,6)=v(228)*v(519)+v(229)*v(543)+v(227)*v(561)+v(240)*v(609)+v(238)*v(615)+v(236)*v(621)
  dsdc(6,1)=v(239)*v(514)+v(241)*v(538)+v(237)*v(556)+v(240)*v(568)+v(238)*v(574)+v(236)*v(580)
  dsdc(6,2)=v(239)*v(515)+v(241)*v(539)+v(237)*v(557)+v(240)*v(569)+v(238)*v(575)+v(236)*v(581)
  dsdc(6,3)=v(239)*v(516)+v(241)*v(540)+v(237)*v(558)+v(240)*v(570)+v(238)*v(576)+v(236)*v(582)
  dsdc(6,4)=v(239)*v(517)+v(241)*v(541)+v(237)*v(559)+v(240)*v(571)+v(238)*v(577)+v(236)*v(583)
  dsdc(6,5)=v(239)*v(518)+v(241)*v(542)+v(237)*v(560)+v(240)*v(572)+v(238)*v(578)+v(236)*v(584)
  dsdc(6,6)=v(239)*v(519)+v(241)*v(543)+v(237)*v(561)+v(240)*v(573)+v(238)*v(579)+v(236)*v(585)
END SUBROUTINE tangentmodulus


SUBROUTINE dettotal(QOLD,QNEW,SE,DSEDCE,DCEDC,S,DSDC)
  IMPLICIT NONE
  DOUBLE PRECISION V(325),QOLD(3,3),QNEW(3,3),SE(6),DSEDCE(6,6),DCEDC(6,6),S(6),DSDC(6,6)
  V(320)=QOLD(3,3)**2
  V(319)=QOLD(3,2)**2
  V(317)=2D0*QOLD(3,1)*QOLD(3,3)
  V(315)=QOLD(3,1)**2
  V(314)=DSEDCE(3,1)*V(320)
  V(313)=DSEDCE(5,1)*V(317)
  V(312)=DSEDCE(2,1)*V(319)
  V(311)=DSEDCE(1,1)*V(315)
  V(310)=2D0*QOLD(3,2)
  V(318)=QOLD(3,1)*V(310)
  V(316)=QOLD(3,3)*V(310)
  V(309)=DSEDCE(5,6)*QOLD(1,1)+DSEDCE(6,6)*QOLD(1,2)+DSEDCE(3,6)*QOLD(1,3)
  V(308)=DSEDCE(5,5)*QOLD(1,1)+DSEDCE(6,5)*QOLD(1,2)+DSEDCE(3,5)*QOLD(1,3)
  V(307)=DSEDCE(5,4)*QOLD(1,1)+DSEDCE(6,4)*QOLD(1,2)+DSEDCE(3,4)*QOLD(1,3)
  V(306)=DSEDCE(5,3)*QOLD(1,1)+DSEDCE(6,3)*QOLD(1,2)+DSEDCE(3,3)*QOLD(1,3)
  V(305)=DSEDCE(5,2)*QOLD(1,1)+DSEDCE(6,2)*QOLD(1,2)+DSEDCE(3,2)*QOLD(1,3)
  V(304)=DSEDCE(5,1)*QOLD(1,1)+DSEDCE(6,1)*QOLD(1,2)+DSEDCE(3,1)*QOLD(1,3)
  V(303)=DSEDCE(4,6)*QOLD(1,1)+DSEDCE(2,6)*QOLD(1,2)+DSEDCE(6,6)*QOLD(1,3)
  V(302)=DSEDCE(4,5)*QOLD(1,1)+DSEDCE(2,5)*QOLD(1,2)+DSEDCE(6,5)*QOLD(1,3)
  V(301)=DSEDCE(4,4)*QOLD(1,1)+DSEDCE(2,4)*QOLD(1,2)+DSEDCE(6,4)*QOLD(1,3)
  V(300)=DSEDCE(4,3)*QOLD(1,1)+DSEDCE(2,3)*QOLD(1,2)+DSEDCE(6,3)*QOLD(1,3)
  V(299)=DSEDCE(4,2)*QOLD(1,1)+DSEDCE(2,2)*QOLD(1,2)+DSEDCE(6,2)*QOLD(1,3)
  V(298)=DSEDCE(4,1)*QOLD(1,1)+DSEDCE(2,1)*QOLD(1,2)+DSEDCE(6,1)*QOLD(1,3)
  V(297)=DSEDCE(5,6)*QOLD(2,1)+DSEDCE(6,6)*QOLD(2,2)+DSEDCE(3,6)*QOLD(2,3)
  V(296)=DSEDCE(5,5)*QOLD(2,1)+DSEDCE(6,5)*QOLD(2,2)+DSEDCE(3,5)*QOLD(2,3)
  V(295)=DSEDCE(5,4)*QOLD(2,1)+DSEDCE(6,4)*QOLD(2,2)+DSEDCE(3,4)*QOLD(2,3)
  V(294)=DSEDCE(5,3)*QOLD(2,1)+DSEDCE(6,3)*QOLD(2,2)+DSEDCE(3,3)*QOLD(2,3)
  V(293)=DSEDCE(5,2)*QOLD(2,1)+DSEDCE(6,2)*QOLD(2,2)+DSEDCE(3,2)*QOLD(2,3)
  V(292)=DSEDCE(5,1)*QOLD(2,1)+DSEDCE(6,1)*QOLD(2,2)+DSEDCE(3,1)*QOLD(2,3)
  V(291)=DSEDCE(4,6)*QOLD(2,1)+DSEDCE(2,6)*QOLD(2,2)+DSEDCE(6,6)*QOLD(2,3)
  V(290)=DSEDCE(4,5)*QOLD(2,1)+DSEDCE(2,5)*QOLD(2,2)+DSEDCE(6,5)*QOLD(2,3)
  V(289)=DSEDCE(4,4)*QOLD(2,1)+DSEDCE(2,4)*QOLD(2,2)+DSEDCE(6,4)*QOLD(2,3)
  V(288)=DSEDCE(4,3)*QOLD(2,1)+DSEDCE(2,3)*QOLD(2,2)+DSEDCE(6,3)*QOLD(2,3)
  V(287)=DSEDCE(4,2)*QOLD(2,1)+DSEDCE(2,2)*QOLD(2,2)+DSEDCE(6,2)*QOLD(2,3)
  V(286)=DSEDCE(4,1)*QOLD(2,1)+DSEDCE(2,1)*QOLD(2,2)+DSEDCE(6,1)*QOLD(2,3)
  V(285)=DSEDCE(1,6)*QOLD(1,1)+DSEDCE(4,6)*QOLD(1,2)+DSEDCE(5,6)*QOLD(1,3)
  V(284)=DSEDCE(1,5)*QOLD(1,1)+DSEDCE(4,5)*QOLD(1,2)+DSEDCE(5,5)*QOLD(1,3)
  V(283)=DSEDCE(1,4)*QOLD(1,1)+DSEDCE(4,4)*QOLD(1,2)+DSEDCE(5,4)*QOLD(1,3)
  V(282)=DSEDCE(1,3)*QOLD(1,1)+DSEDCE(4,3)*QOLD(1,2)+DSEDCE(5,3)*QOLD(1,3)
  V(281)=DSEDCE(1,2)*QOLD(1,1)+DSEDCE(4,2)*QOLD(1,2)+DSEDCE(5,2)*QOLD(1,3)
  V(280)=DSEDCE(1,1)*QOLD(1,1)+DSEDCE(4,1)*QOLD(1,2)+DSEDCE(5,1)*QOLD(1,3)
  V(279)=DSEDCE(1,6)*QOLD(2,1)+DSEDCE(4,6)*QOLD(2,2)+DSEDCE(5,6)*QOLD(2,3)
  V(278)=DSEDCE(1,5)*QOLD(2,1)+DSEDCE(4,5)*QOLD(2,2)+DSEDCE(5,5)*QOLD(2,3)
  V(277)=DSEDCE(1,4)*QOLD(2,1)+DSEDCE(4,4)*QOLD(2,2)+DSEDCE(5,4)*QOLD(2,3)
  V(276)=DSEDCE(1,3)*QOLD(2,1)+DSEDCE(4,3)*QOLD(2,2)+DSEDCE(5,3)*QOLD(2,3)
  V(275)=DSEDCE(1,2)*QOLD(2,1)+DSEDCE(4,2)*QOLD(2,2)+DSEDCE(5,2)*QOLD(2,3)
  V(274)=DSEDCE(1,1)*QOLD(2,1)+DSEDCE(4,1)*QOLD(2,2)+DSEDCE(5,1)*QOLD(2,3)
  V(273)=QOLD(1,3)*SE(3)+QOLD(1,1)*SE(5)+QOLD(1,2)*SE(6)
  V(272)=QOLD(1,2)*SE(2)+QOLD(1,1)*SE(4)+QOLD(1,3)*SE(6)
  V(271)=QOLD(2,3)*SE(3)+QOLD(2,1)*SE(5)+QOLD(2,2)*SE(6)
  V(270)=QOLD(2,2)*SE(2)+QOLD(2,1)*SE(4)+QOLD(2,3)*SE(6)
  V(269)=QOLD(1,1)*SE(1)+QOLD(1,2)*SE(4)+QOLD(1,3)*SE(5)
  V(268)=QOLD(2,1)*SE(1)+QOLD(2,2)*SE(4)+QOLD(2,3)*SE(5)
  V(193)=QOLD(1,1)*V(280)+QOLD(1,2)*V(298)+QOLD(1,3)*V(304)
  V(194)=QOLD(2,1)*V(274)+QOLD(2,2)*V(286)+QOLD(2,3)*V(292)
  V(195)=V(311)+V(312)+V(313)+V(314)+DSEDCE(6,1)*V(316)+DSEDCE(4,1)*V(318)
  V(196)=QOLD(2,1)*V(280)+QOLD(2,2)*V(298)+QOLD(2,3)*V(304)
  V(197)=QOLD(3,1)*V(280)+QOLD(3,2)*V(298)+QOLD(3,3)*V(304)
  V(198)=QOLD(3,1)*V(274)+QOLD(3,2)*V(286)+QOLD(3,3)*V(292)
  V(199)=QOLD(1,1)*V(281)+QOLD(1,2)*V(299)+QOLD(1,3)*V(305)
  V(200)=QOLD(2,1)*V(275)+QOLD(2,2)*V(287)+QOLD(2,3)*V(293)
  V(201)=DSEDCE(1,2)*V(315)+DSEDCE(6,2)*V(316)+DSEDCE(5,2)*V(317)+DSEDCE(4,2)*V(318)+DSEDCE(2,2)*V(319)+DSEDCE(3,2)*V(320&
       &)
  V(202)=QOLD(2,1)*V(281)+QOLD(2,2)*V(299)+QOLD(2,3)*V(305)
  V(203)=QOLD(3,1)*V(281)+QOLD(3,2)*V(299)+QOLD(3,3)*V(305)
  V(204)=QOLD(3,1)*V(275)+QOLD(3,2)*V(287)+QOLD(3,3)*V(293)
  V(205)=QOLD(1,1)*V(282)+QOLD(1,2)*V(300)+QOLD(1,3)*V(306)
  V(206)=QOLD(2,1)*V(276)+QOLD(2,2)*V(288)+QOLD(2,3)*V(294)
  V(207)=DSEDCE(1,3)*V(315)+DSEDCE(6,3)*V(316)+DSEDCE(5,3)*V(317)+DSEDCE(4,3)*V(318)+DSEDCE(2,3)*V(319)+DSEDCE(3,3)*V(320&
       &)
  V(208)=QOLD(2,1)*V(282)+QOLD(2,2)*V(300)+QOLD(2,3)*V(306)
  V(209)=QOLD(3,1)*V(282)+QOLD(3,2)*V(300)+QOLD(3,3)*V(306)
  V(210)=QOLD(3,1)*V(276)+QOLD(3,2)*V(288)+QOLD(3,3)*V(294)
  V(211)=QOLD(1,1)*V(283)+QOLD(1,2)*V(301)+QOLD(1,3)*V(307)
  V(212)=QOLD(2,1)*V(277)+QOLD(2,2)*V(289)+QOLD(2,3)*V(295)
  V(213)=DSEDCE(1,4)*V(315)+DSEDCE(6,4)*V(316)+DSEDCE(5,4)*V(317)+DSEDCE(4,4)*V(318)+DSEDCE(2,4)*V(319)+DSEDCE(3,4)*V(320&
       &)
  V(214)=QOLD(2,1)*V(283)+QOLD(2,2)*V(301)+QOLD(2,3)*V(307)
  V(215)=QOLD(3,1)*V(283)+QOLD(3,2)*V(301)+QOLD(3,3)*V(307)
  V(216)=QOLD(3,1)*V(277)+QOLD(3,2)*V(289)+QOLD(3,3)*V(295)
  V(217)=QOLD(1,1)*V(284)+QOLD(1,2)*V(302)+QOLD(1,3)*V(308)
  V(218)=QOLD(2,1)*V(278)+QOLD(2,2)*V(290)+QOLD(2,3)*V(296)
  V(219)=DSEDCE(1,5)*V(315)+DSEDCE(6,5)*V(316)+DSEDCE(5,5)*V(317)+DSEDCE(4,5)*V(318)+DSEDCE(2,5)*V(319)+DSEDCE(3,5)*V(320&
       &)
  V(220)=QOLD(2,1)*V(284)+QOLD(2,2)*V(302)+QOLD(2,3)*V(308)
  V(221)=QOLD(3,1)*V(284)+QOLD(3,2)*V(302)+QOLD(3,3)*V(308)
  V(222)=QOLD(3,1)*V(278)+QOLD(3,2)*V(290)+QOLD(3,3)*V(296)
  V(223)=QOLD(1,1)*V(285)+QOLD(1,2)*V(303)+QOLD(1,3)*V(309)
  V(224)=QOLD(2,1)*V(279)+QOLD(2,2)*V(291)+QOLD(2,3)*V(297)
  V(225)=DSEDCE(1,6)*V(315)+DSEDCE(6,6)*V(316)+DSEDCE(5,6)*V(317)+DSEDCE(4,6)*V(318)+DSEDCE(2,6)*V(319)+DSEDCE(3,6)*V(320&
       &)
  V(226)=QOLD(2,1)*V(285)+QOLD(2,2)*V(303)+QOLD(2,3)*V(309)
  V(227)=QOLD(3,1)*V(285)+QOLD(3,2)*V(303)+QOLD(3,3)*V(309)
  V(228)=QOLD(3,1)*V(279)+QOLD(3,2)*V(291)+QOLD(3,3)*V(297)
  QNEW(1,1)=QOLD(1,1)
  QNEW(1,2)=QOLD(1,2)
  QNEW(1,3)=QOLD(1,3)
  QNEW(2,1)=QOLD(2,1)
  QNEW(2,2)=QOLD(2,2)
  QNEW(2,3)=QOLD(2,3)
  QNEW(3,1)=QOLD(3,1)
  QNEW(3,2)=QOLD(3,2)
  QNEW(3,3)=QOLD(3,3)
  S(1)=QOLD(1,1)*V(269)+QOLD(1,2)*V(272)+QOLD(1,3)*V(273)
  S(2)=QOLD(2,1)*V(268)+QOLD(2,2)*V(270)+QOLD(2,3)*V(271)
  S(3)=QOLD(3,1)*(QOLD(3,1)*SE(1)+QOLD(3,2)*SE(4)+QOLD(3,3)*SE(5))+QOLD(3,3)*(QOLD(3,3)*SE(3)+QOLD(3,1)*SE(5)+QOLD(3,2&
       &)*SE(6))+QOLD(3,2)*(QOLD(3,2)*SE(2)+QOLD(3,1)*SE(4)+QOLD(3,3)*SE(6))
  S(4)=QOLD(2,1)*V(269)+QOLD(2,2)*V(272)+QOLD(2,3)*V(273)
  S(5)=QOLD(3,1)*V(269)+QOLD(3,2)*V(272)+QOLD(3,3)*V(273)
  S(6)=QOLD(3,1)*V(268)+QOLD(3,2)*V(270)+QOLD(3,3)*V(271)
  DSDC(1,1)=DCEDC(1,1)*V(193)+DCEDC(2,1)*V(199)+DCEDC(3,1)*V(205)+DCEDC(4,1)*V(211)+DCEDC(5,1)*V(217)+DCEDC(6,1)*V(223)
  DSDC(1,2)=DCEDC(1,2)*V(193)+DCEDC(2,2)*V(199)+DCEDC(3,2)*V(205)+DCEDC(4,2)*V(211)+DCEDC(5,2)*V(217)+DCEDC(6,2)*V(223)
  DSDC(1,3)=DCEDC(1,3)*V(193)+DCEDC(2,3)*V(199)+DCEDC(3,3)*V(205)+DCEDC(4,3)*V(211)+DCEDC(5,3)*V(217)+DCEDC(6,3)*V(223)
  DSDC(1,4)=DCEDC(1,4)*V(193)+DCEDC(2,4)*V(199)+DCEDC(3,4)*V(205)+DCEDC(4,4)*V(211)+DCEDC(5,4)*V(217)+DCEDC(6,4)*V(223)
  DSDC(1,5)=DCEDC(1,5)*V(193)+DCEDC(2,5)*V(199)+DCEDC(3,5)*V(205)+DCEDC(4,5)*V(211)+DCEDC(5,5)*V(217)+DCEDC(6,5)*V(223)
  DSDC(1,6)=DCEDC(1,6)*V(193)+DCEDC(2,6)*V(199)+DCEDC(3,6)*V(205)+DCEDC(4,6)*V(211)+DCEDC(5,6)*V(217)+DCEDC(6,6)*V(223)
  DSDC(2,1)=DCEDC(1,1)*V(194)+DCEDC(2,1)*V(200)+DCEDC(3,1)*V(206)+DCEDC(4,1)*V(212)+DCEDC(5,1)*V(218)+DCEDC(6,1)*V(224)
  DSDC(2,2)=DCEDC(1,2)*V(194)+DCEDC(2,2)*V(200)+DCEDC(3,2)*V(206)+DCEDC(4,2)*V(212)+DCEDC(5,2)*V(218)+DCEDC(6,2)*V(224)
  DSDC(2,3)=DCEDC(1,3)*V(194)+DCEDC(2,3)*V(200)+DCEDC(3,3)*V(206)+DCEDC(4,3)*V(212)+DCEDC(5,3)*V(218)+DCEDC(6,3)*V(224)
  DSDC(2,4)=DCEDC(1,4)*V(194)+DCEDC(2,4)*V(200)+DCEDC(3,4)*V(206)+DCEDC(4,4)*V(212)+DCEDC(5,4)*V(218)+DCEDC(6,4)*V(224)
  DSDC(2,5)=DCEDC(1,5)*V(194)+DCEDC(2,5)*V(200)+DCEDC(3,5)*V(206)+DCEDC(4,5)*V(212)+DCEDC(5,5)*V(218)+DCEDC(6,5)*V(224)
  DSDC(2,6)=DCEDC(1,6)*V(194)+DCEDC(2,6)*V(200)+DCEDC(3,6)*V(206)+DCEDC(4,6)*V(212)+DCEDC(5,6)*V(218)+DCEDC(6,6)*V(224)
  DSDC(3,1)=DCEDC(1,1)*V(195)+DCEDC(2,1)*V(201)+DCEDC(3,1)*V(207)+DCEDC(4,1)*V(213)+DCEDC(5,1)*V(219)+DCEDC(6,1)*V(225)
  DSDC(3,2)=DCEDC(1,2)*V(195)+DCEDC(2,2)*V(201)+DCEDC(3,2)*V(207)+DCEDC(4,2)*V(213)+DCEDC(5,2)*V(219)+DCEDC(6,2)*V(225)
  DSDC(3,3)=DCEDC(1,3)*V(195)+DCEDC(2,3)*V(201)+DCEDC(3,3)*V(207)+DCEDC(4,3)*V(213)+DCEDC(5,3)*V(219)+DCEDC(6,3)*V(225)
  DSDC(3,4)=DCEDC(1,4)*V(195)+DCEDC(2,4)*V(201)+DCEDC(3,4)*V(207)+DCEDC(4,4)*V(213)+DCEDC(5,4)*V(219)+DCEDC(6,4)*V(225)
  DSDC(3,5)=DCEDC(1,5)*V(195)+DCEDC(2,5)*V(201)+DCEDC(3,5)*V(207)+DCEDC(4,5)*V(213)+DCEDC(5,5)*V(219)+DCEDC(6,5)*V(225)
  DSDC(3,6)=DCEDC(1,6)*V(195)+DCEDC(2,6)*V(201)+DCEDC(3,6)*V(207)+DCEDC(4,6)*V(213)+DCEDC(5,6)*V(219)+DCEDC(6,6)*V(225)
  DSDC(4,1)=DCEDC(1,1)*V(196)+DCEDC(2,1)*V(202)+DCEDC(3,1)*V(208)+DCEDC(4,1)*V(214)+DCEDC(5,1)*V(220)+DCEDC(6,1)*V(226)
  DSDC(4,2)=DCEDC(1,2)*V(196)+DCEDC(2,2)*V(202)+DCEDC(3,2)*V(208)+DCEDC(4,2)*V(214)+DCEDC(5,2)*V(220)+DCEDC(6,2)*V(226)
  DSDC(4,3)=DCEDC(1,3)*V(196)+DCEDC(2,3)*V(202)+DCEDC(3,3)*V(208)+DCEDC(4,3)*V(214)+DCEDC(5,3)*V(220)+DCEDC(6,3)*V(226)
  DSDC(4,4)=DCEDC(1,4)*V(196)+DCEDC(2,4)*V(202)+DCEDC(3,4)*V(208)+DCEDC(4,4)*V(214)+DCEDC(5,4)*V(220)+DCEDC(6,4)*V(226)
  DSDC(4,5)=DCEDC(1,5)*V(196)+DCEDC(2,5)*V(202)+DCEDC(3,5)*V(208)+DCEDC(4,5)*V(214)+DCEDC(5,5)*V(220)+DCEDC(6,5)*V(226)
  DSDC(4,6)=DCEDC(1,6)*V(196)+DCEDC(2,6)*V(202)+DCEDC(3,6)*V(208)+DCEDC(4,6)*V(214)+DCEDC(5,6)*V(220)+DCEDC(6,6)*V(226)
  DSDC(5,1)=DCEDC(1,1)*V(197)+DCEDC(2,1)*V(203)+DCEDC(3,1)*V(209)+DCEDC(4,1)*V(215)+DCEDC(5,1)*V(221)+DCEDC(6,1)*V(227)
  DSDC(5,2)=DCEDC(1,2)*V(197)+DCEDC(2,2)*V(203)+DCEDC(3,2)*V(209)+DCEDC(4,2)*V(215)+DCEDC(5,2)*V(221)+DCEDC(6,2)*V(227)
  DSDC(5,3)=DCEDC(1,3)*V(197)+DCEDC(2,3)*V(203)+DCEDC(3,3)*V(209)+DCEDC(4,3)*V(215)+DCEDC(5,3)*V(221)+DCEDC(6,3)*V(227)
  DSDC(5,4)=DCEDC(1,4)*V(197)+DCEDC(2,4)*V(203)+DCEDC(3,4)*V(209)+DCEDC(4,4)*V(215)+DCEDC(5,4)*V(221)+DCEDC(6,4)*V(227)
  DSDC(5,5)=DCEDC(1,5)*V(197)+DCEDC(2,5)*V(203)+DCEDC(3,5)*V(209)+DCEDC(4,5)*V(215)+DCEDC(5,5)*V(221)+DCEDC(6,5)*V(227)
  DSDC(5,6)=DCEDC(1,6)*V(197)+DCEDC(2,6)*V(203)+DCEDC(3,6)*V(209)+DCEDC(4,6)*V(215)+DCEDC(5,6)*V(221)+DCEDC(6,6)*V(227)
  DSDC(6,1)=DCEDC(1,1)*V(198)+DCEDC(2,1)*V(204)+DCEDC(3,1)*V(210)+DCEDC(4,1)*V(216)+DCEDC(5,1)*V(222)+DCEDC(6,1)*V(228)
  DSDC(6,2)=DCEDC(1,2)*V(198)+DCEDC(2,2)*V(204)+DCEDC(3,2)*V(210)+DCEDC(4,2)*V(216)+DCEDC(5,2)*V(222)+DCEDC(6,2)*V(228)
  DSDC(6,3)=DCEDC(1,3)*V(198)+DCEDC(2,3)*V(204)+DCEDC(3,3)*V(210)+DCEDC(4,3)*V(216)+DCEDC(5,3)*V(222)+DCEDC(6,3)*V(228)
  DSDC(6,4)=DCEDC(1,4)*V(198)+DCEDC(2,4)*V(204)+DCEDC(3,4)*V(210)+DCEDC(4,4)*V(216)+DCEDC(5,4)*V(222)+DCEDC(6,4)*V(228)
  DSDC(6,5)=DCEDC(1,5)*V(198)+DCEDC(2,5)*V(204)+DCEDC(3,5)*V(210)+DCEDC(4,5)*V(216)+DCEDC(5,5)*V(222)+DCEDC(6,5)*V(228)
  DSDC(6,6)=DCEDC(1,6)*V(198)+DCEDC(2,6)*V(204)+DCEDC(3,6)*V(210)+DCEDC(4,6)*V(216)+DCEDC(5,6)*V(222)+DCEDC(6,6)*V(228)
END SUBROUTINE dettotal



subroutine rootg(xmin, xmax, fu, root, status)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  implicit none
  integer, parameter :: dp = kind(1.0d0)

!> Input
  real(dp), intent(in)  :: xmin, xmax
  external :: fu   ! keep implicit interface: call fu(x, f, df)

!> Output
  real(dp), intent(out) :: root
  integer,  intent(out) :: status

! ---- Status codes
  integer, parameter :: S_OK=0, S_TOL=1, S_MAXIT=2, S_NOBRACK=3, S_BADIN=4, S_NONFIN=5, S_DEGEN=6

! ---- Tunables
  real(dp), parameter :: ftol = 1.0e-8_dp      ! function tolerance
  real(dp), parameter :: atol = 1.0e-13_dp     ! absolute interval tol (~sqrt(eps) for DP scales)
  real(dp), parameter :: rtol = 1.0e-7_dp      ! relative interval tol
  integer,  parameter :: maxit = 100
  real(dp), parameter :: progress_guard = 0.9_dp    ! max fractional step of current width
  real(dp), parameter :: bracket_narrowing_factor = 0.98_dp
! scale-aware denominator guards (heuristics)
  real(dp), parameter :: k_newton = 10.0_dp
  real(dp), parameter :: k_secant = 10.0_dp

! ---- Locals
  real(dp) :: a, b, fa, fb, dfa, dfb
  real(dp) :: x, fx, dfx
  real(dp) :: xprev, fxprev
  real(dp) :: xbest, fbest
  real(dp) :: w, wprev, scale, maxstep, stepscale
  real(dp) :: xtrial, ftrial, dtrial, denom
  logical  :: have_bracket, same_sign, have_newton, have_secant
  logical  :: prefer_bisect
  integer  :: iter
  real(dp) :: mid

! ---- Validate inputs
  if (.not. ieee_is_finite(xmin) .or. .not. ieee_is_finite(xmax)) then
     root = 0.0_dp
     status = S_BADIN
     return
  end if

! Order interval
  a = xmin; b = xmax
  if (b < a) then
     x = a; a = b; b = x
  end if

  w = b - a
  if (w <= atol) then
     root = 0.5_dp*(a+b)
     status = S_DEGEN
     return
  end if

! Endpoint evaluations
  call fu(a, fa, dfa)
  if (.not. ieee_is_finite(fa)) then
     mid = 0.5_dp*(a+b)
     call fu(mid, ftrial, dtrial)
     if (.not. ieee_is_finite(ftrial)) then
        root = mid
        status = S_NONFIN
        return
     else
! replace bad endpoint with finite midpoint value
        a   = mid; fa = ftrial; dfa = dtrial
     end if
  end if

  if (abs(fa) <= ftol) then
     root = a
     status = S_OK
     return
  end if

  call fu(b, fb, dfb)
  if (.not. ieee_is_finite(fb)) then
     mid = 0.5_dp*(a+b)
     call fu(mid, ftrial, dtrial)
     if (.not. ieee_is_finite(ftrial)) then
        root = mid
        status = S_NONFIN
        return
     else
        b   = mid; fb = ftrial; dfb = dtrial
     end if
  end if

  if (abs(fb) <= ftol) then
     root = b
     status = S_OK
     return
  end if

! Check initial bracket
  same_sign   = (fa > 0.0_dp .and. fb > 0.0_dp) .or. (fa < 0.0_dp .and. fb < 0.0_dp)
  have_bracket = .not. same_sign
  if (.not. have_bracket) then
! no sign change and neither endpoint already satisfied ftol
     if (abs(fa) <= abs(fb)) then
        root = a
     else
        root = b
     end if
     status = S_NOBRACK
     return
  end if

! Initial guess: guarded secant from endpoints if safe; else midpoint
  denom = fb - fa
  if (ieee_is_finite(denom) .and. abs(denom) >= k_secant*max(abs(fa),abs(fb))*epsilon(1.0_dp)) then
     x = (a*fb - b*fa) / denom
! enforce bracket and step cap from the midpoint
     mid = 0.5_dp*(a+b)
     x   = max(a, min(b, x))
     if (abs(x - mid) > progress_guard*(b - a)) x = mid + sign(progress_guard*(b - a), x - mid)
     if (x <= a .or. x >= b) x = mid
  else
     x = 0.5_dp*(a+b)
  end if

  call fu(x, fx, dfx)
  if (.not. ieee_is_finite(fx)) then
! fallback to midpoint
     x = 0.5_dp*(a+b)
     call fu(x, fx, dfx)
     if (.not. ieee_is_finite(fx)) then
        root = x
        status = S_NONFIN
        return
     end if
  end if

! Best-so-far tracker
  if (abs(fa) <= abs(fx) .and. abs(fa) <= abs(fb)) then
     xbest = a; fbest = fa
  else if (abs(fb) <= abs(fx)) then
     xbest = b; fbest = fb
  else
     xbest = x; fbest = fx
  end if

  xprev  = x
  fxprev = fx
  wprev  = b - a

  do iter = 1, maxit

! Convergence checks
     scale  = max( 1.0_dp, abs(a), abs(b), abs(x), abs(xbest) )
     if (abs(fx) <= ftol) then
        root = x
        status = S_OK
        return
     end if
     w = b - a
     if (w <= atol + rtol*scale) then
        root = xbest
        status = S_TOL
        return
     end if

! Force bisection if bracket hasn't narrowed enough
     prefer_bisect = ( w >= max( wprev*bracket_narrowing_factor, wprev - (atol + rtol*scale) ) )

! Step caps & scale-aware guards
     maxstep   = progress_guard * w
     stepscale = max( w, spacing( max( abs(a), abs(b), abs(x) ) ) )

     have_newton = ieee_is_finite(dfx) .and. (abs(dfx) >= k_newton * abs(fx) / stepscale)
     have_secant = ieee_is_finite(fxprev) .and. (x /= xprev)
     if (have_secant) then
        denom = fx - fxprev
        have_secant = ieee_is_finite(denom) .and. &
             (abs(denom) >= k_secant*max(abs(fx),abs(fxprev))*epsilon(1.0_dp))
     end if

! Trial selection (Newton → Secant → Bisection)
     if (.not. prefer_bisect .and. have_newton) then
        xtrial = x - fx/dfx
        if (xtrial <= a .or. xtrial >= b .or. abs(xtrial - x) > maxstep) then
! reject Newton if outside or too large; try secant next
           if (.not. prefer_bisect .and. have_secant) then
              xtrial = x - fx * (x - xprev) / (fx - fxprev)
              if (xtrial <= a .or. xtrial >= b .or. abs(xtrial - x) > maxstep) then
                 xtrial = 0.5_dp*(a+b)
              end if
           else
              xtrial = 0.5_dp*(a+b)
           end if
        end if
     else if (.not. prefer_bisect .and. have_secant) then
        xtrial = x - fx * (x - xprev) / (fx - fxprev)
        if (xtrial <= a .or. xtrial >= b .or. abs(xtrial - x) > maxstep) then
           xtrial = 0.5_dp*(a+b)
        end if
     else
        xtrial = 0.5_dp*(a+b)
     end if

! Evaluate trial; on non-finite, try midpoint once; else give up
     call fu(xtrial, ftrial, dtrial)
     if (.not. ieee_is_finite(ftrial)) then
        xtrial = 0.5_dp*(a+b)
        call fu(xtrial, ftrial, dtrial)
        if (.not. ieee_is_finite(ftrial)) then
           root = xbest
           status = S_NONFIN
           return
        end if
     end if

! Update best-so-far
     if (abs(ftrial) < abs(fbest)) then
        xbest = xtrial
        fbest = ftrial
     end if

! Update bracket (maintain sign change)
     if ((fa > 0.0_dp .and. ftrial > 0.0_dp) .or. (fa < 0.0_dp .and. ftrial < 0.0_dp)) then
        a  = xtrial
        fa = ftrial
        dfa = dtrial
     else
        b  = xtrial
        fb = ftrial
        dfb = dtrial
     end if

! Shift iterates
     xprev  = x
     fxprev = fx
     x   = xtrial
     fx  = ftrial
     dfx = dtrial
     wprev = w

  end do

! Max iterations: return best-so-far
  root = xbest
  status = S_MAXIT
end subroutine rootg


