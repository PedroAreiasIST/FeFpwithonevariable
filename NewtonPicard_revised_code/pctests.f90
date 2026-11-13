!=============================================================================
! pctests -- the tension and shear tests of the polycarbonate (Section 8.1.1), with Algorithm 1 and adaptive
! substepping (np_update_sub).
!   tension: E11 prescribed, E22 = E33 found so that S22 = S33 = 0 (the shear components vanish): from the undeformed
!            state to E11 = -0.2 in Ns increments; again from the undeformed state to E11 = 1.2 in Ns increments,
!            followed by unloading to zero stress in Ns equal steps of M11;
!   shear:   the engineering shear 2 E12 prescribed, every other strain component zero: to -0.6, then from the
!            undeformed state to 0.6 and unloading to M12 = 0 in Ns steps.
! usage: pctests <tension|shear> <Ns> <Nn> <out>        (the paper: Ns = 10 with Nn = 1, 5, 10, 0; Ns = 1000 with Nn = 0)
!   Nn > 0: exactly Nn Picard updates per increment (eps_theta = 0), substepped where the scalar iteration fails
!   Nn = 0: the Picard iteration converged to eps_theta = 1e-10 within 50 updates, with substepping
! Output per increment: strain (E11 or 2 E12), Mandel stress (M11 or M12) [Pa], kappa, number of parts, and
! (columns 5-7) the single increment from the committed state to the new strain with the settings of the run, as
! np_update_alg1 returns it before any substepping: its trial ratio seq^tr/Y, its status (0 a stress with the
! Picard iteration converged or, for Nn > 0, after its Nn updates; 5 a stress with the Picard iteration not converged;
! 2, 3, 4 no stress) and the ratio of the last two changes of the flow vector (0 when there are fewer than two updates).
!=============================================================================
program pctests
   use npicard_plasticity
   implicit none
   type(np_material) :: mat
   type(np_state) :: st
   character(len=16) :: test
   character(len=256) :: outf, arg
   integer :: ns, nn, info, uout, i, nsub, is, it
   real(dp) :: eold(6), cp(6), kap, e(6), s(6), cpn(6), kn, x, sm, target, ea, eb, fa, fb, ec, fc, smax
   real(dp) :: dratio, dfac
   integer :: dinfo
   logical :: uni
   call get_command_argument(1, test)
   call get_command_argument(2, arg); read (arg, *) ns
   call get_command_argument(3, arg); read (arg, *) nn
   call get_command_argument(4, outf)
   uni = trim(test) == 'tension'
   mat%mu = 848.0d6; mat%kmod = 3.25d9; mat%iyield = NP_YLD_GHORBELF; mat%rsig = 1.387d0; mat%st0 = 62.0d6
   mat%y0 = 62.0d6; mat%ysat = 78.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   if (nn > 0) then
      mat%nn = nn; mat%eps_theta = 0.0d0
   else
      mat%nn = 50; mat%eps_theta = 1.0d-10
   end if
   call np_check_material(mat, info, arg)
   if (info /= NP_OK) then; write (*, *) 'material: ', trim(arg); stop 1; end if
   open (newunit=uout, file=trim(outf), status='replace')
   write (uout, '(A,A,A,I0,A,I0)') '# pctests ', trim(test), ' Ns ', ns, ' Nn ', nn
   write (uout, '(A)') '# strain  Mandel_stress[Pa]  kappa  parts  trial_ratio  status_single  factor_single'
   ! ---- negative branch
   call reset()
   do i = 1, ns
      call advance(merge(-0.2d0, -0.6d0, uni)*real(i, dp)/real(ns, dp), sm, nsub)
      write (uout, '(3ES18.9,I6,F10.4,I4,F10.4)') strain(), sm, kap, nsub, dratio, dinfo, dfac
   end do
   write (uout, '(A)') ' '
   ! ---- positive branch
   call reset()
   do i = 1, ns
      call advance(merge(1.2d0, 0.6d0, uni)*real(i, dp)/real(ns, dp), sm, nsub)
      write (uout, '(3ES18.9,I6,F10.4,I4,F10.4)') strain(), sm, kap, nsub, dratio, dinfo, dfac
   end do
   ! ---- unloading to zero stress in ns equal steps of the Mandel stress (secant on the prescribed strain)
   smax = sm
   do is = 1, ns
      target = smax*real(ns - is, dp)/real(ns, dp)
      ea = strain(); fa = sm - target
      eb = ea - fa/merge(2.2d9, 8.48d8, uni)          ! elastic estimate
      call trial(eb, fb, nsub); fb = fb - target
      do it = 1, 60
         if (abs(fb) <= 1.0d-9*abs(smax)) exit
         ec = eb - fb*(eb - ea)/(fb - fa)
         call trial(ec, fc, nsub); fc = fc - target
         ea = eb; fa = fb; eb = ec; fb = fc
      end do
      call advance(eb, sm, nsub)
      write (uout, '(3ES18.9,I6,F10.4,I4,F10.4)') strain(), sm, kap, nsub, dratio, dinfo, dfac
   end do
   close (uout)
   write (*, '(A,A,A,I0,A,I0,A)') ' pctests ', trim(test), ' Ns ', ns, ' Nn ', nn, ': done'
contains
   subroutine reset()
      eold = 0.0d0; cp = [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]; kap = 0.0d0
   end subroutine reset
   real(dp) function strain()
      strain = merge(eold(1), eold(4), uni)
   end function strain
   ! one increment from the committed state to the strain en; lev: substep level kept non-decreasing by the caller
   subroutine update(en, lev, so, cpo, ko, mm, nsb, inf)
      real(dp), intent(in) :: en(6)
      integer, intent(inout) :: lev
      real(dp), intent(out) :: so(6), cpo(6), ko, mm
      integer, intent(out) :: nsb, inf
      real(dp) :: d(6, 6), psi
      call np_update_sub(mat, eold, en, cp, kap, so, d, cpo, ko, psi, inf, nsb, 10, st=st, minlev=lev)
      lev = max(lev, nint(log(real(nsb, dp))/log(2.0d0)))
      mm = merge(st%m6(1), st%m6(4)/sqrt(2.0d0), uni)
   end subroutine update
   ! the state at the prescribed strain ev (not committed): Mandel stress and parts; solved lateral strain in xlat
   subroutine solve(ev, so, cpo, ko, mm, nsb, en)
      real(dp), intent(in) :: ev
      real(dp), intent(out) :: so(6), cpo(6), ko, mm, en(6)
      integer, intent(out) :: nsb
      real(dp) :: xa, xb, ga, gb, xc, gc, dx
      integer :: lev, inf, k
      lev = 0
      if (.not. uni) then
         en = 0.0d0; en(4) = ev
         call update(en, lev, so, cpo, ko, mm, nsb, inf)
         if (inf /= NP_OK) then; write (*, *) 'shear: update failed, status ', inf, ' at ', ev; stop 2; end if
         return
      end if
      ! uniaxial stress: g(x) = S22 with E22 = E33 = x; bracket around the elastic estimate, then regula falsi/bisection
      xa = eold(2) - 0.35d0*(ev - eold(1))
      en = [ev, xa, xa, 0.0d0, 0.0d0, 0.0d0]
      call update(en, lev, so, cpo, ko, mm, nsb, inf); ga = so(2)
      if (inf /= NP_OK) then; write (*, *) 'tension: update failed, status ', inf, ' at ', ev, xa; stop 2; end if
      dx = 1.0d-3*max(abs(ev - eold(1)), 1.0d-4)
      xb = xa - sign(dx, ga); gb = ga
      do k = 1, 200
         en = [ev, xb, xb, 0.0d0, 0.0d0, 0.0d0]
         call update(en, lev, so, cpo, ko, mm, nsb, inf); gb = so(2)
         if (inf /= NP_OK) then; write (*, *) 'tension: update failed, status ', inf, ' at ', ev, xb; stop 2; end if
         if (ga*gb <= 0.0d0) exit
         xa = xb; ga = gb; dx = 1.6d0*dx; xb = xb - sign(dx, gb)
      end do
      if (ga*gb > 0.0d0) then; write (*, *) 'tension: no bracket for the lateral strain at ', ev; stop 2; end if
      do k = 1, 200
         if (abs(gb) <= 1.0d-10*mat%y0 .or. abs(xb - xa) <= 1.0d-15) exit
         xc = xb - gb*(xb - xa)/(gb - ga)
         if (mod(k, 3) == 0 .or. (xc - xa)*(xc - xb) >= 0.0d0) xc = 0.5d0*(xa + xb)
         en = [ev, xc, xc, 0.0d0, 0.0d0, 0.0d0]
         call update(en, lev, so, cpo, ko, mm, nsb, inf); gc = so(2)
         if (inf /= NP_OK) then; write (*, *) 'tension: update failed, status ', inf, ' at ', ev, xc; stop 2; end if
         if (ga*gc <= 0.0d0) then
            xb = xc; gb = gc
         else
            xa = xc; ga = gc
         end if
         if (abs(ga) < abs(gb)) then          ! keep the better end in b
            xc = xa; xa = xb; xb = xc; gc = ga; ga = gb; gb = gc
         end if
      end do
      en = [ev, xb, xb, 0.0d0, 0.0d0, 0.0d0]
      call update(en, lev, so, cpo, ko, mm, nsb, inf)
   end subroutine solve
   subroutine trial(ev, mm, nsb)
      real(dp), intent(in) :: ev
      real(dp), intent(out) :: mm
      integer, intent(out) :: nsb
      real(dp) :: so(6), cpo(6), ko, en(6)
      call solve(ev, so, cpo, ko, mm, nsb, en)
   end subroutine trial
   subroutine advance(ev, mm, nsb)
      real(dp), intent(in) :: ev
      real(dp), intent(out) :: mm
      integer, intent(out) :: nsb
      real(dp) :: so(6), cpo(6), ko, en(6)
      real(dp) :: s1(6), d1(6, 6), c1(6), k1, p1
      type(np_state) :: t1
      call solve(ev, so, cpo, ko, mm, nsb, en)
      call np_update_alg1(mat, en, cp, kap, s1, d1, c1, k1, p1, dinfo, st=t1)
      dratio = t1%ratio; dfac = 0.0d0
      if (t1%plastic .and. t1%its >= 2) dfac = t1%dn(t1%its)/max(t1%dn(t1%its - 1), tiny(1.0d0))
      eold = en; cp = cpo; kap = ko
   end subroutine advance
end program pctests
