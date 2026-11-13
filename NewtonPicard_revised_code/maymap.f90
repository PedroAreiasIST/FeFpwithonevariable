!=============================================================================
! maymap -- the iso-error maps of the figures of Section 8.2 (Figures 9 to 13) and their numbers in the text.
! The proportional path from the start state to every point of the grid is divided into N equal increments
! (N = 1, 2, 5, 10, 20, 50), each with one update of the flow direction (np_update_alg1 with N_n = 1,
! eps_theta = 0; an increment whose scalar iteration fails is substepped, np_update_sub); the reference is the same
! update in nref increments; the error is taken over the normal components of the stress. The converged update
! (eps_theta 1e-10, cap 50) in one increment is compared with the converged update in nref increments.
! Materials, maps, plane stress iteration and routines as isomap.f90.
!
! usage: maymap <law> <map> <n> <amax> <nref> <out> [base]
!   base: the staggered-hardening Newton baseline (np_update_base, converged in every increment) in place of the
!         update with one Picard update per increment
!   The maps of the paper: pc virgin 101 0.25 90, pla virgin 101 0.25 90, pla virgin 101 0.25 90 base,
!   barlat point 101 5 90, iso point 101 5 90, barlat virgin 101 0.1 90.
! Columns: 1 a, 2 b, 3 E11, 4 E22, 5 E33 (reference), 6 kappa (reference), 7 failed increments of the reference
!   (>= 1000: no state), 8-13 err_N1 err_N2 err_N5 err_N10 err_N20 err_N50 (N equal increments), 14-19 their failed
!   increments, 20 err_conv (one increment of the converged update, substepped where it fails, against the
!   converged update in nref increments), 21 nsub_conv (its parts), 22 its status, 23 failed increments of its
!   reference, 24 ratio_trial (of the single increment). Error = |S(1:3) - Sref(1:3)|/|Sref(1:3)|, -1: no state.
! The figures show the error fields after a Gaussian filter of radius 4 grid points; the text quotes unsmoothed values.
!=============================================================================
program maymap
   use npicard_plasticity
   use omp_lib
   implicit none
   type(np_material) :: mat
   character(len=16) :: law, map, variant
   character(len=256) :: outf, arg
   integer :: n, nsref, ip, npt, info, uout, i, j
   logical :: pstress, base
   real(dp) :: amax, eb, e33b, a, b, e0(6), de(6), t0
   integer, parameter :: NNS(6) = [1, 2, 5, 10, 20, 50], NCOL = 24
   real(dp), allocatable :: rows(:, :)
   real(dp) :: row(NCOL)
   call get_command_argument(1, law)
   call get_command_argument(2, map)
   call get_command_argument(3, arg); read (arg, *) n
   call get_command_argument(4, arg); read (arg, *) amax
   call get_command_argument(5, arg); read (arg, *) nsref
   call get_command_argument(6, outf)
   variant = ''
   if (command_argument_count() >= 7) call get_command_argument(7, variant)
   base = (trim(variant) == 'base')
   call set_material(trim(law), mat, pstress)
   call np_check_material(mat, info, arg)
   if (info /= NP_OK) then; write (*, *) 'material: ', trim(arg); stop 1; end if
   eb = 0.0d0; e33b = 0.0d0
   if (pstress) call first_yield(eb, e33b)
   if (trim(map) == 'point' .and. .not. pstress) stop 'point maps are plane stress (barlat, iso)'
   npt = n*n
   allocate (rows(NCOL, npt))
   t0 = omp_get_wtime()
!$omp parallel do schedule(dynamic, 1) private(ip, i, j, a, b, e0, de, row)
   do ip = 1, npt
      i = (ip - 1)/n + 1; j = mod(ip - 1, n) + 1
      a = -amax + 2.0d0*amax*real(j - 1, dp)/real(n - 1, dp)
      b = -amax + 2.0d0*amax*real(i - 1, dp)/real(n - 1, dp)
      e0 = 0.0d0; de = 0.0d0
      if (trim(map) == 'point') then
         e0(1) = eb; e0(2) = eb; e0(3) = e33b
         de(1) = a*eb; de(2) = b*eb
      else
         de(1) = a; de(2) = b
      end if
      call one_point(e0, de, row)
      row(1) = a; row(2) = b
      rows(:, ip) = row
   end do
!$omp end parallel do
   open (newunit=uout, file=trim(outf), status='replace')
   write (uout, '(A,A,A,A,A,A,A,I0,A,ES12.4,A,I0,A,ES16.8,A,L1,A,F9.1,A)') '# maymap ', trim(law), ' ', trim(map), ' ', &
      trim(variant), ' n ', n, ' amax ', amax, ' nref ', nsref, ' e_b ', eb, ' planestress ', pstress, ' wall ', &
      omp_get_wtime() - t0, ' s; N equal increments, one update each; reference in nref increments'
   write (uout, '(A)') '# a b E11 E22 E33ref kappa_ref ref_flags err_N1 err_N2 err_N5 err_N10 err_N20 err_N50 '// &
      'fl_N1 fl_N2 fl_N5 fl_N10 fl_N20 fl_N50 err_conv nsub_conv cls_conv ref_conv_flags ratio_trial'
   do ip = 1, npt
      write (uout, '(24ES14.6)') rows(:, ip)
   end do
   close (uout)
   write (*, '(A,A,A,A,A,A,A,F9.1,A)') ' maymap ', trim(law), ' ', trim(map), ' ', trim(variant), ' done in ', &
      omp_get_wtime() - t0, ' s'
contains
   subroutine set_material(lw, m, ps)
      character(len=*), intent(in) :: lw
      type(np_material), intent(out) :: m
      logical, intent(out) :: ps
      integer :: inf
      m%ielast = NP_ELAST_PAPER
      select case (lw)
      case ('barlat', 'iso')
         ! AA2090-T3: Yld2004-18p a = 8 (Barlat et al. 2005); E = 70 GPa, nu = 0.33; Voce + linear fitted to the Swift
         ! law 646 (0.025 + k)^0.227 MPa (Yoon et al. 2000) on [0, 0.3], max deviation 0.62 %
         m%mu = 70.0d9/(2.0d0*1.33d0); m%kmod = 70.0d9/(3.0d0*(1.0d0 - 0.66d0))
         m%iyield = NP_YLD_BARLAT; m%abar = 8
         m%y0 = 279.6178550820906d6; m%ysat = 379.958d6; m%bvoce = 17.7001d0; m%nh = 1.0d0; m%hlin = 411.286d6
         call np_barlat_set(m, merge(1, 0, lw == 'barlat'), inf)
         ps = .true.
      case ('pla')   ! PLA: the transversely isotropic matrix energy of Section 8.1.2, Hosford n = 8
         m%mu = 1.1508d9; m%kmod = 41.924d9; m%iyield = NP_YLD_HOSFORD
         m%ielast = NP_ELAST_GUO; m%a0 = [1.0d0, 0.0d0, 0.0d0]; m%vv = 0.0d0   ! a_0 = e_1
         m%y0 = 52.28d6; m%ysat = 62.0d6; m%bvoce = 14.0d0; m%nh = 1.3d0
         ps = .false.
      case ('pc')    ! polycarbonate: the criterion of Section 8.1.1 (notebook yld_ghorbel.nb), Voce 62 -> 78 MPa,
                     ! b = 14, n_h = 1.3
         m%mu = 848.0d6; m%kmod = 3.25d9; m%iyield = NP_YLD_GHORBELF; m%rsig = 1.387d0; m%st0 = 62.0d6
         m%y0 = 62.0d6; m%ysat = 78.0d6; m%bvoce = 14.0d0; m%nh = 1.3d0
         ps = .false.
      case default
         stop 'law: barlat | iso | pla | pc'
      end select
   end subroutine set_material

   ! the update: one call, optionally substepped (np_update_sub from strain es)
   subroutine law_call(m, es, e, cp, kn, sub, s, d, cpo, ko, st, inf, nsub, minlev, maxlev)
      type(np_material), intent(in) :: m
      real(dp), intent(in) :: es(6), e(6), cp(6), kn
      logical, intent(in) :: sub
      real(dp), intent(out) :: s(6), d(6, 6), cpo(6), ko
      type(np_state), intent(out) :: st
      integer, intent(out) :: inf, nsub
      integer, intent(in), optional :: minlev, maxlev     ! minlev = maxlev: exactly 2**minlev parts
      real(dp) :: psi
      integer :: l0, l1
      l0 = 0; if (present(minlev)) l0 = minlev
      l1 = 10; if (present(maxlev)) l1 = maxlev
      if (sub) then
         call np_update_sub(m, es, e, cp, kn, s, d, cpo, ko, psi, inf, nsub, l1, st=st, minlev=l0)
      else
         nsub = 1
         call np_update_alg1(m, e, cp, kn, s, d, cpo, ko, psi, inf, st=st)
      end if
   end subroutine law_call

   ! plane stress: S33 = 0 by Newton on E33 with backtracking from e(3); best accepted iterate kept if |S33| <= 1e-6 Y0
   ! fixlev (optional, with sub): every evaluation with exactly 2**fixlev parts; S33 is then a smooth function of E33
   ! and its derivative is the tangent of the substepped update
   subroutine ps_update(m, es, e, cp, kn, sub, s, cpo, ko, st, info, nsub, ok, fixlev)
      type(np_material), intent(in) :: m
      real(dp), intent(in) :: es(6), cp(6), kn
      real(dp), intent(inout) :: e(6)
      logical, intent(in) :: sub
      real(dp), intent(out) :: s(6), cpo(6), ko
      type(np_state), intent(out) :: st
      integer, intent(out) :: info, nsub
      logical, intent(out) :: ok
      integer, intent(in), optional :: fixlev
      real(dp) :: d(6, 6), e3good, step, guesses(3), sb(6), cpb(6), kb, e3b, sbest
      type(np_state) :: stb
      integer :: it, ig, infb, nsb, lev
      logical :: good
      ok = .false.; good = .false.; sbest = huge(1.0d0); infb = NP_ERR_INPUT; e3good = e(3); step = 0.0d0; nsb = 1
      lev = 0                    ! the substep level, non-decreasing along the iteration on E33
      guesses = [e(3), e(3) - 0.5d0*abs(e(1) + e(2) - e(3)), e(3) + 0.5d0*abs(e(1) + e(2) - e(3))]
      ig = 1
      do it = 1, merge(25, 60, present(fixlev))
         if (present(fixlev)) then
            call law_call(m, es, e, cp, kn, sub, s, d, cpo, ko, st, info, nsub, minlev=fixlev, maxlev=fixlev)
         else
            call law_call(m, es, e, cp, kn, sub, s, d, cpo, ko, st, info, nsub, minlev=lev)
            if (sub .and. info == NP_OK) lev = max(lev, nint(log(real(nsub, dp))/log(2.0d0)))
         end if
         if (info /= NP_OK .and. info /= NP_UNCONVERGED) then
            if (good) then
               step = 0.5d0*step; e(3) = e3good + step
               if (abs(step) < 1.0d-14) exit
            else
               ig = ig + 1
               if (ig > 3) exit
               e(3) = guesses(ig)
            end if
            cycle
         end if
         good = .true.; e3good = e(3)
         if (abs(s(3)) < sbest) then
            sbest = abs(s(3)); sb = s; cpb = cpo; kb = ko; e3b = e(3); stb = st; infb = info; nsb = nsub
         end if
         if (abs(s(3)) <= 1.0d-8*m%y0) then
            ok = .true.; return
         end if
         if (abs(d(3, 3)) <= tiny(1.0d0)) exit
         step = -s(3)/d(3, 3)
         e(3) = e3good + step
      end do
      if (good .and. sbest <= 1.0d-6*m%y0) then
         s = sb; cpo = cpb; ko = kb; e(3) = e3b; st = stb; info = infb; nsub = nsb; ok = .true.
      else if (good) then
         info = 9
      end if
   end subroutine ps_update

   ! S33 = 0 by an E33 scan (201 samples, e(3) +- w) and bisection on a sign change between accepted updates
   ! (cls 1 found, 3 / 11 none, 12 jump); single steps without substepping
   ! lev (optional, with sub): every evaluation with exactly 2**lev parts. For a fixed number of parts S33 is a
   ! continuous function of E33, so a sign change between two accepted samples holds a root.
   subroutine scan_solve(m, e, cp, kn, w, s, cls, es, sub, cpout, kout, lev, nsubout, nsamp)
      type(np_material), intent(in) :: m
      real(dp), intent(inout) :: e(6)
      real(dp), intent(in) :: cp(6), kn, w
      real(dp), intent(out) :: s(6)
      integer, intent(out) :: cls
      real(dp), intent(in), optional :: es(6)
      logical, intent(in), optional :: sub
      real(dp), intent(out), optional :: cpout(6), kout
      integer, intent(in), optional :: lev, nsamp          ! nsamp: samples of the scan (201 when absent)
      integer, intent(out), optional :: nsubout
      type(np_state) :: st
      real(dp) :: d(6, 6), x(201), f(201), lo, hi, flo, fmid, xm, c, cpo(6), ko, e0s(6)
      logical :: acc(201), anyref, usesub
      integer :: i, inf, it, nsub, la, lb, ns
      usesub = .false.; if (present(sub)) usesub = sub
      ns = 201; if (present(nsamp)) ns = max(3, min(201, nsamp))
      la = 0; lb = 10
      if (present(lev)) then
         la = lev; lb = lev
      end if
      e0s = e; if (present(es)) e0s = es
      c = e(3); anyref = .false.
      do i = 1, ns
         x(i) = c - w + 2.0d0*w*real(i - 1, dp)/real(ns - 1, dp)
         e(3) = x(i)
         call law_call(m, e0s, e, cp, kn, usesub, s, d, cpo, ko, st, inf, nsub, minlev=la, maxlev=lb)
         acc(i) = (inf == NP_OK .or. inf == NP_UNCONVERGED); f(i) = s(3)
         if (.not. acc(i)) anyref = .true.
      end do
      cls = merge(3, 11, anyref)
      do i = 1, ns - 1
         if (.not. (acc(i) .and. acc(i + 1))) cycle
         if (f(i)*f(i + 1) > 0.0d0) cycle
         lo = x(i); hi = x(i + 1); flo = f(i)
         do it = 1, 200
            xm = 0.5d0*(lo + hi)
            e(3) = xm
            call law_call(m, e0s, e, cp, kn, usesub, s, d, cpo, ko, st, inf, nsub, minlev=la, maxlev=lb)
            if (.not. (inf == NP_OK .or. inf == NP_UNCONVERGED)) exit
            fmid = s(3)
            if (abs(fmid) <= 1.0d-8*m%y0) then
               cls = 1
               if (present(cpout)) cpout = cpo
               if (present(kout)) kout = ko
               if (present(nsubout)) nsubout = nsub
               return
            end if
            if (flo*fmid > 0.0d0) then
               lo = xm; flo = fmid
            else
               hi = xm
            end if
            if (hi - lo < 1.0d-15*max(1.0d0, abs(xm))) exit
         end do
         if (.not. (inf == NP_OK .or. inf == NP_UNCONVERGED)) then
            cls = 3
         else
            cls = 12
         end if
      end do
   end subroutine scan_solve

   subroutine first_yield(ebo, e33)
      real(dp), intent(out) :: ebo, e33
      type(np_material) :: m
      type(np_state) :: st
      real(dp) :: lo, hi, e(6), s(6), cpo(6), ko
      integer :: it, inf, ns1
      logical :: ok
      m = mat; m%nn = 1; m%eps_theta = 0.0d0
      lo = 0.0d0; hi = 0.05d0; e33 = 0.0d0
      do it = 1, 200
         ebo = 0.5d0*(lo + hi)
         e = 0.0d0; e(1) = ebo; e(2) = ebo; e(3) = -ebo
         call ps_update(m, [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], &
                        0.0d0, .false., s, cpo, ko, st, inf, ns1, ok)
         if (st%ratio < 1.0d0) then
            lo = ebo; e33 = e(3)
         else
            hi = ebo
         end if
         if (hi - lo < 1.0d-15) exit
      end do
      ebo = lo
   end subroutine first_yield

   ! one increment e0 -> e0 + de (undeformed history), the update called once (substepped inside if sub); plane stress or 3D strain
   subroutine one_step(m, e0, de, sub, hint, s, kap, e33, inf, nsub, st)
      type(np_material), intent(in) :: m
      real(dp), intent(in) :: e0(6), de(6), hint
      logical, intent(in) :: sub
      real(dp), intent(out) :: s(6), kap, e33
      integer, intent(out) :: inf, nsub
      type(np_state), intent(out) :: st
      real(dp) :: e(6), cpo(6), d(6, 6)
      logical :: ok
      integer :: lv, lv0, nsub0
      e = e0 + de
      if (pstress) then
         e(3) = hint
         call ps_update(m, e0, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, sub, s, cpo, kap, st, inf, nsub, ok)
         nsub0 = nsub
         if (.not. ok .and. sub) then
            ! the substepped map jumps where the substep count changes: Newton on E33 may cycle over a jump;
            ! scan E33 and bisect a sign change of S33 between accepted updates
            e = e0 + de; e(3) = hint
            call scan_solve(m, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, max(0.02d0, 2.0d0*abs(hint - e0(3))), &
                            s, inf, es=e0, sub=.true., kout=kap)
            ok = (inf == 1)
            if (.not. ok) then
               ! still no state: the evaluations had different numbers of parts and S33 jumps where the number
               ! changes. With a fixed number of parts S33 is smooth in E33: Newton on E33 with 1, 2, 4, ... parts
               ! in turn, and the first number of parts with a root is taken. This is the rule of Section 6.5 (the
               ! number of parts of a point is not decreased along the iteration that enforces plane stress).
               ! The numbers of parts tried start at the one the failed iteration had reached (it is not decreased).
               lv0 = max(0, min(10, nint(log(real(max(nsub0, 1), dp))/log(2.0d0))))
               do lv = lv0, 10
                  e = e0 + de; e(3) = hint
                  call ps_update(m, e0, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, .true., s, cpo, kap, st, &
                                 inf, nsub, ok, fixlev=lv)
                  if (ok) exit
                  ! a coarse scan with the same number of parts, bisection between accepted samples
                  e = e0 + de; e(3) = hint
                  call scan_solve(m, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, &
                                  max(0.02d0, 2.0d0*abs(hint - e0(3))), s, inf, es=e0, sub=.true., kout=kap, &
                                  lev=lv, nsubout=nsub, nsamp=21)
                  ok = (inf == 1)
                  if (ok) exit
               end do
            end if
            inf = merge(NP_OK, 9, ok)
         end if
         if (.not. ok .and. inf == NP_OK) inf = 9
      else
         call law_call(m, e0, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, sub, s, d, cpo, kap, st, inf, nsub)
      end if
      e33 = e(3)
   end subroutine one_step

   ! the reference: ns converged increments of the proportional path (plane stress at every increment)
   subroutine path(m, e0, de, ns, s, kap, e33, flags)
      type(np_material), intent(in) :: m
      real(dp), intent(in) :: e0(6), de(6)
      integer, intent(in) :: ns
      real(dp), intent(out) :: s(6), kap, e33
      integer, intent(out) :: flags
      real(dp) :: cp(6), kn, e(6), es(6), cpo(6), ko, d(6, 6), psi
      type(np_state) :: st
      integer :: k, inf, nsub
      logical :: ok
      cp = [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]; kn = 0.0d0; flags = 0
      e = e0
      do k = 1, ns
         es = e
         e(1:2) = e0(1:2) + de(1:2)*real(k, dp)/real(ns, dp)
         ! every reference increment through the adaptive substepping (an increment the Picard iteration cannot take is split)
         if (pstress) then
            e(3) = e(3) - (0.33d0/0.67d0)*(de(1) + de(2))/real(ns, dp)
            call ps_update(m, es, e, cp, kn, .true., s, cpo, ko, st, inf, nsub, ok)
            if (.not. ok) then
               call scan_solve(m, e, cp, kn, 0.02d0, s, inf, es=es, sub=.true., cpout=cpo, kout=ko)
               ok = (inf == 1); inf = merge(NP_OK, 9, ok)
            end if
         else
            call np_update_sub(m, es, e, cp, kn, s, d, cpo, ko, psi, inf, nsub, 10, st=st)
            ok = (inf == NP_OK)
         end if
         if (.not. ok .or. inf /= NP_OK) then
            flags = flags + 1
            if (.not. ok) then
               flags = flags + 1000; exit
            end if
         end if
         cp = cpo; kn = ko
      end do
      kap = kn; e33 = e(3)
   end subroutine path

   subroutine one_point(e0, de, row)
      real(dp), intent(in) :: e0(6), de(6)
      real(dp), intent(out) :: row(NCOL)
      type(np_material) :: m1, mc
      type(np_state) :: st
      real(dp) :: sref(6), s(6), kap, e33, kd, e33d, e(6)
      integer :: fl, flr, inf, k, nsub
      row = 0.0d0
      ! the scheme of the maps: one update per increment (N_n = 1), or the baseline (converged in every increment)
      m1 = mat
      if (base) then
         m1%classical = .true.; m1%nn = 50; m1%eps_theta = 1.0d-10
      else
         m1%nn = 1; m1%eps_theta = 0.0d0
      end if
      call path(m1, e0, de, nsref, sref, kap, e33, flr)
      row(3) = e0(1) + de(1); row(4) = e0(2) + de(2); row(5) = e33; row(6) = kap; row(7) = real(flr, dp)
      do k = 1, 6
         call path(m1, e0, de, NNS(k), s, kd, e33d, fl)
         row(7 + k) = errn(s, sref, max(fl, flr))
         row(13 + k) = real(fl, dp)
      end do
      ! the converged update: one increment (substepped where it fails) against the converged update in 90 increments
      mc = mat; mc%nn = 50; mc%eps_theta = 1.0d-10
      if (base) then
         mc%classical = .true.
      end if
      call path(mc, e0, de, nsref, sref, kap, e33, flr)
      call one_step(mc, e0, de, .true., e33, s, kd, e33d, inf, nsub, st)
      row(20) = errn(s, sref, max(merge(0, 1000, inf == NP_OK), flr))
      row(21) = real(nsub, dp); row(22) = real(inf, dp); row(23) = real(flr, dp)
      ! the trial ratio of the single increment
      e = e0 + de; e(3) = e33d
      block
         real(dp) :: sv(6), d(6, 6), cpo(6), ko, psi
         integer :: inf2
         call np_update_alg1(mc, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, sv, d, cpo, ko, psi, inf2, st=st, &
                            noexit=.true.)
         row(24) = st%ratio
      end block
   end subroutine one_point

   ! the error of the maps: normal components only
   pure function errn(s, sref, fl) result(r)
      real(dp), intent(in) :: s(6), sref(6)
      integer, intent(in) :: fl
      real(dp) :: r, nr
      nr = sqrt(sum(sref(1:3)**2))
      if (fl >= 1000) then
         r = -1.0d0
      else if (nr > 0.0d0) then
         r = sqrt(sum((s(1:3) - sref(1:3))**2))/nr
      else
         r = 0.0d0
      end if
   end function errn

   pure function fro(s) result(r)
      real(dp), intent(in) :: s(6)
      real(dp) :: r
      r = sqrt(sum(s(1:3)**2) + 2.0d0*sum(s(4:6)**2))
   end function fro

end program maymap
