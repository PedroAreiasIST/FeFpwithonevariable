!=============================================================================
! basemap -- the baseline panel of the PLA map (Section 8.2): the Newton iteration with staggered hardening
! (np_update_base) against Algorithm 1, on the grid, reference and error measure of isomap.f90 (same routines).
!
! usage: basemap <law> <map> <n> <amax> <nsref> <out>      (the paper: pla virgin 81 0.25 1000)
! Columns: 1 a, 2 b, 3-7 as isomap, 8 estimate of the reference error, 10 err_base (baseline, one increment),
!   16 its status, 24 its Newton iterations, 11 err_np1 (Algorithm 1 converged, one increment, no substepping),
!   17 its status, 25 its Picard updates, 22 err_sub, 23 nsub (Algorithm 1 with substepping, as isomap)
! Error = ||S - S_ref||/||S_ref|| (Frobenius, second Piola-Kirchhoff), -1 when no state was produced.
!=============================================================================
program basemap
   use npicard_plasticity
   use omp_lib
   implicit none
   type(np_material) :: mat
   character(len=16) :: law, map
   character(len=256) :: outf, arg
   integer :: n, nsref, ip, npt, info, uout, i, j
   logical :: pstress
   real(dp) :: amax, eb, e33b, a, b, e0(6), de(6), t0
   integer, parameter :: NNS(6) = [1, 2, 5, 10, 20, 50], NCOL = 26
   real(dp), allocatable :: rows(:, :)
   real(dp) :: row(NCOL)
   call get_command_argument(1, law)
   call get_command_argument(2, map)
   call get_command_argument(3, arg); read (arg, *) n
   call get_command_argument(4, arg); read (arg, *) amax
   call get_command_argument(5, arg); read (arg, *) nsref
   call get_command_argument(6, outf)
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
   write (uout, '(A,A,A,A,A,I0,A,ES12.4,A,I0,A,ES16.8,A,L1,A,F9.1,A)') '# basemap ', trim(law), ' ', trim(map), ' n ', n, &
      ' amax ', amax, ' nsref ', nsref, ' e_b ', eb, ' planestress ', pstress, ' wall ', omp_get_wtime() - t0, &
      ' s; Algorithm 1 (np_update_alg1) + adaptive substepping; baseline np_update_base'
   write (uout, '(A)') '# a b E11 E22 E33ref kappa_ref ref_flags ref_err_est ratio_trial '// &
      'err_Nn1 err_Nn2 err_bsub nsub_bsub err_Nn20 err_Nn50 cls_Nn1 cls_Nn2 cls_bsub cls_Nn10 cls_Nn20 cls_Nn50 '// &
      'err_sub nsub its1 cls_sub E33_sub'
   do ip = 1, npt
      write (uout, '(26ES14.6)') rows(:, ip)
   end do
   close (uout)
   write (*, '(A,A,A,A,A,ES14.6,A,F9.1,A)') ' basemap ', trim(law), ' ', trim(map), ': e_b ', eb, ', done in ', &
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
      case ('pc')    ! polycarbonate: the criterion of Section 8.1.1 (notebook yld_ghorbel.nb), Voce 62 -> 78 MPa
         m%mu = 848.0d6; m%kmod = 3.25d9; m%iyield = NP_YLD_GHORBELF; m%rsig = 1.387d0; m%st0 = 62.0d6
         m%y0 = 62.0d6; m%ysat = 78.0d6; m%bvoce = 14.0d0; m%nh = 1.3d0
         ps = .false.
      case default
         stop 'law: barlat | iso | pla | pc'
      end select
   end subroutine set_material

   ! the update: one call, optionally substepped (np_update_sub from strain es)
   subroutine law_call(m, es, e, cp, kn, sub, s, d, cpo, ko, st, inf, nsub)
      type(np_material), intent(in) :: m
      real(dp), intent(in) :: es(6), e(6), cp(6), kn
      logical, intent(in) :: sub
      real(dp), intent(out) :: s(6), d(6, 6), cpo(6), ko
      type(np_state), intent(out) :: st
      integer, intent(out) :: inf, nsub
      real(dp) :: psi
      ! Algorithm 1 (np_update_alg1); the baseline has its own routine
      if (sub) then
         call np_update_sub(m, es, e, cp, kn, s, d, cpo, ko, psi, inf, nsub, 10, st=st)
      else if (m%classical) then
         nsub = 1
         call np_update_base(m, e, cp, kn, s, d, cpo, ko, psi, inf, st=st)
      else
         nsub = 1
         call np_update_alg1(m, e, cp, kn, s, d, cpo, ko, psi, inf, st=st, noexit=.true.)
      end if
   end subroutine law_call

   ! plane stress: S33 = 0 by Newton on E33 with backtracking from e(3); best accepted iterate kept if |S33| <= 1e-6 Y0
   subroutine ps_update(m, es, e, cp, kn, sub, s, cpo, ko, st, info, nsub, ok)
      type(np_material), intent(in) :: m
      real(dp), intent(in) :: es(6), cp(6), kn
      real(dp), intent(inout) :: e(6)
      logical, intent(in) :: sub
      real(dp), intent(out) :: s(6), cpo(6), ko
      type(np_state), intent(out) :: st
      integer, intent(out) :: info, nsub
      logical, intent(out) :: ok
      real(dp) :: d(6, 6), e3good, step, guesses(3), sb(6), cpb(6), kb, e3b, sbest
      type(np_state) :: stb
      integer :: it, ig, infb, nsb
      logical :: good
      ok = .false.; good = .false.; sbest = huge(1.0d0); infb = NP_ERR_INPUT; e3good = e(3); step = 0.0d0; nsb = 1
      guesses = [e(3), e(3) - 0.5d0*abs(e(1) + e(2) - e(3)), e(3) + 0.5d0*abs(e(1) + e(2) - e(3))]
      ig = 1
      do it = 1, 60
         call law_call(m, es, e, cp, kn, sub, s, d, cpo, ko, st, info, nsub)
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
   subroutine scan_solve(m, e, cp, kn, w, s, cls, es, sub, cpout, kout)
      type(np_material), intent(in) :: m
      real(dp), intent(inout) :: e(6)
      real(dp), intent(in) :: cp(6), kn, w
      real(dp), intent(out) :: s(6)
      integer, intent(out) :: cls
      real(dp), intent(in), optional :: es(6)
      logical, intent(in), optional :: sub
      real(dp), intent(out), optional :: cpout(6), kout
      type(np_state) :: st
      real(dp) :: d(6, 6), x(201), f(201), lo, hi, flo, fmid, xm, c, cpo(6), ko, e0s(6)
      logical :: acc(201), anyref, usesub
      integer :: i, inf, it, nsub
      usesub = .false.; if (present(sub)) usesub = sub
      e0s = e; if (present(es)) e0s = es
      c = e(3); anyref = .false.
      do i = 1, 201
         x(i) = c - w + 2.0d0*w*real(i - 1, dp)/200.0d0
         e(3) = x(i)
         call law_call(m, e0s, e, cp, kn, usesub, s, d, cpo, ko, st, inf, nsub)
         acc(i) = (inf == NP_OK .or. inf == NP_UNCONVERGED); f(i) = s(3)
         if (.not. acc(i)) anyref = .true.
      end do
      cls = merge(3, 11, anyref)
      do i = 1, 200
         if (.not. (acc(i) .and. acc(i + 1))) cycle
         if (f(i)*f(i + 1) > 0.0d0) cycle
         lo = x(i); hi = x(i + 1); flo = f(i)
         do it = 1, 200
            xm = 0.5d0*(lo + hi)
            e(3) = xm
            call law_call(m, e0s, e, cp, kn, usesub, s, d, cpo, ko, st, inf, nsub)
            if (.not. (inf == NP_OK .or. inf == NP_UNCONVERGED)) exit
            fmid = s(3)
            if (abs(fmid) <= 1.0d-8*m%y0) then
               cls = 1
               if (present(cpout)) cpout = cpo
               if (present(kout)) kout = ko
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
      e = e0 + de
      if (pstress) then
         e(3) = hint
         call ps_update(m, e0, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, sub, s, cpo, kap, st, inf, nsub, ok)
         if (.not. ok .and. sub) then
            ! the substepped map jumps where the substep count changes: Newton on E33 may cycle over a jump;
            ! scan E33 and bisect a sign change of S33 between accepted updates
            e = e0 + de; e(3) = hint
            call scan_solve(m, e, [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0], 0.0d0, max(0.02d0, 2.0d0*abs(hint - e0(3))), &
                            s, inf, es=e0, sub=.true., kout=kap)
            ok = (inf == 1)
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
      ! the baseline map (columns in the header)
      ! 12 err_bsub, 13 nsub_bsub, 18 cls_bsub (the baseline with adaptive substepping)
      real(dp), intent(in) :: e0(6), de(6)
      real(dp), intent(out) :: row(NCOL)
      type(np_material) :: mref, mb, mc
      type(np_state) :: st
      real(dp) :: sref(6), s(6), kap, e33, kd, e33d, nref
      integer :: fl, fl2, inf, nsub
      row = 0.0d0
      mref = mat; mref%nn = 50; mref%eps_theta = 1.0d-12
      call path(mref, e0, de, nsref, s, kd, e33d, fl)
      call path(mref, e0, de, 2*nsref, sref, kap, e33, fl2)
      nref = fro(sref)
      row(8) = relerr(s, sref, nref, max(fl, fl2))
      sref = 2.0d0*sref - s
      nref = fro(sref)
      row(3) = e0(1) + de(1); row(4) = e0(2) + de(2); row(5) = e33; row(6) = kap; row(7) = real(fl + fl2, dp)
      mc = mat; mc%nn = 50; mc%eps_theta = 1.0d-10
      call one_step(mc, e0, de, .true., e33, s, kd, e33d, inf, nsub, st)
      row(22) = relerr(s, sref, nref, merge(0, 1000, inf == NP_OK)); row(23) = real(nsub, dp)
      call one_step(mc, e0, de, .false., e33, s, kd, e33d, inf, nsub, st)
      row(11) = relerr(s, sref, nref, merge(0, 1000, inf == NP_OK)); row(17) = real(inf, dp); row(25) = real(st%its, dp)
      mb = mat; mb%classical = .true.; mb%nn = 50; mb%eps_theta = 1.0d-10
      call one_step(mb, e0, de, .false., e33, s, kd, e33d, inf, nsub, st)
      row(10) = relerr(s, sref, nref, merge(0, 1000, inf == NP_OK)); row(16) = real(inf, dp); row(24) = real(st%its, dp)
      row(9) = st%ratio
      ! the baseline with the adaptive substepping of the Newton/Picard update: an increment whose Newton iteration
      ! fails is divided into 2, 4, ... parts
      call one_step(mb, e0, de, .true., e33, s, kd, e33d, inf, nsub, st)
      row(12) = relerr(s, sref, nref, merge(0, 1000, inf == NP_OK)); row(13) = real(nsub, dp); row(18) = real(inf, dp)
   end subroutine one_point

   pure function fro(s) result(r)
      real(dp), intent(in) :: s(6)
      real(dp) :: r
      r = sqrt(sum(s(1:3)**2) + 2.0d0*sum(s(4:6)**2))
   end function fro

   pure function relerr(s, sref, nref, fl) result(r)
      real(dp), intent(in) :: s(6), sref(6), nref
      integer, intent(in) :: fl
      real(dp) :: r
      if (fl >= 1000) then
         r = -1.0d0
      else if (nref > 0.0d0) then
         r = fro(s - sref)/nref
      else
         r = 0.0d0
      end if
   end function relerr
end program basemap
