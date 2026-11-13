! symcheck -- how symmetric is the tangent dS/dE (Section 5.2)?
!   (a) single increments of the converged update (Algorithm 1, eps_theta = 1e-10, cap 50, no substepping) at random
!       general strains (six components, trial ratios up to about 1.6): the maps hold coaxial states only;
!   (b) substepped increments (np_update_sub, exact tangent) on the grids of the maps (E11, E22 prescribed on
!       81 x 81 points, other components zero; range 0.25 for the polycarbonate and PLA, 0.02 for AA2090-T3 and the
!       isotropic set, as dencheck.f90);
!   (c) substepped increments at random general strains up to 6 yield strains.
! Measure: |D - D^T|_F / |D|_F of the returned 6 x 6 tangent (engineering shear strains: symmetric for a symmetric
! fourth-order tangent; the elastic tangent is the control, printed first).
program symcheck
   use npicard_plasticity
   implicit none
   type(np_material) :: m
   type(np_state) :: st
   real(dp) :: e(6), e0(6), cp(6), s(6), d(6, 6), cpo(6), ko, psi, amax, r6(6), amp, ey, a
   real(dp), allocatable :: v(:)
   integer :: info, il, i, j, n, k, nsub, inf2, nv, nsubbed
   character(len=22) :: nm(4)
   nm = ['polycarbonate         ', 'PLA matrix energy     ', 'AA2090-T3 Barlat      ', 'isotropic Hosford a=8 ']
   cp = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]; e0 = 0.0_dp
   n = 81
   allocate (v(n*n + 4000))
   do il = 1, 4
      m = np_material(); amax = 0.25_dp
      select case (il)
      case (1)
         m%mu = 848.0e6_dp; m%kmod = 3.25e9_dp; m%iyield = NP_YLD_GHORBELF; m%rsig = 1.387_dp; m%st0 = 62.0e6_dp
         m%y0 = 62.0e6_dp; m%ysat = 78.0e6_dp; m%bvoce = 14.0_dp; m%nh = 1.3_dp
      case (2)
         m%mu = 1.1508e9_dp; m%kmod = 41.924e9_dp; m%iyield = NP_YLD_HOSFORD
         m%y0 = 52.28e6_dp; m%ysat = 62.0e6_dp; m%bvoce = 14.0_dp; m%nh = 1.3_dp
         m%ielast = NP_ELAST_GUO; m%a0 = [1.0_dp, 0.0_dp, 0.0_dp]; m%vv = 0.0_dp
      case (3, 4)
         m%mu = 70.0e9_dp/(2.0_dp*1.33_dp); m%kmod = 70.0e9_dp/(3.0_dp*(1.0_dp - 0.66_dp)); m%iyield = NP_YLD_BARLAT
         m%abar = 8
         m%y0 = 279.6178550820906e6_dp; m%ysat = 379.958e6_dp; m%bvoce = 17.7001_dp; m%nh = 1.0_dp; m%hlin = 411.286e6_dp
         call np_barlat_set(m, merge(1, 0, il == 3), inf2)
         amax = 0.02_dp
      end select
      m%nn = 50; m%eps_theta = 1.0e-10_dp
      ey = m%y0/(2.0_dp*m%mu)
      write (*, '(1X,A)') nm(il)
      ! control: an elastic state
      e = 0.2_dp*ey*[1.0_dp, -0.3_dp, 0.2_dp, 0.4_dp, -0.1_dp, 0.25_dp]
      call np_update_alg1(m, e, cp, 0.0_dp, s, d, cpo, ko, psi, info, st=st)
      write (*, '(4X,A,ES9.2,A,L1)') 'control, elastic general state: ', asy_of(d), '  plastic ', st%plastic
      ! (a) single increments at random general strains
      call random_seed(put=[(4711 + 13*i, i=1, 64)])
      nv = 0
      do k = 1, 4000
         call random_number(r6); call random_number(amp)
         e = (r6 - 0.5_dp)*2.0_dp*ey*(0.6_dp + 1.6_dp*amp)
         call np_update_alg1(m, e, cp, 0.0_dp, s, d, cpo, ko, psi, info, st=st)
         if (info /= NP_OK .or. .not. st%plastic) cycle
         nv = nv + 1; v(nv) = asy_of(d)
      end do
      call report('(a) single increments, random general strains (4000 drawn)', v, nv, -1)
      ! (b) substepped increments on the grid of the map
      nv = 0; nsubbed = 0
      do i = 1, n
         do j = 1, n
            e = 0.0_dp
            e(1) = -amax + 2.0_dp*amax*real(i - 1, dp)/real(n - 1, dp)
            e(2) = -amax + 2.0_dp*amax*real(j - 1, dp)/real(n - 1, dp)
            call np_update_sub(m, e0, e, cp, 0.0_dp, s, d, cpo, ko, psi, info, nsub, 10, st=st)
            if (info /= NP_OK .or. nsub < 2) cycle
            nv = nv + 1; v(nv) = asy_of(d)
         end do
      end do
      call report('(b) substepped increments on the 81 x 81 grid (coaxial)', v, nv, -1)
      ! (c) substepped increments at random general strains up to 6 yield strains
      nv = 0
      do k = 1, 2000
         call random_number(r6); call random_number(amp)
         e = (r6 - 0.5_dp)*2.0_dp*ey*(0.5_dp + 5.5_dp*amp)
         call np_update_sub(m, e0, e, cp, 0.0_dp, s, d, cpo, ko, psi, info, nsub, 10, st=st)
         if (info /= NP_OK .or. nsub < 2) cycle
         nv = nv + 1; v(nv) = asy_of(d)
      end do
      call report('(c) substepped increments, random general strains (2000 drawn)', v, nv, -1)
   end do
contains
   function asy_of(dm) result(r)
      real(dp), intent(in) :: dm(6, 6)
      real(dp) :: r
      r = sqrt(sum((dm - transpose(dm))**2))/sqrt(sum(dm**2))
   end function asy_of
   subroutine report(what, a, na, dummy)
      character(len=*), intent(in) :: what
      real(dp), intent(inout) :: a(:)
      integer, intent(in) :: na, dummy
      real(dp) :: t
      integer :: i1, j1
      if (na < 1) then
         write (*, '(4X,A,A)') what, ': none'; return
      end if
      do i1 = 2, na
         t = a(i1); j1 = i1 - 1
         do while (j1 >= 1)
            if (a(j1) <= t) exit
            a(j1 + 1) = a(j1); j1 = j1 - 1
         end do
         a(j1 + 1) = t
      end do
      write (*, '(4X,A,A,I6,A,ES9.2,A,ES9.2,A,ES9.2,A,F5.1,A)') what, ': n ', na, '  median ', a(max(1, na/2)), '  95th ', &
         a(max(1, (95*na)/100)), '  max ', a(na), '  above 1e-3: ', 100.0_dp*real(count(a(1:na) > 1.0e-3_dp), dp)/real(na, dp), ' %'
      if (dummy > 0) continue
   end subroutine report
end program symcheck
