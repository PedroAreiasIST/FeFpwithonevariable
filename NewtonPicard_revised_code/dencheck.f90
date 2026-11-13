! dencheck -- sign of the denominator of the multiplier increment and of det(I - R_N:H:T) at the converged single
! increments of the grids with prescribed strains (E11, E22 on 81 x 81 points, other components zero), and the
! unsymmetric part of the tangent dS/dE (Section 5.2). Polycarbonate, PLA, PLA with the isotropic energy,
! AA2090-T3, isotropic set (Hosford a = 8). Optional argument: points per side of the grid.
program dencheck
   use npicard_plasticity
   implicit none
   type(np_material) :: m
   type(np_state) :: st
   real(dp) :: e(6), cp(6), s(6), d(6, 6), cpo(6), ko, psi, amax, asy, asymax, rmin, rmax, dsym(6, 6), dmin, hmin, detmin, &
               asymed(7000)
   integer :: info, il, i, j, n, nconv, nneg, ndet, nsame, nas, inf2, k
   integer :: nnear = 0
   real(dp) :: asynear = 0.0_dp, asyfar = 0.0_dp
   character(len=22) :: nm(5)
   real(dp), parameter :: SF(6) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
   nm = ['polycarbonate         ', 'PLA matrix energy     ', 'PLA isotropic energy  ', 'AA2090-T3 Barlat      ', &
         'isotropic Hosford a=8 ']
   cp = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   n = 81
   block
      character(len=32) :: arg
      call get_command_argument(1, arg)      ! optional: points per side of the grid (default 81, the maps)
      if (len_trim(arg) > 0) read (arg, *) n
   end block
   do il = 1, 5
      m = np_material(); amax = 0.25_dp
      select case (il)
      case (1)
         m%mu = 848.0e6_dp; m%kmod = 3.25e9_dp; m%iyield = NP_YLD_GHORBELF; m%rsig = 1.387_dp; m%st0 = 62.0e6_dp
         m%y0 = 62.0e6_dp; m%ysat = 78.0e6_dp; m%bvoce = 14.0_dp; m%nh = 1.3_dp
      case (2, 3)
         m%mu = 1.1508e9_dp; m%kmod = 41.924e9_dp; m%iyield = NP_YLD_HOSFORD
         m%y0 = 52.28e6_dp; m%ysat = 62.0e6_dp; m%bvoce = 14.0_dp; m%nh = 1.3_dp
         if (il == 2) then
            m%ielast = NP_ELAST_GUO; m%a0 = [1.0_dp, 0.0_dp, 0.0_dp]; m%vv = 0.0_dp
         end if
      case (4, 5)
         m%mu = 70.0e9_dp/(2.0_dp*1.33_dp); m%kmod = 70.0e9_dp/(3.0_dp*(1.0_dp - 0.66_dp)); m%iyield = NP_YLD_BARLAT; m%abar = 8
         m%y0 = 279.6178550820906e6_dp; m%ysat = 379.958e6_dp; m%bvoce = 17.7001_dp; m%nh = 1.0_dp; m%hlin = 411.286e6_dp
         call np_barlat_set(m, merge(1, 0, il == 4), inf2)
         amax = 0.02_dp
      end select
      m%nn = 50; m%eps_theta = 1.0e-10_dp
      nconv = 0; nneg = 0; ndet = 0; nsame = 0; nas = 0; asymax = 0.0_dp; rmin = huge(1.0_dp); rmax = 0.0_dp
      dmin = huge(1.0_dp); hmin = huge(1.0_dp); detmin = huge(1.0_dp)
      do i = 1, n
         do j = 1, n
            e = 0.0_dp
            e(1) = -amax + 2.0_dp*amax*real(i - 1, dp)/real(n - 1, dp)
            e(2) = -amax + 2.0_dp*amax*real(j - 1, dp)/real(n - 1, dp)
            call np_update_alg1(m, e, cp, 0.0_dp, s, d, cpo, ko, psi, info, st=st)
            if (info /= NP_OK .or. .not. st%plastic) cycle
            nconv = nconv + 1
            if (st%dg_den < 0.0_dp) then
               nneg = nneg + 1; rmin = min(rmin, st%ratio); rmax = max(rmax, st%ratio)
            end if
            if (il == 1) then          ! polycarbonate: next to the band of the filter (|I1| <= 20 l) and away from it
               if (abs(sum(st%m6(1:3))) <= 20.0_dp*0.075_dp*m%st0) then
                  nnear = nnear + 1; asynear = max(asynear, asy_of(d))
               else
                  asyfar = max(asyfar, asy_of(d))
               end if
            end if
            if (st%dg_det < 0.0_dp) ndet = ndet + 1
            if ((st%dg_den < 0.0_dp) .eqv. (st%dg_det < 0.0_dp)) nsame = nsame + 1
            dmin = min(dmin, st%dg_den/st%dg_h0); detmin = min(detmin, st%dg_det)
            ! unsymmetric part of dS/dE in tensor components (the update returns engineering shear columns: undo the factor)
            dsym = d
            asy = sqrt(sum((dsym - transpose(dsym))**2))/sqrt(sum(dsym**2))
            asymax = max(asymax, asy)
            if (nconv <= size(asymed)) asymed(nconv) = asy
         end do
      end do
      ! median of the unsymmetric part
      k = min(nconv, size(asymed))
      call sort(asymed(1:k))
      write (*, '(1X,A,A,I5,A,I5,A,F5.1,A,I5,A,I5,A,ES10.2,A,ES10.2,A,ES9.2,A,ES9.2)') nm(il), ': converged single increments ', nconv, &
         '; denominator < 0 at ', nneg, ' (', 100.0_dp*real(nneg, dp)/real(max(nconv, 1), dp), ' %), det < 0 at ', ndet, &
         ', same sign at ', nsame, '; min den/(H*+alpha eta) ', dmin, ', min det ', detmin, '; |D - D^T|/|D| median ', &
         asymed(max(1, k/2)), ' max ', asymax
      if (nneg > 0) write (*, '(4X,A,F6.2,A,F6.2)') 'trial ratios of the states with a negative denominator: ', rmin, ' to ', rmax
      if (il == 1) write (*, '(4X,A,I4,A,ES9.2,A,ES9.2)') 'unsymmetric part: ', nnear, ' increments end with |I1| <= 20 l, max ', asynear, &
         '; away from the band max ', asyfar
   end do
contains
   function asy_of(dm) result(r)
      real(dp), intent(in) :: dm(6, 6)
      real(dp) :: r
      r = sqrt(sum((dm - transpose(dm))**2))/sqrt(sum(dm**2))
   end function asy_of
   subroutine sort(a)
      real(dp), intent(inout) :: a(:)
      real(dp) :: t
      integer :: i1, j1
      do i1 = 2, size(a)
         t = a(i1); j1 = i1 - 1
         do while (j1 >= 1)
            if (a(j1) <= t) exit
            a(j1 + 1) = a(j1); j1 = j1 - 1
         end do
         a(j1 + 1) = t
      end do
   end subroutine sort
end program dencheck
