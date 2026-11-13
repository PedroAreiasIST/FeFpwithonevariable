!=============================================================================
! raydiag -- why a single increment fails (Section 8.2). For the plastic states of the two maps with prescribed
! strains (pc, pla: E11, E22 in [-0.25, 0.25] on 81 x 81 points, other components zero, from the undeformed state):
! what happens on the first frozen ray of Algorithm 1, whose direction is the trial normal.
!   ratio          trial ratio seq_tr/Y
!   g0             the Newton estimate f(0)/s(0) from which the bracket is doubled
!   fmin           the minimum of f/Y_tr on the ray (> 0: the ray has no root)
!   g1, g2         the first root and the end of the interval f <= 0 (0, 0: no root)
!   kdbl           >= 0: the doubling finds f <= 0 after kdbl doublings; -1: it finds no sign change (it passes the
!                  interval [g1, g2] when the ray has a root)
!   status1        status of np_update_alg1 with exactly one update (0 = stress returned, 3 = bracketing failed)
! usage: raydiag <pc|pla> <out>
!=============================================================================
program raydiag
   use npicard_plasticity
   use omp_lib
   implicit none
   type(np_material) :: mat, m1
   character(len=16) :: law
   character(len=256) :: outf, arg
   integer, parameter :: n = 81, NCOL = 11
   integer :: ip, info, uout, i, j, kdbl, inf1
   real(dp), parameter :: amax = 0.25d0
   real(dp) :: a, b, e(6), cp0(6), s(6), d(6, 6), cpo(6), ko, psi, ratio, g0, gmin, fmin, g1, g2
   type(np_state) :: st
   real(dp) :: rows(NCOL, n*n), row(NCOL)
   call get_command_argument(1, law)
   call get_command_argument(2, outf)
   mat%ielast = NP_ELAST_PAPER
   select case (trim(law))
   case ('pla')   ! as isomap.f90
      mat%mu = 1.1508d9; mat%kmod = 41.924d9; mat%iyield = NP_YLD_HOSFORD
      mat%ielast = NP_ELAST_GUO; mat%a0 = [1.0d0, 0.0d0, 0.0d0]; mat%vv = 0.0d0
      mat%y0 = 52.28d6; mat%ysat = 62.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   case ('pc')    ! as isomap.f90
      mat%mu = 848.0d6; mat%kmod = 3.25d9; mat%iyield = NP_YLD_GHORBELF; mat%rsig = 1.387d0; mat%st0 = 62.0d6
      mat%y0 = 62.0d6; mat%ysat = 78.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   case default
      stop 'law: pc | pla'
   end select
   call np_check_material(mat, info, arg)
   if (info /= NP_OK) then; write (*, *) 'material: ', trim(arg); stop 1; end if
   m1 = mat; m1%nn = 1; m1%eps_theta = 0.0d0
   cp0 = [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]
!$omp parallel do schedule(dynamic, 1) private(ip, i, j, a, b, e, row, s, d, cpo, ko, psi, ratio, g0, gmin, fmin, g1, g2, &
!$omp&   kdbl, info, inf1, st)
   do ip = 1, n*n
      i = (ip - 1)/n + 1; j = mod(ip - 1, n) + 1
      a = -amax + 2.0d0*amax*real(j - 1, dp)/real(n - 1, dp)
      b = -amax + 2.0d0*amax*real(i - 1, dp)/real(n - 1, dp)
      e = 0.0d0; e(1) = a; e(2) = b
      row = 0.0d0; row(1) = a; row(2) = b
      call np_raydiag(mat, e, cp0, 0.0d0, ratio, g0, gmin, fmin, g1, g2, kdbl, info)
      row(3) = ratio; row(4) = g0; row(5) = gmin; row(6) = fmin; row(7) = g1; row(8) = g2; row(9) = real(kdbl, dp)
      row(10) = real(info, dp); row(11) = -1.0d0
      if (info == 0) then
         call np_update_alg1(m1, e, cp0, 0.0d0, s, d, cpo, ko, psi, inf1, st)
         row(11) = real(inf1, dp)
      end if
      rows(:, ip) = row
   end do
!$omp end parallel do
   open (newunit=uout, file=trim(outf), status='replace')
   write (uout, '(A,A,A,I0,A,ES12.4)') '# raydiag ', trim(law), ' n ', n, ' amax ', amax
   write (uout, '(A)') '# E11 E22 ratio g0 gmin fmin g1 g2 kdbl info status1'
   do ip = 1, n*n
      write (uout, '(11ES14.6)') rows(:, ip)
   end do
   close (uout)
   write (*, '(A,A,A)') ' raydiag ', trim(law), ': done'
end program raydiag
