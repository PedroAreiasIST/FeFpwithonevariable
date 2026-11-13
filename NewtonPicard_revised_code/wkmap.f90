!=============================================================================
! wkmap -- the stored elastic energy W and the hardening variable kappa as functions of E11 and E22 (Section 8.1):
! the states of the iso-error maps of the polycarbonate and of PLA (E11, E22 in [-0.25, 0.25] on 81 x 81 points,
! the other components zero), each reached from the undeformed state along the proportional path in nstep
! converged increments (Algorithm 1 with substepping, eps_theta 1e-10). Materials as isomap.f90.
! usage: wkmap <pc|pla> <nstep> <out>            columns: a b W[Pa] kappa      (the paper: nstep = 200)
!=============================================================================
program wkmap
   use npicard_plasticity
   use omp_lib
   implicit none
   type(np_material) :: mat
   type(np_state) :: st
   character(len=16) :: law
   character(len=256) :: outf, arg
   integer, parameter :: n = 81
   integer :: nstep, ip, i, j, k, info, nsub, uout, nbad
   real(dp), parameter :: amax = 0.25d0
   real(dp) :: a, b, e0(6), e1(6), cp(6), cpn(6), kap, kn, s(6), d(6, 6), psi, rows(4, n*n)
   call get_command_argument(1, law)
   call get_command_argument(2, arg); read (arg, *) nstep
   call get_command_argument(3, outf)
   mat%ielast = NP_ELAST_PAPER
   select case (trim(law))
   case ('pla')
      mat%mu = 1.1508d9; mat%kmod = 41.924d9; mat%iyield = NP_YLD_HOSFORD
      mat%ielast = NP_ELAST_GUO; mat%a0 = [1.0d0, 0.0d0, 0.0d0]; mat%vv = 0.0d0
      mat%y0 = 52.28d6; mat%ysat = 62.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   case ('pc')
      mat%mu = 848.0d6; mat%kmod = 3.25d9; mat%iyield = NP_YLD_GHORBELF; mat%rsig = 1.387d0; mat%st0 = 62.0d6
      mat%y0 = 62.0d6; mat%ysat = 78.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   case default
      stop 'law: pc | pla'
   end select
   mat%nn = 50; mat%eps_theta = 1.0d-10
   call np_check_material(mat, info, arg)
   if (info /= NP_OK) then; write (*, *) 'material: ', trim(arg); stop 1; end if
   nbad = 0
!$omp parallel do schedule(dynamic, 1) private(ip, i, j, a, b, e0, e1, cp, cpn, kap, kn, s, d, psi, k, info, nsub, st) &
!$omp&   reduction(+:nbad)
   do ip = 1, n*n
      i = (ip - 1)/n + 1; j = mod(ip - 1, n) + 1
      a = -amax + 2.0d0*amax*real(j - 1, dp)/real(n - 1, dp)
      b = -amax + 2.0d0*amax*real(i - 1, dp)/real(n - 1, dp)
      e0 = 0.0d0; cp = [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]; kap = 0.0d0; psi = 0.0d0
      do k = 1, nstep
         e1 = 0.0d0; e1(1) = a*real(k, dp)/real(nstep, dp); e1(2) = b*real(k, dp)/real(nstep, dp)
         call np_update_sub(mat, e0, e1, cp, kap, s, d, cpn, kn, psi, info, nsub, 10, st=st)
         if (info /= NP_OK) nbad = nbad + 1
         e0 = e1; cp = cpn; kap = kn
      end do
      rows(:, ip) = [a, b, psi, kap]
   end do
!$omp end parallel do
   open (newunit=uout, file=trim(outf), status='replace')
   write (uout, '(A,A,A,I0,A,I0)') '# wkmap ', trim(law), ' nstep ', nstep, ' increments without a converged update ', nbad
   write (uout, '(A)') '# a b W kappa'
   do ip = 1, n*n
      write (uout, '(4ES18.9)') rows(:, ip)
   end do
   close (uout)
   write (*, '(A,A,A,I0)') ' wkmap ', trim(law), ': done, failed increments ', nbad
end program wkmap
