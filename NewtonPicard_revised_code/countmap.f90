!=============================================================================
! countmap -- iteration counts of Algorithm 1 (Section 6.4, "Hardening subproblem"), measured on the states of the
! iso-error maps that converge in one increment (eps_theta = 1e-10, cap 50, no substepping), with the material data
! of isomap.f90 (pla | pc).
!
! usage: countmap <law> <n> <amax>        E11, E22 in [-amax, amax], other components zero (the paper: 81 0.25)
! Prints, over the plastic states with a converged single increment: Picard updates per increment (st%its), scalar
! Newton iterations on gamma per Picard update (st%newton/st%its), evaluations of f_gamma per increment and per
! update (st%nfeval: the bracketing and the Newton iterations; each needs the hardening root kappa(gamma)), and the
! Newton iterations of the hardening equation per solution (nkap/nksolve of np_update_alg1).
!=============================================================================
program countmap
   use npicard_plasticity
   implicit none
   type(np_material) :: mat
   type(np_state) :: st
   character(len=16) :: law
   character(len=64) :: arg
   integer :: n, i, j, info, m, nconv, nplast, nkap, nks
   real(dp), allocatable :: kpu(:)
   real(dp) :: amax, e(6), cp(6), s(6), d(6, 6), cpo(6), ko, psi
   integer, allocatable :: its(:), nwt(:), nfe(:), kit(:)
   real(dp), allocatable :: nwpu(:), nfpu(:)
   call get_command_argument(1, law)
   call get_command_argument(2, arg); read (arg, *) n
   call get_command_argument(3, arg); read (arg, *) amax
   mat%ielast = NP_ELAST_PAPER
   select case (trim(law))
   case ('pla')
      mat%mu = 1.1508d9; mat%kmod = 41.924d9; mat%iyield = NP_YLD_HOSFORD
      mat%ielast = NP_ELAST_GUO; mat%a0 = [1.0d0, 0.0d0, 0.0d0]; mat%vv = 0.0d0   ! the matrix energy of Section 8.1.2, a_0 = e_1
      mat%y0 = 52.28d6; mat%ysat = 62.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   case ('pc')   ! the criterion of Section 8.1.1 (notebook yld_ghorbel.nb), Voce 62 -> 78 MPa
      mat%mu = 848.0d6; mat%kmod = 3.25d9; mat%iyield = NP_YLD_GHORBELF; mat%rsig = 1.387d0; mat%st0 = 62.0d6
      mat%y0 = 62.0d6; mat%ysat = 78.0d6; mat%bvoce = 14.0d0; mat%nh = 1.3d0
   case default
      stop 'law: pla | pc'
   end select
   mat%nn = 50; mat%eps_theta = 1.0d-10
   allocate (its(n*n), nwt(n*n), nfe(n*n), kit(n*n), nwpu(n*n), nfpu(n*n), kpu(n*n))
   nconv = 0; nplast = 0
   do i = 1, n
      do j = 1, n
         e = 0.0d0
         e(1) = -amax + 2.0d0*amax*real(j - 1, dp)/real(n - 1, dp)
         e(2) = -amax + 2.0d0*amax*real(i - 1, dp)/real(n - 1, dp)
         cp = [1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]
         call np_update_alg1(mat, e, cp, 0.0d0, s, d, cpo, ko, psi, info, st=st, nkap=nkap, nksolve=nks, noexit=.true.)
         if (.not. st%plastic) cycle
         nplast = nplast + 1
         if (info /= NP_OK .or. st%its >= mat%nn) cycle
         nconv = nconv + 1
         its(nconv) = st%its; nwt(nconv) = st%newton; nfe(nconv) = st%nfeval
         nwpu(nconv) = real(st%newton, dp)/real(st%its, dp); nfpu(nconv) = real(st%nfeval, dp)/real(st%its, dp)
         kpu(nconv) = real(nkap, dp)/real(max(nks, 1), dp)
      end do
   end do
   write (*, '(A,A,A,I0,A,I0,A,I0,A)') ' countmap ', trim(law), ': ', n*n, ' states, ', nplast, ' plastic, ', nconv, &
      ' converge in one step'
   m = nconv
   call report('Picard updates per step       ', real(its(1:m), dp))
   call report('Newton iterations per update  ', nwpu(1:m))
   call report('f_gamma evaluations per step  ', real(nfe(1:m), dp))
   call report('f_gamma evaluations per update', nfpu(1:m))
   call report('hardening Newton its per solve', kpu(1:m))
contains
   subroutine report(name, x)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: v(:)
      integer :: k
      if (size(x) == 0) return
      v = x
      call sort(v)
      k = size(v)
      write (*, '(2X,A,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2)') name, ': min ', v(1), '  median ', v((k + 1)/2), &
         '  mean ', sum(v)/k, '  95th ', v(max(1, nint(0.95d0*k))), '  max ', v(k)
   end subroutine report

   subroutine sort(v)
      real(dp), intent(inout) :: v(:)
      integer :: a, b
      real(dp) :: t
      do a = 2, size(v)
         t = v(a); b = a - 1
         do while (b >= 1)
            if (v(b) <= t) exit
            v(b + 1) = v(b); b = b - 1
         end do
         v(b + 1) = t
      end do
   end subroutine sort
end program countmap
