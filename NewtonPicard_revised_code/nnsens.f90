! nnsens -- sensitivity to the strain of the update with exactly N_n Picard updates (eps_theta = 0) (Section 8.3).
! PLA of the framed beam (matrix energy, a_0 = e_1, Hosford n = 8, Voce 40.92 -> 48 MPa, b = 20), single
! increments from the undeformed state along fixed strain directions scaled to a trial ratio r. For a strain
! perturbation dE = delta |E| p (central differences):
!     amp = |S(E + dE) - S(E - dE)| / |2 D dE|,   D = the tangent returned by the update at E
! (1 = the tangent predicts the change of the stress). Printed per direction and ratio for N_n = 1, 5, 10 with the
! Picard residual |N_new - N| after the last update, and for the converged update (eps_theta = 1e-10, cap 50, no
! substepping; -1 = one of the three updates did not converge).
!   NNS_PERT = 1 (default), 2: generic perturbation directions; 3: p = (0, 1, -1, 0, 0, 0)/sqrt(2), which breaks the
!   transverse symmetry of the axial direction 1;  NNS_DELTA = delta (default 1e-7; the result is the same for 1e-6
!   to 1e-10: the amplification is that of the derivative of the map with N_n updates).
! The paper: direction (1, -0.45, -0.45), NNS_PERT = 3, trial ratio 2.
program nnsens
   use npicard_plasticity
   implicit none
   type(np_material) :: m0, m
   type(np_state) :: st, st1, st2
   real(dp) :: e(6), ep(6), em(6), de(6), s(6), sp(6), sm(6), d(6,6), dd(6,6), cp(6), cpo(6), ko, psi, dirs(6,4), pert(6)
   real(dp) :: t, r, lo, hi, amp(4), res(4), rat
   integer :: id, ir, k, info, i1, i2, in, nns(4)
   real(dp), parameter :: RR(16) = [1.05d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.8d0, 2.0d0, 2.5d0, 3.0d0, 4.0d0, 6.0d0, 8.0d0, 12.0d0, 20.0d0]
   logical :: ok
   character(len=16) :: nm(4)
   character(len=32) :: arg
   real(dp) :: delta = 1.0d-7
   nm = [character(len=16) :: 'uniaxial 1', 'plane strain', 'shear 12', 'biaxial']
   dirs(:,1) = [1.0d0, -0.45d0, -0.45d0, 0d0, 0d0, 0d0]
   dirs(:,2) = [1.0d0, 0.0d0, -0.9d0, 0d0, 0d0, 0d0]
   dirs(:,3) = [0.02d0, -0.01d0, -0.01d0, 1.0d0, 0d0, 0d0]
   dirs(:,4) = [1.0d0, 0.7d0, -1.5d0, 0.1d0, 0d0, 0d0]
   pert = [0.3d0, -0.5d0, 0.1d0, 0.6d0, -0.4d0, 0.2d0]
   call get_environment_variable('NNS_PERT', arg)
   if (trim(arg) == '2') pert = [-0.2d0, 0.1d0, 0.7d0, -0.3d0, 0.5d0, 0.4d0]
   if (trim(arg) == '3') pert = [0.0d0, 1.0d0, -1.0d0, 0.0d0, 0.0d0, 0.0d0]
   pert = pert/sqrt(sum(pert**2))
   call get_environment_variable('NNS_DELTA', arg)
   if (len_trim(arg) > 0) read (arg, *) delta
   nns = [1, 5, 10, 50]
   cp = [1d0,1d0,1d0,0d0,0d0,0d0]
   m0 = np_material()
   m0%mu = 1.1508d9; m0%kmod = 41.924d9; m0%iyield = NP_YLD_HOSFORD
   m0%ielast = NP_ELAST_GUO; m0%a0 = [1.0d0, 0.0d0, 0.0d0]; m0%vv = 0.0d0
   m0%y0 = 40.92d6; m0%ysat = 48.0d6; m0%bvoce = 20.0d0; m0%nh = 1.0d0
   write (*, '(A)') ' direction        ratio |  N_n = 1: amp, |dN| |  N_n = 5: amp, |dN| | N_n = 10: amp, |dN| | converged (50, 1e-10): amp, its, info'
   do id = 1, 4
      do ir = 1, 16
         r = RR(ir)
         ! scale of the direction with trial ratio r (bisection)
         lo = 0.0d0; hi = 1.0d0
         do k = 1, 60
            t = 0.5d0*(lo + hi)
            call np_trial(m0, t*dirs(:, id), cp, rat, ok)
            if (.not. ok .or. rat > r) then; hi = t; else; lo = t; end if
         end do
         e = t*dirs(:, id); de = delta*sqrt(sum(e**2))*pert
         ep = e + de; em = e - de
         do in = 1, 4
            m = m0; m%nn = nns(in); m%eps_theta = merge(1.0d-10, 0.0d0, in == 4)
            call np_update_alg1(m, e, cp, 0.0d0, s, d, cpo, ko, psi, info, st=st, noexit=.true.)
            call np_update_alg1(m, ep, cp, 0.0d0, sp, dd, cpo, ko, psi, i1, st=st1, noexit=.true.)
            call np_update_alg1(m, em, cp, 0.0d0, sm, dd, cpo, ko, psi, i2, st=st2, noexit=.true.)
            amp(in) = -1.0d0; res(in) = st%dn(max(1, min(st%its, size(st%dn))))
            if (info == NP_OK .and. i1 == NP_OK .and. i2 == NP_OK) amp(in) = sqrt(sum((sp - sm)**2))/sqrt(sum((2.0d0*matmul(d, de))**2))
            if (in == 4) write (*, '(1X,A14,F7.2,3(A,ES9.2,ES9.1),A,ES9.2,I4,I3)') nm(id), r, ' |', amp(1), res(1), ' |', amp(2), res(2), &
               ' |', amp(3), res(3), ' |', amp(4), st%its, info
         end do
      end do
   end do
contains
   subroutine np_trial(mm, ee, cc, ratio, okk)
      type(np_material), intent(in) :: mm
      real(dp), intent(in) :: ee(6), cc(6)
      real(dp), intent(out) :: ratio
      logical, intent(out) :: okk
      type(np_material) :: m1
      type(np_state) :: s1
      real(dp) :: ss(6), d1(6,6), c1(6), k1, p1
      integer :: inf
      m1 = mm; m1%nn = 1; m1%eps_theta = 0.0d0
      call np_update_alg1(m1, ee, cc, 0.0d0, ss, d1, c1, k1, p1, inf, st=s1, noexit=.true.)
      ratio = s1%ratio; okk = (inf == NP_OK) .or. s1%ratio > 0.0d0
      if (s1%ratio <= 0.0d0) ratio = 0.5d0     ! elastic: below yield
   end subroutine np_trial
end program nnsens
