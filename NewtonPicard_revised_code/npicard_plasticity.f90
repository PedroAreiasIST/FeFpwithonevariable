!=============================================================================
! NPICARD_PLASTICITY -- the Newton/Picard return mapping of the paper
!   "A Newton/Picard semi-implicit algorithm for finite-strain plasticity"
!   (P. Areias, T. Rabczuk, N. Silvestre, M.F. Vaz), Computational Mechanics, revised version:
! finite-strain plasticity with a separable yield function f = seq(Ms) - Y(kappa), multiplicative split,
! exponential integration of the flow rule, symmetrized Mandel stress Ms, power-balance hardening.
!
! MODEL
!   kinematics   Ce^tr = U C U with U = sqrt(Cp_n^-1) (the stored history is the symmetric Cp^-1);
!                Ce = exp(-gamma N) Ce^tr exp(-gamma N); Cp_{n+1}^-1 = U exp(-2 gamma N) U;
!                S = U exp(-gamma N) Se exp(-gamma N) U (second Piola-Kirchhoff stress).
!   elastic law  ielast = 1: W = mu/2 (tr Ce - 3 - ln det Ce) + K/2 (sqrt(det Ce) - 1)^2
!                ielast = 2: the transversely isotropic energy of the PLA model of Section 8.1.2 with the
!                matrix constants mu, kmod, the direction a0 and the void fraction vv, which is constant
!                during the update (generated code, npicard_guo_kernel.f90).
!                Se = 2 dW/dCe, Ms = (Ce Se + Se Ce)/2, N = d seq/d Ms (not normalized).
!   yield        1 von Mises sqrt(3 J2);
!                2 Hosford with exponent 8 in the polynomial form (129 J2^4 - 324 J2 J3^2)^(1/8);
!                3 Barlat Yld2004-18p: X' = C' dev Ms, X'' = C'' dev Ms, phi = sum_ij |X'_i - X''_j|^a over the
!                principal values, seq = (phi/4)^(1/a). For an even exponent a the double sum is the trace
!                polynomial phi = sum_k binom(a, k) (-1)^k tr(X'^k) tr(X''^(a-k)) (tr X^0 = 3), so seq, its gradient
!                and its Hessian need no eigenvalues and stay smooth where principal values coincide (the
!                equibiaxial state). Coefficients: np_barlat_set (0 isotropic, 1 AA2090-T3);
!                4 the criterion of the polycarbonate of Section 8.1.1: (q/r)(1 + (r - 1) h(I1)), q = sqrt(3 J2),
!                h = (1 + I1/sqrt(I1^2 + l^2))/2, l = 0.075 st0 (generated code, notebook yld_ghorbel.nb).
!   hardening    Y(k) = Ysat - (Ysat - Y0) exp(-b k^nh) + H k; (k - kn) Y(k) = gamma Ms:N (power balance).
!
! ROUTINES
!   np_update_alg1   Algorithm 1 of the paper: the scalar iteration on gamma at a frozen direction N (bracket by
!                    doubling, safeguarded Newton, the hardening equation solved at every evaluation) inside the
!                    Picard update N <- d seq/d Ms, and the consistent tangent at the returned state (Section 5).
!   np_update_sub    the same update with adaptive substepping of the strain increment (Section 6.5) and the
!                    exact tangent of the substepped increment (np_update_dh: derivatives with respect to the
!                    history).
!   np_update_base   the baseline of Section 8.2: Newton iteration on (Ce, gamma) with the yield stress frozen
!                    and kappa updated afterwards (mat%classical = .true.).
!   np_raydiag       diagnostic of Section 8.2: the residual of the scalar iteration along the first ray.
!   np_yield, np_elastic, np_barlat_set, np_check_material.
!
! CONVENTIONS
!   Interface Voigt order {11,22,33,12,13,23}. Strain in: engineering Green-Lagrange {E11,E22,E33,2E12,2E13,2E23},
!   so C = I + 2E reads C_ii = 1 + 2 E_i, C_ij = E_(ij); stress out: tensor components of the second
!   Piola-Kirchhoff stress; dsde is dS/dE with respect to the engineering strain; Cp^-1 in tensor components.
!   Internally every symmetric tensor is a Mandel 6-vector (shears times sqrt 2: an orthonormal basis of Sym) and
!   every linear map on Sym a 6x6 matrix in that basis, so transposes are adjoints and dot products are double
!   contractions.
!
! STATUS (info): 0 a stress with the Picard iteration converged (or after its nn updates when eps_theta = 0);
!   1 invalid input (material, C or Cp^-1 not positive definite); 2 the scalar iteration did not converge;
!   3 no bracket for the scalar iteration; 4 a non-finite number; 5 a stress, with the Picard iteration not
!   converged within nn updates (eps_theta > 0). On info = 1..4 every output is zero except cpinv = cpinv_n,
!   kappa = kappa_n: the increment is substepped (np_update_sub) or the caller cuts the load step.
!
! No dependencies besides npicard_guo_kernel; no saved state (thread-safe).
!=============================================================================
module npicard_plasticity
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use npicard_guo_kernel, only: npg_modifiedse
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: NP_YLD_VONMISES = 1, NP_YLD_HOSFORD = 2, NP_YLD_BARLAT = 3, NP_YLD_GHORBELF = 4
   integer, parameter, public :: NP_BARLAT_MAXA = 16
   integer, parameter, public :: NP_ELAST_PAPER = 1, NP_ELAST_GUO = 2
   integer, parameter, public :: NP_OK = 0, NP_ERR_INPUT = 1, NP_ERR_DIVERGED = 2, NP_ERR_NOROOT = 3, &
                                 NP_ERR_NONFINITE = 4, NP_UNCONVERGED = 5
   integer, parameter, public :: NP_MAXNN = 400
   real(dp), parameter :: SQ2 = 1.4142135623730951_dp
   real(dp), parameter :: HOS_A = 129.0_dp, HOS_B = 324.0_dp   ! (129 J2^4 - 324 J2 J3^2)^(1/8) = Hosford n = 8
   real(dp), parameter :: GOLD = 0.6180339887498949_dp         ! (sqrt 5 - 1)/2
   real(dp), parameter :: EXPARG = 700.0_dp                    ! largest |argument| handed to exp (exp(709.78) = huge)
   real(dp), parameter :: I6(6) = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   real(dp), parameter :: SFAC(6) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/SQ2, 1.0_dp/SQ2, 1.0_dp/SQ2]  ! Mandel -> Voigt stress
   real(dp), parameter :: CFAC(6) = [2.0_dp, 2.0_dp, 2.0_dp, SQ2, SQ2, SQ2]                      ! dC(Mandel)/dE(engineering)

   type, public :: np_material
      real(dp) :: mu = 0.0_dp, kmod = 0.0_dp          ! shear and bulk moduli [stress]
      integer  :: ielast = NP_ELAST_PAPER
      integer  :: iyield = NP_YLD_VONMISES
      real(dp) :: nhos = 8.0_dp                       ! Hosford exponent (8 is the only one implemented)
      real(dp) :: rsig = 1.0_dp, st0 = 0.0_dp         ! polycarbonate criterion: r and sigma_t0
      real(dp) :: y0 = 0.0_dp, ysat = 0.0_dp, bvoce = 0.0_dp, nh = 1.0_dp   ! Voce
      real(dp) :: hlin = 0.0_dp                       ! linear hardening modulus H [stress] (0 = Voce only)
      integer  :: nn = 1                              ! cap of the Picard updates N_n
      real(dp) :: eps_theta = 1.0e-12_dp              ! Picard tolerance on ||N_new - N|| (0 = exactly nn updates)
      real(dp) :: eps_f = 1.0e-10_dp                  ! tolerance of the scalar iteration, relative to Y(kappa_n)
      logical  :: classical = .false.                 ! the baseline of Section 8.2 (np_update_base): Newton iteration on
                                                      ! (Ce, gamma) with N(Ce), the yield stress frozen at Y(kappa_n)
                                                      ! and kappa updated explicitly (staggered hardening)
      real(dp) :: vv = 0.0_dp                         ! PLA energy (ielast 2): void fraction, constant in the update
      real(dp) :: a0(3) = [1.0_dp, 0.0_dp, 0.0_dp]    ! PLA energy: unit deposition direction
      integer  :: abar = 8                            ! Barlat Yld2004-18p: even exponent a (2..NP_BARLAT_MAXA)
      real(dp) :: bl1(6, 6) = 0.0_dp, bl2(6, 6) = 0.0_dp   ! C' P and C'' P in the Mandel basis (P = deviatoric projector)
   end type np_material

   type, public :: np_state
      logical  :: plastic = .false.
      real(dp) :: seqtr = 0.0_dp, ytr = 0.0_dp, ratio = 0.0_dp, ftr = 0.0_dp    ! trial state; ratio = seq^tr/Y
      real(dp) :: gamma = 0.0_dp, kappa = 0.0_dp, seq = 0.0_dp, f = 0.0_dp
      integer  :: its = 0, newton = 0, nfeval = 0     ! Picard updates, scalar Newton iterations, evaluations of f
      integer  :: noroot = 0, noroot_last = 0         ! rays without a bracket
      logical  :: unconverged = .false.               ! cap nn reached with eps_theta > 0
      real(dp) :: dn(NP_MAXNN) = 0.0_dp               ! ||N_new - N|| of each Picard update
      real(dp) :: n6(6) = 0.0_dp, m6(6) = 0.0_dp, ce6(6) = 0.0_dp, cetr6(6) = 0.0_dp, se6(6) = 0.0_dp
      ! diagnostics of the tangent (Section 5.2): the denominator of the multiplier increment, its value at
      ! gamma = 0, and det(I - R_N : H : T)
      real(dp) :: dg_den = 0.0_dp, dg_det = 0.0_dp, dg_h0 = 0.0_dp
   end type np_state

   ! the frozen direction of a ray: N and its eigendecomposition (exp(-gamma N) = Q exp(-gamma lam) Q^T)
   type :: ray_t
      real(dp) :: n3(3, 3), n6(6), lam(3), q(3, 3)
   end type ray_t
   ! the state along a ray at gamma, with the operators of the tangent
   type :: eval_t
      real(dp) :: gam, e3(3, 3), ce3(3, 3), se3(3, 3), m3(3, 3), w, seq, kappa, y, dy, f, df
      real(dp) :: n6(6), m6(6), phat6(6), dk, alpha, eta, prod
      real(dp) :: c66(6, 6), t66(6, 6), h66(6, 6)
   end type eval_t

   ! Derivatives of the single update (C, A = Cp_n^-1, kappa_n) -> (S, A' = Cp_n+1^-1, kappa') in the Mandel basis
   ! (np_update_dh): the blocks of the consistent tangent of a substepped increment (np_update_sub)
   type, public :: np_dhist
      logical  :: ok = .false.
      real(dp) :: s_c(6, 6) = 0.0_dp, s_a(6, 6) = 0.0_dp, s_k(6) = 0.0_dp   ! dS/dC, dS/dA, dS/dkappa_n
      real(dp) :: a_c(6, 6) = 0.0_dp, a_a(6, 6) = 0.0_dp, a_k(6) = 0.0_dp   ! dA'/dC, dA'/dA, dA'/dkappa_n
      real(dp) :: k_c(6) = 0.0_dp, k_a(6) = 0.0_dp, k_k = 0.0_dp            ! dkappa'/dC, dkappa'/dA, dkappa'/dkappa_n
   end type np_dhist

   public :: np_update_alg1, np_update_sub, np_update_dh, np_update_base
   public :: np_check_material, np_yield, np_elastic, np_barlat_set, np_raydiag
contains

   !==========================================================================
   ! Mandel helpers
   !==========================================================================
   pure function to6(a) result(v)
      real(dp), intent(in) :: a(3, 3)
      real(dp) :: v(6)
      v(1) = a(1, 1); v(2) = a(2, 2); v(3) = a(3, 3)
      v(4) = 0.5_dp*SQ2*(a(1, 2) + a(2, 1))
      v(5) = 0.5_dp*SQ2*(a(1, 3) + a(3, 1))
      v(6) = 0.5_dp*SQ2*(a(2, 3) + a(3, 2))
   end function to6

   pure function from6(v) result(a)
      real(dp), intent(in) :: v(6)
      real(dp) :: a(3, 3)
      a(1, 1) = v(1); a(2, 2) = v(2); a(3, 3) = v(3)
      a(1, 2) = v(4)/SQ2; a(2, 1) = a(1, 2)
      a(1, 3) = v(5)/SQ2; a(3, 1) = a(1, 3)
      a(2, 3) = v(6)/SQ2; a(3, 2) = a(2, 3)
   end function from6

   pure function basis(k) result(b)
      integer, intent(in) :: k
      real(dp) :: b(3, 3), v(6)
      v = 0.0_dp; v(k) = 1.0_dp
      b = from6(v)
   end function basis

   pure function eye3() result(a)
      real(dp) :: a(3, 3)
      a = 0.0_dp; a(1, 1) = 1.0_dp; a(2, 2) = 1.0_dp; a(3, 3) = 1.0_dp
   end function eye3

   pure function eye6() result(a)
      real(dp) :: a(6, 6)
      integer :: i
      a = 0.0_dp
      do i = 1, 6
         a(i, i) = 1.0_dp
      end do
   end function eye6

   ! X -> A X A (A symmetric)
   pure function sandwich66(a) result(s)
      real(dp), intent(in) :: a(3, 3)
      real(dp) :: s(6, 6)
      integer :: k
      do k = 1, 6
         s(:, k) = to6(matmul(a, matmul(basis(k), a)))
      end do
   end function sandwich66

   ! X -> X A B + B A X (A, B symmetric; the result is symmetric)
   pure function sandwich2_66(a, b) result(s)
      real(dp), intent(in) :: a(3, 3), b(3, 3)
      real(dp) :: s(6, 6), x(3, 3), y(3, 3)
      integer :: k
      do k = 1, 6
         x = basis(k)
         y = matmul(x, matmul(a, b))
         s(:, k) = to6(y + transpose(y))
      end do
   end function sandwich2_66

   ! X -> A X + X A (A symmetric)
   pure function prodsym66(a) result(s)
      real(dp), intent(in) :: a(3, 3)
      real(dp) :: s(6, 6), x(3, 3)
      integer :: k
      do k = 1, 6
         x = basis(k)
         s(:, k) = to6(matmul(a, x) + matmul(x, a))
      end do
   end function prodsym66

   pure function outer6(a, b) result(s)
      real(dp), intent(in) :: a(6), b(6)
      real(dp) :: s(6, 6)
      integer :: i, j
      do j = 1, 6
         do i = 1, 6
            s(i, j) = a(i)*b(j)
         end do
      end do
   end function outer6

   pure function det3(a) result(d)
      real(dp), intent(in) :: a(3, 3)
      real(dp) :: d
      d = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
          + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
   end function det3

   pure function inv3(a, d) result(b)
      real(dp), intent(in) :: a(3, 3), d
      real(dp) :: b(3, 3)
      b(1, 1) = (a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2))/d
      b(1, 2) = (a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3))/d
      b(1, 3) = (a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))/d
      b(2, 1) = (a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3))/d
      b(2, 2) = (a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1))/d
      b(2, 3) = (a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3))/d
      b(3, 1) = (a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))/d
      b(3, 2) = (a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2))/d
      b(3, 3) = (a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))/d
   end function inv3

   ! 6x6 inverse by Gauss-Jordan with partial pivoting; ok = .false. on a zero pivot
   pure subroutine inv6(a, b, ok)
      real(dp), intent(in) :: a(6, 6)
      real(dp), intent(out) :: b(6, 6)
      logical, intent(out) :: ok
      real(dp) :: w(6, 12), t, piv
      integer :: i, j, k, ip
      w(:, 1:6) = a; w(:, 7:12) = eye6()
      ok = .true.
      do k = 1, 6
         ip = k; piv = abs(w(k, k))
         do i = k + 1, 6
            if (abs(w(i, k)) > piv) then
               piv = abs(w(i, k)); ip = i
            end if
         end do
         if (.not. (piv > 0.0_dp)) then
            ok = .false.; b = 0.0_dp; return
         end if
         if (ip /= k) then
            do j = 1, 12
               t = w(k, j); w(k, j) = w(ip, j); w(ip, j) = t
            end do
         end if
         t = 1.0_dp/w(k, k)
         w(k, :) = w(k, :)*t
         do i = 1, 6
            if (i /= k) then
               t = w(i, k)
               if (t /= 0.0_dp) w(i, :) = w(i, :) - t*w(k, :)
            end if
         end do
      end do
      b = w(:, 7:12)
   end subroutine inv6

   !==========================================================================
   ! Symmetric 3x3 eigendecomposition (cyclic Jacobi, to round-off): a = q diag(lam) q^T
   !==========================================================================
   pure subroutine sym_eig3(a, lam, q)
      real(dp), intent(in) :: a(3, 3)
      real(dp), intent(out) :: lam(3), q(3, 3)
      real(dp) :: b(3, 3), off, dia, theta, t, c, s, tau, g, h
      integer :: sweep, p, r, m
      b = 0.5_dp*(a + transpose(a))
      q = eye3()
      do sweep = 1, 60
         off = b(1, 2)**2 + b(1, 3)**2 + b(2, 3)**2
         dia = b(1, 1)**2 + b(2, 2)**2 + b(3, 3)**2
         if (off == 0.0_dp .or. off <= 1.0e-36_dp*dia) exit
         do p = 1, 2
            do m = p + 1, 3
               if (b(p, m) == 0.0_dp) cycle
               ! |b_pm| below the half diagonal gap over sqrt(huge): theta^2 overflows and the rotation is the identity
               ! to round-off (t = 1/Inf = 0 in a fast-math build, a trap under -fpe0), so b_pm is dropped
               if (abs(b(p, m)) < abs(0.5_dp*(b(m, m) - b(p, p)))/sqrt(huge(1.0_dp))) then
                  b(p, m) = 0.0_dp; b(m, p) = 0.0_dp; cycle
               end if
               theta = (b(m, m) - b(p, p))/(2.0_dp*b(p, m))
               t = sign(1.0_dp, theta)/(abs(theta) + sqrt(theta*theta + 1.0_dp))
               c = 1.0_dp/sqrt(t*t + 1.0_dp)
               s = t*c
               tau = s/(1.0_dp + c)
               h = t*b(p, m)
               b(p, p) = b(p, p) - h
               b(m, m) = b(m, m) + h
               b(p, m) = 0.0_dp; b(m, p) = 0.0_dp
               do r = 1, 3
                  if (r == p .or. r == m) cycle
                  g = b(r, p); h = b(r, m)
                  b(r, p) = g - s*(h + g*tau); b(p, r) = b(r, p)
                  b(r, m) = h + s*(g - h*tau); b(m, r) = b(r, m)
               end do
               do r = 1, 3
                  g = q(r, p); h = q(r, m)
                  q(r, p) = g - s*(h + g*tau)
                  q(r, m) = h + s*(g - h*tau)
               end do
            end do
         end do
      end do
      lam(1) = b(1, 1); lam(2) = b(2, 2); lam(3) = b(3, 3)
   end subroutine sym_eig3

   ! exp(A) for symmetric A, from its eigendecomposition, at the scale s: exp(s A)
   pure function expm_eig(lam, q, s) result(e)
      real(dp), intent(in) :: lam(3), q(3, 3), s
      real(dp) :: e(3, 3), d(3, 3)
      integer :: i
      d = 0.0_dp
      do i = 1, 3
         d(i, i) = exp(s*lam(i))
      end do
      e = matmul(q, matmul(d, transpose(q)))
   end function expm_eig

   ! dense linear solve by Gaussian elimination with partial pivoting (the 7 x 7 system of the classical baseline)
   pure subroutine solve_dense(n, a, b, x, ok)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n), b(n)
      real(dp), intent(out) :: x(n)
      logical, intent(out) :: ok
      real(dp) :: m(n, n + 1), t(n + 1), piv, sc
      integer :: i, k, p
      m(:, 1:n) = a; m(:, n + 1) = b
      sc = maxval(abs(a))
      ok = sc > 0.0_dp
      x = 0.0_dp
      if (.not. ok) return
      do k = 1, n
         p = k - 1 + maxloc(abs(m(k:n, k)), 1)
         piv = m(p, k)
         if (abs(piv) <= 1.0e-14_dp*sc) then
            ok = .false.; return
         end if
         if (p /= k) then
            t = m(k, :); m(k, :) = m(p, :); m(p, :) = t
         end if
         do i = k + 1, n
            m(i, k:n + 1) = m(i, k:n + 1) - (m(i, k)/m(k, k))*m(k, k:n + 1)
         end do
      end do
      do i = n, 1, -1
         x(i) = (m(i, n + 1) - dot_product(m(i, i + 1:n), x(i + 1:n)))/m(i, i)
      end do
   end subroutine solve_dense

   ! Frechet derivative of exp at the symmetric A = q diag(alam) q^T, as a Mandel 6x6:
   ! L[X] = q ((q^T X q) o Phi) q^T, Phi_ij = (e^ai - e^aj)/(ai - aj) (Daleckii-Krein)
   pure function expm_frechet66(alam, q) result(l)
      real(dp), intent(in) :: alam(3), q(3, 3)
      real(dp) :: l(6, 6), phi(3, 3), y(3, 3), d, m, sc
      integer :: i, j, k
      do i = 1, 3
         do j = 1, 3
            if (i == j) then
               phi(i, i) = exp(alam(i))
            else
               d = 0.5_dp*(alam(i) - alam(j)); m = 0.5_dp*(alam(i) + alam(j))
               if (abs(d) > 1.0e-3_dp) then
                  sc = sinh(d)/d
               else
                  sc = 1.0_dp + d*d/6.0_dp + d**4/120.0_dp
               end if
               phi(i, j) = exp(m)*sc
            end if
         end do
      end do
      do k = 1, 6
         y = matmul(transpose(q), matmul(basis(k), q))
         y = y*phi
         l(:, k) = to6(matmul(q, matmul(y, transpose(q))))
      end do
   end function expm_frechet66

   ! symmetric square root of an SPD matrix; ok = .false. if not positive definite
   pure subroutine sqrtm_spd(a, u, ok)
      real(dp), intent(in) :: a(3, 3)
      real(dp), intent(out) :: u(3, 3)
      logical, intent(out) :: ok
      real(dp) :: lam(3), q(3, 3), d(3, 3)
      integer :: i
      call sym_eig3(a, lam, q)
      ok = all(lam > 0.0_dp)
      d = 0.0_dp
      if (ok) then
         do i = 1, 3
            d(i, i) = sqrt(lam(i))
         end do
      end if
      u = matmul(q, matmul(d, transpose(q)))
   end subroutine sqrtm_spd

   !==========================================================================
   ! Elastic law: energy, Se = 2 dW/dCe and C_e = dSe/dCe (Mandel 6x6, if wanted)
   !==========================================================================
   pure subroutine np_elastic(mat, c3, w, se3, c66, wantc, ok)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: c3(3, 3)
      real(dp), intent(out) :: w, se3(3, 3), c66(6, 6)
      logical, intent(in) :: wantc
      logical, intent(out) :: ok
      real(dp) :: d, rj, cinv(3, 3), cinv6(6), phi, trc, cev(6), sev(6), kgv(6, 6), sv(6)
      integer :: i, j
      d = det3(c3)
      ok = d > 0.0_dp
      w = 0.0_dp; se3 = 0.0_dp; c66 = 0.0_dp
      if (.not. ok) return
      rj = sqrt(d)
      cinv = inv3(c3, d)
      cinv6 = to6(cinv)
      trc = c3(1, 1) + c3(2, 2) + c3(3, 3)
      select case (mat%ielast)
      case (NP_ELAST_GUO)
         ! the energy of the PLA model (Section 8.1.2): mu = matrix shear modulus, kmod = matrix bulk modulus
         cev = [c3(1, 1), c3(2, 2), c3(3, 3), c3(1, 2), c3(1, 3), c3(2, 3)]
         ! (npicard_guo_kernel, generated code): CE, SE tensor components {11,22,33,12,13,23}, DSEDCE(i, j) = dSE_i/dCE_j
         ! with CE(4..6) each standing for the symmetric pair of shear components
         call npg_modifiedse(mat%vv, mat%mu, mat%kmod, mat%a0, cev, w, sev, kgv)
         se3 = reshape([sev(1), sev(4), sev(5), sev(4), sev(2), sev(6), sev(5), sev(6), sev(3)], [3, 3])
         if (wantc) then
            sv = [1.0_dp, 1.0_dp, 1.0_dp, SQ2, SQ2, SQ2]       ! tensor -> Mandel: Se_M = s Se, Ce_M = s Ce
            do j = 1, 6
               do i = 1, 6
                  c66(i, j) = sv(i)*kgv(i, j)/sv(j)
               end do
            end do
         end if
      case default   ! NP_ELAST_PAPER
         w = 0.5_dp*mat%mu*(trc - 3.0_dp - log(d)) + 0.5_dp*mat%kmod*(rj - 1.0_dp)**2
         phi = mat%kmod*rj*(rj - 1.0_dp) - mat%mu
         se3 = mat%mu*eye3() + phi*cinv
         if (wantc) then
            c66 = (0.5_dp*mat%kmod*(2.0_dp*rj - 1.0_dp)*rj)*outer6(cinv6, cinv6) - phi*sandwich66(cinv)
         end if
      end select
   end subroutine np_elastic

   ! Ms = (Ce Se + Se Ce)/2 and T = dMs/dCe (if wanted)
   pure subroutine mandel_stress(ce3, se3, c66, m3, t66, wantt)
      real(dp), intent(in) :: ce3(3, 3), se3(3, 3), c66(6, 6)
      real(dp), intent(out) :: m3(3, 3), t66(6, 6)
      logical, intent(in) :: wantt
      m3 = 0.5_dp*(matmul(ce3, se3) + matmul(se3, ce3))
      t66 = 0.0_dp
      if (wantt) t66 = 0.5_dp*prodsym66(se3) + 0.5_dp*matmul(prodsym66(ce3), c66)
   end subroutine mandel_stress

   !==========================================================================
   ! Hardening: Y(k) and Y'(k)
   !==========================================================================
   pure subroutine hardening(mat, k, y, dy)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: k
      real(dp), intent(out) :: y, dy
      real(dp) :: ex, kp
      kp = max(k, 0.0_dp)
      ex = exp(-mat%bvoce*kp**mat%nh)
      y = mat%ysat - (mat%ysat - mat%y0)*ex + mat%hlin*kp
      dy = (mat%ysat - mat%y0)*mat%bvoce*mat%nh*kp**(mat%nh - 1.0_dp)*ex + mat%hlin
   end subroutine hardening

   ! (k - kn) Y(k) = prod by a safeguarded Newton iteration on [kn, kn + prod/Y(kn)]
   pure function kappa_solve(mat, kn, prod) result(k)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: kn, prod
      real(dp) :: k, hlo, hhi, y, dy, res, dres, knew
      integer :: it
      k = kn
      if (.not. (prod > 0.0_dp)) return
      call hardening(mat, kn, y, dy)
      hlo = kn; hhi = kn + prod/y
      k = hhi
      do it = 1, 100
         call hardening(mat, k, y, dy)
         res = (k - kn)*y - prod
         dres = y + (k - kn)*dy
         if (res > 0.0_dp) then
            hhi = k
         else
            hlo = k
         end if
         if (abs(res) <= 1.0e-12_dp*prod) exit
         if (dres > 0.0_dp) then
            knew = k - res/dres
         else
            knew = hlo - 1.0_dp
         end if
         if (.not. (knew > hlo .and. knew < hhi)) knew = 0.5_dp*(hlo + hhi)
         if (abs(knew - k) <= 4.0e-16_dp*abs(knew)) then
            k = knew; exit
         end if
         k = knew
      end do
   end function kappa_solve

   !==========================================================================
   ! Yield function: seq, n = d seq/d Ms (Mandel 6), H = d2 seq/d Ms2 (Mandel 6x6, if wanted)
   !==========================================================================
   pure subroutine np_yield(mat, m3, seq, n6, h66, wanth)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: m3(3, 3)
      real(dp), intent(out) :: seq, n6(6), h66(6, 6)
      logical, intent(in) :: wanth
      real(dp) :: trm, s3(3, 3), s6(6), rj2, rj3, p66(6, 6), g, g2, g3, g22, g23, g33, dj3(6), dg(6), &
                  ss6(6), k66(6, 6), d2j3(6, 6), mt(6), dst(6), nt(6), dnt(6, 6)
      seq = 0.0_dp; n6 = 0.0_dp; h66 = 0.0_dp
      trm = m3(1, 1) + m3(2, 2) + m3(3, 3)
      s3 = m3 - (trm/3.0_dp)*eye3()
      s6 = to6(s3)
      rj2 = 0.5_dp*dot_product(s6, s6)
      p66 = eye6() - outer6(I6, I6)/3.0_dp
      select case (mat%iyield)
      case (NP_YLD_VONMISES)
         seq = sqrt(3.0_dp*rj2)
         if (seq > 0.0_dp) then
            n6 = 1.5_dp*s6/seq
            if (wanth) h66 = (1.5_dp/seq)*p66 - outer6(n6, n6)/seq
         end if
      case (NP_YLD_HOSFORD)
         rj3 = det3(s3)
         g = HOS_A*rj2**4 - HOS_B*rj2*rj3**2
         if (g > 0.0_dp) then
            seq = g**0.125_dp
            ss6 = to6(matmul(s3, s3))
            dj3 = matmul(p66, ss6)                                  ! dJ3/dM = dev(s^2)
            g2 = 4.0_dp*HOS_A*rj2**3 - HOS_B*rj3**2
            g3 = -2.0_dp*HOS_B*rj2*rj3
            dg = g2*s6 + g3*dj3
            n6 = 0.125_dp*g**(-0.875_dp)*dg
            if (wanth) then
               g22 = 12.0_dp*HOS_A*rj2**2
               g23 = -2.0_dp*HOS_B*rj3
               g33 = -2.0_dp*HOS_B*rj2
               k66 = prodsym66(s3)
               d2j3 = matmul(p66, matmul(k66, p66))
               h66 = g2*p66 + g3*d2j3 + g22*outer6(s6, s6) + g23*(outer6(s6, dj3) + outer6(dj3, s6)) &
                     + g33*outer6(dj3, dj3)
               h66 = 0.125_dp*g**(-0.875_dp)*h66 - (7.0_dp/64.0_dp)*g**(-1.875_dp)*outer6(dg, dg)
            end if
         end if
      case (NP_YLD_BARLAT)
         call barlat_yield(mat, to6(m3), seq, n6, h66, wanth)
      case (NP_YLD_GHORBELF)
         ! the criterion of the polycarbonate, notebook yld_ghorbel.nb (generated code below):
         ! (q/r)(1 + (r - 1) h(I1)), h the smoothed step of I1, l = 0.075 st0. The kernel works with tensor
         ! components {11, 22, 33, 12, 13, 23}; n and H go to the Mandel basis. On the hydrostatic axis (q = 0)
         ! seq = 0 and the gradient is left zero, as for von Mises.
         if (rj2 > 0.0_dp) then
            mt = [m3(1, 1), m3(2, 2), m3(3, 3), 0.5_dp*(m3(1, 2) + m3(2, 1)), 0.5_dp*(m3(1, 3) + m3(3, 1)), &
                  0.5_dp*(m3(2, 3) + m3(3, 2))]
            call np_ghorbelf_kernel(mat%rsig, mat%st0, mt, seq, dst, nt, dnt)
            n6(1:3) = nt(1:3); n6(4:6) = SQ2*nt(4:6)
            if (wanth) then
               h66(1:3, 1:3) = dnt(1:3, 1:3)
               h66(4:6, 1:3) = SQ2*dnt(4:6, 1:3)
               h66(1:3, 4:6) = dnt(1:3, 4:6)/SQ2
               h66(4:6, 4:6) = dnt(4:6, 4:6)
            end if
         end if
      end select
   end subroutine np_yield

   !==========================================================================
   ! np_ghorbelf_kernel -- the equivalent stress of the polycarbonate (Section 8.1.1, yield 4):
   !     seq(M) = (q/r) (1 + (r - 1) h(I1)),  q = sqrt(3 J2(M)),  h = (1 + I1/sqrt(I1^2 + l^2))/2,  l = 0.075 st0
   ! (r = rsigma = sigma_c/sigma_t, st0 = sigma_t0). seq -> q for I1 >> l (tension) and seq -> q/r for I1 << -l
   ! (compression); h is a smoothed step of the first invariant, which makes seq C-infinity away from q = 0.
   ! The statements below are the AceGen output of the notebook yld_ghorbel.nb. mandel = {M11, M22, M33, M12, M13,
   ! M23} (tensor components); dseqdmandel = d seq/d mandel (the shear entries are derivatives with respect to the
   ! symmetric pair); norm = the flow direction in tensor components (shear entries dseqdmandel/2); dnormdmandel =
   ! d norm/d mandel. Do not edit by hand: regenerate from the notebook. Changes to the generated file: the work
   ! vector v is local and the routine is PURE with declared intents.
   !==========================================================================
   pure subroutine np_ghorbelf_kernel(rsigma, st0, mandel, seq, dseqdmandel, norm, dnormdmandel)
      double precision, intent(in) :: rsigma, st0, mandel(6)
      double precision, intent(out) :: seq, dseqdmandel(6), norm(6), dnormdmandel(6, 6)
      double precision :: v(181)
  v(167)=2d0*mandel(6)**2
  v(166)=2d0*mandel(5)**2
  v(165)=2d0*mandel(4)**2
  v(164)=0.5625d-2*st0**2
  v(163)=mandel(1)+mandel(2)+mandel(3)
  v(162)=1d0/rsigma
  v(161)=(-1d0)+rsigma
  v(173)=0.5d0*v(161)
  v(97)=v(162)*v(173)
  v(73)=2d0*v(163)
  v(57)=(v(163)*v(163))
  v(56)=v(164)+v(57)
  v(74)=1d0/sqrt(v(56))
  v(171)=v(73)*v(74)
  v(75)=(-0.5d0)*v(171)/v(56)
  v(86)=v(163)*v(75)
  v(87)=v(162)*(0.15d1+v(161)*(0.75d0+0.75d0*v(163)*v(74)))
  v(55)=1d0/sqrt(v(56))
  v(100)=v(55)+v(86)
  v(96)=v(100)*v(97)
  v(49)=(-1d0/3d0)*v(163)
  v(65)=mandel(3)+v(49)
  v(63)=mandel(2)+v(49)
  v(60)=mandel(1)+v(49)
  v(50)=(v(165)+v(166)+v(167)+(v(60)*v(60))+(v(63)*v(63))+(v(65)*v(65)))/2d0
  v(170)=sqrt(3d0*v(50))
  v(76)=1d0/v(170)
  v(169)=0.15d1*v(76)
  v(168)=3d0*v(76)
  v(85)=mandel(6)*v(168)
  v(105)=v(85)*v(96)
  v(83)=mandel(5)*v(168)
  v(104)=v(83)*v(96)
  v(81)=mandel(4)*v(168)
  v(103)=v(81)*v(96)
  v(79)=v(169)*v(65)
  v(78)=v(169)*v(63)
  v(77)=v(169)*v(60)
  v(99)=v(170)*((v(171)*v(57))/v(56)**2+(3d0-v(57)/v(56))*v(75))*v(97)
  v(102)=v(79)*v(96)+v(99)
  v(101)=v(78)*v(96)+v(99)
  v(98)=v(77)*v(96)+v(99)
  v(90)=1d0/v(170)**2
  v(172)=-(v(87)*v(90))
  v(95)=v(172)*v(85)
  v(94)=v(172)*v(83)
  v(93)=v(172)*v(81)
  v(89)=v(169)*(v(74)+v(86))*v(97)
  v(92)=v(172)*v(79)+v(89)
  v(91)=v(172)*v(78)+v(89)
  v(88)=v(172)*v(77)+v(89)
  v(62)=v(76)*v(87)
  v(130)=(2d0/3d0)*v(62)
  v(124)=(-1d0/3d0)*v(62)
  v(176)=v(124)+v(98)
  v(175)=v(102)+v(124)
  v(174)=v(101)+v(124)
  v(115)=2d0*v(62)
  v(58)=v(170)*v(96)
  v(61)=v(58)+v(60)*v(62)
  v(64)=v(58)+v(62)*v(63)
  v(66)=v(58)+v(62)*v(65)
  v(67)=mandel(4)*v(115)
  v(68)=mandel(5)*v(115)
  v(69)=mandel(6)*v(115)
  v(146)=mandel(4)*v(94)
  v(147)=mandel(4)*v(95)
  v(152)=mandel(5)*v(95)
  seq=v(162)*v(170)*(1d0+v(173)*(1d0+v(163)*v(55)))
  dseqdmandel(1)=v(61)
  dseqdmandel(2)=v(64)
  dseqdmandel(3)=v(66)
  dseqdmandel(4)=v(67)
  dseqdmandel(5)=v(68)
  dseqdmandel(6)=v(69)
  norm(1)=v(61)
  norm(2)=v(64)
  norm(3)=v(66)
  norm(4)=v(67)/2d0
  norm(5)=v(68)/2d0
  norm(6)=v(69)/2d0
  dnormdmandel(1,1)=v(130)+v(60)*v(88)+v(98)
  dnormdmandel(1,2)=v(174)+v(60)*v(91)
  dnormdmandel(1,3)=v(175)+v(60)*v(92)
  dnormdmandel(1,4)=v(103)+v(60)*v(93)
  dnormdmandel(1,5)=v(104)+v(60)*v(94)
  dnormdmandel(1,6)=v(105)+v(60)*v(95)
  dnormdmandel(2,1)=v(176)+v(63)*v(88)
  dnormdmandel(2,2)=v(101)+v(130)+v(63)*v(91)
  dnormdmandel(2,3)=v(175)+v(63)*v(92)
  dnormdmandel(2,4)=v(103)+v(63)*v(93)
  dnormdmandel(2,5)=v(104)+v(63)*v(94)
  dnormdmandel(2,6)=v(105)+v(63)*v(95)
  dnormdmandel(3,1)=v(176)+v(65)*v(88)
  dnormdmandel(3,2)=v(174)+v(65)*v(91)
  dnormdmandel(3,3)=v(102)+v(130)+v(65)*v(92)
  dnormdmandel(3,4)=v(103)+v(65)*v(93)
  dnormdmandel(3,5)=v(104)+v(65)*v(94)
  dnormdmandel(3,6)=v(105)+v(65)*v(95)
  dnormdmandel(4,1)=mandel(4)*v(88)
  dnormdmandel(4,2)=mandel(4)*v(91)
  dnormdmandel(4,3)=mandel(4)*v(92)
  dnormdmandel(4,4)=(v(115)+2d0*mandel(4)*v(93))/2d0
  dnormdmandel(4,5)=v(146)
  dnormdmandel(4,6)=v(147)
  dnormdmandel(5,1)=mandel(5)*v(88)
  dnormdmandel(5,2)=mandel(5)*v(91)
  dnormdmandel(5,3)=mandel(5)*v(92)
  dnormdmandel(5,4)=v(146)
  dnormdmandel(5,5)=(v(115)+2d0*mandel(5)*v(94))/2d0
  dnormdmandel(5,6)=v(152)
  dnormdmandel(6,1)=mandel(6)*v(88)
  dnormdmandel(6,2)=mandel(6)*v(91)
  dnormdmandel(6,3)=mandel(6)*v(92)
  dnormdmandel(6,4)=v(147)
  dnormdmandel(6,5)=v(152)
  dnormdmandel(6,6)=(v(115)+2d0*mandel(6)*v(95))/2d0
   end subroutine np_ghorbelf_kernel

   !==========================================================================
   ! Barlat Yld2004-18p for an even exponent a (yield 3). In units of Y0:
   ! x1 = L1 m/Y0, x2 = L2 m/Y0 (L = C P, Mandel), X1, X2 their tensors, t1(k) = tr X1^k,
   ! phi = sum_k c_k t1(k) t2(a-k), c_k = binom(a,k)(-1)^k, seq = Y0 (phi/4)^(1/a).
   ! dt(k)/dx = k to6(X^(k-1)); d to6(X^j)/dx = the map dX -> sum_{m<j} X^m dX X^(j-1-m).
   !==========================================================================
   pure subroutine barlat_yield(mat, m6, seq, n6, h66, wanth)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: m6(6)
      real(dp), intent(out) :: seq, n6(6), h66(6, 6)
      logical, intent(in) :: wanth
      integer :: a, k, j
      real(dp) :: y, x1(6), x2(6), p1(3, 3, 0:NP_BARLAT_MAXA), p2(3, 3, 0:NP_BARLAT_MAXA), t1(0:NP_BARLAT_MAXA), &
                  t2(0:NP_BARLAT_MAXA), c(0:NP_BARLAT_MAXA), phi, g1(6), g2(6), gphi(6), w1(0:NP_BARLAT_MAXA), &
                  w2(0:NP_BARLAT_MAXA), h11(6, 6), h22(6, 6), h12(6, 6), hphi(6, 6), fac
      seq = 0.0_dp; n6 = 0.0_dp; h66 = 0.0_dp
      a = mat%abar; y = mat%y0
      x1 = matmul(mat%bl1, m6)/y; x2 = matmul(mat%bl2, m6)/y
      p1(:, :, 0) = eye3(); p2(:, :, 0) = eye3()
      p1(:, :, 1) = from6(x1); p2(:, :, 1) = from6(x2)
      do k = 2, a
         p1(:, :, k) = matmul(p1(:, :, k - 1), p1(:, :, 1))
         p2(:, :, k) = matmul(p2(:, :, k - 1), p2(:, :, 1))
      end do
      do k = 0, a
         t1(k) = p1(1, 1, k) + p1(2, 2, k) + p1(3, 3, k)
         t2(k) = p2(1, 1, k) + p2(2, 2, k) + p2(3, 3, k)
      end do
      c(0) = 1.0_dp
      do k = 1, a
         c(k) = -c(k - 1)*real(a - k + 1, dp)/real(k, dp)
      end do
      phi = 0.0_dp
      do k = 0, a
         phi = phi + c(k)*t1(k)*t2(a - k)
      end do
      if (.not. (phi > 0.0_dp)) return
      ! gradients with respect to x1, x2: w1(j) multiplies to6(X1^j), w2(j) multiplies to6(X2^j)
      w1 = 0.0_dp; w2 = 0.0_dp
      do k = 1, a
         w1(k - 1) = w1(k - 1) + c(k)*real(k, dp)*t2(a - k)
      end do
      do k = 0, a - 1
         w2(a - k - 1) = w2(a - k - 1) + c(k)*real(a - k, dp)*t1(k)
      end do
      g1 = 0.0_dp; g2 = 0.0_dp
      do j = 0, a - 1
         g1 = g1 + w1(j)*to6(p1(:, :, j))
         g2 = g2 + w2(j)*to6(p2(:, :, j))
      end do
      gphi = matmul(transpose(mat%bl1), g1) + matmul(transpose(mat%bl2), g2)   ! d phi / d(m/Y0)
      seq = y*(0.25_dp*phi)**(1.0_dp/real(a, dp))
      fac = seq/(real(a, dp)*phi)                                             ! d seq/d phi, in stress units
      n6 = (fac/y)*gphi
      if (.not. wanth) return
      h11 = 0.0_dp; h22 = 0.0_dp; h12 = 0.0_dp
      do j = 1, a - 1
         h11 = h11 + w1(j)*powsand66(p1, j)
         h22 = h22 + w2(j)*powsand66(p2, j)
      end do
      do k = 1, a - 1
         h12 = h12 + c(k)*real(k, dp)*real(a - k, dp)*outer6(to6(p1(:, :, k - 1)), to6(p2(:, :, a - k - 1)))
      end do
      hphi = matmul(transpose(mat%bl1), matmul(h11, mat%bl1)) + matmul(transpose(mat%bl2), matmul(h22, mat%bl2)) &
             + matmul(transpose(mat%bl1), matmul(h12, mat%bl2)) + matmul(transpose(mat%bl2), matmul(transpose(h12), mat%bl1))
      h66 = (fac/(y*y))*(hphi + (1.0_dp/real(a, dp) - 1.0_dp)*outer6(gphi, gphi)/phi)
   end subroutine barlat_yield

   ! d to6(X^j)/dx: the map dX -> sum_{m=0}^{j-1} X^m dX X^(j-1-m), powers of X in p(:,:,0:j)
   pure function powsand66(p, j) result(s)
      real(dp), intent(in) :: p(3, 3, 0:NP_BARLAT_MAXA)
      integer, intent(in) :: j
      real(dp) :: s(6, 6), y(3, 3), b(3, 3)
      integer :: k, m
      do k = 1, 6
         b = basis(k); y = 0.0_dp
         do m = 0, j - 1
            y = y + matmul(p(:, :, m), matmul(b, p(:, :, j - 1 - m)))
         end do
         s(:, k) = to6(y)
      end do
   end function powsand66

   !==========================================================================
   ! Barlat Yld2004-18p coefficient sets. Barlat et al. (2005), Int J Plast 21:1009, eq. (7):
   ! C' = [[0,-c12,-c13],[-c21,0,-c23],[-c31,-c32,0]] (+) diag(c44, c55, c66) on the deviator
   ! s = (xx, yy, zz, yz, zx, xy); here mapped to the Mandel order (11, 22, 33, 12, 13, 23):
   ! c66 (xy) -> 4, c55 (zx) -> 5, c44 (yz) -> 6 (a diagonal shear factor is the same on tensor
   ! and Mandel components). bl = C P with P the deviatoric projector, so bl acts on Ms itself.
   !   iset 0: isotropic (every c = 1): seq = Hosford with exponent a
   !   iset 1: AA2090-T3, Barlat et al. (2005), Table 2, a = 8
   !==========================================================================
   subroutine np_barlat_set(mat, iset, info)
      type(np_material), intent(inout) :: mat
      integer, intent(in) :: iset
      integer, intent(out) :: info
      real(dp) :: cp(9), cpp(9)
      ! order: c12 c13 c21 c23 c31 c32 c44 c55 c66
      info = NP_OK
      select case (iset)
      case (0)
         cp = 1.0_dp; cpp = 1.0_dp
      case (1)
         cp = [-0.069888_dp, 0.936408_dp, 0.079143_dp, 1.003060_dp, 0.524741_dp, 1.363180_dp, &
               1.023770_dp, 1.069060_dp, 0.954322_dp]
         cpp = [0.981171_dp, 0.476741_dp, 0.575316_dp, 0.866827_dp, 1.145010_dp, -0.079294_dp, &
                1.051660_dp, 1.147100_dp, 1.404620_dp]
      case default
         info = NP_ERR_INPUT; return
      end select
      mat%bl1 = barlat_mandel(cp); mat%bl2 = barlat_mandel(cpp)
   end subroutine np_barlat_set

   pure function barlat_mandel(c) result(b)
      real(dp), intent(in) :: c(9)
      real(dp) :: b(6, 6), cm(6, 6), p66(6, 6)
      cm = 0.0_dp
      cm(1, 2) = -c(1); cm(1, 3) = -c(2); cm(2, 1) = -c(3); cm(2, 3) = -c(4); cm(3, 1) = -c(5); cm(3, 2) = -c(6)
      cm(6, 6) = c(7); cm(5, 5) = c(8); cm(4, 4) = c(9)
      p66 = eye6() - outer6(I6, I6)/3.0_dp
      b = matmul(cm, p66)
   end function barlat_mandel

   !==========================================================================
   ! Material validation
   !==========================================================================
   subroutine np_check_material(mat, info, msg)
      type(np_material), intent(in) :: mat
      integer, intent(out) :: info
      character(len=*), intent(out) :: msg
      info = NP_ERR_INPUT
      msg = ''
      if (.not. (mat%mu > 0.0_dp .and. mat%kmod > 0.0_dp)) then
         msg = 'mu > 0 and K > 0 are required'; return
      end if
      if (mat%ielast /= NP_ELAST_PAPER .and. mat%ielast /= NP_ELAST_GUO) then
         msg = 'elastic law choice must be 1 or 2'; return
      end if
      if (mat%ielast == NP_ELAST_GUO .and. .not. (mat%vv >= 0.0_dp .and. mat%vv < 1.0_dp .and. &
          abs(sqrt(sum(mat%a0**2)) - 1.0_dp) < 1.0e-10_dp)) then
         msg = 'PLA energy: 0 <= vv < 1 and a unit a0 are required'; return
      end if
      if (mat%iyield < 1 .or. mat%iyield > 4) then
         msg = 'yield choice must be 1..4'; return
      end if
      if (mat%iyield == NP_YLD_BARLAT) then
         if (mat%abar < 2 .or. mat%abar > NP_BARLAT_MAXA .or. mod(mat%abar, 2) /= 0) then
            msg = 'Barlat Yld2004-18p: the exponent a must be even, 2..16 (trace-polynomial form)'; return
         end if
         if (.not. (maxval(abs(mat%bl1)) > 0.0_dp .and. maxval(abs(mat%bl2)) > 0.0_dp)) then
            msg = 'Barlat Yld2004-18p: coefficients not set (np_barlat_set)'; return
         end if
      end if
      if (mat%iyield == NP_YLD_HOSFORD .and. abs(mat%nhos - 8.0_dp) > 0.0_dp) then
         msg = 'Hosford: only the exponent n = 8 (polynomial form) is implemented'; return
      end if
      if (mat%iyield == NP_YLD_GHORBELF .and. .not. (mat%st0 > 0.0_dp .and. mat%rsig > 0.0_dp)) then
         msg = 'polycarbonate criterion: r > 0 and st0 > 0 are required'; return
      end if
      if (.not. (mat%y0 > 0.0_dp)) then
         msg = 'Y0 > 0 is required'; return
      end if
      if (.not. (mat%ysat >= mat%y0 .and. mat%bvoce >= 0.0_dp .and. mat%nh >= 1.0_dp .and. mat%hlin >= 0.0_dp)) then
         msg = 'Hardening: Ysat >= Y0, b >= 0, nh >= 1, H >= 0 are required'; return
      end if
      if (mat%nn < 1 .or. mat%nn > NP_MAXNN) then
         msg = 'nn must be in 1..NP_MAXNN'; return
      end if
      if (.not. (mat%eps_theta >= 0.0_dp .and. mat%eps_f > 0.0_dp)) then
         msg = 'eps_theta >= 0 and eps_f > 0 are required'; return
      end if
      info = NP_OK
   end subroutine np_check_material

   !==========================================================================
   ! The state along the ray gamma -> exp(-gamma N) Ce^tr exp(-gamma N) at frozen N.
   ! level 0: f only; 1: + the exact frozen-path slope df/dgamma; 2: + the Hessian H
   !==========================================================================
   subroutine ray_eval(mat, ray, cetr3, kn, gam, level, ev, ok)
      type(np_material), intent(in) :: mat
      type(ray_t), intent(in) :: ray
      real(dp), intent(in) :: cetr3(3, 3), kn, gam
      integer, intent(in) :: level
      type(eval_t), intent(out) :: ev
      logical, intent(out) :: ok
      real(dp) :: dm(6)
      ev%gam = gam
      ! an exponent beyond +-EXPARG is not evaluated (the state is refused; the negated test also refuses a NaN)
      if (.not. all(abs(gam*ray%lam) <= EXPARG)) then
         ok = .false.; ev%f = huge(1.0_dp); ev%df = 0.0_dp
         return
      end if
      ev%e3 = expm_eig(ray%lam, ray%q, -gam)
      ev%ce3 = matmul(ev%e3, matmul(cetr3, ev%e3))
      call np_elastic(mat, ev%ce3, ev%w, ev%se3, ev%c66, level >= 1, ok)
      if (.not. ok) return
      call mandel_stress(ev%ce3, ev%se3, ev%c66, ev%m3, ev%t66, level >= 1)
      call np_yield(mat, ev%m3, ev%seq, ev%n6, ev%h66, level >= 2)
      ev%m6 = to6(ev%m3)
      ev%eta = dot_product(ev%m6, ray%n6)
      ev%prod = gam*ev%eta
      ev%kappa = kappa_solve(mat, kn, ev%prod)
      call hardening(mat, ev%kappa, ev%y, ev%dy)
      ev%f = ev%seq - ev%y
      ev%dk = ev%y + ev%dy*(ev%kappa - kn)
      ev%alpha = ev%dy/ev%dk
      ev%df = 0.0_dp; ev%phat6 = 0.0_dp
      if (level >= 1) then
         ev%phat6 = to6(matmul(ray%n3, ev%ce3) + matmul(ev%ce3, ray%n3))   ! dCe/dgamma = -phat
         dm = -matmul(ev%t66, ev%phat6)
         ev%df = dot_product(ev%n6, dm) - ev%alpha*(ev%eta + gam*dot_product(ray%n6, dm))
      end if
      ok = ieee_is_finite(ev%f) .and. ieee_is_finite(ev%df)
   end subroutine ray_eval

   pure function c_from_e(eps) result(c3)
      real(dp), intent(in) :: eps(6)
      real(dp) :: c3(3, 3)
      c3(1, 1) = 1.0_dp + 2.0_dp*eps(1)
      c3(2, 2) = 1.0_dp + 2.0_dp*eps(2)
      c3(3, 3) = 1.0_dp + 2.0_dp*eps(3)
      c3(1, 2) = eps(4); c3(2, 1) = eps(4)
      c3(1, 3) = eps(5); c3(3, 1) = eps(5)
      c3(2, 3) = eps(6); c3(3, 2) = eps(6)
   end function c_from_e

   pure function voigt_to_3(v) result(a)
      real(dp), intent(in) :: v(6)
      real(dp) :: a(3, 3)
      a(1, 1) = v(1); a(2, 2) = v(2); a(3, 3) = v(3)
      a(1, 2) = v(4); a(2, 1) = v(4)
      a(1, 3) = v(5); a(3, 1) = v(5)
      a(2, 3) = v(6); a(3, 2) = v(6)
   end function voigt_to_3

   pure function voigt_from_3(a) result(v)
      real(dp), intent(in) :: a(3, 3)
      real(dp) :: v(6)
      v(1) = a(1, 1); v(2) = a(2, 2); v(3) = a(3, 3)
      v(4) = 0.5_dp*(a(1, 2) + a(2, 1)); v(5) = 0.5_dp*(a(1, 3) + a(3, 1)); v(6) = 0.5_dp*(a(2, 3) + a(3, 2))
   end function voigt_from_3
   !==========================================================================
   ! The baseline of Section 8.2 (Newton iteration with staggered hardening): Newton iteration on x = (Ce [Mandel],
   ! gamma) for
   !    R1 = Ce - exp(-gamma N(Ce)) Ce^tr exp(-gamma N(Ce)) = 0,   R2 = (seq(Ms(Ce)) - Y_n)/Y_n = 0,
   ! N = d seq/d Ms re-evaluated at every iterate, the yield stress frozen at Y_n = Y(kappa_n); after convergence
   ! kappa = kappa_n + gamma Ms:N / Y_n (explicit, power balance at Y_n). Analytic Jacobian (derivative of the
   ! exponential, the yield Hessian, dMs/dCe); full Newton steps from the trial state (gamma = 0), at most
   ! 50 iterations to max(|R1|, |R2|) < 1e-12; no convergence or a singular Jacobian -> NP_ERR_DIVERGED.
   ! Tangent: the exact derivative of this map (Y_n does not depend on the strain).
   ! Arguments and conventions as np_update_alg1.
   !==========================================================================
   subroutine np_update_base(mat, eps, cpinv_n, kappa_n, stress, dsde, cpinv, kappa, psi, info, st)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: eps(6), cpinv_n(6), kappa_n
      real(dp), intent(out) :: stress(6), dsde(6, 6), cpinv(6), kappa, psi
      integer, intent(out) :: info
      type(np_state), intent(out), optional :: st
      type(np_state) :: s
      real(dp) :: c3(3, 3), u(3, 3), cetr3(3, 3), uhat(6, 6), w, se3(3, 3), c66(6, 6), m3(3, 3), t66(6, 6), &
                  ntr6(6), htr66(6, 6), ytr, dytr, ds66(6, 6), s3(3, 3)
      logical :: ok
      info = NP_OK
      stress = 0.0_dp; dsde = 0.0_dp; cpinv = cpinv_n; kappa = kappa_n; psi = 0.0_dp
      !--- kinematics: Ce^tr = U C U, U = sqrt(Cp_n^-1)
      c3 = c_from_e(eps)
      call sqrtm_spd(voigt_to_3(cpinv_n), u, ok)
      if (.not. ok) then
         info = NP_ERR_INPUT; call finish(); return
      end if
      cetr3 = matmul(u, matmul(c3, u))
      uhat = sandwich66(u)
      s%cetr6 = to6(cetr3)
      !--- Ce^tr must be positive definite (leading principal minors > 0); the negated form also rejects a NaN
      if (.not. (cetr3(1, 1) > 0.0_dp .and. cetr3(1, 1)*cetr3(2, 2) - cetr3(1, 2)*cetr3(2, 1) > 0.0_dp &
                 .and. det3(cetr3) > 0.0_dp)) then
         info = NP_ERR_INPUT; call finish(); return
      end if
      !--- trial state
      call np_elastic(mat, cetr3, w, se3, c66, .true., ok)
      if (.not. ok) then
         info = NP_ERR_INPUT; call finish(); return
      end if
      call mandel_stress(cetr3, se3, c66, m3, t66, .true.)
      call np_yield(mat, m3, s%seqtr, ntr6, htr66, .false.)
      call hardening(mat, kappa_n, ytr, dytr)
      s%ytr = ytr; s%ftr = s%seqtr - ytr; s%ratio = s%seqtr/ytr
      if (.not. (ieee_is_finite(s%ftr))) then
         info = NP_ERR_NONFINITE; call finish(); return
      end if
      if (s%ftr <= 0.0_dp) then
         !--- elastic: S = U Se U, dS/dC = Uhat C_e Uhat, history unchanged
         s%plastic = .false.
         s%n6 = ntr6; s%m6 = to6(m3); s%ce6 = s%cetr6; s%se6 = to6(se3); s%seq = s%seqtr; s%kappa = kappa_n
         s3 = matmul(u, matmul(se3, u))
         ds66 = matmul(uhat, matmul(c66, uhat))
         psi = w
         call pack_out(s3, ds66)
         call finish(); return
      end if
      !--- plastic
      s%plastic = .true.
      call classical_return()
   contains
      ! the Newton iteration, the state and the tangent
      subroutine classical_return()
         real(dp) :: ce6(6), gam, yn, ce3c(3, 3), wc, sec(3, 3), c66c(6, 6), mc(3, 3), t66c(6, 6), seqc, nc6(6), &
                     hc66(6, 6), lamc(3), qc(3, 3), ea(3, 3), r(7), jac(7, 7), dx(7), lc66(6, 6), dedA(6, 6), &
                     rhs(7, 6), dxdc(7, 6), dadc(6, 6), dpdc(6, 6), p3(3, 3), e2(3, 3), rn
         integer :: it, k
         logical :: okc
         yn = ytr
         ce6 = s%cetr6; gam = 0.0_dp
         do it = 1, 51
            ce3c = from6(ce6)
            call np_elastic(mat, ce3c, wc, sec, c66c, .true., okc)
            if (.not. okc .or. .not. all(ieee_is_finite(ce6))) exit
            call mandel_stress(ce3c, sec, c66c, mc, t66c, .true.)
            call np_yield(mat, mc, seqc, nc6, hc66, .true.)
            call sym_eig3(from6(nc6), lamc, qc)
            ea = expm_eig(lamc, qc, -gam)
            r(1:6) = ce6 - to6(matmul(ea, matmul(cetr3, ea)))
            r(7) = (seqc - yn)/yn
            rn = maxval(abs(r))
            if (rn < 1.0e-12_dp) then
               okc = .true.; exit
            end if
            okc = .false.
            if (it == 51) exit
            lc66 = expm_frechet66(-gam*lamc, qc)
            dedA = matmul(sandwich2_66(cetr3, ea), lc66)            ! dE/dA, E = e^A Ce^tr e^A, A = -gamma N
            jac = 0.0_dp
            jac(1:6, 1:6) = eye6() + gam*matmul(dedA, matmul(hc66, t66c))
            jac(1:6, 7) = matmul(dedA, nc6)
            jac(7, 1:6) = matmul(nc6, t66c)/yn
            call solve_dense(7, jac, -r, dx, okc)
            if (.not. okc) exit
            ce6 = ce6 + dx(1:6); gam = gam + dx(7)
            s%its = it
         end do
         if (.not. okc) then
            info = NP_ERR_DIVERGED; call finish(); return
         end if
         s%its = it - 1
         ! state, history, stress
         s%gamma = gam; s%seq = seqc; s%f = seqc - yn
         s%kappa = kappa_n + gam*dot_product(to6(mc), nc6)/yn
         s%n6 = nc6; s%m6 = to6(mc); s%ce6 = ce6; s%se6 = to6(sec)
         kappa = s%kappa
         e2 = expm_eig(lamc, qc, -2.0_dp*gam)
         cpinv = voigt_from_3(matmul(u, matmul(e2, u)))
         p3 = matmul(ea, matmul(sec, ea))
         ! tangent: dx/dCe^tr = J^-1 [sandwich66(e^A); 0] at the converged state
         lc66 = expm_frechet66(-gam*lamc, qc)
         dedA = matmul(sandwich2_66(cetr3, ea), lc66)
         jac = 0.0_dp
         jac(1:6, 1:6) = eye6() + gam*matmul(dedA, matmul(hc66, t66c))
         jac(1:6, 7) = matmul(dedA, nc6)
         jac(7, 1:6) = matmul(nc6, t66c)/yn
         rhs = 0.0_dp; rhs(1:6, :) = sandwich66(ea)
         do k = 1, 6
            call solve_dense(7, jac, rhs(:, k), dxdc(:, k), okc)
            if (.not. okc) then
               info = NP_ERR_DIVERGED; call finish(); return
            end if
         end do
         ! P = e^A Se e^A: dP = (dE^A Se e^A + e^A Se dE^A) + e^A dSe e^A, dA = -dgamma N - gamma dN, dN = H T dCe
         dadc = -outer6(nc6, dxdc(7, :)) - gam*matmul(hc66, matmul(t66c, dxdc(1:6, :)))
         dpdc = matmul(sandwich2_66(sec, ea), matmul(lc66, dadc)) + matmul(sandwich66(ea), matmul(c66c, dxdc(1:6, :)))
         psi = wc
         call pack_out(matmul(u, matmul(p3, u)), matmul(uhat, matmul(dpdc, uhat)))
         call finish()
      end subroutine classical_return
      subroutine pack_out(s3, ds66)
         real(dp), intent(in) :: s3(3, 3), ds66(6, 6)
         integer :: i, j
         stress = voigt_from_3(s3)
         do j = 1, 6
            do i = 1, 6
               dsde(i, j) = SFAC(i)*CFAC(j)*ds66(i, j)
            end do
         end do
      end subroutine pack_out
      subroutine finish()
         if (present(st)) st = s
      end subroutine finish
   end subroutine np_update_base

   !==========================================================================
   ! Algorithm 1 of the paper: the return mapping with the inner scalar iteration on gamma at a frozen direction N
   ! and the outer Picard update N <- d seq/d Ms.
   !   elastic predictor; N <- d seq/d Ms at the trial state;
   !   outer loop (Picard), at each frozen N:
   !     gamma_L = 0; f(0) <= eps_f Y: gamma = 0 (small correction);
   !     gamma_R by doubling, from the Newton estimate f(0)/s(0), until f(gamma_R) <= 0; a ray on which no
   !     sign change is found up to GAMMA_CAP, or whose states the elastic law cannot evaluate, fails the update
   !     ("bracketing fails"): status NP_ERR_NOROOT, the caller substeps;
   !     safeguarded Newton from gamma_L with the slope s = H* + dY/dgamma built from the frozen direction,
   !     H* = N : T : (N Ce + Ce N), dY/dgamma = alpha (Ms:N - gamma H*) (Sections 3 and 4), bisection when the
   !     step leaves the bracket, at most NMAX_SCALAR iterations (then NP_ERR_DIVERGED);
   !     at every evaluation the hardening equation is solved by Newton from kappa_n with the safeguard
   !     kappa <- max(kappa, kappa_n), to 1e-12 gamma Ms:N (Section 6.4);
   !     state at the root, N_new = d seq/d Ms; stop when ||N_new - N|| < eps_theta or after nn updates.
   !   With eps_theta = 0 exactly nn updates are made. With eps_theta > 0 a cap reached returns the last update
   !   with status NP_UNCONVERGED.
   !   Early exit (eps_theta > 0 only): the contraction factor is proportional to the increment (Section 7), so an
   !   iteration that will not reach eps_theta within the cap is recognized from the measured two-update ratio
   !   q2 = ||r_k||/||r_k-2|| (k >= 4): it is abandoned (NP_UNCONVERGED) when q2 >= 1 at two consecutive updates,
   !   or when k + 2 ln(eps_theta/||r_k||)/ln(q2) exceeds 1.5 nn; noexit = .true. disables it.
   !   Tangent: the closed form of Section 5 at the returned state.
   ! In: mat; eps engineering Green-Lagrange strain; cpinv_n, kappa_n the history. Out: stress (second
   ! Piola-Kirchhoff), dsde = dS/dE, cpinv, kappa, psi (elastic energy), info (status, header), st (optional: the
   ! state and the counts), nkap, nksolve (optional: hardening Newton iterations and hardening solutions of the
   ! call, the counts of Section 6.4).
   !==========================================================================
   subroutine np_update_alg1(mat, eps, cpinv_n, kappa_n, stress, dsde, cpinv, kappa, psi, info, st, nkap, nksolve, noexit)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: eps(6), cpinv_n(6), kappa_n
      real(dp), intent(out) :: stress(6), dsde(6, 6), cpinv(6), kappa, psi
      integer, intent(out) :: info
      type(np_state), intent(out), optional :: st
      integer, intent(out), optional :: nkap, nksolve
      logical, intent(in), optional :: noexit
      real(dp), parameter :: GAMMA_CAP = 8.0_dp
      integer, parameter :: NMAX_SCALAR = 100
      type(np_state) :: s
      type(ray_t) :: ray
      type(eval_t) :: ev
      real(dp) :: c3(3, 3), u(3, 3), cetr3(3, 3), uhat(6, 6), w, se3(3, 3), c66(6, 6), m3(3, 3), t66(6, 6), ntr6(6), &
                  htr66(6, 6), ytr, dytr, r6(6), dn, eps_theta, gam, d66(6, 6), rl66(6, 6), rn66(6, 6), s66(6, 6), &
                  tnt(6), pbar(6), jbar(6), v(6), denom, gvec(6), sp(6), cfull(6, 6), nmat(6, 6), da66(6, 6), &
                  de66(6, 6), dz66(6, 6), ds66(6, 6), s3(3, 3), e2(3, 3), q2, q2prev
      logical :: ok, converged, early
      integer :: nn, k, j, istat, nit, nfe, kapit, kapn
      info = NP_OK
      stress = 0.0_dp; dsde = 0.0_dp; cpinv = cpinv_n; kappa = kappa_n; psi = 0.0_dp
      kapit = 0; kapn = 0
      nn = max(1, min(NP_MAXNN, mat%nn))
      eps_theta = mat%eps_theta
      early = eps_theta > 0.0_dp
      if (present(noexit)) early = early .and. .not. noexit
      !--- kinematics: Ce^tr = U C U, U = sqrt(Cp_n^-1)
      c3 = c_from_e(eps)
      call sqrtm_spd(voigt_to_3(cpinv_n), u, ok)
      if (.not. ok) then
         info = NP_ERR_INPUT; call finish(); return
      end if
      cetr3 = matmul(u, matmul(c3, u))
      uhat = sandwich66(u)
      s%cetr6 = to6(cetr3)
      if (.not. (cetr3(1, 1) > 0.0_dp .and. cetr3(1, 1)*cetr3(2, 2) - cetr3(1, 2)*cetr3(2, 1) > 0.0_dp &
                 .and. det3(cetr3) > 0.0_dp)) then
         info = NP_ERR_INPUT; call finish(); return
      end if
      !--- elastic predictor
      call np_elastic(mat, cetr3, w, se3, c66, .true., ok)
      if (.not. ok) then
         info = NP_ERR_INPUT; call finish(); return
      end if
      call mandel_stress(cetr3, se3, c66, m3, t66, .true.)
      call np_yield(mat, m3, s%seqtr, ntr6, htr66, .false.)
      call hardening(mat, kappa_n, ytr, dytr)
      s%ytr = ytr; s%ftr = s%seqtr - ytr; s%ratio = s%seqtr/ytr
      if (.not. (ieee_is_finite(s%ftr))) then
         info = NP_ERR_NONFINITE; call finish(); return
      end if
      if (s%ftr <= 0.0_dp) then
         s%plastic = .false.
         s%n6 = ntr6; s%m6 = to6(m3); s%ce6 = s%cetr6; s%se6 = to6(se3); s%seq = s%seqtr; s%kappa = kappa_n
         s3 = matmul(u, matmul(se3, u))
         ds66 = matmul(uhat, matmul(c66, uhat))
         psi = w
         call pack_out(s3, ds66)
         call finish(); return
      end if
      !--- N <- d seq/d Ms at the trial state
      s%plastic = .true.
      ray%n6 = ntr6
      converged = .false.; q2prev = 0.0_dp
      !--- outer Picard update of N
      do k = 1, nn
         ray%n3 = from6(ray%n6)
         call sym_eig3(ray%n3, ray%lam, ray%q)
         call scalar_alg1(ev, istat, nit, nfe)
         s%newton = s%newton + nit; s%nfeval = s%nfeval + nfe
         if (istat /= 0) then
            if (istat == 1) then
               info = NP_ERR_NOROOT; s%noroot = s%noroot + 1; s%noroot_last = 1
            else if (istat == 2) then
               info = NP_ERR_DIVERGED
            else
               info = NP_ERR_NONFINITE
            end if
            s%its = k
            call finish(); return
         end if
         r6 = ev%n6 - ray%n6
         dn = sqrt(dot_product(r6, r6))
         s%dn(k) = dn; s%its = k
         s%gamma = ev%gam; s%kappa = ev%kappa; s%seq = ev%seq; s%f = ev%f
         s%n6 = ray%n6; s%m6 = ev%m6; s%ce6 = to6(ev%ce3); s%se6 = to6(ev%se3)
         if (.not. ieee_is_finite(dn)) then
            info = NP_ERR_NONFINITE; call finish(); return
         end if
         converged = dn < eps_theta
         if (converged .or. k == nn) exit
         if (early .and. k >= 4) then
            ! the measured contraction: abandon an iteration that will not reach eps_theta within the cap
            q2 = dn/max(s%dn(k - 2), tiny(1.0_dp))
            if (q2 >= 1.0_dp) then
               if (q2prev >= 1.0_dp) exit
            else if (q2 > 0.0_dp) then
               if (real(k, dp) + 2.0_dp*log(eps_theta/dn)/log(q2) > 1.5_dp*real(nn, dp)) exit
            end if
            q2prev = q2
         else if (early .and. k == 3) then
            q2prev = dn/max(s%dn(1), tiny(1.0_dp))
         end if
         ray%n6 = ev%n6                                  ! refresh direction: N <- d seq/d Ms at Ms(gamma)
      end do
      if (.not. converged .and. eps_theta > 0.0_dp) then
         s%unconverged = .true.; info = NP_UNCONVERGED
      end if
      gam = ev%gam
      !--- the consistent tangent (Section 5) at the returned state
      d66 = sandwich66(ev%e3)
      rl66 = expm_frechet66(-gam*ray%lam, ray%q)
      rn66 = -gam*matmul(sandwich2_66(cetr3, ev%e3), rl66)
      call inv6(eye6() - matmul(rn66, matmul(ev%h66, ev%t66)), s66, ok)
      if (.not. ok) then
         info = NP_ERR_NONFINITE; call finish(); return
      end if
      tnt = matmul(transpose(ev%t66), ev%n6)
      pbar = matmul(transpose(s66), tnt)
      jbar = matmul(transpose(s66), matmul(transpose(ev%t66), ev%n6 + matmul(transpose(ev%h66), ev%m6)))
      v = pbar - ev%alpha*gam*jbar
      denom = dot_product(pbar, ev%phat6) + ev%alpha*ev%eta - ev%alpha*gam*dot_product(jbar, ev%phat6)
      s%dg_den = denom; s%dg_h0 = dot_product(tnt, ev%phat6) + ev%alpha*ev%eta; s%dg_det = det6x(eye6() - matmul(rn66, matmul(ev%h66, ev%t66)))
      gvec = matmul(transpose(d66), v)/denom
      sp = matmul(s66, ev%phat6)
      cfull = matmul(s66, d66) - outer6(sp, gvec)
      nmat = matmul(ev%h66, matmul(ev%t66, cfull))
      !--- S = U E Se E U and dS/dC through E = exp(-gamma N), Se(Ce)
      do j = 1, 6
         da66(:, j) = -(ray%n6*gvec(j) + gam*nmat(:, j))
      end do
      de66 = matmul(rl66, da66)
      dz66 = matmul(sandwich2_66(ev%se3, ev%e3), de66) + matmul(sandwich66(ev%e3), matmul(ev%c66, cfull))
      ds66 = matmul(uhat, matmul(dz66, uhat))
      s3 = matmul(u, matmul(ev%e3, matmul(ev%se3, matmul(ev%e3, u))))
      e2 = expm_eig(ray%lam, ray%q, -2.0_dp*gam)
      cpinv = voigt_from_3(matmul(u, matmul(e2, u)))
      kappa = ev%kappa
      psi = ev%w
      call pack_out(s3, ds66)
      if (.not. (all(ieee_is_finite(stress)) .and. all(ieee_is_finite(dsde)))) then
         info = NP_ERR_NONFINITE
         stress = 0.0_dp; dsde = 0.0_dp; cpinv = cpinv_n; kappa = kappa_n; psi = 0.0_dp
      end if
      call finish()
   contains
      ! the state on the frozen ray at gamma with the hardening Newton iteration of Algorithm 1; level as ray_eval
      subroutine eval_alg1(gam_in, level, e, oke)
         real(dp), intent(in) :: gam_in
         integer, intent(in) :: level
         type(eval_t), intent(out) :: e
         logical, intent(out) :: oke
         real(dp) :: dm(6), kk, yy, dyy, res, dres
         integer :: it
         e%gam = gam_in
         if (.not. all(abs(gam_in*ray%lam) <= EXPARG)) then
            oke = .false.; e%f = huge(1.0_dp); e%df = 0.0_dp
            return
         end if
         e%e3 = expm_eig(ray%lam, ray%q, -gam_in)
         e%ce3 = matmul(e%e3, matmul(cetr3, e%e3))
         call np_elastic(mat, e%ce3, e%w, e%se3, e%c66, level >= 1, oke)
         if (.not. oke) return
         call mandel_stress(e%ce3, e%se3, e%c66, e%m3, e%t66, level >= 1)
         call np_yield(mat, e%m3, e%seq, e%n6, e%h66, level >= 2)
         e%m6 = to6(e%m3)
         e%eta = dot_product(e%m6, ray%n6)
         e%prod = gam_in*e%eta
         ! hardening subproblem: Newton from kappa_n, clamped at kappa_n
         kk = kappa_n
         kapn = kapn + 1
         if (e%prod > 0.0_dp) then
            do it = 1, 100
               call hardening(mat, kk, yy, dyy)
               res = yy*(kk - kappa_n) - e%prod
               if (abs(res) <= 1.0e-12_dp*e%prod) exit
               dres = dyy*(kk - kappa_n) + yy
               kk = max(kk - res/dres, kappa_n)
               kapit = kapit + 1
            end do
         end if
         e%kappa = kk
         call hardening(mat, e%kappa, e%y, e%dy)
         e%f = e%seq - e%y
         e%dk = e%y + e%dy*(e%kappa - kappa_n)
         e%alpha = e%dy/e%dk
         e%df = 0.0_dp; e%phat6 = 0.0_dp
         if (level >= 1) then
            e%phat6 = to6(matmul(ray%n3, e%ce3) + matmul(e%ce3, ray%n3))
            dm = -matmul(e%t66, e%phat6)
            ! the slope of the scalar iteration, with the frozen N in both slots: -(H* + dY/dgamma)
            e%df = dot_product(ray%n6, dm) - e%alpha*(e%eta + gam_in*dot_product(ray%n6, dm))
         end if
         oke = ieee_is_finite(e%f) .and. ieee_is_finite(e%df)
      end subroutine eval_alg1
      ! the scalar iteration of Algorithm 1 at the frozen direction of ray. istat: 0 root, 1 bracketing failed, 2 no convergence
      ! in NMAX_SCALAR iterations, 3 the state at gamma = 0 cannot be evaluated
      subroutine scalar_alg1(e, istat, nit, nfe)
         type(eval_t), intent(out) :: e
         integer, intent(out) :: istat, nit, nfe
         type(eval_t) :: e1
         real(dp) :: gl, gr, g, gnew, sl, fr
         logical :: oke
         integer :: it, jd
         istat = 0; nit = 0; nfe = 0
         call eval_alg1(0.0_dp, 1, e1, oke); nfe = nfe + 1
         if (.not. oke) then
            istat = 3; return
         end if
         if (e1%f <= mat%eps_f*ytr) then
            ! small correction, gamma = 0
            call eval_alg1(0.0_dp, 2, e, oke); nfe = nfe + 1
            if (.not. oke) istat = 3
            return
         end if
         ! gamma_L = 0, gamma_R by doubling from the Newton estimate
         sl = -e1%df
         if (.not. (sl > 0.0_dp)) then
            istat = 1; return
         end if
         gl = 0.0_dp; gr = e1%f/sl
         fr = huge(1.0_dp)
         do jd = 1, 60
            if (gr > GAMMA_CAP) exit
            call eval_alg1(gr, 0, e, oke); nfe = nfe + 1
            if (.not. oke) exit                     ! a state the elastic law cannot evaluate: no bracket
            fr = e%f
            if (fr <= 0.0_dp) exit
            gr = 2.0_dp*gr
         end do
         if (.not. (fr <= 0.0_dp)) then
            istat = 1; return                       ! bracketing fails
         end if
         ! safeguarded Newton from gamma_L
         g = gl
         do it = 1, NMAX_SCALAR
            nit = it
            if (it > 1) then
               call eval_alg1(g, 1, e1, oke); nfe = nfe + 1
               if (.not. oke) then
                  istat = 3; return
               end if
            end if
            if (abs(e1%f) <= mat%eps_f*ytr) exit
            if (e1%f > 0.0_dp) then
               gl = g
            else
               gr = g
            end if
            if (gr - gl <= 1.0e-16_dp*max(1.0_dp, gr)) exit
            sl = -e1%df
            if (abs(sl) < tiny(1.0_dp)) sl = sign(tiny(1.0_dp), sl)
            gnew = g + e1%f/sl
            if (gnew <= gl .or. gnew >= gr) gnew = 0.5_dp*(gl + gr)
            g = min(max(gnew, gl), gr)
            if (it == NMAX_SCALAR) then
               istat = 2; return
            end if
         end do
         call eval_alg1(g, 2, e, oke); nfe = nfe + 1
         if (.not. oke) istat = 3
      end subroutine scalar_alg1
      subroutine pack_out(sx, dsx)
         real(dp), intent(in) :: sx(3, 3), dsx(6, 6)
         integer :: i, jj
         stress = voigt_from_3(sx)
         do jj = 1, 6
            do i = 1, 6
               dsde(i, jj) = SFAC(i)*CFAC(jj)*dsx(i, jj)
            end do
         end do
      end subroutine pack_out
      subroutine finish()
         if (present(st)) st = s
         if (present(nkap)) nkap = kapit
         if (present(nksolve)) nksolve = kapn
      end subroutine finish
   end subroutine np_update_alg1

   !==========================================================================
   ! Derivatives of the single update with respect to the strain and the history: with U = sqrt(A), A = Cp_n^-1,
   ! the update is
   !     Ce^tr = U C U,   S = U P U,   A' = U B U,   P = E Se E,   B = E^2,   E = exp(-gamma N),
   ! and the core map (Ce^tr, kappa_n) -> (P, B, kappa') is differentiated at the returned state st of
   ! np_update_alg1 (the same call arguments): the closed form of Section 5 for the strain part, and for kappa_n
   ! the same consistency equations with the extra term of the hardening residual,
   !     D_kappa dkappa' = eta dgamma + gamma d(Ms:N) + Y(kappa') dkappa_n,
   ! which gives dgamma = gvec.dCe^tr - (alpha Y/denom) dkappa_n. dU = (U (.) + (.) U)^-1 dA (Sylvester).
   ! An elastic update has P = Se(Ce^tr), B = I, kappa' = kappa_n. The baseline (mat%classical) is differentiated
   ! through its own 7 x 7 Newton system (Y_n = Y(kappa_n) enters the yield equation and the explicit kappa).
   ! All blocks are in the Mandel basis (np_dhist). dh%ok = .false. if the state cannot be differentiated.
   ! For exactly nn updates (eps_theta = 0) the blocks are those of the same formulas evaluated at the returned
   ! state.
   !==========================================================================
   subroutine np_update_dh(mat, eps, cpinv_n, kappa_n, st, dh)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: eps(6), cpinv_n(6), kappa_n
      type(np_state), intent(in) :: st
      type(np_dhist), intent(out) :: dh
      type(ray_t) :: ray
      type(eval_t) :: ev
      real(dp) :: c3(3, 3), u(3, 3), cetr3(3, 3), uhat(6, 6), w, se3(3, 3), c66(6, 6), z6(6), z66(6, 6), gam, &
                  d66(6, 6), rl66(6, 6), rn66(6, 6), s66(6, 6), tnt(6), pbar(6), jbar(6), v(6), denom, gvec(6), &
                  sp(6), cfull(6, 6), nmat(6, 6), da66(6, 6), dz66(6, 6), gk, ck(6), jv(6), nk(6), dak(6), pk(6), &
                  rl2(6, 6), bc(6, 6), bk(6), kc(6), kk, e2(3, 3), p3(3, 3)
      logical :: ok
      integer :: j
      z6 = 0.0_dp; z66 = 0.0_dp
      c3 = c_from_e(eps)
      call sqrtm_spd(voigt_to_3(cpinv_n), u, ok)
      if (.not. ok) return
      cetr3 = matmul(u, matmul(c3, u))
      uhat = sandwich66(u)
      if (.not. st%plastic) then
         !--- elastic: P = Se(Ce^tr), B = I, kappa' = kappa_n
         call np_elastic(mat, cetr3, w, se3, c66, .true., ok)
         if (.not. ok) return
         call blocks(se3, eye3(), c66, z6, z66, z6, z6, 1.0_dp)
         return
      end if
      if (mat%classical) then
         call classical_blocks()
         return
      end if
      !--- the returned state on its last ray
      ray%n6 = st%n6; ray%n3 = from6(ray%n6)
      call sym_eig3(ray%n3, ray%lam, ray%q)
      gam = st%gamma
      call ray_eval(mat, ray, cetr3, kappa_n, gam, 2, ev, ok)
      if (.not. ok) return
      !--- the closed form of Section 5, as in np_update_alg1
      d66 = sandwich66(ev%e3)
      rl66 = expm_frechet66(-gam*ray%lam, ray%q)
      rn66 = -gam*matmul(sandwich2_66(cetr3, ev%e3), rl66)
      call inv6(eye6() - matmul(rn66, matmul(ev%h66, ev%t66)), s66, ok)
      if (.not. ok) return
      tnt = matmul(transpose(ev%t66), ev%n6)
      pbar = matmul(transpose(s66), tnt)
      jv = matmul(transpose(ev%t66), ev%n6 + matmul(transpose(ev%h66), ev%m6))     ! d(Ms:N) = jv.dCe
      jbar = matmul(transpose(s66), jv)
      v = pbar - ev%alpha*gam*jbar
      denom = dot_product(pbar, ev%phat6) + ev%alpha*ev%eta - ev%alpha*gam*dot_product(jbar, ev%phat6)
      if (.not. (abs(denom) > 0.0_dp)) return
      gvec = matmul(transpose(d66), v)/denom
      sp = matmul(s66, ev%phat6)
      cfull = matmul(s66, d66) - outer6(sp, gvec)
      nmat = matmul(ev%h66, matmul(ev%t66, cfull))
      do j = 1, 6
         da66(:, j) = -(ray%n6*gvec(j) + gam*nmat(:, j))
      end do
      dz66 = matmul(sandwich2_66(ev%se3, ev%e3), matmul(rl66, da66)) + matmul(sandwich66(ev%e3), matmul(ev%c66, cfull))
      !--- kappa_n: dgamma, dCe, dN, dA = d(-gamma N), dP
      gk = -ev%alpha*ev%y/denom
      ck = -sp*gk
      nk = matmul(ev%h66, matmul(ev%t66, ck))
      dak = -(ray%n6*gk + gam*nk)
      pk = matmul(sandwich2_66(ev%se3, ev%e3), matmul(rl66, dak)) + matmul(sandwich66(ev%e3), matmul(ev%c66, ck))
      !--- B = E^2 = exp(-2 gamma N)
      rl2 = expm_frechet66(-2.0_dp*gam*ray%lam, ray%q)
      bc = 2.0_dp*matmul(rl2, da66)
      bk = 2.0_dp*matmul(rl2, dak)
      !--- kappa'
      kc = (ev%eta*gvec + gam*matmul(transpose(cfull), jv))/ev%dk
      kk = (ev%eta*gk + gam*dot_product(jv, ck) + ev%y)/ev%dk
      e2 = expm_eig(ray%lam, ray%q, -2.0_dp*gam)
      p3 = matmul(ev%e3, matmul(ev%se3, ev%e3))
      call blocks(p3, e2, dz66, pk, bc, bk, kc, kk)
   contains
      ! from the core derivatives (with respect to Ce^tr and kappa_n) to those with respect to C, A and kappa_n
      subroutine blocks(pp, bb, dpc, dpk, dbc, dbk, dkc, dkk)
         real(dp), intent(in) :: pp(3, 3), bb(3, 3), dpc(6, 6), dpk(6), dbc(6, 6), dbk(6), dkc(6), dkk
         real(dp) :: dsq(6, 6), wc(6, 6), wp(6, 6), wb(6, 6)
         logical :: okb
         call inv6(prodsym66(u), dsq, okb)                   ! dU/dA: U dU + dU U = dA
         if (.not. okb) return
         wc = sandwich2_66(c3, u)                            ! dCe^tr = wc dU + uhat dC
         wp = sandwich2_66(pp, u)                            ! dS = wp dU + uhat dP
         wb = sandwich2_66(bb, u)                            ! dA' = wb dU + uhat dB
         dh%s_c = matmul(uhat, matmul(dpc, uhat))
         dh%s_a = matmul(wp + matmul(uhat, matmul(dpc, wc)), dsq)
         dh%s_k = matmul(uhat, dpk)
         dh%a_c = matmul(uhat, matmul(dbc, uhat))
         dh%a_a = matmul(wb + matmul(uhat, matmul(dbc, wc)), dsq)
         dh%a_k = matmul(uhat, dbk)
         dh%k_c = matmul(transpose(uhat), dkc)
         dh%k_a = matmul(transpose(dsq), matmul(transpose(wc), dkc))
         dh%k_k = dkk
         dh%ok = all(ieee_is_finite(dh%s_a)) .and. all(ieee_is_finite(dh%a_a)) .and. all(ieee_is_finite(dh%s_c)) &
                 .and. all(ieee_is_finite(dh%a_c)) .and. all(ieee_is_finite(dh%k_a)) .and. ieee_is_finite(dh%k_k)
      end subroutine blocks
      ! the baseline: x = (Ce, gamma) solves R1 = Ce - e^A Ce^tr e^A = 0 (A = -gamma N(Ce)), R2 = (seq - Y_n)/Y_n = 0
      subroutine classical_blocks()
         real(dp) :: ce3c(3, 3), wc, sec(3, 3), c66c(6, 6), mc(3, 3), t66c(6, 6), seqc, nc6(6), hc66(6, 6), &
                     lamc(3), qc(3, 3), ea(3, 3), jac(7, 7), lc66(6, 6), deda(6, 6), rhs(7), dxdc(7, 6), dxk(7), &
                     dadc(6, 6), dpdc(6, 6), yn, dyn, etac, sd(6, 6)
         logical :: okc
         integer :: k
         call hardening(mat, kappa_n, yn, dyn)
         gam = st%gamma
         ce3c = from6(st%ce6)
         call np_elastic(mat, ce3c, wc, sec, c66c, .true., okc)
         if (.not. okc) return
         call mandel_stress(ce3c, sec, c66c, mc, t66c, .true.)
         call np_yield(mat, mc, seqc, nc6, hc66, .true.)
         call sym_eig3(from6(nc6), lamc, qc)
         ea = expm_eig(lamc, qc, -gam)
         lc66 = expm_frechet66(-gam*lamc, qc)
         deda = matmul(sandwich2_66(cetr3, ea), lc66)
         jac = 0.0_dp
         jac(1:6, 1:6) = eye6() + gam*matmul(deda, matmul(hc66, t66c))
         jac(1:6, 7) = matmul(deda, nc6)
         jac(7, 1:6) = matmul(nc6, t66c)/yn
         sd = sandwich66(ea)
         do k = 1, 6
            rhs = 0.0_dp; rhs(1:6) = sd(:, k)
            call solve_dense(7, jac, rhs, dxdc(:, k), okc)
            if (.not. okc) return
         end do
         rhs = 0.0_dp; rhs(7) = dyn/yn                          ! -dR2/dkappa_n at seq = Y_n
         call solve_dense(7, jac, rhs, dxk, okc)
         if (.not. okc) return
         dadc = -outer6(nc6, dxdc(7, :)) - gam*matmul(hc66, matmul(t66c, dxdc(1:6, :)))
         dpdc = matmul(sandwich2_66(sec, ea), matmul(lc66, dadc)) + matmul(sandwich66(ea), matmul(c66c, dxdc(1:6, :)))
         dak = -nc6*dxk(7) - gam*matmul(hc66, matmul(t66c, dxk(1:6)))
         pk = matmul(sandwich2_66(sec, ea), matmul(lc66, dak)) + matmul(sandwich66(ea), matmul(c66c, dxk(1:6)))
         rl2 = expm_frechet66(-2.0_dp*gam*lamc, qc)
         bc = 2.0_dp*matmul(rl2, dadc)
         bk = 2.0_dp*matmul(rl2, dak)
         jv = matmul(transpose(t66c), nc6 + matmul(transpose(hc66), to6(mc)))
         etac = dot_product(to6(mc), nc6)
         kc = (etac*dxdc(7, :) + gam*matmul(transpose(dxdc(1:6, :)), jv))/yn
         kk = 1.0_dp + (etac*dxk(7) + gam*dot_product(jv, dxk(1:6)))/yn - gam*etac*dyn/yn**2
         e2 = expm_eig(lamc, qc, -2.0_dp*gam)
         p3 = matmul(ea, matmul(sec, ea))
         call blocks(p3, e2, dpdc, pk, bc, bk, kc, kk)
      end subroutine classical_blocks
   end subroutine np_update_dh

   !==========================================================================
   ! Adaptive substepping of the strain increment (Section 6.5). The Picard iteration contracts only for
   ! sufficiently small increments (Section 7), so an increment whose update does not succeed (info /= NP_OK:
   ! no stress, or the Picard cap reached with eps_theta > 0) is split into m = 2, 4, ..., 2**maxlev equal parts
   ! of the engineering Green-Lagrange strain from eps0 (the converged strain of cpinv_n, kappa_n) to eps, each
   ! integrated from the history of the previous one; the first m at which every part succeeds is kept
   ! (nsub = m; nsub = 1: the single update, its own tangent).
   ! Tangent for m > 1: the exact derivative of the substepped map at fixed m, by the chain rule through the
   ! parts (consistent tangent of a substepping scheme, Perez-Foguet, Rodriguez-Ferran and Huerta 2001). With
   ! E_k = E_0 + (k/m)(E - E_0) and the blocks of np_update_dh at part k,
   !     dA_k/dE = a_c (k/m) dC/dE + a_a dA_k-1/dE + a_k dkappa_k-1/dE,
   !     dkappa_k/dE = k_c (k/m) dC/dE + k_a dA_k-1/dE + k_k dkappa_k-1/dE,
   !     dS/dE = [dS/dE of the last update] + s_a dA_m-1/dE + s_k dkappa_m-1/dE.
   ! fdtan = .true. (optional): central differences of the substepped map at fixed m (step 1e-7), the numerical
   ! tangent of Perez-Foguet et al. (2000), 12 more substepped evaluations; a probe that fails takes the last
   ! part's analytical column. All levels failed: the single update's outputs and status are returned
   ! (nsub = 1). minlev (optional): start at m = 2**minlev - a caller keeps a point's level non-decreasing within
   ! a load step, because a number of parts that changes between global iterations makes the map discontinuous
   ! and stalls the global Newton iteration. noexit (optional): as np_update_alg1.
   ! Thread-safe.
   !==========================================================================
   subroutine np_update_sub(mat, eps0, eps, cpinv_n, kappa_n, stress, dsde, cpinv, kappa, psi, info, nsub, maxlev, st, &
                            minlev, fdtan, noexit)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: eps0(6), eps(6), cpinv_n(6), kappa_n
      real(dp), intent(out) :: stress(6), dsde(6, 6), cpinv(6), kappa, psi
      integer, intent(out) :: info, nsub
      integer, intent(in) :: maxlev
      type(np_state), intent(out), optional :: st
      integer, intent(in), optional :: minlev      ! start at m = 2**minlev substeps (the caller's sticky level)
      logical, intent(in), optional :: fdtan       ! .true.: central differences instead of the exact chain
      logical, intent(in), optional :: noexit      ! .true.: no early exit of a non-contracting iteration
      type(np_state) :: sl, s1st
      real(dp) :: dlast(6, 6), ep(6), sp(6), sm(6), cdum(6), kdum, pdum, ddum(6, 6), h
      real(dp) :: s1(6), d1(6, 6), cp1(6), k1, p1
      integer :: lev, m, j, ip, im, inf1, lev0
      logical :: exact, alg1, nox
      exact = .true.
      if (present(fdtan)) exact = .not. fdtan
      alg1 = .not. mat%classical
      nox = .false.
      if (present(noexit)) nox = noexit
      lev0 = 0
      if (present(minlev)) lev0 = max(0, min(maxlev, minlev))
      nsub = 1
      if (lev0 == 0) then
         call one(eps, cpinv_n, kappa_n, stress, dsde, cpinv, kappa, psi, info, sl)
         if (present(st)) st = sl
         if (info == NP_OK) return
         s1 = stress; d1 = dsde; cp1 = cpinv; k1 = kappa; p1 = psi; inf1 = info; s1st = sl
      end if
      do lev = max(1, lev0), maxlev
         m = 2**lev
         call substepped(m, eps, stress, dlast, cpinv, kappa, psi, info, sl, exact)
         if (info /= NP_OK) cycle
         nsub = m
         if (present(st)) st = sl
         if (exact) then
            dsde = dlast
            return
         end if
         h = 1.0e-7_dp
         do j = 1, 6
            ep = eps; ep(j) = ep(j) + h
            call substepped(m, ep, sp, ddum, cdum, kdum, pdum, ip, sl, .false.)
            ep = eps; ep(j) = ep(j) - h
            call substepped(m, ep, sm, ddum, cdum, kdum, pdum, im, sl, .false.)
            if (ip == NP_OK .and. im == NP_OK) then
               dsde(:, j) = (sp - sm)/(2.0_dp*h)
            else
               dsde(:, j) = dlast(:, j)
            end if
         end do
         return
      end do
      ! every level failed: the single update and its status (with NP_UNCONVERGED the outputs are a valid state)
      if (lev0 > 0) then
         call one(eps, cpinv_n, kappa_n, s1, d1, cp1, k1, p1, inf1, s1st)
      end if
      stress = s1; dsde = d1; cpinv = cp1; kappa = k1; psi = p1; info = inf1
      if (present(st)) st = s1st
   contains
      ! one update: Algorithm 1 or the baseline
      subroutine one(ex, cpx, kx, sx, dx, cpo, ko, po, infx, stx)
         real(dp), intent(in) :: ex(6), cpx(6), kx
         real(dp), intent(out) :: sx(6), dx(6, 6), cpo(6), ko, po
         integer, intent(out) :: infx
         type(np_state), intent(out) :: stx
         if (alg1) then
            call np_update_alg1(mat, ex, cpx, kx, sx, dx, cpo, ko, po, infx, st=stx, noexit=nox)
         else
            call np_update_base(mat, ex, cpx, kx, sx, dx, cpo, ko, po, infx, st=stx)
         end if
      end subroutine one
      subroutine substepped(mm, e1, s1, d1, cp1, k1, p1, inf, stl, chainrule)
         integer, intent(in) :: mm
         real(dp), intent(in) :: e1(6)
         real(dp), intent(out) :: s1(6), d1(6, 6), cp1(6), k1, p1
         integer, intent(out) :: inf
         type(np_state), intent(out) :: stl
         logical, intent(in) :: chainrule
         type(np_dhist) :: dhk
         real(dp) :: cpk(6), kk, ek(6), th, ga(6, 6), gkap(6), gan(6, 6), gkn(6), col(6)
         integer :: k, i, jj
         cpk = cpinv_n; kk = kappa_n
         ga = 0.0_dp; gkap = 0.0_dp                    ! dA_k/dE (A Mandel, E engineering) and dkappa_k/dE
         do k = 1, mm
            th = real(k, dp)/real(mm, dp)
            ek = eps0 + (e1 - eps0)*th
            call one(ek, cpk, kk, s1, d1, cp1, k1, p1, inf, stl)
            if (inf /= NP_OK) return
            if (chainrule) then
               call np_update_dh(mat, ek, cpk, kk, stl, dhk)
               if (.not. dhk%ok) then
                  inf = NP_ERR_NONFINITE; return
               end if
               if (k == mm) then
                  ! dS/dE = the last update's own tangent + the history terms
                  do jj = 1, 6
                     col = matmul(dhk%s_a, ga(:, jj)) + dhk%s_k*gkap(jj)
                     do i = 1, 6
                        d1(i, jj) = d1(i, jj) + SFAC(i)*col(i)
                     end do
                  end do
               else
                  do jj = 1, 6
                     gan(:, jj) = dhk%a_c(:, jj)*(CFAC(jj)*th) + matmul(dhk%a_a, ga(:, jj)) + dhk%a_k*gkap(jj)
                     gkn(jj) = dhk%k_c(jj)*(CFAC(jj)*th) + dot_product(dhk%k_a, ga(:, jj)) + dhk%k_k*gkap(jj)
                  end do
                  ga = gan; gkap = gkn
               end if
            end if
            cpk = cp1; kk = k1
         end do
      end subroutine substepped
   end subroutine np_update_sub

   !==========================================================================
   ! np_raydiag -- diagnostic of Section 8.2: the residual of the scalar iteration of Algorithm 1 on the ray of
   ! direction ndir (Mandel 6; the trial normal when absent) from the trial state.
   !   ratio  seq_tr/Y(kappa_n);  g0 = f(0)/s(0), the Newton estimate from which the bracket is doubled
   !   fmin, gmin   the minimum of f/Y_tr on the ray (geometric march from g0/100 by 10^(1/200), then golden section)
   !   g1, g2       the first root of f and the end of the interval f <= 0 (0, 0 when the ray has no root; g2 = the
   !                last evaluable gamma <= 8 when f stays <= 0)
   !   kdbl   >= 0: the doubling from g0 finds f <= 0 after kdbl doublings; -1: it finds no sign change up to 8
   !   info   0 plastic state; -1 elastic; 1 input / evaluation error
   !==========================================================================
   subroutine np_raydiag(mat, eps, cpinv_n, kappa_n, ratio, g0, gmin, fmin, g1, g2, kdbl, info, ndir, ntrial)
      type(np_material), intent(in) :: mat
      real(dp), intent(in) :: eps(6), cpinv_n(6), kappa_n
      real(dp), intent(out) :: ratio, g0, gmin, fmin, g1, g2
      integer, intent(out) :: kdbl, info
      real(dp), intent(in), optional :: ndir(6)
      real(dp), intent(out), optional :: ntrial(6)
      real(dp), parameter :: GMAXD = 8.0_dp
      type(ray_t) :: ray
      type(eval_t) :: ev
      real(dp) :: c3(3, 3), u(3, 3), cetr3(3, 3), w, se3(3, 3), c66(6, 6), m3(3, 3), t66(6, 6), ntr6(6), htr66(6, 6), &
                  seqtr, ytr, dytr, f0, s0, g, fac, fprev, gprev, glo, ghi, a, b, fa, x1, x2, f1, f2, gpp, dm(6)
      logical :: ok, in
      integer :: j, k
      ratio = 0.0_dp; g0 = 0.0_dp; gmin = 0.0_dp; fmin = 0.0_dp; g1 = 0.0_dp; g2 = 0.0_dp; kdbl = -1; info = 1
      c3 = c_from_e(eps)
      call sqrtm_spd(voigt_to_3(cpinv_n), u, ok)
      if (.not. ok) return
      cetr3 = matmul(u, matmul(c3, u))
      call np_elastic(mat, cetr3, w, se3, c66, .true., ok)
      if (.not. ok) return
      call mandel_stress(cetr3, se3, c66, m3, t66, .true.)
      call np_yield(mat, m3, seqtr, ntr6, htr66, .false.)
      call hardening(mat, kappa_n, ytr, dytr)
      if (present(ntrial)) ntrial = ntr6
      ratio = seqtr/ytr
      if (seqtr - ytr <= 0.0_dp) then
         info = -1; return
      end if
      ray%n6 = ntr6
      if (present(ndir)) ray%n6 = ndir
      ray%n3 = from6(ray%n6)
      call sym_eig3(ray%n3, ray%lam, ray%q)
      call ray_eval(mat, ray, cetr3, kappa_n, 0.0_dp, 1, ev, ok)
      if (.not. ok) return
      f0 = ev%f
      dm = -matmul(ev%t66, ev%phat6)
      s0 = -(dot_product(ray%n6, dm) - ev%alpha*ev%eta)          ! the slope of the paper at gamma = 0
      if (.not. (s0 > 0.0_dp)) return
      g0 = f0/s0
      info = 0
      ! march
      fac = 10.0_dp**(1.0_dp/200.0_dp)
      g = 0.01_dp*g0; gprev = 0.0_dp; fprev = f0; gpp = 0.0_dp
      fmin = f0; gmin = 0.0_dp; glo = 0.0_dp; ghi = g
      in = .false.
      do j = 1, 100000
         if (g > GMAXD) exit
         call ray_eval(mat, ray, cetr3, kappa_n, g, 0, ev, ok)
         if (.not. ok) exit
         if (ev%f < fmin) then
            fmin = ev%f; gmin = g; glo = gprev; ghi = g*fac
         end if
         if (.not. in .and. ev%f <= 0.0_dp) then
            ! first root in (gprev, g]
            a = gprev; b = g
            do k = 1, 200
               x1 = 0.5_dp*(a + b)
               call ray_eval(mat, ray, cetr3, kappa_n, x1, 0, ev, ok)
               if (ev%f > 0.0_dp) then
                  a = x1
               else
                  b = x1
               end if
               if (b - a <= 1.0e-14_dp*b) exit
            end do
            g1 = b; in = .true.; g2 = g
         else if (in .and. ev%f > 0.0_dp) then
            a = gprev; b = g
            do k = 1, 200
               x1 = 0.5_dp*(a + b)
               call ray_eval(mat, ray, cetr3, kappa_n, x1, 0, ev, ok)
               if (ev%f <= 0.0_dp) then
                  a = x1
               else
                  b = x1
               end if
               if (b - a <= 1.0e-14_dp*b) exit
            end do
            g2 = a
            exit
         else if (in) then
            g2 = g
         end if
         gpp = gprev; gprev = g; fprev = ev%f
         g = g*fac
      end do
      if (.not. in) then
         ! refine the minimum by golden section on (glo, ghi)
         a = glo; b = min(ghi, GMAXD)
         do k = 1, 200
            x1 = b - GOLD*(b - a); x2 = a + GOLD*(b - a)
            call ray_eval(mat, ray, cetr3, kappa_n, x1, 0, ev, ok); f1 = ev%f
            if (.not. ok) f1 = huge(1.0_dp)
            call ray_eval(mat, ray, cetr3, kappa_n, x2, 0, ev, ok); f2 = ev%f
            if (.not. ok) f2 = huge(1.0_dp)
            if (f1 < f2) then
               b = x2
            else
               a = x1
            end if
            if (min(f1, f2) < fmin) then
               fmin = min(f1, f2); gmin = merge(x1, x2, f1 < f2)
            end if
            if (b - a <= 1.0e-13_dp*max(b, 1.0e-300_dp)) exit
         end do
         if (fmin <= 0.0_dp) then
            g1 = gmin; g2 = gmin
         end if
      else
         ! the minimum inside the interval f <= 0 is not refined (not needed)
         continue
      end if
      fmin = fmin/ytr
      ! the doubling of Algorithm 1 from the Newton estimate
      g = g0
      do k = 0, 60
         if (g > GMAXD) exit
         call ray_eval(mat, ray, cetr3, kappa_n, g, 0, ev, ok)
         if (.not. ok) exit
         if (ev%f <= 0.0_dp) then
            kdbl = k; exit
         end if
         g = 2.0_dp*g
      end do
      fa = fprev; fa = gpp     ! (unused)
   end subroutine np_raydiag

   ! determinant of a 6 x 6 matrix by elimination with partial pivoting (diagnostic dg_det)
   pure function det6x(a) result(d)
      real(dp), intent(in) :: a(6, 6)
      real(dp) :: d, w(6, 6), t(6)
      integer :: k, p, i1
      w = a; d = 1.0_dp
      do k = 1, 6
         p = k - 1 + maxloc(abs(w(k:6, k)), 1)
         if (w(p, k) == 0.0_dp) then
            d = 0.0_dp; return
         end if
         if (p /= k) then
            t = w(k, :); w(k, :) = w(p, :); w(p, :) = t; d = -d
         end if
         d = d*w(k, k)
         do i1 = k + 1, 6
            w(i1, k:6) = w(i1, k:6) - (w(i1, k)/w(k, k))*w(k, k:6)
         end do
      end do
   end function det6x

end module npicard_plasticity
