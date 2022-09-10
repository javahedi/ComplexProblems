!-----------------------------------------------------------------------------------------------------------------------------------
!                                                       Entanglement
!-----------------------------------------------------------------------------------------------------------------------------------
!real(8) function eentropy(psi, d, da, db, ssys)  ! Returns the ENTANGLEMENT ENTROPY of a bipartite pure state
!implicit none
!integer :: d, da, db  ! Dimension of whole state space and of the marginal spaces
!complex(8) :: psi(1:d)  ! The bipartite state vector
!complex(8), allocatable :: rho(:,:)  ! The density operator corresponding to psi, i.e., rho = |psi><psi|
!complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! Reduced density operators
!character(1) :: ssys  ! The sub-system whose entropy is to be computed (ssys is 'a' or 'b')
!real(8) :: neumann  ! For the von Neumann's entropy function

!allocate( rho(1:d,1:d) ) ;   call projector(psi, d, rho)
!if ( ssys == 'a' ) then
!  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a) ;   deallocate( rho )
!  eentropy = neumann(da, rho_a) ;   deallocate( rho_a )
!else if ( ssys == 'b' ) then
!  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b) ;   deallocate( rho )
!  eentropy = neumann(db, rho_b) ;   deallocate( rho_b )
!endif

!end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine schmidt_coefficients(psi, d, da, db, schcoeff, eigvec_a, eigvec_b)  ! Returns the Schmidt coeff. for a bipartite pure state
! The matrix evecs has dimension ds x ds and contains in its columns the eigenvectors of rho_s
integer :: d, da, db  ! Dimension of whole state space and of the marginal spaces
complex(8) :: psi(1:d)  ! The bipartite state vector
complex(8), allocatable :: rho(:,:)  ! The density operator corresponding to psi, i.e., rho = |psi><psi|
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! Reduced density operators
real(8), allocatable :: Wa(:), Wb(:)  ! For the eigenvalues of rho_a and rho_b
real(8) :: schcoeff(1:min(da,db))  ! For the Schmidt coefficients
complex(8) :: eigvec_a(1:da,1:da), eigvec_b(1:db,1:db)  ! For the eigenvectors of the reduced density matrices

allocate( rho(1:d,1:d) ) ;   call projector(psi, d, rho)
allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
allocate( rho_b(1:da,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
eigvec_a = rho_a ;   eigvec_b = rho_b ;   deallocate( rho, rho_a, rho_b )
allocate( Wa(1:da), Wb(1:db) ) ;   call lapack_zheevd('V', da, eigvec_a, Wa) ;   call lapack_zheevd('V', db, eigvec_b, Wb)
if ( da <= db ) then !;   allocate( schcoeff(1:da) ) ;
  forall(j=1:da) schcoeff(j) = sqrt(Wa(j))
else if ( da > db ) then !;  allocate( schcoeff(1:db) ) ;
  forall(j=1:db) schcoeff(j) = sqrt(Wb(j))
endif
deallocate( Wa, Wb )

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function concurrence_2qb(rho)  ! Returns the entanglement measure concurrence, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(4,4)  ! Density matrix we want to compute the concurrence
complex(8) :: R(4,4), rho_tilde(4,4), s2_kp_s2(4,4)  ! Auxiliary matrices
complex(8) :: egv(4) ! Eigenvalues of R = rho*rho^tilde
real(8) :: egv_max  ! The greater eigenvalue of R
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, s2_kp_s2) ;   rho_tilde = matmul( matmul(s2_kp_s2,conjg(rho)) , s2_kp_s2 )
R = matmul(rho,rho_tilde) ;   call lapack_zgeev('N', 4, R, egv)

egv_max = max( real(egv(1)), real(egv(2)), real(egv(3)), real(egv(4)) )
 concurrence_2qb = max( 0.d0, (2.d0*sqrt(egv_max)-sqrt(real(egv(1)))-sqrt(real(egv(2)))-sqrt(real(egv(3)))-sqrt(real(egv(4)))))

end
!-----------------------------------------------------------------------------------------------------------------------------------
!(8) function EoF_2qb(rho)  ! Returns the entanglement of formation, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
!implicit none
!complex(8) :: rho(1:4,1:4)  ! Density matrix we want to compute the concurrence
!real(8) :: concurrence_2qb  ! For the concurrence function
!real(8) :: pv(1:2), shannon  ! Probability vector and Shannon's entropy

!pv(1) = (1.d0 + sqrt(1.d0 - concurrence_2qb(rho)**2.d0))/2.d0 ;   pv(2) = 1.d0 - pv(1) ;   EoF_2qb = shannon(2, pv)

!end
!-----------------------------------------------------------------------------------------------------------------------------------
!real(8) function negativity(d, rho_pt)
  ! Returns the entanglement negativity of a "bipartite" system
! This is an simplified version of the subroutine below.
! Here only the partial transposed matrix is given as input.
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement,
! Phys. Rev. A 65, 032314 (2002).
!implicit none
!integer :: d ! Dimension of the state space
!complex(8) :: rho_pt(1:d,1:d)  ! Partial transposed of a state
!real(8) :: norm_tr  ! For the trace norm function

!negativity = 0.5d0*(norm_tr(d, rho_pt) - 1.d0)
! It's equal to the sum of the negative eigenvalues of rho_pt

!--------------------------------------------------
!real(8) function negativity(da, db, ssys, rho)  ! Returns the entanglement negativity of a bipartite system
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
!implicit none
!character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
!integer :: da, db ! Dimensions of the subsystems
!complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
!real(8) :: norm_tr  ! For the trace norm function
!if (ssys == 'a') then
!  call partial_transpose_a(da, db, rho, rho_pt)
!else if (ssys == 'b') then
!  call partial_transpose_b(da, db, rho, rho_pt)
!endif
!negativity = 0.5d0*(norm_tr(da*db, rho_pt) - 1.d0)
!end
!--------------------------------------------------
!end
!-----------------------------------------------------------------------------------------------------------------------------------
!real(8) function log_negativity(d, rho_pt)  ! Returns the entanglement logaritmic negativity of a "bipartite" system
! This is an simplified version of the code below. Here only the partial transposed matrix is given as input.
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
!implicit none
!integer :: d ! Dimension of the state space
!complex(8) :: rho_pt(1:d,1:d)  ! Partial transposed of a state
!real(8) :: norm_tr  ! For the trace norm function
!real(8) :: log2  ! For the log base two
!real(8) :: negativity  ! For the negativity of entanglement

!log_negativity = log2(2.d0*negativity(d, rho_pt)+1.d0)
!log_negativity = log2( norm_tr(d, rho_pt) )

!--------------------------------------------------
!real(8) function log_negativity(da, db, ssys, rho)  ! Returns the entanglement logaritmic negativity of a bipartite system
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
!implicit none
!character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
!integer :: da, db ! Dimensions of the subsystems
!complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
!real(8) :: norm_tr  ! For the trace norm function
!real(8) :: log2  ! For the log base two
!if (ssys == 'a') then
!  call partial_transpose_a(da, db, rho, rho_pt)
!else if (ssys == 'b') then
!  call partial_transpose_b(da, db, rho, rho_pt)
!endif
!log_negativity = log2( norm_tr(da*db, rho_pt) )
!end
!--------------------------------------------------
!end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine entanglement_hs(d, rho_pt, Ehs, css)  ! Returns the Hilbert-Schmidt entanglement of two qudits
! Ref: J. Maziero, Computing partial transposes and related entanglement measures, Braz. J. Phys. 46, 605 (2016),  arXiv:1609.00323
implicit none
integer :: d, dm, dp, dpp  ! For the dimensions (d is the whole system dimension)
complex(8) :: rho_pt(1:d,1:d)  ! The partial transpose of the state under analysis (input). On exit, if css = 'y' and
                               ! Ehs > 0 then the closest separable state (CSS) is returned via this variable
real(8) :: Ehs  ! For the Hilbert-Schmidt entanglement
character(1) :: css  ! If css = 'y' the CSS is computed and returned in rho_pt, if css = 'n' and/or Ehs = 0 the CSS is not
                     ! computed and rho_pt is not modified
complex(8), allocatable :: A(:,:)  ! Auxiliary variable for sending the PT to Lapack (it may returns the eigenvectors on exit)
real(8), allocatable :: W(:), Wd(:)  ! For the eigenvalues of the PT, in ascending and descending order, respectively
real(8) :: sw, sn1, sn2, sp1, sp2, xi  ! Auxiliary variable for the sums of eigenvalues and for xi
integer :: j  ! Auxiliary variable for counters
complex(8), allocatable :: proj(:,:)  ! Auxiliary variable for projectors

allocate( A(1:d,1:d), W(1:d) ) ;   A = rho_pt  ! Computes the eigenvalues and eigenvectors of the PT
if ( css == 'n' ) then ; call lapack_zheevd('N', d, A, W) ; else if ( css == 'y' ) then ; call lapack_zheevd('V', d, A, W) ; endif

j = 0 ;   dm = 0 ;   do ;   j = j + 1 ;  if ( W(j) >= 0.d0 ) exit ;   dm = j ;   enddo  ! Computes d-

if ( dm == 0 ) then  ! In this case Ehs is null
  Ehs = 0.d0
else if ( dm > 0 ) then  ! Computes the auxiliary dimension and Ehs
  sn1 = 0.d0 ;   sn2 = 0.d0 ;   do j = 1, dm ;   sn1 = sn1 + dabs(W(j)) ;   sn2 = sn2 + (W(j))**2.d0 ;   enddo
  sp1 = 0.d0 ;   sp2 = 0.d0 ;  dp = d - dm ;   allocate( Wd(1:d) ) ;   forall ( j = 1:d ) Wd(j) = W(d-j+1)
  sw = 0.d0 ;   j = 0 ;   do ;   j = j + 1 ;   sw = sw + Wd(j)   ;  if ( (sw > 1.d0) .or. (j >= dp) ) exit ; enddo ;   dpp = j
  if ( dp > dpp ) then ;   do j = dpp+1 , dp ;   sp1 = sp1 + Wd(j) ;   sp2 = sp2 + (Wd(j))**2.d0 ;   enddo ;    endif
  Ehs = dsqrt( (sn1 - sp1)**2.d0 + sp2 + sn2 )
  if ( css == 'y' ) then  ! For computing the closest separable state
    allocate( proj(1:d,1:d) ) ;   rho_pt = 0.d0 ;   xi = 0.d0
    do j = 1, dpp-1 ;   call projector(A(:,d-j+1), d, proj) ;   rho_pt = rho_pt + Wd(j)*proj ;   xi = xi + Wd(j) ;   enddo
    j = dpp ;   xi = 1.d0 - xi ;   call projector(A(:,d-j+1), d, proj) ;   rho_pt = rho_pt + xi*proj ;   deallocate( proj )
  endif ;   deallocate( Wd )
endif

deallocate( A, W )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_c(M1, nr1, nc1, M2, nr2, nc2, M1_kp_M2)  ! Returns the tensor product of two general complex matrices
implicit none
integer :: nr1, nc1, nr2, nc2  ! Number of rows and columns of the two matrices
complex(8) :: M1(1:nr1,1:nc1), M2(1:nr2,1:nc2)  ! Matrices to take the tensor product of
complex(8) :: M1_kp_M2(1:nr1*nr2,1:nc1*nc2)  ! Matrix containing the tensor product of M1 and M2
integer :: i, j  ! Auxiliary variables for counters

M1_kp_M2 = 0.d0 ;   forall ( i = 1:nr1 , j = 1:nc1 ) M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2

end
!-----------------------------------------------------------------------------------------------------------------------------------

subroutine identity_c(d, identity)  ! Returns the complex dxd identity matrix
implicit none
integer :: d  ! Dimension of the identity matrix
complex(8) :: identity(1:d,1:d)  ! Identity matrix
integer :: j, k  ! Auxiliary variable for counters

forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = (0.d0,0.d0) ;   forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = (1.d0,0.d0)

end


!-----------------------------------------------------------------------------------------------------------------------------------
subroutine pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)  ! Defines the three Pauli's matrices and the identity matrix
implicit none
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

sigma_0 = 0.d0 ;   sigma_0(1,1) = 1.d0 ;   sigma_0(2,2) = 1.d0
sigma_1 = 0.d0 ;   sigma_1(1,2) = 1.d0 ;   sigma_1(2,1) = 1.d0
sigma_2 = 0.d0 ;   sigma_2(1,2) = -(0.d0,1.d0) ;   sigma_2(2,1) = (0.d0,1.d0)
sigma_3 = 0.d0 ;   sigma_3(1,1) = 1.d0 ;   sigma_3(2,2) = -1.d0

end



!-----------------------------------------------------------------------------------------------------------------------------------
subroutine projector(vec, d, proj)  ! Returns a PROJECTOR on the provided COMPLEX vector
implicit none
integer :: d  ! Dimension of the vector
complex(8) :: vec(1:d)  ! Vector we want the projector on
complex(8) :: proj(1:d,1:d)  ! Projector on vec
integer :: j, k  ! Auxiliary variables for counters

forall ( j=1:d, k=1:d, j <= k ) proj(j,k) = vec(j)*conjg(vec(k))  ! Elements in the diagonal and above
forall ( j=1:d, k=1:d, j > k ) proj(j,k) = conjg(proj(k,j))  ! Elements below the diagonal

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine projector_re(vec, d, proj)  ! Returns a PROJECTOR on the provided REAL vector
implicit none
integer :: d  ! Dimension of the vector
real(8) :: vec(1:d)  ! Vector we want the projector on
real(8) :: proj(1:d,1:d)  ! Projector on vec
integer :: j, k  ! Auxiliary variables for counters

forall ( j=1:d, k=1:d, j <= k ) proj(j,k) = vec(j)*vec(k)  ! Elements in the diagonal and above
forall ( j=1:d, k=1:d, j > k ) proj(j,k) = proj(k,j)  ! Elements below the diagonal

end



!-----------------------------------------------------------------------------------------------------------------------------------
!                                                     LAPACK callers
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zgeev(JOBVR, N, A, Wc)  ! Calls LAPACK's eigensolver for GENERAL complex matrices
! ZGEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
!character(1) :: JOBVL  !  JOBVL is CHARACTER*1; = 'N': left eigenvectors of A are not computed; = 'V': left eigenvectors of are computed.
character(1) :: JOBVR  !  JOBVR is CHARACTER*1; = 'N': right eigenvectors of A are not computed;  = 'V': right eigenvectors of A are computed.
character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA = N  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(1:N,1:N)  ! A is COMPLEX*16 array, dimension (LDA,N); On entry, the N-by-N matrix A. On exit, A has been overwritten.
complex(8) :: Wc(1:N)  ! Wc is COMPLEX*16 array, dimension (N). Wc contains the computed eigenvalues.
!integer :: LDVL = N  !LDVL is INTEGER. The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
complex(8) :: VL(1:N,1:N)  ! VL is COMPLEX*16 array, dimension (LDVL,N); If JOBVL = 'V', the left eigenvectors u(j) are stored one
! after another in the columns of VL, in the same order as their eigenvalues.
! If JOBVL = 'N', VL is not referenced. u(j) = VL(:,j), the j-th column of VL.
!integer :: LDVR = N  ! LDVR is INTEGER. The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
complex(8) :: VR(1:N,1:N)  ! VR is COMPLEX*16 array, dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one
! after another in the columns of VR, in the same order their eigenvalues.
! If JOBVR = 'N', VR is not referenced. v(j) = VR(:,j), the j-th column of VR.

!integer :: LWORK = 2*N  ! LWORK is INTEGER; The dimension of the array WORK.  LWORK >= max(1,2*N). For good performance, LWORK must generally be larger.
! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
! this value as the first entry of the WORK array, and no error related to LWORK is issued by XERBLA.
complex(8) :: WORK(1:2*N)  ! WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
real(8) :: RWORK(1:2*N)  ! RWORK is DOUBLE PRECISION array, dimension (2*N)
integer :: INFO   ! INFO is INTEGER
! = 0:  successful exit
! < 0:  if INFO = -i, the i-th argument had an illegal value.
! > 0:  if INFO = i, the QR algorithm failed to compute all the
! eigenvalues, and no eigenvectors have been computed; elements and i+1:N of W contain eigenvalues which have converged.

call zgeev ('N',   JOBVR, N, A, N,   Wc, VL, N,    VR, N,    WORK, 2*N,   RWORK, INFO)
!call zgeev (JOBVL, JOBVR, N, A, LDA, Wc, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zheevd(JOBZ, N, A, W)  ! Calls LAPACK's eigensolver for HERMITIAN complex matrices
! ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
character(1) :: JOBZ  ! JOBZ is CHARACTER*1; = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(N,N)  ! A is COMPLEX array, dimension (LDA, N). On entry, the Hermitian matrix A.  If UPLO = 'U', the
! leading N-by-N upper triangular part of A contains the upper triangular part of the matrix A.
! If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower triangular part of the
! matrix A. On exit, if JOBZ = 'V', then if INFO = 0, A contains the orthonormal eigenvectors of the
! matrix A. If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') or the upper triangle
! (if UPLO='U') of A, including the diagonal, is destroyed.
real(8) :: W(N)  ! W is REAL array, dimension (N). If INFO = 0, the eigenvalues in ascending order.
integer :: LWORK  ! LWORK is INTEGER
!The length of the array WORK.
!If N <= 1,                LWORK must be at least 1.
!If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
!If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.

!If LWORK = -1, then a workspace query is assumed; the routine
!only calculates the optimal sizes of the WORK, RWORK and
!IWORK arrays, returns these values as the first entries of
!the WORK, RWORK and IWORK arrays, and no error message
!related to LWORK or LRWORK or LIWORK is issued by XERBLA.
complex(8), allocatable :: WORK(:)  !WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: LRWORK   ! LRWORK is INTEGER
!The dimension of the array RWORK.
!If N <= 1,                LRWORK must be at least 1.
!If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
!If JOBZ  = 'V' and N > 1, LRWORK must be at least
!              1 + 5*N + 2*N**2.

!If LRWORK = -1, then a workspace query is assumed; the
!routine only calculates the optimal sizes of the WORK, RWORK
!and IWORK arrays, returns these values as the first entries
!of the WORK, RWORK and IWORK arrays, and no error message
!related to LWORK or LRWORK or LIWORK is issued by XERBLA.
real(8), allocatable :: RWORK(:)    ! RWORK is DOUBLE PRECISION array,
! dimension (LRWORK)
! On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
integer :: LIWORK   ! LIWORK is INTEGER
!The dimension of the array IWORK.
!If N <= 1,                LIWORK must be at least 1.
!If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
!If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.

!If LIWORK = -1, then a workspace query is assumed; the
!routine only calculates the optimal sizes of the WORK, RWORK
!and IWORK arrays, returns these values as the first entries
!of the WORK, RWORK and IWORK arrays, and no error message
!related to LWORK or LRWORK or LIWORK is issued by XERBLA.
integer, allocatable :: IWORK(:)   ! IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
integer :: INFO   ! INFO is INTEGER; = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value
! > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed to converge; i off-diagonal elements
! of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = 'V', then the
! algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns
! INFO/(N+1) through mod(INFO,N+1).

if (JOBZ == 'N') then ;   LWORK = N + 1 ;   LRWORK = N ;   LIWORK = 1
else if (JOBZ == 'V') then ;   LWORK = 2*N + N**2 ;   LRWORK = 1 + 5*N + 2*N**2 ;   LIWORK = 3 + 5*N
endif
allocate( WORK(1:LWORK), RWORK(1:LRWORK), IWORK(1:LIWORK) )

call zheevd(JOBZ,'U',N,A,N,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
!call zheevd(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)

deallocate( WORK, RWORK, IWORK )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_dsyevd(JOBZ, N, A, W)  ! Calls LAPACK's eigensolver for SYMMETRIC real matrices
! DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
! real symmetric matrix A. If eigenvectors are desired, it uses a
! divide and conquer algorithm.!

! The divide and conquer algorithm makes very mild assumptions about
! floating point arithmetic. It will work on machines with a guard
! digit in add/subtract, or on those binary machines without guard
! digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
! Cray-2. It could conceivably fail on hexadecimal or decimal machines
! without guard digits, but we know of none.

! Because of large use of BLAS of level 3, DSYEVD needs N**2 more
! workspace than DSYEVX.
character(1) :: JOBZ  ! JOBZ is CHARACTER*1; = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
real(8) :: A(1:N,1:N)  ! A is DOUBLE PRECISION array, dimension (LDA, N)
!On entry, the symmetric matrix A.  If UPLO = 'U', the
!leading N-by-N upper triangular part of A contains the
!upper triangular part of the matrix A.  If UPLO = 'L',
!the leading N-by-N lower triangular part of A contains
!the lower triangular part of the matrix A.
!On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!orthonormal eigenvectors of the matrix A.
!If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!or the upper triangle (if UPLO='U') of A, including the
!diagonal, is destroyed.
real(8) :: W(1:N)  ! W is DOUBLE PRECISION array, dimension (N)
!If INFO = 0, the eigenvalues in ascending order.
integer :: LWORK  ! LWORK is INTEGER
!The dimension of the array WORK.
!If N <= 1,               LWORK must be at least 1.
!If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
!If JOBZ = 'V' and N > 1, LWORK must be at least
!                                      1 + 6*N + 2*N**2.

!If LWORK = -1, then a workspace query is assumed; the routine
!only calculates the optimal sizes of the WORK and IWORK
!arrays, returns these values as the first entries of the WORK
!and IWORK arrays, and no error message related to LWORK or
!LIWORK is issued by XERBLA.
real(8), allocatable :: WORK(:)  !WORK is DOUBLE PRECISION array,
!                              dimension (LWORK)
!On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: LIWORK   ! LIWORK is INTEGER
!The dimension of the array IWORK.
!If N <= 1,                LIWORK must be at least 1.
!If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
!If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.!

! If LIWORK = -1, then a workspace query is assumed; the
! routine only calculates the optimal sizes of the WORK and
! IWORK arrays, returns these values as the first entries of
! the WORK and IWORK arrays, and no error message related to
! LWORK or LIWORK is issued by XERBLA.
integer, allocatable :: IWORK(:)   ! IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
integer :: INFO   ! INFO is INTEGER; = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value
! > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed to converge; i off-diagonal elements
! of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = 'V', then the
! algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns
! INFO/(N+1) through mod(INFO,N+1).

if (JOBZ == 'N') then ;   LWORK = 2*N+1 ;   LIWORK = 1
else if (JOBZ == 'V') then ;   LWORK =  1 + 6*N + 2*N**2 ;   LIWORK = 3 + 5*N
endif
allocate( WORK(1:LWORK), IWORK(1:LIWORK) )

call dsyevd(JOBZ,'U',N,A,N,W,WORK,LWORK,IWORK,LIWORK,INFO)
!call dsyevd(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO)

deallocate(WORK, IWORK)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zgeqrfp(N, Z, Q, R)  ! Calls LAPACK QR factorization
! For a complex square matrix A, it returns the orthogonalized Q and upper triangular R matrices
implicit none
!ZGEQRFP computes a QR factorization of a complex M-by-N matrix A:
! A = Q * R. The diagonal entries of R are real and nonnegative.

!integer :: M !is INTEGER
!The number of rows of the matrix A.  M >= 0.
integer :: N !is INTEGER
!The number of columns of the matrix A.  N >= 0.
!integer :: LDA !is INTEGER
!          !The leading dimension of the array A.  LDA >= max(1,M).
complex(8) :: A(1:N,1:N) !is COMPLEX*16 array, dimension (LDA,N)
!On entry, the M-by-N matrix A.
!On exit, the elements on and above the diagonal of the array
!contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!upper triangular if m >= n); the elements below the diagonal,
!with the array TAU, represent the orthogonal matrix Q as a
!product of min(m,n) elementary reflectors (see Further
!Details).
complex(8) :: TAU(1:N) !is COMPLEX*16 array, dimension (min(M,N))
!The scalar factors of the elementary reflectors
!Further Details
!The matrix Q is represented as a product of elementary reflectors
!  Q = H(1) H(2) . . . H(k), where k = min(m,n).
!Each H(i) has the form
!  H(i) = I - tau * v * v**H
!where tau is a real scalar, and v is a real vector with
!v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!and tau in TAU(i).
!integer :: LWORK !is INTEGER
!The dimension of the array WORK. The dimension can be divided into three parts.
! 1) The part for the triangular factor T. If the very last T is not bigger
!    than any of the rest, then this part is NB x ceiling(K/NB), otherwise,
!   NB x (K-NT), where K = min(M,N) and NT is the dimension of the very last T
!2) The part for the very last T when T is bigger than any of the rest T.
!   The size of this part is NT x NT, where NT = K - ceiling ((K-NX)/NB) x NB,
!   where K = min(M,N), NX is calculated by
!        NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
! 3) The part for dlarfb is of size max((N-M)*K, (N-M)*NB, K*NB, NB*NB)
! So LWORK = part1 + part2 + part3
! If LWORK = -1, then a workspace query is assumed; the routine
! only calculates the optimal size of the WORK array, returns
! this value as the first entry of the WORK array, and no error
! message related to LWORK is issued by XERBLA.
complex(8) :: WORK(1:3*N*N) !is COMPLEX*16 array, dimension (MAX(1,LWORK))
!On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: INFO !is INTEGER
!= 0:  successful exit
!< 0:  if INFO = -i, the i-th argument had an illegal value
integer :: j, k  ! Auxiliary variables
complex(8) :: identity(1:N,1:N), proj(1:N,1:N), v(1:N)
complex(8) :: Z(1:N,1:N), Q(1:N,1:N), R(1:N,1:N)  ! Input and output matrices

A = Z
call zgeqrf(N, N, A, N, TAU, WORK, 3*N*N, INFO) ! It did not find zgeqrfp in the last compilation
!call zgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO)

R = 0.d0 ;   forall (j = 1:N, k = 1:N , k >= j ) R(j,k) = A(j,k)

! Computing the unitary
call identity_c(N, identity) ;   Q = identity
do j = 1, N
if (j > 1) v(1:j-1) = 0.d0 ;   v(j) = 1.d0 ;   v(j+1:N) = A(j+1:N,j) ;   call projector(v, N, proj)
Q = matmul(Q,(identity-TAU(j)*proj))
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------


