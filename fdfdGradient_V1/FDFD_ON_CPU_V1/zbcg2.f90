subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
                  mxmatvec,work,info)

! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)
!
! Improved "vanilla" BiCGStab(2) iterative method
!
! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements 
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on original work of D.R.Fokkema:
!
! subroutine zbistbl v1.1 1998    
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted, 
! provided that the copies are not made or distributed 
! for resale, and that the copyright notice and this 
! notice are retained.
!
! }}
!
! Your bug reports, comments, etc. are welcome: 
! m.a.botchev@math.utwente.nl
!
! ------------------------------
! Description of the parameters:
! ------------------------------
!
! print_resid (input) LOGICAL. If print_resid=.true. the number of 
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve 
!
! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not. 
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are 
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of 
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix 
!            vector multiplications allowed to be done.  On output: 
!            if info>=0, mxmatvec is set to the actual number of matrix 
!            vector multiplications done
!
! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and 
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value 
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
implicit none

! Parameters:
logical,   intent(in)   :: print_resid,nonzero_x
integer,   intent(in)   :: l, n
integer,   intent(inout):: mxmatvec
integer,   intent(out)  :: info
complex*16,intent(inout):: x(n)
complex*16,intent(in)   :: rhs(n)
real*8,    intent(inout):: toler
complex*16,intent(out)  :: work(n,3+2*(l+1))
external                   matvec,precond

! Local variables:
complex*16 :: matrix_z(l+1,l+1),y0(l+1),yl(l+1),zy0(l+1),zyl(l+1)
logical    :: rcmp, xpdt
integer    :: i, j, k, nmatvec
complex*16 :: alpha, beta, omega, rho0, rho1, sigma
complex*16 :: varrho, hatgamma
real*8     :: rnrm0, rnrm
real*8     :: mxnrmx, mxnrmr
complex*16 :: kappa0, kappal

! Aliases for the parts of the work array:
integer          :: rr, r, u, xp, bp

! Constants:
real*8,    parameter :: delta = 1d-2
complex*16,parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

! Functions:
real*8     :: dnorm2_bcg
complex*16 :: zdot_bcg


info = 0

if (l<1 .or. l>2) info = -2
if (n<1) info = -3
if (toler<=0d0) info = -9
if (mxmatvec<0) info = -10

rr = 1
r = rr+1
u = r+(l+1)
xp = u+(l+1)
bp = xp+1

if (info/=0) return 

! Initialize first residual

if (nonzero_x) then
   call matvec (n, x, work(1:n,r))
   work(1:n,r) = rhs - work(1:n,r)
   nmatvec = 1
else
   work(1:n,r) = rhs
   nmatvec = 0
end if
call precond (n, work(1:n,r))

! Initialize iteration loop

work(1:n,rr) = work(1:n,r)
work(1:n,bp) = work(1:n,r)
work(1:n,xp) = x
x = zzero

rnrm0 = dnorm2_bcg (n, work(1:n,r))
rnrm = rnrm0

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zzero
omega = zone
sigma = zone
rho0  = zone

! Iterate

do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)

! =====================
! The BiCG part ---
! =====================

   rho0 = -omega*rho0
   do k=1,l
      rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
      if (rho0.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
      enddo
      call matvec (n, work(1:n,u+k-1), work(1:n,u+k))
      call precond (n, work(1:n,u+k))
      nmatvec = nmatvec+1

      sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
      if (sigma.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      alpha = rho1/sigma
      x(1:n) = alpha*work(1:n,u) + x(1:n)
      do j=0,k-1
         work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
      enddo
      call matvec (n, work(1:n,r+k-1), work(1:n,r+k))
      call precond (n, work(1:n,r+k))
      nmatvec = nmatvec+1
      rnrm = dnorm2_bcg (n, work(1,r))
      mxnrmx = max (mxnrmx, rnrm)
      mxnrmr = max (mxnrmr, rnrm)
   enddo

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

   do i=1,l+1
      do j=1,i
         matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
      end do
   end do

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
   do j=2,l+1
      matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
   end do

!  small vectors y0 and yl

y0(1) = -zone
y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
y0(l+1) = zzero

yl(1) = zzero
yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
yl(l+1) = -zone

!  --- Convex combination

! compute Z*y0 and Z*yl
zy0 = zzero
zyl = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
   zyl = zyl + matrix_z(:,j)*yl(j)
end do

kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)

hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 ) 

y0 = y0 - (hatgamma*kappa0/kappal)*yl


!  --- Update

omega = y0(l+1)

do j=1,l
   work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
   x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
   work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
enddo

! y0 has changed; compute Z*y0 once more
zy0 = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
end do

rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )

! ================================
! The reliable update part ---
! ================================

mxnrmx = max (mxnrmx, rnrm)
mxnrmr = max (mxnrmr, rnrm)
xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
if (rcmp) then
   call matvec (n, x, work(1:n,r))
   call precond (n, work(1:n,r))
   nmatvec = nmatvec + 1
   work(1:n,r) =  work(1:n,bp) - work(1:n,r)
   mxnrmr = rnrm
   if (xpdt) then

      work(1:n,xp) = x(1:n) + work(1:n,xp)
      x = zzero
      work(1:n,bp) = work(1:n,r)

      mxnrmx = rnrm
   endif
endif

if (print_resid) print *,nmatvec,' ',rnrm

enddo
print *,nmatvec,' ',rnrm
! =========================
! End of iterations ---
! =========================

x(1:n) = x(1:n) + work(1:n,xp)

if (rnrm>toler*rnrm0) info = 1

! compute the true residual:

! --------------------- One matvec can be saved by commenting out this:
call matvec (n, x, work(1:n,r) )
work(1:n,r) = rhs(1:n) - work(1:n,r)   
call precond (n, work(1:n,r))
rnrm = dnorm2_bcg(n,work(1:n,r))
nmatvec = nmatvec+1
! --------------------- One matvec can be saved by commenting out this^

toler = rnrm/rnrm0
mxmatvec = nmatvec

end subroutine zbcg2


complex*16 function zdot_bcg(n,zx,zy)

! complex inner product function

implicit none
integer,       intent(in):: n
complex*16,intent(in):: zx(n),zy(n)

zdot_bcg = sum( conjg(zx) * zy )

end function zdot_bcg


real*8 function dnorm2_bcg(n,zx)

! l2 norm function

implicit none
integer,       intent(in):: n
complex*16,intent(in):: zx(n)
complex*16,external  :: zdot_bcg

dnorm2_bcg = sqrt( abs( zdot_bcg(n, zx, zx) ) )

end function dnorm2_bcg
