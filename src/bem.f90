module bem
  implicit none
  real, parameter:: pi = 4.0*atan(1.0)
  contains
    subroutine dostuff(x)
     real:: x
     write(*,*) "Hello from Fortran land"
     write(*,*) "x =",x
    end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine CPF(xi, eta, xk, yk, nkx, nky, L, PF1, PF2)
   real, intent(in):: xi, eta, xk, yk, nkx, nky, L
   real:: A, B, E, D, BA, EA
   real, intent(out):: PF1, PF2

   A = L**2.0
   B = 2.0*L * (-nky * (xk-xi) * nkx * (yk-eta))
   E = (xk-xi)**2.0 + (yk-eta)**2.0
   D = sqrt(abs(4.0*A * E-B**2.0))
   BA = B/A
   EA = E/A

   if (D < 1.0e-12) then
    PF1 = 0.5*L * (log(L) + (1.0+0.5*BA) * log(abs(1.0+0.5*BA)) &
        - 0.5*BA*log(abs(0.5*BA)) - 1.0)
    PF2 = 0.0
   else
    PF1 = 0.25*L * (2.0*(log(L)-1.0)-0.5*BA*log(abs(EA)) &
        + (1.0+0.5*BA)*log(abs(1.0+BA+EA))+(D/A)*(atan((2.0*A+B)/D)-atan(B/D)))
    PF2 = L*(nkx*(xk-xi)+nky*(yk-eta))/D*(atan((2.0*A+B)/D)-atan(B/D))
   end if

  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine CELAP1(N, xm, ym, xb, yb, nx, ny, lg, BCT, BCV, phi, dphi)
    integer:: m, k, N, BCT(N)
    real:: xm(N), ym(N), xb(N+1), yb(N+1)
    real:: nx(N), ny(N), lg(N), BCV(N), Z(N)
    real:: A(N, N), B(N), phi(N), dphi(N)
    real:: PF1, PF2, del, F1, F2
    
    do m = 1, N
      B(m) = 0.0
      do k = 1, N
        call CPF(xm(m), ym(m), xb(k), yb(k), nx(k), ny(k), lg(k), PF1, PF2)
        F1 = PF1/pi
        F2 = PF2/pi
        if (k == m) then
          del = 1.0
        else
          del = 0.0
        end if
        if (BCT(k) == 0) then
          A(m, k) = -F1
          B(m) = B(m) + BCV(k)*(-F2+0.5*del)
        else
          A(m, k) = F2-0.5*del
          B(m) = B(m) + BCV(k) * F1
        end if
      end do
    end do

    call solver(A, B, N, 1, Z)

    do m = 1, N
      if (BCT(m) == 0) then
        phi(m) = BCV(m)
        dphi(m) = Z(m)
      else
        phi(m) = Z(m)
        dphi(m) = BCV(m)
      end if
    end do

  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine CELAP2(N, xi, eta, xb, yb, nx, ny, lg, phi, dphi, pint)
    integer:: N, i
    real, intent(in):: xi, eta, xb(N+1), yb(N+1), nx(N), ny(N)
    real, intent(in):: lg( N), phi(N), dphi(N)
    real:: sum, pint, PF1, PF2

    sum = 0.0

    do i = 1,N
      call CPF(xi, eta, xb(i), yb(i), nx(i), ny(i), lg(i), PF1, PF2)
      sum = sum + phi(i) * PF2 - dphi(i) * PF1
    end do

    pint = sum / pi
  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine solver(A, B, N, lud, Z)
    integer, intent(in):: N
    integer:: lda, info, lud, IDAMAX, i, j, k, kp1, l, nm1, kb
    integer:: ipvt(N)
    real:: A(N, N), AMD(N, N), B(N), Z(N), t

!    common/ludcmp/ipvt, AMD
    
    nm1 = N - 1

    do i = 1, N
      Z(i) = B(i)
    end do

    if (lud == 0) goto 99

    do i = 1, N
      do j = 1, N
        AMD(i, j) = A(i, j)
      end do
    end do

    info = 0
    
    if (nm1 < 1) goto 70

    do k = 1, nm1
      kp1 = k + 1
      l = IDAMAX(N - k + 1, AMD(k, k), 1) + k-1
      ipvt(k) = l
      if (AMD(l, k) == 0.0) goto 40
      if (l == k) goto 10
      t = AMD(l, k)
      AMD(l,k) = AMD(k, k)
      AMD(k,k) = t
  10  continue
      t = -1.0/AMD(k, k)
      call DSCAL(N-k, t, AMD(k+1, k), 1)
      do j = kp1, N
        t = AMD(l, j)
        if (l == k) goto 20
        AMD(l, j) = AMD(k, j)
        AMD(k, j) = t
  20  continue
        call DAXPY(N-k, t, AMD(k+1, k), 1, AMD(k+1, j), 1)
    end do
    goto 50
  40 continue
    info = k
  50 continue
    end do
  70 continue
    ipvt(N) = N

    if (AMD(N, N) == 0.0) info = N
    if (info /= 0) write(*,*) "Division by zero in SOLVER"
  99 continue

    if (nm1 < 1) goto 130

    do k = 1, nm1
      l = ipvt(k)
      t = Z(l)
      if (l /= k) then
        Z(l) = Z(k)
        Z(k) = t
      end if

      call DAXPY(N-k, t, AMD(k+1, k), 1, Z(k+1), 1)
  
    end do

  130 continue

    do kb = 1, N
      k = N+1 - kb
      Z(k) = Z(k) / AMD(k, k)
      t = -Z(k)
      call DAXPY(k-1, t, AMD(1, k), 1, Z(1), 1)
    end do

  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine chap1ex1(NO)
    integer:: NO, N, i, ians
    integer:: BCT(4*NO)
    real:: xb(4*NO+1), yb(4*NO+1), xm(4*NO), ym(4*NO), nx(4*NO), ny(4*NO), lg(4*NO), BCV(4*NO)
    real:: phi(4*NO), dphi(4*NO), pint, dl, xi, eta

    N = 4 * NO
    dl = 1.0/real(NO)

    do i = 1, NO
      xb(i)= real(i-1)*dl
      yb(i)=0.0
      xb(NO+i)=1.0
      yb(NO+i)=xb(i)
      xb(2*NO+i)=1.0-xb(i)
      yb(2*NO+i)=1.0
      xb(3*NO+i)=0.0
      yb(3*NO+i)=1.0-xb(i)
    end do

    xb(N+1)=xb(1)
    yb(N+1)=yb(1)

    do i = 1, N
      xm(i) = 0.5 * (xb(i) + xb(i+1))
      ym(i) = 0.5 * (yb(i) + yb(i+1))
      lg(i)=sqrt((xb(i+1)-xb(i))**2.0 + (yb(i+1)-yb(i))**2.0)
      nx(i) = (yb(i+1) - yb(i)) / lg(i)
      ny(i) = (xb(i) - xb(i+1)) / lg(i)
    end do

    do i = 1, N
      if (i <= NO) then
        BCT(i) = 1
        BCV(i)=0.0
      else if ((i > NO) .and. (i <= (2.0*NO))) then
        BCT(i) = 0
        BCV(i) = cos(pi*ym(i))
      else if ((i > (2.0*NO)) .and. (i <= (3.0*NO))) then 
        BCT(i) = 1
        BCV(i) = 0.0
      else
        BCT(i) = 0
        BCV(i) = 0.0
      end if
    end do

    call CELAP1(N, xm, ym, xb, yb, nx, ny, lg, BCT, BCV, phi, dphi)

    xi = 0.1
    eta = 0.2

    call CELAP2(N, xi, eta, xb, yb, nx, ny, lg, phi, dphi, pint)

    write(*,*) pint

  end subroutine

end module
