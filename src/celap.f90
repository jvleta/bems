module celap
  use iso_c_binding
  use scifor
  implicit none
  private
  public:: CELAP1
  public:: CELAP2
  contains
  subroutine run_celap_ex1() bind(c, name='run_celap_ex1')
    use iso_c_binding
    call celap_ex1(20)
  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine CPF(xi, eta, xk, yk, nkx, nky, L, PF1, PF2)
   real(8), intent(in):: xi, eta, xk, yk, nkx, nky, L
   real(8):: A, B, E, D, BA, EA
   real(8), intent(out):: PF1, PF2

   A = L**2.0
   B = 2.0*L * (-nky * (xk-xi) + nkx * (yk-eta))
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
    real(8):: xm(N), ym(N), xb(N+1), yb(N+1)
    real(8):: nx(N), ny(N), lg(N), BCV(N), Z(N)
    real(8):: A(N, N), B(N), phi(N), dphi(N)
    real(8):: PF1, PF2, del, F1, F2
    
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
    
    Z = B
    call solve(A, Z)

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
    real(8), intent(in):: xi, eta, xb(N+1), yb(N+1), nx(N), ny(N)
    real(8), intent(in):: lg( N), phi(N), dphi(N)
    real(8):: sum, pint, PF1, PF2

    sum = 0.0

    do i = 1, N
      call CPF(xi, eta, xb(i), yb(i), nx(i), ny(i), lg(i), PF1, PF2)
      sum = sum+phi(i) * PF2-dphi(i) * PF1
    end do

    pint = sum/pi
  
  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine celap_ex1(NO)
    integer:: NO, N, i, ians
    integer:: BCT(4*NO)
    real(8):: xb(4*NO+1), yb(4*NO+1), xm(4*NO), ym(4*NO), nx(4*NO), ny(4*NO), lg(4*NO), BCV(4*NO)
    real(8):: phi(4*NO), dphi(4*NO), pint, dl, xi, eta
    real(8):: xi_pts(9), eta_pts(9)

    open(unit = 123, file="../tests/celap/celap_example1.out")

    N = 4*NO
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

    xi_pts = (/0.10, 0.10, 0.10, 0.50, 0.50, 0.50, 0.90, 0.90, 0.90/)
    eta_pts = (/0.20, 0.30, 0.40, 0.20, 0.30, 0.40, 0.20, 0.30, 0.40/)
    
    write(123, 100) "number of elements:", N
100 format(A, I4)
    do i = 1, 9
      xi = xi_pts(i)
      eta = eta_pts(i)
      call CELAP2(N, xi, eta, xb, yb, nx, ny, lg, phi, dphi, pint)
      write(123, 200) xi, eta, pint
200   format(F6.3, F6.3, F10.6)
    end do
    write(123, *) ""
    close(123)

  end subroutine

end module
