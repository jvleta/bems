module bem
  implicit none
  real, parameter:: pi = 4.0*atan(1.0)
  integer, parameter:: max_size = 1000
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
  subroutine CLAP1(N, xm, ym, xb, yb, nx, ny, lg, BCT, BCV, phi, dphi)
    integer:: m, k, N, BCT(max_size)
    real:: xm(max_size), ym(max_size), xb(max_size), yb(max_size)
    real:: nx(max_size), ny(max_size), lg(max_size), BCV(max_size), Z(max_size)
    real:: A(max_size, max_size), B(max_size), phi(max_size), dphi(max_size)
    real:: PF1, PF2, del, F1, F2
  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine CLAP2(N, xi, eta, xb, yb, nx, ny, lg, phi, dphi, pint)
    integer:: N, i
    real, intent(in):: xi, eta, xb(max_size), yb(max_size), nx(max_size), ny(max_size)
    real, intent(in):: lg(max_size), phi(max_size), dphi(max_size), pint
  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine solver(A, B, N, lud, Z)
    integer:: lda, N, info, lud, IDAMAX, j, k, kp1, l, nm1, kb
    integer:: ipvt(max_size)
    real:: A(max_size, max_size), B(max_size), Z(max_size)
  end subroutine
  !
  !
  !
  subroutine run_chap1_ex1()
  end subroutine
  !
  !----------------------------------------------------------------------------
  !
  subroutine run_chap1_ex2()
  end subroutine 
  !
  !----------------------------------------------------------------------------
  !
end module
