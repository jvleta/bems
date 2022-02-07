module linalg
  implicit none
  contains

  subroutine solve(A, b, n)
    integer, intent(in) :: n
    real, intent(inout) :: A(n,n)
    real, intent(inout) :: b(n)
  end subroutine

end module