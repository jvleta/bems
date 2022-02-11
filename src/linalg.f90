module linalg
   implicit none
   type :: PermutationResults
      integer :: count = 0
      integer, allocatable :: indices(:)
   end type

   type :: LUDecompositionResults
      real, allocatable :: LU(:,:)
      type(PermutationResults) :: permutation
   end type


contains


   function lu_decompose(A, b, n) result(r)
      integer, intent(in) :: n
      real, intent(in) :: A(n, n)
      real, intent(in) :: b(n)
      type(LUDecompositionResults) :: r
      integer :: row_index_of_max_value = 0
      real :: max_value_in_column = 0.0
      integer :: i, col_index, row_index, temp_int
      real :: val = 0.0, temp_real

      allocate(r%LU(n, n))
      allocate(r%permutation%indices(n))

      r%LU = A

      do i=1,n
         row_index_of_max_value = i
         max_value_in_column = 0.0
         do row_index=i,n
            val = abs(r%LU(row_index, i))
            if (val > max_value_in_column) then
               row_index_of_max_value = row_index
               max_value_in_column = val
            end if
         end do

         if (max_value_in_column < 0.001) then
            write(*,*) "You get nothing!"
         end if

         if (row_index_of_max_value /= i) then
            temp_int = r%permutation%indices(i)
            r%permutation%indices(i) = r%permutation%indices(row_index_of_max_value)
            r%permutation%indices(row_index_of_max_value) = temp_int

            do col_index=1,n
               temp_real = r%LU(i, col_index)
               r%LU(i, col_index) = r%LU(row_index_of_max_value, col_index)
               r%LU(row_index_of_max_value, col_index) = temp_real
            end do
            
            r%permutation%count = r%permutation%count + 1

         end if

         do row_index=i + 1,n
            r%LU(row_index, i) = r%LU(row_index, i) / r%LU(i, i)
            do col_index=i + 1,n
               r%LU(row_index, col_index) = r%LU(row_index, col_index) - r%LU(row_index, i) * r%LU(i, col_index)
            end do
         end do
      end do


   end function

   subroutine solve(A, b, n)
      integer, intent(in) :: n
      real, intent(inout) :: A(n,n)
      real, intent(inout) :: b(n)

   end subroutine

end module
