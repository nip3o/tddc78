program laplsolv
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem
! on a square using the Jacobi method.
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
<<<<<<< HEAD
  integer, parameter                  :: n = 100, maxiter = 10
=======
  integer, parameter                  :: n = 1000, maxiter = 1000
>>>>>>> ed9413660b7639863292ea3ddbc8397dc70dfeb1
  double precision, parameter          :: tol = 1.0E-3
  double precision, dimension(0:n+1, 0:n+1) :: T
  double precision, dimension(n)       :: tmp1, tmp2, tmp3
  double precision                    :: error,x
  real                                :: t1,t0
  integer                             :: i, j, k, last_row
  character(len = 20)                   :: str
  integer :: omp_get_num_threads, omp_get_thread_num

  ! Set boundary conditions and initial values for the unknowns
  T = 0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0


  ! Solve the linear system of equations using the Jacobi method
  call cpu_time(t0)

  do k = 1, maxiter
<<<<<<< HEAD

!     error = 0.0D0
=======
     error = 0.0D0
>>>>>>> ed9413660b7639863292ea3ddbc8397dc70dfeb1

    !$omp parallel private(j, tmp3, last_row) shared(T) reduction (max:error)
    last_row = (omp_get_thread_num() + 1) * (n / omp_get_num_threads())

    tmp1 = T(1:n, max(0, omp_get_thread_num() * (n / omp_get_num_threads()) - 1))
    tmp3 = T(1:n, last_row)

     !$omp do
     do j = 1, n
        tmp2 = T(1:n, j)

        if (j == last_row) then
            T(1:n, j) = ( T(0:n-1, j) + T(2:n+1, j) + tmp3 + tmp1) / 4.0D0
        else
            T(1:n, j) = ( T(0:n-1, j) + T(2:n+1, j) + T(1:n, j+1) + tmp1) / 4.0D0
        end if
        !$omp critical
        error = max(error, maxval(abs(tmp2 - T(1:n, j))))
        !$omp end critical
        tmp1 = tmp2

     end do
     !$omp end do
     !$omp end parallel

     if (error < tol) then
        exit
     end if

  end do

  call cpu_time(t1)

  write(unit=*, fmt=*) 'Time:', t1 - t0, 'Number of Iterations:', k
  write(unit=*, fmt=*) 'Temperature of element T(1,1)  = ', T(1, 1)

  ! Uncomment the next part if you want to write the whole solution
  ! to a file. Useful for plotting.

  !open(unit=7,action='write',file='result.dat',status='unknown')
  !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
  !do i=0,n+1
  !   write (unit=7,fmt=str) T(i,0:n+1)
  !end do
  !close(unit=7)

end program laplsolv
