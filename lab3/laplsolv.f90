program laplsolv
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem
! on a square using the Jacobi method.
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
  use omp_lib
  integer, parameter                  :: n = 1000, maxiter = 1000
  double precision, parameter          :: tol = 1.0E-3
  double precision, dimension(0:n+1, 0:n+1) :: T
  double precision, dimension(n)       :: tmp1, tmp2, tmp3
  double precision                    :: error,x,t0,t1
  integer                             :: i, j, k, kEnd
  character(len = 20)                   :: str
  integer :: iterations, exiting, last_row, first_row, interval

  ! Set boundary conditions and initial values for the unknowns
  T = 0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0


  ! Solve the linear system of equations using the Jacobi method
  t0 = omp_get_wtime() !call cpu_time(t0) !omp_time
  !$omp parallel private(j, k, tmp1, tmp2, tmp3, first_row, last_row, interval) shared(T, error)
  interval = n / omp_get_num_threads()
  last_row = (omp_get_thread_num() + 1) * interval
  first_row = max(0, omp_get_thread_num() * interval)

  do k = 1, maxiter

    error = 0.0D0

    tmp1 = T(1:n, first_row)
    tmp3 = T(1:n, last_row+1)

    !$omp barrier
    !$omp do reduction(max : error)
    do j = 1, n
      tmp2 = T(1:n, j)

      if (j == last_row) then
        T(1:n, j) = ( T(0:n-1, j) + T(2:n+1, j) + tmp3 + tmp1) / 4.0D0
      else
        T(1:n, j) = ( T(0:n-1, j) + T(2:n+1, j) + T(1:n, j+1) + tmp1) / 4.0D0
      end if

      error = max(error, maxval(abs(tmp2 - T(1:n, j))))
      tmp1 = tmp2
    end do
    !$omp end do

    if (error < tol) then
      if (omp_get_thread_num() == 0) then
        kEnd = k
      end if
      exit
    end if
    !$omp barrier
  end do
  !$omp end parallel
  
  t1 = omp_get_wtime() !call cpu_time(t1)

  write(unit=*, fmt=*) 'Time:', t1 - t0, 'Number of Iterations:', kEnd
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
