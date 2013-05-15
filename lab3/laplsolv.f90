program laplsolv
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem
! on a square using the Jacobi method.
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
  integer, parameter                  :: n = 100, maxiter = 100
  double precision, parameter          :: tol = 1.0E-3
  double precision, dimension(0:n+1, 0:n+1) :: T
  double precision, dimension(n)       :: tmp1, tmp2, tmp3
  double precision                    :: error,x
  real                                :: t1,t0
  integer                             :: i, j, k
  character(len = 20)                   :: str
  integer omp_get_num_threads, omp_get_thread_num, intval

  ! Set boundary conditions and initial values for the unknowns
  T = 0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0

  ! Solve the linear system of equations using the Jacobi method
  call cpu_time(t0)

  do k = 1, maxiter

    error = 0.0D0
    
    !$omp parallel private(tmp1, tmp2, tmp3, intval)
    !tmp1 = T(1:n, 0)
    intval = n / omp_get_num_threads()
    tmp1 = T(1:n, omp_get_thread_num() * intval)
    !$omp barrier
    do j = omp_get_thread_num() * intval, (omp_get_thread_num()+1) * intval
      print *, j
      tmp2 = T(1:n, j)
      tmp3 = T(1:n, j+1)
      T(1:n, j) = ( T(0:n-1, j) + T(2:n+1, j) + tmp3 + tmp1) / 4.0D0

      !$omp critical
        error = max(error, maxval(abs(tmp2 - T(1:n, j))))
      !$omp end critical

      tmp1 = tmp2
      !$omp barrier
    end do
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
