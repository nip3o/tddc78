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

  integer :: nthreads, rank, start, end, interval, omp_get_num_threads, omp_get_thread_num

  ! Set boundary conditions and initial values for the unknowns
  T = 0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0

  ! Solve the linear system of equations using the Jacobi method
  call cpu_time(t0)

  !$omp parallel private(nthreads, rank, start, end, interval, tmp1, tmp2, tmp3, error, j, k)


  rank = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  interval = maxiter / nthreads
  start = rank * interval 
  end = (rank + 1) * interval
  
  if (rank .eq. 0) then
    start = 1
  end if

  print *, 'start: ', start, ', end: ', end, ', from: ', rank

  ! Solve the linear system of equations using the Jacobi method
  call cpu_time(t0)

  do k = start, end

     tmp1 = T(1:n, 0)
     error = 0.0D0

!New(i,j) = f ( Old(i-1,j), Old(i,j-1), Old(i,j+1), Old(i+1,j))
!while (...) do
!for all tasks i,j do in parallel
  !send Old(i,j) to each neighbor
  !receive Old(i-1,j), Old(i,j-1), Old(i,j+1), Old(i+1,j) from neighbors 
  !compute New(i,j) f ( Old(i-1,j), Old(i,j-1), Old(i,j+1), Old(i+1,j) )
!od
!(then copy back Old(i,j) New(i,j) i, j to prepare for next iteration) 
!od

     do j = 1, n
        tmp2 = T(1:n, j)
        tmp3 = ( T(0:n-1, j) + T(2:n+1, j) + T(1:n, j+1) + tmp1) / 4.0D0

        !omp critical
          T(1:n, j) = tmp3
        !omp end critcal

        error = max(error, maxval(abs(tmp2 - T(1:n, j))))
        tmp1 = tmp2
     end do

     if (error < tol) then
        exit
     end if

  end do

  !$omp end parallel

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
