program matmulSkel
!---------------------------------------------------------
!   This is a skeleton file for matrix - matrix 
!   multiplication on a 2D processor grid between  
!   matrices of the same size using the MPI interface!
!
!   Berkant Savas March 2006 
!   E-mail besav@mai.liu.se if you have any questions.
!---------------------------------------------------------

include 'mpif.h'
integer,parameter:: DIMM     = 1000     ! Dimension of matrix

integer,parameter:: cartDim  = 2        ! Dimension of map
logical,parameter:: reorder  = .true.   ! Reorder nodes 
logical,parameter:: WRAP     = .false.  ! Wrap around corners

integer,parameter:: leftRight= 0
integer,parameter:: upDown   = 1
integer,parameter:: up       = 1
integer,parameter:: down     = 2
integer,parameter:: left     = 3
integer,parameter:: right    = 4
integer,parameter:: X        = 1
integer,parameter:: Y        = 2

integer:: numprocs, myId, myX, myY, ierr, cartComm
integer:: cartSize(1:2), cartDimX, cartDimY, temp, myM, myN
integer:: rowPos(1:20), colPos(1:20), myCoord(1:2), link(1:4)
logical:: period(1:2)
double precision,dimension(:,:),allocatable:: myA, myB, myC
double precision,dimension(:,:),allocatable:: workA, workB
integer:: status(MPI_STATUS_SIZE)
!---------------------------------------------------------
!--- Additional variables: 
!---------------------------------------------------------


!---------------------------------------------------------
!--- Initiate the MPI interface
!---------------------------------------------------------
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myId,ierr)

!---------------------------------------------------------
!--- Build a Cartesian grid with optimal size 
!---------------------------------------------------------
cartSize=(/0,0/)
period=(/ WRAP, WRAP /)

call MPI_DIMS_CREATE(numprocs,cartDim,cartSize,ierr)
call MPI_CART_CREATE(MPI_COMM_WORLD,cartDim,cartSize,&
     period,reorder,cartComm,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myId,ierr)
call MPI_CART_COORDS(cartComm,myId,cartDim,myCoord,ierr)
call MPI_CART_SHIFT(cartComm,leftRight,1,link(left),& 
     link(right),ierr)
call MPI_CART_SHIFT(cartComm,upDown,1,link(up),&
     link(down) ,ierr)

cartDimX=cartSize(X)  ! X and Y dimensions of the   
cartDimY=cartSize(Y)  ! Cartesian grid. 

myX=myCoord(X)        ! X and Y coordinates of the
myY=myCoord(Y)        ! processors in the grid. 

!---------------------------------------------------------
!--- Determine row and column indices of A for the local 
!--- block partitioning myA contained in each processor.
!---------------------------------------------------------
rowPos(1)=1
temp=floor(real(DIMM)/real(cartDimY))
do i=2,cartDimY
   rowPos(i)=rowPos(i-1)+temp;
end do
rowPos(cartDimY+1)=DIMM+1

colPos(1)=1
temp=floor((real(DIMM))/real(cartDimX));
do i=2,cartDimX
   colPos(i)=colPos(i-1)+temp
end do
colPos(cartDimX+1)=DIMM+1
  
!---------------------------------------------------------
!--- Determine the dimensions of the local matrix blocks.
!---------------------------------------------------------
myM = rowPos(myY+2) - rowPos(myY+1)
myN = colPos(myX+2) - colPos(myX+1)

!---------------------------------------------------------
!--- Let each processor allocate and initiate their own 
!--- block partitionings of A, B and C.  
!---------------------------------------------------------
allocate(myA(1:myM,1:myN))
allocate(myB(1:myM,1:myN))
allocate(myC(1:myM,1:myN))

myA=1.5D0
myB=1.7D0
myC=0.0D0

call MPI_BARRIER(cartComm,ierr)

!-------------------------------------------------------------------
!---------------------  Description of variables -------------------
!-------------------------------------------------------------------
!
!  The processors are arranged in a 2D grid and the links are 
!  specified in the vector link, se below. The matrices is 
!  partitioned by the same structure as the processors and 
!  distributed to corresponding processor. myA, myB and myC are 
!  the local blocks on each processor. 
!  If you want a torus instead of a grid change 
!  WRAP=.false. to WRAP=.true. in the beginning.
!  
!    DIMM        : Number of rows (and also columns) in A
!    myM         : Number of rows in my partition of A
!    myN         : Number of columns in my partition of A
!                  myM and myN depend on the partitioning of 
!                  the processors. They are usually not the same
!    myX,myY     : X and Y coordinates of the processors in the grid
!    myA,myB,myC : Local partitionings of the matrices A, B and C. 
!                  Their dimensions are (myM x myN). Note that no
!                  processor contains the whole matrices, each 
!                  one contains its own block partitioning only. 
!    cartComm    : Communicator of the grid topology
!    rowPos      : Array (1:cartDimY+1), rowPos(myX) is the row index 
!                  in the whole A which corresponds to the first row 
!                  in the local block partitioning myA. 
!    colPos      : Array (1:cartDimX+1), here we have the indices 
!                  for the columns instead.
!    link        : links to the neighborhood, use up, down, left and 
!                  right as argument 
!    link(left)  : Backward in x 
!    link(right) : Forward  in x
!    link(up)    : Backward in y
!    link(down)  : Forward  in y
!-------------------------------------------------------------------
!  Send and receive calls have the following syntax:
!    call MPI_SEND(workA,size,MPI_DOUBLE_PRECISION,link(direction),& 
!         tag,cartComm,ierr)
!    call MPI_RECV(workA,size,MPI_DOUBLE_PRECISION,link(direction),&
!         tag,cartComm,status,ierr)
!  Change size, direction and tag to the correct values.  
!-------------------------------------------------------------------

!---------------------------------------------------------
!--- Examine the different variables and ensure yourself 
!--- of their meaning by writing them out in a series of 
!--- scenarios. Below are a few examples, just uncomment.
!---------------------------------------------------------
!if (myId==0) then
!   print *,"Dimension of matrix:",DIMM
!   print *,"Number of proc:",numprocs
!   print *,"Cart dim X:",cartDimX
!   print *,"Cart dim Y:",cartDimY
!end if
!print *,"myId and coord in grid", myId, myX, myY
!print *,"myId, myM and myN:", myID, myM, myN
!print *,"myId and link", myId, link 

!---------------------------------------------------------
!--- Write your matrix multiplication routine here.
!--- The answer should be returned in the matrix myC.
!--- Use workA and workB as work arrays.
!---------------------------------------------------------
allocate(workA(1:myM,1))
allocate(workB(1,1:myN))
!=========================================================
! My code here.....


!=========================================================
deallocate(workA)
deallocate(workB)
deallocate(myA)
deallocate(myB)
deallocate(myC)
call MPI_FINALIZE(ierr);
end program matmulSkel
