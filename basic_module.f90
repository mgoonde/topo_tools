module basic_module

  implicit none
  public

contains


  function cross(a,b) result(res)
    implicit none
    real, dimension(3) :: a, b
    real, dimension(3) :: res

    res(1) = a(2)*b(3) - a(3)*b(2)
    res(2) = a(3)*b(1) - a(1)*b(3)
    res(3) = a(1)*b(2) - a(2)*b(1)
  end function cross


  real function norm(r)
    implicit none
    real, dimension(:) :: r

    norm = sqrt( dot_product( r,r ) )
  end function norm


  function lower(s1)  result (s2)
    character(*)       :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: duc = ichar('A') - ichar('a')
    integer            :: i

    do i = 1,len(s1)
       ch = s1(i:i)
       if (ch >= 'A'.and. ch <= 'Z') ch = char(ichar(ch)-duc)
       s2(i:i) = ch
    end do
  end function lower


  subroutine periodic(c)
    !--------------------------------
    ! general periodic boundary condition, for 1, 2 or 3 dimensional vector input.
    !--------------------------------
    implicit none
    real, dimension(:),intent(inout) :: c
    integer :: n
    n=size(c)
    if(c(1) < (-0.5)) c(1) = c(1) + 1.0
    if(c(1) >= 0.5) c(1) = c(1) - 1.0
    if ( n .gt. 1 ) then
       if(c(2) < (-0.5)) c(2) = c(2) + 1.0
       if(c(2) >= 0.5) c(2) = c(2) - 1.0
    endif
    if (n .gt. 2 ) then
       if(c(3) < (-0.5)) c(3) = c(3) + 1.0
       if(c(3) >= 0.5) c(3) = c(3) - 1.0
    endif
  end subroutine periodic


  subroutine cart_to_crist(xpp,ct)
  !----------------------------
  ! cartesian to crystallographic coordinates transform, in 3-dimension
  ! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
  ! --------
  ! xpp(3)      ==> input vector of position in cartesian
  ! ct(3,3)     ==> conversion matrix, vectors of the Bravais lattice
  !----------------------------
  ! bt(3,3) ==> inverse matrix of ct, used locally
  ! xc(3)   ==> copy of xpp, used locally
  ! detct   ==> determinant of ct, used locally
  !
   implicit none
   real, dimension(3), intent(inout) :: xpp
   real, dimension(3,3), intent(in) :: ct

   real,dimension(3) :: xc
   real :: detct
   real, dimension(3,3) :: bt

       bt(:,:)=0.0
       xc(:) = 0.0

  ! -----------------------------------------------
  !  inverse matrix of ct(:,:)
  !------------------------------------------------
       detct=ct(1,1)*ct(2,2)*ct(3,3)+&
             ct(1,2)*ct(2,3)*ct(3,1)+&
             ct(2,1)*ct(3,2)*ct(1,3)&
            -ct(1,3)*ct(2,2)*ct(3,1)&
            -ct(3,2)*ct(2,3)*ct(1,1)&
            -ct(1,2)*ct(2,1)*ct(3,3)

       bt(1,1)= ct(2,2)*ct(3,3)-ct(2,3)*ct(3,2)
       bt(1,2)=-(ct(1,2)*ct(3,3)-ct(1,3)*ct(3,2))
       bt(1,3)= ct(1,2)*ct(2,3)-ct(1,3)*ct(2,2)
       bt(2,1)=-(ct(2,1)*ct(3,3)-ct(2,3)*ct(3,1))
       bt(2,2)= ct(1,1)*ct(3,3)-ct(3,1)*ct(1,3)
       bt(2,3)=-(ct(1,1)*ct(2,3)-ct(1,3)*ct(2,1))
       bt(3,1)= ct(2,1)*ct(3,2)-ct(2,2)*ct(3,1)
       bt(3,2)=-(ct(1,1)*ct(3,2)-ct(1,2)*ct(3,1))
       bt(3,3)= ct(1,1)*ct(2,2)-ct(2,1)*ct(1,2)
  !------------------------------------------------

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))/detct
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))/detct
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))/detct

       xpp(:) = xc(:)

       return
  end subroutine cart_to_crist


  subroutine crist_to_cart(xpp,bt)
  !--------------------------------
  ! crystallographic to cartesian transformation in 3-dimensions
  ! R_cart = B * v_crist; where B is the matrix formed by cell vectors horizontally
  ! -----------
  ! xpp(3)    ==> input vector in crystallographic, output vector in cartesian
  ! bt(3,3)   ==> input conversion matrix, vectors of the Bravais lattice
  !-----
  ! xc(3)   ==> local vector
  !
   implicit none
   real, dimension(3), intent(inout) :: xpp
   real, dimension(3,3), intent(in) :: bt
   real, dimension(3) :: xc

       xc(:) = 0.0

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))

       xpp(:)=0.0

       xpp(:) = xc(:)

       return
  end subroutine crist_to_cart


  subroutine find_loc(n,array,val,direction,loc)
    ! Find location of value 'val' in array 'array', location given in 'loc'.
    implicit none
    integer, intent(in) :: n !! dimension of array
    integer, dimension(n), intent(in) :: array !! array to look in
    integer, intent(in) :: val !! value that is searched
    integer, intent(in) :: direction !! direction to search (1:frwd, -1:bckwd)
    integer, intent(out) :: loc !! location found

    integer :: i

    !! at failure to find val
    loc = 0

    do i = 1, n
       if( array( i ) == val ) then
          loc = i
          !! found first occurence of val
          if( direction == 1 ) exit
          !! continue searching for last occurence of val
       endif
    end do

  end subroutine find_loc


  subroutine read_line(fd, line, end_of_file)
    !--------------------
    ! read a line, makes possible to use # for comment lines, skips empty lines,
    !  is pretty much a copy from QE.
    !---------
    ! fd ==> file descriptor
    ! line ==> what it reads
    implicit none
    integer, intent(in) :: fd
    integer             :: ios
    character(len=100), intent(out) :: line
    logical, optional, intent(out) :: end_of_file
    logical :: tend

    tend = .false.
101 read(fd,fmt='(A256)',END=111, iostat=ios) line
    !print*, line 
    if (ios /= 0)then
       print*, " Reading Problem..."; stop; endif
    if(line == ' ' .or. line(1:1) == '#') go to 101
    go to 105
111 continue
    tend = .true.
    go to 105
105 continue

    if( present(end_of_file)) then
       end_of_file = tend
    endif
  end subroutine read_line


  subroutine number_boxes( lat, rcut, boxes )
    !! Stefano's function to check how many boxes you need for
    !! the given cutoff in order to have all atoms inside.
    !! Provided by Ruggero.
    implicit none
    real, dimension(3,3), intent(in) :: lat
    real, intent(in) :: rcut
    integer, dimension(3), intent(out) :: boxes

    real :: detlat
    real, dimension(3,3) :: invlat
    integer :: i

    detlat = lat(1,1)*lat(2,2)*lat(3,3)+&
        lat(1,2)*lat(2,3)*lat(3,1)+&
        lat(2,1)*lat(3,2)*lat(1,3)&
        -lat(1,3)*lat(2,2)*lat(3,1)&
        -lat(3,2)*lat(2,3)*lat(1,1)&
        -lat(1,2)*lat(2,1)*lat(3,3)

    invlat(1,1) = lat(2,2)*lat(3,3) - lat(2,3)*lat(3,2)
    invlat(1,2) = -(lat(1,2)*lat(3,3) - lat(1,3)*lat(3,2))
    invlat(1,3) = lat(1,2)*lat(2,3) - lat(1,3)*lat(2,2)
    invlat(2,1) = -(lat(2,1)*lat(3,3) - lat(2,3)*lat(3,1))
    invlat(2,2) = lat(1,1)*lat(3,3) - lat(3,1)*lat(1,3)
    invlat(2,3) = -(lat(1,1)*lat(2,3) - lat(1,3)*lat(2,1))
    invlat(3,1) = lat(2,1)*lat(3,2) - lat(2,2)*lat(3,1)
    invlat(3,2) = -(lat(1,1)*lat(3,2) - lat(1,2)*lat(3,1))
    invlat(3,3) = lat(1,1)*lat(2,2) - lat(2,1)*lat(1,2)

    invlat(:,:) = invlat(:,:)/detlat

    do i = 1, 3
!       boxes(i) = 2*nint( rcut * norm(invlat(i,:))) !! from - to + direction
       boxes(i) = nint( rcut * norm(invlat(i,:)))  !! only + direction
    end do

    return
  end subroutine number_boxes


  subroutine set_random_seed()
    ! ----- setting a random seed, based on current time -----
    integer :: i_seed
    integer, dimension(:), allocatable :: a_seed
    integer, dimension(1:8) :: dt_seed

    ! ----- Set up random seed portably -----
    call random_seed(size=i_seed)
    allocate(a_seed(1:i_seed))
    call random_seed(get=a_seed)
    call date_and_time(values=dt_seed)
    a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    call random_seed(put=a_seed)
    deallocate(a_seed)
  end subroutine set_random_seed




end module basic_module
