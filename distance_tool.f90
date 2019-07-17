program distance_tool

  use organize_module, only: point_set_registry
  implicit none

  integer :: nat_1, nat_2
  integer, allocatable :: typ_1(:), typ_2(:)
  real, allocatable :: coords_1(:,:), coords_2(:,:)
  integer :: i
  integer, allocatable :: found(:)
  real, allocatable :: dists(:)

  integer :: n

  !! read first structure
  read(*,*) nat_1
  allocate( typ_1(1:nat_1) )
  allocate( coords_1(1:3, 1:nat_1) )
  read(*,*)
  do i = 1, nat_1
     read(*,*) typ_1(i), coords_1(1,i), coords_1(2,i), coords_1(3,i)
  end do

  !! read second structure
  read(*,*) nat_2
  allocate( typ_2(1:nat_2) )
  allocate( coords_2(1:3, 1:nat_2) )
  read(*,*)
  do i = 1, nat_2
     read(*,*) typ_2(i), coords_2(1,i), coords_2(2,i), coords_2(3,i)
  end do

  n = max(nat_1, nat_2)
  allocate( found(1:n) )
  allocate( dists(1:n) )
  found(:) = 0
  dists(:) = 0.0

  call point_set_registry( nat_1, typ_1, coords_1, &
                           nat_2, typ_2, coords_2, &
                           found, dists)

  write(*,*) 'set 1 index, set 2 index, distance'
  do i = 1, n
     write(*,*) i, found(i), dists(i)
  end do

  write(*,*)
  write(*,*) 'norm of dists vector:', sqrt( dot_product( dists, dists) )
  write(*,*) 'maximum value of distance:', maxval( dists, 1 )

  deallocate( typ_1, coords_1 )
  deallocate( typ_2, coords_2 )
  deallocate( found, dists )


end program distance_tool
