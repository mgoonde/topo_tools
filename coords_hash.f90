program coords_hash

  !! program generates the connectivity matrix, canonical labeling, and hash
  !! from a configuration given in input file.
  !!
  use f90nautyinterf
  use graph_module
  implicit none


  integer :: ntyp, i, j, k, nat
  real, allocatable :: color_cut(:,:)
  real :: dist
  integer, allocatable :: typ(:)
  real, allocatable :: coords(:,:)
  integer, allocatable :: connect(:,:)
  integer :: hash, h1, h2, h3
  integer, allocatable :: lab(:), color(:)



  !! open the cutoff file
  open( unit = 444, file= 'color_cutoff.dat', status = 'old')

  !! how many types of atoms are there
  read(444,*) ntyp

  allocate( color_cut( 1:ntyp, 1:ntyp) )



  !! read the color cutoffs for each combination of types
  do i = 1, ntyp*(ntyp+1)/2
     read(444,*) j, k, dist
     color_cut(j,k) = dist
     color_cut(k,j) = dist
  end do




  !! read configuration in xyz format:
  !
  ! read number of atoms
  read(*,*) nat
  allocate( typ(1:nat) )
  allocate( coords(1:3, 1:nat) )
  ! read empty line
  read(*,*)
  do i = 1, nat
     !  type_of_atom    x_coord      y_coord      z_coord
     read(*,*) typ(i), coords(1,i), coords(2,i), coords(3,i)
  end do




  !! ------ generate the connectivity matrix -----
  allocate( connect(1:nat, 1:nat) )
  call generate_connect( nat, coords, typ, ntyp, color_cut, connect )

  write(*,*) 'connectivity matrix:'
  write(*,*)
  do i = 1, nat
     write(*,'(*(i2))') connect(:,i)
  end do

  write(*,*)
  !!------------------------------------------------




  !! ----- generate hash of the connectivity ------

  ! allocate the array of labels
  allocate( lab(1:nat) )

  ! allocate array of colors
  allocate( color(1:nat) )


  ! prepare arrays for nauty
  call prepare_nauty( nat, typ, lab, color )
  write(*,*) 'nauty-format color array:'
  write(*,'(*(i3))') color(:)
  write(*,*)

  ! call nauty
  call c_ffnauty( nat, connect, lab, color, typ, h1, h2, h3 )

  write(*,*) 'canonical labeling:'
  write(*,'(*(i3))') lab(:) + 1

  write(*,*)
  ! final hash is a mix of three output hashes from nauty
  hash = modulo (modulo (h1,104729)+ modulo(h2, 15485863)+ modulo(h3, 882377) - 1, 1299709)+1

  write(*,*) 'hash is:', hash
  !!--------------------------------------------------





  deallocate( color_cut )
  deallocate( typ, coords )
  deallocate( connect )




end program coords_hash
