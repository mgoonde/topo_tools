program site_hash

  !! Go over all sites in a configuration, generate small local environment
  !! with a spherical radius cutoff, generate connectivity and hash of
  !! each such small local environment,
  use f90nautyinterf
  use graph_module

  integer :: ntyp
  real, allocatable :: color_cut(:,:)
  real :: dist
  !
  integer :: nat
  integer, allocatable :: typ(:)
  real, allocatable :: coords(:,:)
  !
  real :: Rcut
  integer :: isite, k
  integer :: nat_loc
  integer, allocatable :: typ_loc(:)
  real, allocatable :: coords_loc(:,:)
  integer, allocatable :: connect_loc(:,:)
  integer, allocatable :: vec(:)
  real, dimension(3) ::r_i, r_c
  integer :: hash



  !! open the Rcut file and read Rcut
  open(unit = 333, file = 'rcut.dat', status = 'old')
  read(333,*) Rcut




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





  !! allocate vector for keeping track of indexes
  allocate( vec(1:nat) )
  !! go over each site in the system
  do isite = 1, nat

     !! empty the memory of atom indexes
     vec(:) = 0

     !! calculate the distance from this site to all other sites
     do i = 1, nat

        !! skip the same site
        if( i .eq. isite) cycle

        r_i(:) = coords(:, i) - coords(:, isite)
        dist = sqrt( dot_product( r_i, r_i ) )

        !! if distance is less than Rcut, remember this atom index
        if( dist .le. Rcut ) vec(i) = 1

     end do

     !! now, 'vec' contains atom indexes that are within Rcut sphere of isite
     !! except fot isite. Total number of atoms in the local environment is
     !! then: sum(vec) + 1
     nat_loc = sum( vec ) + 1
     if( nat_loc .lt. 2 ) then
        write(*,*) 'WARNING: local environment contains only central atom'
        write(*,*) 'check Rcut, and color_cut parameters!'
     endif


     !! allocate memory for local environment
     allocate( typ_loc(1:nat_loc) )
     allocate( coords_loc(1:3, 1:nat_loc))

     !! on the first spot in the array of local environment, put isite
     !! and center the local configuration around (0.0, 0.0, 0.0)
     r_c(:) = coords(:,isite)
     typ_loc(1) = typ( isite )
     coords_loc(:,1) = coords(:, isite ) - r_c(:)

     !! fill up the rest of local arrays with info from 'vec'
     k = 2
     do i = 1, nat

        typ_loc(k) = typ(i) * vec(i)
        coords_loc(:,k) = ( coords(:,i) - r_c(:) ) * vec(i)

        !! increase index if atom is present
        k = k + vec(i)

        !! exit this loop when k > nat_loc: all atoms have been filled
        if( k .gt. nat_loc ) exit

     end do

     write(*,*) '> isite:',isite
     write(*,*) '    local environment'
     write(*,*) nat_loc
     write(*,*)
     do i = 1, nat_loc
        write(*,*) typ_loc(i), coords_loc(:,i)
     end do

     !! now we know the local environment around isite, get connect
     allocate( connect_loc(1:nat_loc, 1:nat_loc) )
     call generate_connect( nat_loc, coords_loc, typ_loc, ntyp, color_cut, connect_loc )

     write(*,*) '    local connect'
     do i = 1, nat_loc
        write(*,'(10x,*(i2))') connect_loc(:,i)
     end do

     !! get hash
     call prepare_hash( nat_loc, typ_loc, coords_loc, ntyp, color_cut, hash )

     write(*,*) '    local hash:',hash


     write(*,*)

     !! deallocate from memory the local environment'
     deallocate( typ_loc, coords_loc, connect_loc )

  end do









end program site_hash
