module graph_module

  implicit none
  public

contains

  subroutine generate_connect(nat,coords,types,n_col,color_cutoff,connect,lat)
    !! generate just the connectivity matrix from coords
    !! @NOTE: This does not do PBC by default, have to input lattice.
    !!
    !! nat (in)           : number of atoms
    !! coords (in)        : coordinates
    !! types (in)         : atomic types
    !! n_col (in)         : number of types
    !! color_cutoff (in)  : cutoff for connections
    !! connect (out)      : connectivity matrix
    !! lat (in, optional) : lattice, if we want connect in PBC
    use basic_module, only: cart_to_crist, &
                            periodic, &
                            crist_to_cart
    implicit none
    integer, intent(in) :: nat
    real, dimension(3,nat), intent(in) :: coords
    integer, dimension(nat), intent(in) :: types
    integer, intent(in) :: n_col
    real, dimension(n_col,n_col), intent(in) :: color_cutoff
    integer, dimension(nat,nat), intent(out) :: connect
    real, dimension(3,3), intent(in), optional :: lat

    integer :: i, j
    real :: dij
    real, dimension(3) :: rij

    connect(:,:) = 0

    do i = 1, nat
       do j = i+1, nat
          dij = 0.0
          rij(:) = coords(:,j) - coords(:,i)

          !! if lat present, do PBX
          if( present(lat) ) then
             call cart_to_crist( rij, lat )
             call periodic( rij )
             call crist_to_cart( rij, lat )
          endif


          dij = sqrt( dot_product( rij, rij) )
          connect(i,j) = NINT( 0.5*erfc(dij-color_cutoff( types(i),types(j) )))
          connect(j,i) = connect(i,j)
       end do
    end do
  end subroutine generate_connect


  subroutine prepare_nauty(nat,typ,lab,color)
    !! prepare the order vector and color vector for nauty.
    !! lab(:) is the order vector for C input
    !! color(:) is the (0,1)-format vector specific of nauty
    implicit none
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ
    integer, dimension(nat), intent(out) :: color
    integer, dimension(nat), intent(out) :: lab

    integer :: i, j, dum

    !! assume starting order
    do i = 1, nat
       lab(i) = i
    end do
    !! reorder in ascending values of typ
    do j = 1, nat
       do i = 1, nat-1
          if( typ( lab(i) ) > typ( lab(i+1) ) ) then
             dum = lab(i)
             lab(i) = lab(i+1)
             lab(i+1) = dum
          endif
       end do
    end do
    !! now generate color vector
    do i = 1, nat
       color(i) = typ( lab(i) )
    end do
    do i = 1, nat-1
       if ( color(i) .ne. color(i+1) ) color(i) = 0
    end do
    color(nat) = 0
    !! is for C, so -1
    do i = 1, nat
       lab(i) = lab(i) - 1
    end do
  end subroutine prepare_nauty


  subroutine identify_cluster(nat,connect,isite,res)
  !! identify a cluster around isite, based on connectivity matrix
  !! matrix-vector multiplication of connect with a vector that is 1 at indices of atoms we wish
  !! to get connections of, gives all connections to those atoms. Like this can start from atom
  !! 'isite', and grow the vector with indices step by step to include further connections, while
  !! also keeping the old ones.
  !! Once the resulting vector is the same as input vector (no new connections have been found), the
  !! algorithm can stop, we have a cluster.
   implicit none
   integer, intent(in) :: nat
   integer, dimension(nat,nat), intent(in) :: connect
   integer, intent(in) :: isite
   integer, dimension(nat), intent(out) :: res !! gives a vector with 1 on index which is in cluster

   integer, dimension(nat) :: vector, res_old
   integer :: i, j

   res(:) = 0
   res_old(:) = 0
   vector(:) = 0
   vector(isite) = 1
   DO j = 1, nat
     !! do matrix-vector product where vector contains all previously found indices
     do i = 1, nat
       res(i) = dot_product( connect(i,:), vector )
       if( res(i) .ne. 0) res(i) = 1   !! keeping values at 0, only interested if there is connection
     end do
!write(*,*) 'res'
!write(*,'(40I3)') (res(i),i=1,nat)
     !! add newly found connections to the vector
     vector = vector + res

     !! if there are no new connections, exit
     if( all(res .eq. res_old )) exit
!write(*,'(64I2)') res(:) - res_old(:)

     res_old = res
   END DO
!write(*,'(40I3)') (res(i),i=1,nat)
  end subroutine identify_cluster


  subroutine prepare_hash( nat, typ, coords, ntyp, color_cut, hash )

    use f90nautyinterf

    implicit none
    integer,                    intent(in) :: nat
    integer, dimension(nat),    intent(in) :: typ
    real, dimension(3,nat),     intent(in) :: coords
    integer,                    intent(in) :: ntyp
    real, dimension(ntyp,ntyp), intent(in) :: color_cut
    integer,                    intent(out) :: hash

    integer, allocatable :: connect(:,:)
    integer, allocatable :: lab(:), color(:)
    integer :: i, h1, h2, h3

!    write(*,*) 'input conf:'
!    write(*,*) nat
!    write(*,*)
!    do i = 1, nat
!       write(*,*) typ(i), coords(i,:)
!    end do



    !! allocate and generate connectivity
    allocate( connect(1:nat,1:nat))
    call generate_connect( nat, coords, typ, ntyp, color_cut, connect )

!    do i = 1, nat
!       write(*,'(64i2)') connect(i,:)
!    end do

    !! hash the connect
    allocate( lab(1:nat) )
    allocate( color(1:nat) )
    call prepare_nauty( nat, typ, lab, color)
    call c_ffnauty( nat, connect, lab, color, typ, h1, h2, h3 )

    hash = modulo (modulo (h1,104729)+ modulo(h2, 15485863)+ modulo(h3, 882377) - 1, 1299709)+1


    !! output final info

!    write(*,*) 'found hashes:',h1,h2,h3, hash



    deallocate( lab, color )
    deallocate( connect )

  end subroutine prepare_hash


   subroutine generate_neighbor_connectivity(nat, c_in, neigh, c_out)
     !! generate a connectivity matrix which is directed such that each vertex has a directed
     !! connection FROM a vertex that is up to 'neigh' connections away, including self.
     !! If the vertex 4 is found from vertex 6, the connectivity element c_6,4 = 1.
     !! Like this, if you multiply thuis connectivity with a vector with 1 at some index i, you
     !! get back all indices up to 'neigh' away from i.
     implicit none
     integer, intent(in) :: nat
     integer, dimension(nat,nat), intent(in) :: c_in
     integer, intent(in) :: neigh
     integer, dimension(nat,nat), intent(out) :: c_out

     integer :: i, j, site, n_count
     integer, dimension(nat) :: vec, res, res_old
     logical, dimension(nat) :: ci, vi

     c_out(:,:) = c_in(:,:)

     do site = 1, nat
       res(:) = 0
       res_old(:) = 0
       vec(:) = 0
       vec(site) = 1
       n_count = 1

       !! for each site in the graph go 'neigh'-deep and find all connected vertices
       do while( n_count .le. neigh )

         !! matrix-vector 'multiplication' with logicals maybe faster for large matrices
         do i = 1, nat
           ci = c_in(i,:)
           vi = vec(:)
           if( any( ci(:) .and. vi(:) ) ) res(i) = 1
         end do

         !! add the newly found directed connections
         c_out(:,site) = c_out(:,site) + res(:)

         !! if nothing new found, exit
         if( all( res(:) .eq. res_old(:) )) exit

         !! increment with keeping all previous vertices
         vec(:) = vec(:) + res(:)
         res_old(:) = res(:)
         n_count = n_count + 1
       end do
     end do

     do i = 1, nat
       do j = 1, nat
         if( c_out(i,j) .ne. 0 ) c_out(i,j) = 1
       end do
     end do

   end subroutine generate_neighbor_connectivity


   subroutine selective_local_connect(n_in, c_in, list, n_out, c_out)
     !! Generate a smaller connect from the global. Pick out only the vertices given in
     !! the input vector 'list' which is a (0,1)-vector.
     implicit none
     integer, intent(in) :: n_in
     integer, dimension(n_in,n_in), intent(in) :: c_in
     integer, dimension(n_in), intent(in) :: list
     integer, intent(in) :: n_out
     integer, dimension(n_out,n_out), intent(out) :: c_out

     integer :: i, j, m, n

     m = 1
     n = 1
     do i = 1, n_in
        n = 1
        do j = 1, n_in
           c_out(m,n) = c_in(i,j)*list(j)
           n = n + list(j)
           if( n .eq. n_out + 1 ) exit
        end do
        m = m + list(i)
        if( m .eq. n_out + 1 ) exit
     end do

   end subroutine selective_local_connect


   subroutine set_new_connections( n, idx_in, list_new, connect )
     !!
     !! using outer product to generate conenctivity matrix parts
     !! from input site which has its connections in 'list_new'
     !!
     !!  ( i )     ) ( output )        ( o )     ) ( input  )
     !!  ( n )     ) (        )        ( u )     ) (        )
     !!  ( p ) --> ) (    ^   )   +    ( t ) --> ) (    ^   )    =    M_1 + M_2
     !!  ( u )     ) (    |   )        ( p )     ) (    |   )
     !!  ( t )     ) (    |   )        ( t )     ) (    |   )
     !!
     !!============================
     !! Here connect is assumed to be proper size already!!
     implicit none
     integer,                 intent(in) :: n
     integer,                 intent(in) :: idx_in
     integer, dimension(n),   intent(in) :: list_new
     integer, dimension(n,n), intent(inout) :: connect
     !!
     integer, allocatable :: list_in(:)
     integer :: i, j
     integer, dimension(n,n) :: M_1
!     integer, dimension(n,n) :: M_2

     allocate( list_in(1:n) )
     list_in(:) = 0
     list_in(idx_in) = 1

     !! "outer product"
     !! must compute both parts, since don't know on which triangular
     !! side of connectivity the connections are;
     !! or just compute one, and then force symmetric connect later
     M_1(:,:) = 0
!     M_2(:,:) = 0
     do i = 1, n
        M_1(:,i) = list_in(:) * list_new(i)
!        M_2(:,i) = list_new(:) * list_in(i)
     end do

     connect(:,:) = connect(:,:) + M_1 !+ M_2

      do i = 1, n
         do j = 1, n
            if( connect(i,j) .ne. 0) then
               connect(i,j) = 1
               connect(j,i) = 1
            endif
         end do
      end do

      deallocate( list_in )

   end subroutine set_new_connections


   subroutine find_new_connections( idx_site, nat, typ, coords, ntyp, color_cut, list_out, lat )
     !! find connections from site idx_site, output them in list_out
     use basic_module, only: cart_to_crist, &
                             periodic, &
                             crist_to_cart
     implicit none
     integer, intent(in) :: idx_site
     integer, intent(in) :: nat
     integer, dimension(nat), intent(in) :: typ
     real, dimension(3,nat), intent(in) :: coords
     integer, intent(in) :: ntyp
     real, dimension(ntyp, ntyp), intent(in) :: color_cut
     integer, dimension(nat), intent(out) :: list_out
     real ,dimension(3,3), optional :: lat

     integer :: i
     real :: dist
     real, dimension(3) :: r_isite, r_i

     list_out(:) = 0

     !! the vector of site
     r_isite(:) = coords(:,idx_site)

     do i = 1, nat
        !! skip himself
        if( i .eq. idx_site ) cycle
        r_i = coords(:, i) - r_isite(:)
        if( present(lat) ) then
           call cart_to_crist( r_i, lat )
           call periodic( r_i )
           call crist_to_cart( r_i, lat )
        endif
        dist = sqrt( dot_product( r_i, r_i) )
        if( dist .le. color_cut( typ( idx_site), typ(i) ) ) then
           list_out(i) = 1
        endif
     end do

   end subroutine find_new_connections


   subroutine clean_connect_matrix( nat, list_in, connect )
     !! remove connections given by list_in from the connect matrix
     implicit none
     integer, intent(in) :: nat
     integer, dimension(nat), intent(in) :: list_in
     integer, dimension(nat,nat), intent(inout) :: connect

     integer, dimension(nat) :: list
     integer :: idx

     list(:) = list_in(:)
     idx = maxloc(list, 1)
     do while( list( idx ) .ne. 0)
        connect(idx, :) = 0
        connect(:, idx) = 0
        list(idx) = 0
        idx = maxloc(list, 1)
     end do

   end subroutine clean_connect_matrix


end module graph_module
