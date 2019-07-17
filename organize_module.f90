module organize_module
  implicit none
  public

  contains

  real function global_function(vec, q_alpha, alpha, mu_k, sigma_k)
    !> @brief Global function at point q_alpha.
    !!
    !! @detail
    !! The value of the function:\n
    !!
    !! exp( | R | / mu_k) *
    !! exp( (q_alpha - <Q, e_alpha>)**2 / (2*sigma_k**2) ).
    !!
    !! at point q_alpha is calculated for single input vector R.
    !! Since the vector should be in orthonormal basis, the norm of R
    !! should be preserved, so no need to convert it to 'crist' vector Q.
    !! The projection <Q, e_alpha> is taken simply as value of Q along
    !! the current axis alpha -> Q( alpha ). That can be done because
    !! q_alpha is going along the axes of the basis.
    !!==========================================
    !! @param vec       input vector Q (or R)
    !! @param q_alpha   point in which the function is evaluated
    !! @param alpha     axis
    !! @param mu_k      'penetration depth' corresponding to input vector typ k
    !! @param sigma_k   'gaussian width' corresponding to input vector typ k
    use basic_module, only: norm
    implicit none
    real, dimension(3) :: vec
    real               :: q_alpha
    integer               :: alpha
    real               :: mu_k
    real               :: sigma_k

    real :: g_h, g, pi

    pi = 4.0*atan(1.0)
    !! the gaussian height
    g_h = exp( -norm( vec(:) ) / mu_k ) * 1.0 / (sigma_k*sqrt(2.0*pi))

    !! the value of gaussian
    g = exp( -( q_alpha - vec( alpha ) )**2 / (2.0 * sigma_k**2) )

    !! multiply
    global_function = g_h * g

  end function global_function


  subroutine set_orthonorm_bas(vec1,vec2,basis,fail)
    !> @brief Set orthonormal basis from input vectors, third is cross of
    !! the first two, fail happens when input is collinear.
    use basic_module, only: norm, cross
    implicit none
    real, dimension(3), intent(in) :: vec1, vec2
    real, dimension(3,3),intent(out) :: basis
    logical, intent(out) :: fail

    real :: prod
    real :: collinearity_thr

    collinearity_thr = 0.9

    fail = .false.
    basis(:,:) = 0.0

    !! first vector, normalize
    basis(1,:) = vec1(:) / norm( vec1(:) )

    !! second vector, normalize
    basis(2,:) = vec2(:) / norm( vec2(:) )

    !! check projection
    prod = dot_product( basis(1,:), basis(2,:) )

    !! if vectors are collinear, return with fail = .true.
    if( abs(prod) .gt. collinearity_thr ) then
       fail = .true.
       return
    endif

    !! take orthogonal component of second vector, normalize
    basis(2,:) = basis(2,:) - prod*basis(1,:)
    basis(2,:) = basis(2,:) / norm( basis(2,:) )

    !! third vector is cross product of first two, normalize
    basis(3,:) = cross( basis(1,:), basis(2,:) )
    basis(3,:) = basis(3,:) / norm( basis(3,:) )

  end subroutine set_orthonorm_bas


  subroutine compute_overlap(nat1, typ1, coords1, g1_ind, &
                   nat2, typ2, coords2, g2_ind, &
                   ntyp, sigma_k, mu_k, Rcut, overlap,step,site)
    !> @brief
    !! calculate the overlap on the fly (NOT YET COMPLETELY TESTED).
    !!
    !! @detail
    !! typ1 and coords1 are the reference geometry, typ2 and coords2 are the
    !! system geometry.
    !! The overlap is normalized to the function of reference.
    !! The function f_alpha^k is evaluated for both geometries at once,\n
    !!
    !! f_alpha^k (q_alpha) = sum_i {
    !! exp( | R_i^k | / mu_k) *
    !! exp( (step - <Q_i^k, e_alpha>)**2 / (2*sigma_k**2) )
    !!                              }.
    !!
    !! alpha is the axis, k is the color, q_alpha is a variable along alpha,
    !! e_alpha is a basis vector along alpha, R_i is the cartesian vector of
    !! atom i, Q_i is the same vector written in the local basis
    !! (if using orthonormal basis they can be the same vector since distances
    !! are preserved). The overlap is computed as\n
    !!
    !! O_alpha^k = sum_{q_alpha} {
    !! min( g_alpha^k( q_alpha ), f_alpha^k( q_alpha ) )
    !!                            }.
    !!
    !! where g_alpha^k( q_alpha ) is the function evaluated for the first
    !! geometry and f_alpha^k( q_alpha ) is the function for the second
    !! geometry. At the same time the total sum of g_alpha^k is computed and
    !! this sum used to normalize the overlap.
    !!================================================
    !! arguments to pass:
    !!
    !! @param nat1       [in]  number of atoms in geometries 1
    !! @param nat2       [in]  number of atoms in geometries 2
    !! @param typ1       [in]  array of atomic types (color) of geometry 1
    !! @param typ2       [in]  array of atomic types (color) of geometry 2
    !! @param coords1    [in]  coords of geometry 1
    !! @param coords2    [in]  coords of geometry 2
    !! @param g1_ind     [in]  indices of atoms used in basis of geometry 1
    !!                         where (1) is central atom; (2), (3) are basis
    !! @param g2_ind     [in]  indices of atoms used in basis of geometry 2
    !! @param ntyp       [in]  total number of different atomic types
    !! @param sigma_k    [in]  array giving gaussian width for each color k
    !! @param mu_k       [in]  array giving 'penetration depth' for each color k
    !! @param Rcut       [in]  cutoff, start evaluating at -Rcut until +Rcut
    !! @param overlap    [out]  overlap value for each color
    !!==========================================
    !! local variables:
    !!
    !! alpha        axis
    !! step         value of the step in evaluation
    !! q_alpha      variable over which to evaluate
    !! normal       normalization
    !! f1, f2       the two functions to evaluate
    !! g            value of the gaussian at the current point
    !! g_h          height of the gaussian
    !! present1     indicator wether a color is present in geometry 1.
    !!              Used in normalization.
    !! present2     indicator wether a color is present in geometry 2.
    use basic_module, only: norm
   implicit none
   integer, intent(in)                       :: nat1, nat2
   integer, dimension(nat1), intent(in)      :: typ1
   integer, dimension(nat2), intent(in)      :: typ2
   real   , dimension(3,nat1), intent(in)    :: coords1
   real   , dimension(3,nat2), intent(in)    :: coords2
   integer, dimension(3), intent(in)         :: g1_ind, g2_ind
   integer, intent(in)                       :: ntyp
   real   , dimension(ntyp), intent(in)      :: sigma_k
   real   , dimension(ntyp), intent(in)      :: mu_k
   real   , intent(in)                       :: Rcut
   real   , dimension(ntyp+1), intent(out)     :: overlap
   real   , intent(in)                       :: step

   integer, optional :: site

   integer :: i, k, alpha
   real    :: q_alpha
   real   , dimension(ntyp) :: normal
   real, dimension(ntyp)    :: f1, f2
   real    :: g
   !real    :: g_h
   logical, dimension(ntyp) :: present1, present2
   !character(len=100) :: fmt
   real :: heights1, heights2

!   step = 1e-4

   overlap(:) = 0.0
   normal(:) = 0.0

   present1(:) = .false.
   present2(:) = .false.


   heights1 = 0.0
   heights2 = 0.0
   !! compute the sum of all gaussian heights present in graphs,
   !! needed for normalization
   do alpha = 1, 3
      do k = 1, ntyp
         do i = 1, nat1
            if ( i .eq. g1_ind(1) ) cycle
            if( typ1(i) .ne. k ) cycle
            if(( alpha .eq. 2) .and. (i .eq. g1_ind(2) )) cycle
            if( &
                 ( alpha .eq. 3) .and. &
                 ( (i .eq. g1_ind(2)) .or. (i .eq. g1_ind(3)) ) &
                 ) cycle

            heights1 = heights1 + exp( -norm(coords1(:,i))/mu_k(k))
         end do
         do i = 1, nat2
            if( i .eq. g2_ind(1) ) cycle
            if( typ2(i) .ne. k ) cycle
            if(( alpha .eq. 2) .and. (i .eq. g2_ind(2) )) cycle
            if( &
                 (alpha .eq. 3) .and. &
                 ( (i .eq. g2_ind(2)) .or. (i .eq. g2_ind(3)) ) &
                 ) cycle

            heights2 = heights2 + exp( -norm(coords2(:,i))/mu_k(k))
         end do
      end do
   end do


   !! each axis separately -- (Q:why?)(A:maybe eventually axis-resolved overlap?)
   do alpha = 1, 3

         !! starting point of evaluation
         q_alpha = -Rcut
         do while( q_alpha .le. Rcut )
           f1(:) = 0.0
           f2(:) = 0.0
           !! each color separately
           do k = 1, ntyp
              !! sum up the contributions to the function of the first geometry
              !! in this q_alpha. Start from index 1, central atom index is
              !! skipped in g_ind(1)
              do i = 1, nat1
                 if ( i .eq. g1_ind(1) ) cycle

                 !! only take atoms of the current color k
                 if( typ1(i) .ne. k ) cycle

                 if(( alpha .eq. 2) .and. (i .eq. g1_ind(2) )) cycle

                 if( &
                      ( alpha .eq. 3) .and. &
                      ( (i .eq. g1_ind(2)) .or. (i .eq. g1_ind(3)) ) &
                      ) cycle

                 g = global_function(coords1(:,i),q_alpha,alpha,mu_k(k),sigma_k(k))

                 f1(k) = f1(k) + g/heights1
              end do
              !! If there is a contribution to f1, this color is present in geometry 1
              if( f1(k) .gt. 1e-8 ) present1(k) = .true.


              do i = 1, nat2
                 if( i .eq. g2_ind(1) ) cycle

                 if( typ2(i) .ne. k ) cycle

                 if(( alpha .eq. 2) .and. (i .eq. g2_ind(2) )) cycle

                 if( &
                      (alpha .eq. 3) .and. &
                      ( (i .eq. g2_ind(2)) .or. (i .eq. g2_ind(3)) ) &
                      ) cycle

                 g = global_function(coords2(:,i),q_alpha,alpha,mu_k(k),sigma_k(k))
                 f2(k) = f2(k) + g/heights2
              end do
              !! If there is a contribution to f2, this color is present in geometry 2
              if( f2(k) .gt. 1e-8 ) present2(k) = .true.

              !! 'integrate' the overlap
              overlap(k) = overlap(k) + min(f1(k),f2(k))*step

              !! normalization relative to the first geometry
              normal(k) = normal(k) + f1(k)*step


           end do
           !! advance to next point
           q_alpha = q_alpha + step

           !! if printout, printout
           if( present(site) ) then
!              if( alpha .eq. 1) write(601,*) q_alpha,sum(f1), sum(f2)
!              if( alpha .eq. 2) write(602,*) q_alpha,sum(f1), sum(f2)
!              if( alpha .eq. 3) write(603,*) q_alpha,sum(f1), sum(f2)
              if( alpha .eq. 1) write(601,*) q_alpha,(f1(i),i=1,ntyp),sum(f1),(f2(i),i=1,ntyp),sum(f2)
              if( alpha .eq. 2) write(602,*) q_alpha,(f1(i),i=1,ntyp),sum(f1),(f2(i),i=1,ntyp),sum(f2)
              if( alpha .eq. 3) write(603,*) q_alpha,(f1(i),i=1,ntyp),sum(f1),(f2(i),i=1,ntyp),sum(f2)

           endif
        end do


  end do

  !! put in the last element of overlap the total overlap normalized per total sum
  overlap(ntyp+1) = sum(overlap)/sum(normal)
  !! normalize
!  write(*,*) 'normal, overlap before norm'
!  write(*,*) normal(:)
!  write(*,*) overlap(:)
  do k = 1, ntyp
     overlap(k) = overlap(k) / normal(k)

     !! If we have divided zero by zero, which gives NaN.
     !! This can happen when ntyp is larger than the actual
     !! number of colors in any of the geometries
     if ( overlap(k) .ne. overlap(k) ) then

        !! Assuming some color(s) NOT present in both geometries at once,
        !! overlap is actually ok.
        overlap(k) = 1.0

        !! If some color present in one but not the other, overlap is not ok.
        if( present1(k) .neqv. present2(k) ) overlap(k) = 0.0

     endif
  end do

!  overlap(ntyp+1) = sum(overlap(1:ntyp))/ntyp
!  write(*,*) 'heights1'
!  write(*,*) heights1
!  write(*,*) 'heights2'
!  write(*,*) heights2
!  write(*,*) 'normal, overlap after norm'
!  write(*,*) normal(:)
!  write(*,*) overlap(:)

  end subroutine compute_overlap


  subroutine point_set_registry( n_s, typ_s, coords_s, &
                                 n_m, typ_m, coords_m, &
                                 found, dists,lat )
    !! Do the point-set-registration.
    !! In two sets of data points find the corresponding (matching) points.
    !! To each point from S find point in M such that the distance
    !! between them is minimal. The minimal distance should be very close to zero,
    !! because there should be no significant spatial transformation from S to M,
    !! and both sets S and M should be written in equivalent frames of reference.
    !! If the minimum distance for some point from S is found significantly
    !! larger than 0, then this point from S is not found in M!
    !! NOTE: Set S must be larger!! (probably?)
    !!===============
    !! -> On output, array 'found' gives the correspondance,
    !!    e.g.
    !!    found(2) = 5 means that point index 2 from S has a corresponding point
    !!    in M at index 5.
    !! -> array 'dists' gives the distances from S to M, correspoding to how the
    !!    points from S were found
    !!    e.g.
    !!    dists(3) gives the distance from point 3 in S to the point in M which
    !!    is at index found(3)
    use basic_module, only: cart_to_crist, &
                            periodic, &
                            crist_to_cart, &
                            find_loc
    implicit none
    integer, intent(in) :: n_s, n_m
    real, dimension(3,n_s), intent(in) :: coords_s
    real, dimension(3,n_m), intent(in) :: coords_m
    integer, dimension(n_s), intent(in) :: typ_s
    integer, dimension(n_m), intent(in) :: typ_m
    integer, dimension(n_s), intent(out) :: found
    real, dimension(n_s),intent(out) :: dists
    real, dimension(3,3), optional :: lat
    integer :: i, j, ind, here, n_search
    real :: dist, dist_old
    real, dimension(3) :: r_i

    found(:) = 0
    ind = 0
    dists(:) = 9999.9

    !!
    !! 'loop' through found dists
    !! i is point in S
    !!
    i = maxloc( dists, 1)
    n_search = 0
    do while( dists( i ) .gt. 1000.0 )
       !!
       !! set initial distance to something large
       !!
       dist_old = 9999.9
       ind = 0
       !!
!      write(*,*)
!      write(*,*) '>>i',i
       !!
       !! calculate distance to all points in M, keep track of minimum
       !!
       do j = 1, n_m
          !!
!          write(*,*) 'j',j
          !!
          !! sensitive to type
          !!
          if( typ_s(i) .ne. typ_m(j)) cycle
          !!
          !! calculate the distance from point in S to point in M
          !!
          r_i(:) = coords_m(:,j) - coords_s(:,i)
          if( present( lat ) ) then
             call cart_to_crist( r_i, lat )
             call periodic( r_i )
             call crist_to_cart( r_i, lat )
          endif
          dist = sqrt( dot_product( r_i, r_i ))
          !!
!          write(*,*) 'dist',dist
          !!
          if( any(found .eq. j) ) then
             !! if index j is already found for some other atom i,
             !! find which atom that is: 'here'
             call find_loc(n_s,found, j,1,here)
!             write(*,*) 'collision with', here, dists(here)
!             write(*,'(20i4)') found(:)
             if( dist .gt. dists(here))then
                !! if the distance from i to j is larger than from 'here' to j,
                !! then do nothing, cycle this atom j
!                write(*,*) 'dist',dist,'larger than prev dist found here', dists(here)
!                write(*,*) 'cycling'
                cycle
             else
                !! if distance i-j is smaller than 'here'-j, set this distance and
                !! overwrite the previously used atom 'here' such that another
                !! search will be performed there
!                write(*,*) 'dist',dist,'smaller than prev dist found here',dists(here)
!                write(*,*) 'overwriting'
                dists(here) = 9999.9
                found(here) = 0
             endif
             !!
          endif
          !!
!          write(*,*) i,j,dist
          !!
          !! if distance smaller than previous lowest, keep track
          !!
          if( dist .le. dist_old ) then
             dist_old = dist
             ind = j
          endif
          !!
!          write(*,*) 'setting found to',ind
!          write(*,*) 'setting dists to',dist_old
          !!
          !! set the found info, always overwriting old with new
          !!
          found(i) = ind
          dists(i) = dist_old

       end do
       !!
!       write(*,*) 'situation now'
!       write(*,'(20i4)') found(:)
!       write(*,*) ind,dist_old
       !!
       !! if nobody is found, set found index to 0, and sit to some large distance
       !! which is smaller than the condition for while loop
       !!
       if( dist_old .gt. 9000.0 ) then
          found(i) = 0
          dists(i) = 500.0
       endif
       !!
       !! set next search
       !!
       i = maxloc( dists, 1)
       !!
       n_search = n_search + 1
       if( n_search .gt. 2.0*n_s*n_m ) then
          write(*,*) 'PROBLEM: in point_set_registry'
          write(*,*) '         cannot find corresponding atoms.'
          exit
       endif
       !print*,"n_search,", n_search
    end do

    !  write(*,*) 'final'
    !  do i = 1, n_s
    !     write(*,*) i, found(i), dists(i)
    !  end do

  end subroutine point_set_registry

  !!!!!--------------------------------------------------------------


  !! experimental-----------------------------------------------
  subroutine set_gyration(nat,coords,gyr)
    !! gyration tensor as per wikipedia
    integer, intent(in) :: nat
    real, dimension(3,nat), intent(in) :: coords
    real, dimension(3,3), intent(out) :: gyr

    integer :: i

    gyr(:,:) = 0.0
    do i = 1, nat
       gyr(1,1) = gyr(1,1) + coords(1,i)**2
       gyr(2,2) = gyr(2,2) + coords(2,i)**2
       gyr(3,3) = gyr(3,3) + coords(3,i)**2
       gyr(1,2) = gyr(1,2) + coords(1,i)*coords(2,i)
       gyr(1,3) = gyr(1,3) + coords(1,i)*coords(3,i)
       gyr(2,3) = gyr(2,3) + coords(2,i)*coords(3,i)
    end do
    gyr(2,1) = gyr(1,2)
    gyr(3,1) = gyr(1,3)
    gyr(3,2) = gyr(2,3)

  end subroutine set_gyration

!  subroutine diag(A,k,eigvalues,vec)
!    ! if vec=1, also compute the eigenvectors
!    ! take care of the Upper or Lower triangle matrix, the
!    ! DSYEV reads only a triangle, and always thinks that
!    ! the matrix is symmetric!
!    !
!    ! A is the input matrix,
!    ! k is the dimension of the matrix (k,k)
!    ! eigvalues is output vector of dimension(k)
!    ! vec is integer(0,1): 0 if not needed eigenvectors, 1 if yes
!    !
!    !calling list
!    integer, intent(in)             :: k,vec
!    real, intent(inout) :: A(k,k), eigvalues(k)
!    !
!    !local
!    real ,allocatable :: work(:)
!    integer                      :: lwork,info
!
!    lwork = max(1,3*k-1)
!    allocate(work(lwork))
!
!    if(vec==0) then
!       call dsyev('N','U',k,A,k,eigvalues,WORK,LWORK,info)
!       ! call zhbgv('N','U',k,A,k,eigvalues,WORK,LWORK,info)
!       ! call zheev
!    else
!       call dsyev('V','U',k,A,k,eigvalues,WORK,LWORK,info)
!       ! call zhbgv('V','U',k,A,k,eigvalues,WORK,LWORK,info)
!    endif
!  end subroutine diag

  subroutine set_inertia(nat,typ,coords,ntyp,mu_k,inertia)
    implicit none

    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ
    real, dimension(3,nat), intent(in) :: coords
    integer, intent(in) :: ntyp
    real, dimension(ntyp) :: mu_k
    real, dimension(3,3), intent(out) :: inertia

    integer :: i
    real :: weight, pi

    pi = 4.0*atan(1.0)
    inertia(:,:) = 0.0
    do i = 1, nat
!       weight = 10*exp( - norm(coords(i,:)) / mu_k(typ(i)) )
!       weight = cos( norm(coords(i,:))*pi / mu_k(typ(i)) )
!       weight = -atan( norm(coords(i,:)) - mu_k(typ(i)) )/atan(mu_k(typ(i)))
!       weight = exp( - norm(coords(:,i)) / mu_k(typ(i)) )
       weight =  (typ(i)*1.0)**2
       inertia(1,1) = inertia(1,1) + (coords(2,i)**2+coords(3,i)**2)*weight
       inertia(2,2) = inertia(2,2) + (coords(1,i)**2+coords(3,i)**2)*weight
       inertia(3,3) = inertia(3,3) + (coords(1,i)**2+coords(2,i)**2)*weight
       inertia(1,2) = inertia(1,2) - coords(1,i)*coords(2,i)*weight
       inertia(1,3) = inertia(1,3) - coords(1,i)*coords(3,i)*weight
       inertia(2,3) = inertia(2,3) - coords(2,i)*coords(3,i)*weight
    end do
    inertia(2,1) = inertia(1,2)
    inertia(3,1) = inertia(1,3)
    inertia(3,2) = inertia(2,3)

  end subroutine set_inertia

 !!----------------------------------------------------------------

end module organize_module
