module topo_params
  use iso_c_binding
  implicit none
  public

  save
  !! distance cutoff for local environment
!  real( c_double )              :: Rcut
  real( c_double )              :: Rcut_common

  !! number of atomic types
  integer( c_int )              :: ntyp

  !! distance cutoff for connections between different atomic types
  real( c_double ), allocatable :: color_cut(:,:)

  !! 'penetration depth' of atomic types for overlap calculation
  real( c_double ), allocatable :: mu_k(:)

  !! gaussian width for each atomic type for overlap calculation
  real( c_double ), allocatable :: sigma_k(:)

  !! file descriptor integer for event library
  integer( c_int ), parameter   :: evlib_fd = 777

  !! threshold for small moves: atomic move smaller than this cannot happen
!  real( c_double )              :: small_move_thr

  !! Rcut mode
  character(len=4)              :: rcut_mode

end module topo_params
