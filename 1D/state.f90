module state
  use data_structures
  use parameters
  use shorthands

  implicit none

  ! state of the system: organisms

  type(org), allocatable :: o_list(:, :)      ! Properties of all organisms
  type(org_sts) :: o_stats(2)                 ! Statistical properties of the organisms
  type(fld) :: field(2)                       ! Field variables

  integer :: py = 2                           ! whether time is odd (1) or even (2)
  integer :: py_N = 1                         ! whether next time is odd (1) or even (2)

  type(org), allocatable :: o_list_prev(:)    ! stated of the system in the parent population
  type(org_sts) :: o_stats_prev

  ! State of the system: collectives

  type(coll) :: c_list(MAX_COLL, 2)           ! Properties of all collectives
  type(coll_sts) :: c_stats(2)                ! Statistical properties of collectives

  integer :: pyC_N = 2                        ! whether rep_coll is odd (1) or even (2)
  integer :: pyC = 1                          ! whether (rep_coll -1) is odd (1) or even (2)

  integer :: origins(MAX_COLL,MAX_COLL) = 0   ! from which ancestoral collective did the particles in each collective originate?


contains

  subroutine prepare_state()
    write(*, '(a)', advance = "no") "# Allocating memory for lists..."
    allocate (o_list(MAX_POP,2)) ! coordinates and phenotype of all organisms
    allocate (o_list_prev(MAX_POP))
    print *, "Done."
  end subroutine prepare_state

  subroutine cleanup_state()
    write(*, '(a)', advance = "no") "# Deallocating memory for lists..."
    deallocate(o_list)
    deallocate(o_list_prev)
    print *, "Done."
  end subroutine cleanup_state

end module state
