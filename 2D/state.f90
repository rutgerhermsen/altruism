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

contains

  subroutine prepare_state()
    write(*, '(a)', advance = "no") "# Allocating memory for lists..."
    allocate (o_list(MAX_POP,2)) ! coordinates and phenotype of all organisms
    print *, "Done."
  end subroutine prepare_state

  subroutine cleanup_state()
    write(*, '(a)', advance = "no") "# Deallocating memory for lists..."
    deallocate(o_list)
    print *, "Done."
  end subroutine cleanup_state

end module state
