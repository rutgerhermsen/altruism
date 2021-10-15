module fft
  use parameters

  implicit none

  ! prepare the memory required for the FFT
  integer(4), parameter :: lensav = 3*N
  integer(4) :: ier = 0
  integer(4), parameter :: lenwrk = 2*N
  real(DP), allocatable :: work(:)
  real(DP), allocatable :: wsave(:)

contains

  subroutine prepare_fft()

    write(*, '(a)', advance = "no") "# Preparing objects for fast fourier transform."
    ! allocate spaces for FFT
    allocate(work(lenwrk))
    allocate(wsave(lensav))
    ! prepare wsave for FFT
    call cfft1i(N, wsave, lensav, ier)
    print *, "Done."

  end subroutine prepare_fft


  subroutine cleanup_fft()
    deallocate(work)
    deallocate(wsave)
  end subroutine cleanup_fft

end module fft
