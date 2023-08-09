module lapack_wrap
  use iso_fortran_env, only: real64

  interface gtsv
     procedure zgtsv
  end interface gtsv

  interface
     pure subroutine zgtsv(n, nrhs, dl, d, du, b, ldb, info)
       import real64
       integer, intent(in) :: n, nrhs
       integer, intent(in) :: ldb
       integer, intent(out) :: info
       complex(kind=real64), dimension(*), intent(in out) :: dl, d, du
       complex(kind=real64), dimension(ldb, *), intent(in out) :: b
     end subroutine zgtsv
  end interface
end module lapack_wrap
