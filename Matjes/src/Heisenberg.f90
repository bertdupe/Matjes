      module m_Heisenberg
      type R_order_param
       real (Kind=DEF_DBL_PREC), dimension(:,:) :: vector
       integer :: nline,ncol
      end type

      type I_order_param
       complex (Kind=DEF_DBL_PREC), dimension(:,:) :: vector
       integer :: nline,ncol
      end type

      contains

      end module m_Heisenberg
