module m_dipolar_fft_internal
!contains different functions to get dipolar effect field or set dipolar internal magnetization
!which are called in M_internal or H_internal
public


abstract interface
    subroutine int_set_M(M,M_in,dim_lat,N_rep,dim_mode)
        integer,intent(in)          :: dim_mode
        integer,intent(in)          :: dim_lat(3)
        integer,intent(in)          :: N_rep(3)
        real(8),intent(inout)       :: M   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
        real(8),intent(in)          :: M_in(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))
    end subroutine
end interface

abstract interface
    subroutine int_get_H(H,H_out,dim_lat,N_rep,dim_mode)
        integer,intent(in)          :: dim_mode
        integer,intent(in)          :: dim_lat(3)
        integer,intent(in)          :: N_rep(3)
        real(8),intent(in)          :: H    (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
        real(8),intent(inout)       :: H_out(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))
    end subroutine
end interface

contains
subroutine set_M_period_TTT(M,M_in,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(inout)       :: M   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(in)          :: M_in(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))

    M=M_in
end subroutine

subroutine set_M_period_TTF(M,M_in,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(inout)       :: M   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(in)          :: M_in(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))

    M=0.0d0
    M(:,:,:,1:dim_lat(3))=M_in
end subroutine

subroutine set_M_period_TFF(M,M_in,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(inout)       :: M   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(in)          :: M_in(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))
    integer     ::  i3

    M=0.0d0
    do i3=1,dim_lat(3)
        M(:,:,1:dim_lat(2),i3)=M_in(:,:,:,i3)
    enddo
end subroutine

subroutine set_M_period_FFF(M,M_in,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(inout)       :: M   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(in)          :: M_in(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))
    integer     ::  i3, i2

    M=0.0d0
    do i2=1,dim_lat(2)
        do i3=1,dim_lat(3)
            M(:,1:dim_lat(1),i2,i3)=M_in(:,:,i2,i3)
        enddo
    enddo
end subroutine

subroutine set_H_period_TTT(H,H_out,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(in)          :: H   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(inout)       :: H_out(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))

    H_out=H
end subroutine

subroutine set_H_period_TTF(H,H_out,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(in)          :: H   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(inout)       :: H_out(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))

    H_out=H(:,:,:,1:dim_lat(3))
end subroutine

subroutine set_H_period_TFF(H,H_out,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(in)          :: H   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(inout)       :: H_out(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))
    integer     ::  i3

    do i3=1,dim_lat(3)
        H_out(:,:,:,i3)=H(:,:,1:dim_lat(2),i3)
    enddo
end subroutine

subroutine set_H_period_FFF(H,H_out,dim_lat,N_rep,dim_mode)
    integer,intent(in)          :: dim_mode
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(in)          :: H   (dim_mode,N_rep  (1),N_rep  (2),N_rep  (3))
    real(8),intent(inout)       :: H_out(dim_mode,dim_lat(1),dim_lat(2),dim_lat(3))
    integer     ::  i3, i2

    do i2=1,dim_lat(2)
        do i3=1,dim_lat(3)
            H_out(:,:,i2,i3)=H(:,1:dim_lat(1),i2,i3)
        enddo
    enddo
end subroutine
end module
