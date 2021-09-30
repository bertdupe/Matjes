#if 0
module m_store_relaxation
       interface store_relaxation
        module procedure store_relaxation_serial_2d
       end interface store_relaxation

       interface write_relaxation
        module procedure write_relaxation_kt
        module procedure write_relaxation_mpi
        module procedure write_relaxation_3d
       end interface write_relaxation
       contains

       subroutine store_relaxation_serial_2d(Relaxation,i_pos,pos,E_total,E,N_cell,kt,Magnetization,rate,cone,qp,qm)
       use m_constants, only : k_b
       implicit none
       real(kind=8), intent(inout) :: Relaxation(:,:)
       integer, intent(in) :: i_pos
       real(kind=8), intent(in) :: pos,E_total,N_cell,kt,Magnetization(:),rate,cone,qp,qm,E(:)


       Relaxation(1,i_pos)=pos
       Relaxation(2,i_pos)=E_total/N_cell
       Relaxation(3,i_pos)=kT/k_b
       Relaxation(4,i_pos)=sqrt(Magnetization(1)**2+Magnetization(2)**2+Magnetization(3)**2)/N_cell
       Relaxation(5,i_pos)=rate*100.0d0
       Relaxation(6,i_pos)=cone
       Relaxation(7,i_pos)=qp+qm
       Relaxation(8:15,i_pos)=E

       end subroutine store_relaxation_serial_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       subroutine write_relaxation_kt(Relax,kT,n_sizerelax,Cor_log)
       use m_Corre
       use m_constants, only : k_b
       implicit none
       real(kind=8), intent(inout) :: Relax(:,:)
       real(kind=8), intent(in) :: kT
       integer, intent(in) :: n_sizerelax
       logical, intent(in) :: Cor_log
! interal
       character(len=30) :: fname,toto
       integer :: i,j

       write(fname,'(f8.4)') kT/k_B
       toto=trim(adjustl(fname))

       write(fname,'(a,14a,a)')'Equilibriumphi-',(toto(i:i),i=1,len_trim(toto)),'.dat'
       OPEN(8,FILE=fname,status='unknown')

       Write(8,'(a)') "#   1:n_MC  2:Etot  3:T  4:M  5:rate  6:cone  7:Q  8:Exch  9:Zeeman  10:Ani &
   &         11:4S  12:DM  13:biq  14:dip  15:stoner  16:Chi(E)(t,T)  17:Chi(M)(t,T)  18:Chi(T)(t)"

       if (Cor_log) Relax(16,:)=correlation(Relax(2,:),n_sizerelax)
       if (Cor_log) Relax(17,:)=correlation(Relax(4,:),n_sizerelax)
       if (Cor_log) Relax(18,:)=correlation(Relax(7,:),n_sizerelax)

       Do i=1,n_sizerelax
           Write(8,'(i10,18(2x,E20.10E3))') int(Relax(1,i)),(Relax(j,i),j=2,18)
       enddo

       close(8)

       write(6,'(/,a,/)') 'Equilibrium files are written.'

       end subroutine write_relaxation_kt

       subroutine write_relaxation_mpi(Relax,kT_all,irank_working,n_sizerelax,size_table,Cor_log)
       use m_Corre
       use m_constants, only : k_b
       implicit none
       real(kind=8), intent(inout) :: Relax(:,:,:)
       real(kind=8), intent(in) :: kT_all(:)
       integer, intent(in) :: n_sizerelax,irank_working,size_table
       logical, intent(in) :: Cor_log
! interal
       character(len=30) :: fname,toto,name_proc
       integer :: i,j,k

       write(toto,'(I5)') irank_working
       name_proc=trim(adjustl(toto))

       do k=1,size_table

              write(fname,'(f8.4)') kT_all(k)/k_B
              toto=trim(adjustl(fname))

              write(fname,'(a,14a,a,14a,a)')'Equilibriumphi-',(toto(i:i),i=1,len_trim(toto)),'proc-',(name_proc(i:i),i=1,len_trim(name_proc)),'.dat'
              OPEN(8,FILE=fname,status='unknown')

              Write(8,'(a)') "#   1:n_MC  2:Etot  3:T  4:M  5:rate  6:cone  7:Q  8:Exch  9:Zeeman  10:Ani &
   &         11:4S  12:DM  13:biq  14:dip  15:stoner  16:Chi(E)(t,T)  17:Chi(M)(t,T)  18:Chi(T)(t)"

             if (Cor_log) Relax(16,:,k)=correlation(Relax(2,:,k),n_sizerelax)
             if (Cor_log) Relax(17,:,k)=correlation(Relax(4,:,k),n_sizerelax)
             if (Cor_log) Relax(18,:,k)=correlation(Relax(7,:,k),n_sizerelax)

             Do i=1,n_sizerelax
                Write(8,'(i10,18(2x,E20.10E3))') int(Relax(1,i,k)),(Relax(j,i,k),j=2,18)
             enddo

             close(8)
       enddo

       write(6,'(/,a,/)') 'Equilibrium files are written.'

       end subroutine write_relaxation_mpi

       subroutine write_relaxation_3d(Relax,kT_all,n_sizerelax,size_table,Cor_log)
       use m_Corre
       use m_constants, only : k_b
       implicit none
       real(kind=8), intent(inout) :: Relax(:,:,:)
       real(kind=8), intent(in) :: kT_all(:)
       integer, intent(in) :: n_sizerelax,size_table
       logical, intent(in) :: Cor_log
! interal
       character(len=30) :: fname,toto
       integer :: i,j,k

       do k=1,size_table

           write(fname,'(f8.4)') kT_all(k)/k_B
           toto=trim(adjustl(fname))

           write(fname,'(a,14a,a)')'Equilibriumphi-',(toto(i:i),i=1,len_trim(toto)),'.dat'
           OPEN(8,FILE=fname,status='unknown')

           Write(8,'(a)') "#   1:n_MC  2:Etot  3:T  4:M  5:rate  6:cone  7:Q  8:Exch  9:Zeeman  10:Ani &
   &     11:4S  12:DM  13:biq  14:dip  15:stoner  16:Chi(E)(t,T)  17:Chi(M)(t,T)  18:Chi(T)(t)"

           if (Cor_log) Relax(16,:,k)=correlation(Relax(2,:,k),n_sizerelax)
           if (Cor_log) Relax(17,:,k)=correlation(Relax(4,:,k),n_sizerelax)
           if (Cor_log) Relax(18,:,k)=correlation(Relax(7,:,k),n_sizerelax)

           Do i=1,n_sizerelax
              Write(8,'(i10,18(2x,E20.10E3))') int(Relax(1,i,k)),(Relax(j,i,k),j=2,18)
           enddo

           close(8)
       enddo

       write(6,'(/,a,/)') 'Equilibrium files are written.'

       end subroutine write_relaxation_3d

       end module m_store_relaxation
 #endif
