module m_fermi
    implicit none

    public :: get_fermi_k, calc_fermi,fermi_distrib
    private
    contains

        subroutine get_fermi_k(eps,N_electrons,N_cell,kt,E_f)
            implicit none
            real(kind=8), intent(out) :: E_f
            real(kind=8), intent(in) :: kt          !k_b*T for fermi dirac distribution
            real(kind=8), intent(in) :: eps(:,:) !assume [Ns,Nk] with Nk=number of k-points, Ns=number of states
            real(kind=8), intent(in) :: N_electrons !number of states to be occupied at each k
            integer,intent(in)       :: N_cell

            integer     ::  Nk,Ns
            real(8)     ::  N_occ !overall aim of number of occupied states

            Ns=size(eps(:,1))
            Nk=size(eps(1,:))
            N_occ=real(Nk,8)*N_electrons*N_cell

            Call calc_fermi(reshape(eps,[Nk*Ns]),N_occ,kt,E_f)

        end subroutine


        subroutine calc_fermi(eig, N_occ, kt, E_f)
            !subroutine which calculates the Fermi energy
            !first guess is zero temperature, difference between lowest unoccupied and higher occupied energy
            !then iterate to achieve the aimed occupation using the fermi-dirac distribution
            use m_sort
            implicit none
            real(kind=8), intent(out) :: E_f
            real(kind=8), intent(in)  :: kt      !k_b*T for fermi dirac distribution
            real(kind=8), intent(in)  :: eig(:)  !all eigen energies
            real(kind=8), intent(in)  :: N_occ !number of states to be occupied

            ! Internal variable
            integer     ::  N_eig
            integer     :: i
            integer, allocatable :: indices(:)
            real(kind=8),allocatable :: tmp_E(:) !sorted energy array

            !fermi-dirac temporary variables
            !!variables determine considered energy range
            real(8),parameter     ::  cutoff_kt=10.0d0 !how many instances of kt are included in the considered energy range
                                      !arbitraily choose 10 times kt, which should be large enough n_f(-10)=0.99995
            real(8)               ::  min_E,max_E !minimal/maximal energy for calculating the fermi energy
            integer               ::  min_ind,max_ind
            !!variables within occupation convergence loop
            integer     ::  weight_ind  
            real(8)     ::  weight_aim !fermi dirac occupation aim ( in numbers of occupied states in tmp_E(min_ind:max_ind)
            real(8)     ::  dE
            !!parameters for convergence loop control
            integer,parameter     ::  fd_loop_max=50 !maximal number of fermi-dirac iterations
            real(8),parameter     ::  occ_cutoff=1.0d-4 !accurancy of the fermi dirac occupation  
            !components for first derivative:
            real(8),parameter     ::  c_deriv(7)=[1.0d0,-9.0d0,45.0d0,0.0d0,-45.0d0,9.0d0,-1.0d0]/60.0d0 
            !temporary iteration parameters for obtaining the fermi energy from the fermi-dirac distribution
            real(8)    :: E_F_ext(2),weight_ext(2)
            real(8)    :: E_F_tmp,weight_tmp

            !get sorted eigenvalues in tmp_E
            N_eig=size(eig)
            allocate(tmp_E(N_eig),source=eig)
            allocate(indices(N_eig),source=0)
            call sort(N_eig, tmp_E, indices, 1.0d-5)

            !trivial guess for fermi energy
            if(nint(N_occ)+1 > size(tmp_E)) STOP 'Too many electrons to calculate fermi energy, reduce N_electrons?'
            E_f=(tmp_E(nint(N_occ))+tmp_E(nint(N_occ)+1))*0.5d0

            !Use Fermi-dirac distribution
            !!prepare considered energy range
            !!will assume all states below min_ind are fully occupied and above max_ind are not occupied
            min_E=E_F-cutoff_kt*kt
            max_E=E_F+cutoff_kt*kt
            min_ind=1
            do i=nint(N_occ)+1,1,-1
                if(tmp_E(i) <= min_E)then !could be done faster with an temporary array
                    min_ind=i
                    exit
                endif
            enddo
            max_ind=N_eig
            do i=nint(N_occ),N_eig
                if(tmp_E(i) >= max_E)then
                    max_ind=i
                    exit
                endif
            enddo
            weight_aim=N_occ-real(min_ind-1,kind=8)

            !Get initial fermi energy guesses and weights 
            weight_ind=floor(weight_aim)
            if (weight_ind-3<1 .or. weight_ind+3> N_eig)then
                write(*,'(A)') "Warning, very few energy values above or below the aimed occupied state"
                write(*,'(A)') "Fermi-Dirac implementation not sensible, will use initial guess Fermi energy"
                write(6,'(/a,2x,f12.6,2x,a/)') 'fermi_level = ',  E_f, '[eV]'
                return
            endif
            dE=max(dot_product(c_deriv,tmp_E(weight_ind-3:weight_ind+3)),1.0e-2)
            E_F_ext=[E_F-dE,E_F+dE]
            weight_ext(1)=fermi_weight_sum(E_F_ext(1),tmp_E(min_ind:max_ind),kt)-weight_aim
            weight_ext(2)=fermi_weight_sum(E_F_ext(2),tmp_E(min_ind:max_ind),kt)-weight_aim
            weight_ext=weight_ext

            !make sure the E_F_ext correspond to weights below and above the aim
            do while (weight_ext(1)>0.0d0 .or. weight_ext(2)<0.0d0)
                if(weight_ext(1)>0.0d0)then
                    weight_tmp=weight_ext(1)
                    E_F_tmp=E_F_ext(1)
                    E_F_ext(1)=E_F_ext(1)-(E_F_ext(2)-E_F_ext(1))*(weight_ext(1)+2.0d0)
                    weight_ext(1)=fermi_weight_sum(E_F_ext(1),tmp_E(min_ind:max_ind),kt)-weight_aim
                    weight_ext(2)=weight_tmp
                    E_F_ext(2)=E_F_tmp
                endif
                if(weight_ext(2)<0.0d0)then
                    weight_tmp=weight_ext(2)
                    E_F_tmp=E_F_ext(2)
                    E_F_ext(2)=E_F_ext(2)+(E_F_ext(2)-E_F_ext(1))*(-weight_ext(2)+2.0d0)
                    weight_ext(2)=fermi_weight_sum(E_F_ext(2),tmp_E(min_ind:max_ind),kt)-weight_aim
                    weight_ext(1)=weight_tmp
                    E_F_ext(1)=E_F_tmp
                endif
            end do

            !sucessively divide E_F and update E_F_ext and the weights until weight threshold is reached
            i=1
            weight_tmp=occ_cutoff+1.0d0
            do while (abs(weight_tmp)>occ_cutoff)
                E_F_tmp=sum(E_F_ext)*0.5d0
                weight_tmp=fermi_weight_sum(E_F_tmp,tmp_E(min_ind:max_ind),kt)-weight_aim
                if(weight_tmp<0.0d0)then
                    weight_ext(1)=weight_tmp
                    E_F_ext(1)=E_F_tmp
                else
                    weight_ext(2)=weight_tmp
                    E_F_ext(2)=E_F_tmp
                endif
                i=i+1
                if(i>fd_loop_max)then
                    write(*,'(A)') "warning, fermi dirac occupation has not reached the wanted cutoff"
                    write(*,'(A,E16.8)') "weight-difference",weight_tmp
                    write(*,'(A)') "Possibly continuing with wrong Fermi energy"
                    exit
                endif
            end do 
            E_f=E_F_tmp
            write(6,'(/a,2x,f12.6,2x,a/)') 'fermi_level = ',  E_f, '[eV]'

        end subroutine


        function fermi_weight_sum(E_F,energy,kt)result(weight)
            real(8)     :: weight
            real(8)     :: E_F
            real(8)     :: kt
            real(8)     :: energy(:)

            integer     ::  i
            weight=0.0d0
            do i=1,size(energy)
                weight=weight+fermi_distrib(E_F, energy(i), kt)
            end do
        end function

        ! Function implementing the FD distribution for a given energy
        real(kind=8) function fermi_distrib(E_F, energy, kt)
            implicit none
            real(kind=8),intent(in) :: E_F
            real(kind=8),intent(in) :: kt
            real(kind=8),intent(in) :: energy

            real(8)                 ::  exp_val
            real(8)                 ::  cutoff=30.0d0

            exp_val=(energy-E_F)/kt
            if(exp_val>cutoff)then
                fermi_distrib=0.0d0
            elseif(exp_val<cutoff)then
                fermi_distrib=1.0d0
            else
                fermi_distrib = 1.0d0/( 1.0d0 + exp(exp_val ) )
            endif
        end function fermi_distrib

end module m_fermi
