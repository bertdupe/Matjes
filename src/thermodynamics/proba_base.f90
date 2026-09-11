module m_proba_base
implicit none

type :: proba_data
   real(8),allocatable :: Pdistrib(:)
   real(8),allocatable :: Energy(:)
   logical             :: io_Pdist=.false.
   logical             :: is_set=.false.
   integer             :: N_bins=10

   procedure(Psampling), pointer :: sampling_type
   procedure(Pplot), pointer     :: plot_type

end type

abstract interface
    subroutine Psampling(this,delta_Energy,kT)
        import proba_data
        class(proba_data), intent(inout) :: this
        real(8),            intent(in)   ::  kT
        real(8),            intent(in)   ::  delta_Energy(:,:)
    end subroutine

    subroutine Pplot(this,tag)
        import proba_data
        class(proba_data), intent(in) :: this
        integer,intent(in)            ::  tag
    end subroutine
end interface

private
public proba_data,alloc_P_distrib

contains

!!!!!!!!!!!!!!!!!
! allocate the probability distribution
!!!!!!!!!!!!!!!!!

subroutine alloc_P_distrib(N,distrib_energy,distrib_proba)
integer, intent(in) :: N
real(8),allocatable, intent(out)  :: distrib_energy(:),distrib_proba(:)

if (allocated(distrib_proba)) stop "probability_distrib already allocated"

allocate(distrib_proba(N),source=0.0d0)
allocate(distrib_energy(N),source=0.0d0)

end subroutine

end module
