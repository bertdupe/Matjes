      module m_skyrmion
      real(kind=8),dimension(:),allocatable :: XSky,YSky,RSky,Tthreshold
      real(kind=8),dimension(:),allocatable :: XantiSky,YantiSky,RantiSky,Tantithreshold
      real(kind=8),dimension(:),allocatable :: XstarSky,YstarSky,RstarSky,Tstarthreshold
      real(kind=8),dimension(:),allocatable :: XuserSky,YuserSky,RuserSky,Tuserthreshold
      real(kind=8),dimension(:),allocatable :: Xtargetsky,Ytargetsky,Rtargetsky,Ttargetthreshold
      Logical, dimension(:),allocatable :: TabPinning
      Logical, dimension(:),allocatable :: TabantiPinning
      Logical, dimension(:),allocatable :: TabstarPinning
      Logical, dimension(:),allocatable :: TabtargetPinning
      logical :: SkX
      real(kind=8) :: qskx,coeffx,coeffy,starx,stary
      integer :: NSkyAdd,NSkyAnnihilation, N_Pinning,NuserSkyAdd
      integer :: NantiSkyAdd,NantiSkyAnnihilation, Nanti_Pinning
      integer :: NstarSkyAdd,NstarSkyAnnihilation, Nstar_Pinning
      integer :: Ntarget_Pinning,NtargetSkyAdd
      end module m_skyrmion

      subroutine rw_skyrmion()
      use m_skyrmion
      use m_SkX_utils
      use m_rw_lattice, only : net,dim_lat
      implicit none
      integer, parameter  :: io=1
      integer :: fin,i
      character(len=100) :: str
      character(len=10) :: dummy
      character(len=1)  :: pinning


      open (io,file='rwskyrmion.in',form='formatted',status='old', &
        action='read')

      SkX=.False.
      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
        if ( str(1:7) == 'NSkyAdd') then
           backspace(io)
           read(io,*) dummy,NSkyAdd
          endif
        if ( str(1:16) == 'NSkyAnnihilation') then
           backspace(io)
           read(io,*) dummy,NSkyAnnihilation
          endif
        if ( str(1:9) == 'N_Pinning') then
           backspace(io)
           read(io,*) dummy,N_Pinning
          endif
! part antiskyrmion
        if ( str(1:11) == 'NantiSkyAdd') then
           backspace(io)
           read(io,*) dummy,NantiSkyAdd
          endif
        if ( str(1:20) == 'NantiSkyAnnihilation') then
           backspace(io)
           read(io,*) dummy,NantiSkyAnnihilation
          endif
        if ( str(1:13) == 'Nanti_Pinning') then
           backspace(io)
           read(io,*) dummy,Nanti_Pinning
          endif
! part starskyrmion
        if ( str(1:11) == 'NstarSkyAdd') then
           backspace(io)
           read(io,*) dummy,NstarSkyAdd
          endif
        if ( str(1:20) == 'NstarSkyAnnihilation') then
           backspace(io)
           read(io,*) dummy,NstarSkyAnnihilation
          endif
        if ( str(1:13) == 'Nstar_Pinning') then
           backspace(io)
           read(io,*) dummy,Nstar_Pinning
          endif
! user define
        if ( str(1:11) == 'NuserSkyAdd') then
           backspace(io)
           read(io,*) dummy,NuserSkyAdd
          endif
!
        if ( str(1:7) == 'lattice') then
           backspace(io)
           read(io,*) dummy,SkX
          endif
        if ( str(1:8) == 'qlattice') then
           backspace(io)
           read(io,*) dummy,qskx
          endif
! target-skyrmion define
        if ( str(1:13) == 'NtargetSkyAdd') then
           backspace(io)
           read(io,*) dummy,NtargetSkyAdd
          endif
        if ( str(1:15) == 'Ntarget_Pinning') then
           backspace(io)
           read(io,*) dummy,Ntarget_Pinning
          endif
!
      enddo
      
      if(NSkyAdd+NSkyAnnihilation+N_Pinning .eq. 0) then

              Allocate(XSky(1))
              Allocate(YSky(1))
              Allocate(RSky(1))
              Allocate(Tthreshold(1))
              Allocate(TabPinning(1))
              Xsky(1) = 0.0d0
              YSky(1) = 0.0d0
              RSky(1) = -1.0d0
              Tthreshold(1) = -1.0d0

      else
       if (SkX) NSkyAdd=nint(qskx*dim_lat(1)*qskx*dim_lat(2))

       Allocate(XSky(NSkyAdd+NSkyAnnihilation+N_Pinning))
       Allocate(YSky(NSkyAdd+NSkyAnnihilation+N_Pinning))
       Allocate(RSky(NSkyAdd+NSkyAnnihilation+N_Pinning))
       Allocate(Tthreshold(NSkyAdd+NSkyAnnihilation+N_Pinning))
       Allocate(TabPinning(NSkyAdd+NSkyAnnihilation+N_Pinning))

      endif

! part antiskyrmion

      if(NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning .eq. 0) then

              Allocate(XantiSky(1))
              Allocate(YantiSky(1))
              Allocate(RantiSky(1))
              Allocate(Tantithreshold(1))
              Allocate(TabantiPinning(1))
              Xantisky(1) = 0.0d0
              YantiSky(1) = 0.0d0
              RantiSky(1) = -1.0d0
              Tantithreshold(1) = -1.0d0

      else
      if (SkX) NantiSkyAdd=nint(qskx*dim_lat(1)*qskx*dim_lat(2))

      Allocate(XantiSky(NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning))
      Allocate(YantiSky(NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning))
      Allocate(RantiSky(NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning))
      Allocate(Tantithreshold(NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning))
      Allocate(TabantiPinning(NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning))

      endif

! part starskyrmion

      if(NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning .eq. 0) then

              Allocate(XstarSky(1))
              Allocate(YstarSky(1))
              Allocate(RstarSky(1))
              Allocate(Tstarthreshold(1))
              Allocate(TabstarPinning(1))
              Xstarsky(1) = 0.0d0
              YstarSky(1) = 0.0d0
              RstarSky(1) = -1.0d0
              Tstarthreshold(1) = -1.0d0

      else
      if (SkX) NstarSkyAdd=nint(qskx*dim_lat(1)*qskx*dim_lat(2))

      Allocate(XstarSky(NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning))
      Allocate(YstarSky(NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning))
      Allocate(RstarSky(NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning))
      Allocate(Tstarthreshold(NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning))
      Allocate(TabstarPinning(NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning))

      endif

! part target skyrmion

      if(NtargetSkyAdd+Ntarget_Pinning .eq. 0) then

              Allocate(XtargetSky(1))
              Allocate(YtargetSky(1))
              Allocate(RtargetSky(1))
              Allocate(Ttargetthreshold(1))
              Allocate(TabtargetPinning(1))
              Xtargetsky(1) = 0.0d0
              YtargetSky(1) = 0.0d0
              RtargetSky(1) = -1.0d0
              Ttargetthreshold(1) = -1.0d0
      else

      if (SkX) NstarSkyAdd=nint(qskx*dim_lat(1)*qskx*dim_lat(2))

      Allocate(XtargetSky(NtargetSkyAdd+Ntarget_Pinning))
      Allocate(YtargetSky(NtargetSkyAdd+Ntarget_Pinning))
      Allocate(RtargetSky(NtargetSkyAdd+Ntarget_Pinning))
      Allocate(Ttargetthreshold(NtargetSkyAdd+Ntarget_Pinning))
      Allocate(TabtargetPinning(NtargetSkyAdd+Ntarget_Pinning))

      endif
! part if user defined skyrmion

      if(NuserSkyAdd .eq. 0) then

              Allocate(XuserSky(1))
              Allocate(YuserSky(1))
              Allocate(RuserSky(1))
              Allocate(Tuserthreshold(1))
              Xusersky(1) = 0.0d0
              YuserSky(1) = 0.0d0
              RuserSky(1) = -1.0d0
              Tuserthreshold(1) = -1.0d0

      else
      if (SkX) NuserSkyAdd=nint(qskx*dim_lat(1)*qskx*dim_lat(2))

      Allocate(XuserSky(NuserSkyAdd))
      Allocate(YuserSky(NuserSkyAdd))
      Allocate(RuserSky(NuserSkyAdd))
      Allocate(Tuserthreshold(NuserSkyAdd))

      endif

       rewind(io)
       do
       read (io,'(a)',iostat=fin) str
       if (fin /= 0) exit

        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
        if ( str(1:8) == 'PXYRTSky') then

         if (.not.Skx) then
          do i=1,NSkyAdd+NSkyAnnihilation+N_Pinning
           read(io,*) TabPinning(i), XSky(i),YSky(i),RSky(i),Tthreshold(i)
          enddo
         else
          read(io,*) TabPinning(1), XSky(1),YSky(1),RSky(1),Tthreshold(1)
          TabPinning=TabPinning(1)
          Tthreshold=Tthreshold(1)
          RSky=RSky(1)
         endif
        endif

! part antiskyrmion
        if ( str(1:12) == 'PXYRTantiSky') then
         if (.not.Skx) then
          do i=1,NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning
           read(io,*) TabPinning(i), XantiSky(i),YantiSky(i),RantiSky(i),Tantithreshold(i)
          enddo
         else
          read(io,*) TabantiPinning(1), XantiSky(1),YantiSky(1),RantiSky(1),Tantithreshold(1)
          TabantiPinning=TabantiPinning(1)
          Tantithreshold=Tantithreshold(1)
          RantiSky=RantiSky(1)
         endif
        endif

! part antiskyrmion
        if ( str(1:12) == 'PXYRTstarSky') then
         if (.not.Skx) then
          do i=1,NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning
           read(io,*) TabPinning(i), XstarSky(i),YstarSky(i),RstarSky(i),Tstarthreshold(i)
          enddo
         else
          read(io,*) TabstarPinning(1), XstarSky(1),YstarSky(1),RstarSky(1),Tstarthreshold(1)
          TabstarPinning=TabstarPinning(1)
          Tstarthreshold=Tstarthreshold(1)
          RstarSky=RstarSky(1)
         endif
        endif

! part user defined skyrmion
        if ( str(1:12) == 'PXYRTuserSky') then
         read(io,*) coeffx,coeffy,starx,stary
         if (.not.Skx) then
          do i=1,NuserSkyAdd
           read(io,*) XuserSky(i),YuserSky(i),RuserSky(i),Tuserthreshold(i)
          enddo
         else
          read(io,*) XuserSky(1),YuserSky(1),RuserSky(1),Tuserthreshold(1)
          Tuserthreshold=Tuserthreshold(1)
          RuserSky=RuserSky(1)
         endif
        endif

! part target skyrmions
        if ( str(1:14) == 'PXYRTtargetSky') then
         if (.not.Skx) then
          do i=1,NtargetSkyAdd
           read(io,*) TabtargetPinning(i),XtargetSky(i),YtargetSky(i),RtargetSky(i),Ttargetthreshold(i)
          enddo
         else
          read(io,*) TabtargetPinning(1),XtargetSky(1),YtargetSky(1),RtargetSky(1),Ttargetthreshold(1)
          Ttargetthreshold=Ttargetthreshold(1)
          RtargetSky=RtargetSky(1)
         endif
        endif
       enddo

       close(io)
       
       if (SkX) then
        if (NSkyAdd.ne.0) call find_XYsky(XSky,YSky,NSkyAdd,qSkX)
        if (NantiSkyAdd.ne.0) call find_XYsky(XantiSky,YantiSky,NantiSkyAdd,qSkX)
        if (NstarSkyAdd.ne.0) call find_XYsky(XstarSky,YstarSky,NstarSkyAdd,qSkX)
        if (NuserSkyAdd.ne.0) call find_XYsky(XuserSky,YuserSky,NuserSkyAdd,qSkX)
        if (NstarSkyAdd.ne.0) call find_XYsky(XtargetSky,YtargetSky,NtargetSkyAdd,qSkX)
        if ((NSkyAdd.eq.0).and.(NantiSkyAdd.eq.0).and.(NstarSkyAdd.eq.0).and.&
     & (NuserSkyAdd.eq.0)) then
         write(6,'(a)') 'In case of SkX, NSkyAdd or NantiSkyAdd or NstarSkyAdd or NtargetSkyAdd should be different from 0'
         stop
        endif
       endif

      end subroutine rw_skyrmion

      subroutine kornevien(kt)
      use m_skyrmion
      use m_constants, only : k_B
      use m_vector, only : cross
      use m_rw_lattice, only : net
      implicit none
      real(kind=8), intent(in) :: kt
      integer :: i,j,k
      real(kind=8) :: rad(3),qvec(3),kv(3,3)

      if (NSkyAdd+NSkyAnnihilation+N_Pinning /= 0) then
                    
      Do i = 1,NSkyAdd+NSkyAnnihilation+N_Pinning

         if((Tthreshold(i).gt. kT/k_B).and.(i.le.NSkyAdd)) then

           Call Creation(XSky(i),YSky(i),RSky(i), TabPinning(i),0.0d0,1.0d0)

! Once created, we do not want the program to create a skyrmion
! at each temperature step. Thus we put a negative radius, which
! basically means that Creation or annihilation does nothing.
              RSky(i) = -1.0d0
         elseif ((Tthreshold(i).gt.kT/k_B).and.(i .le. &
                          NSkyAdd+NSkyAnnihilation)) then
            Call Annihilation(XSky(i),YSky(i),RSky(i))
              RSky(i) = -1.0d0
         elseif (Tthreshold(i).gt.kT/k_B) then
            Call Pinned(XSky(i),YSky(i),RSky(i), TabPinning(i))
         endif
      enddo

      endif

      if (NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning /= 0) then

      Do i = 1,NantiSkyAdd+NantiSkyAnnihilation+Nanti_Pinning

         if((Tantithreshold(i).gt. kT/k_B).and.(i.le.NantiSkyAdd)) then
            Call Creation(XantiSky(i),YantiSky(i),RantiSky(i),TabantiPinning(i),1.0d0,1.0d0)
! Once created, we do not want the program to create a Askyrmion
! at each temperature step. Thus we put a negative radius, which
! basically means that Creation or annihilation does nothing.
              RantiSky(i) = -1.0d0
         elseif ((Tantithreshold(i).gt.kT/k_B).and.(i .le. &
                          NantiSkyAdd+NantiSkyAnnihilation)) then
            Call Annihilation(XantiSky(i),YantiSky(i),RantiSky(i))
              RantiSky(i) = -1.0d0
         elseif (Tantithreshold(i).gt.kT/k_B) then
            Call Pinned(XantiSky(i),YantiSky(i),RantiSky(i), TabantiPinning(i))
         endif
      enddo

      endif

      if (NstarSkyAdd+NstarSkyAnnihilation+Nanti_Pinning /= 0) then

      Do i = 1,NstarSkyAdd+NstarSkyAnnihilation+Nstar_Pinning

         if((Tstarthreshold(i).gt. kT/k_B).and.(i.le.NstarSkyAdd)) then
            Call Creation(XstarSky(i),YstarSky(i),RstarSky(i),TabstarPinning(i),1.0d0,3.0d0)
! Once created, we do not want the program to create a star-skyrmion
! at each temperature step. Thus we put a negative radius, which
! basically means that Creation or annihilation does nothing.
              RstarSky(i) = -1.0d0
         elseif ((Tstarthreshold(i).gt.kT/k_B).and.(i .le. &
                          NstarSkyAdd+NstarSkyAnnihilation)) then
            Call Annihilation(XstarSky(i),YstarSky(i),RstarSky(i))
              RstarSky(i) = -1.0d0
         elseif (Tstarthreshold(i).gt.kT/k_B) then
            Call Pinned(XstarSky(i),YstarSky(i),RstarSky(i), TabstarPinning(i))
         endif
      enddo

      endif

      if (NuserSkyAdd /= 0) then

      Do i = 1,NuserSkyAdd

         if((Tuserthreshold(i).gt. kT/k_B).and.(i.le.NuserSkyAdd)) then
            Call usercreation(XuserSky(i),YuserSky(i),RuserSky(i),.False.)
! Once created, we do not want the program to create a star-skyrmion
! at each temperature step. Thus we put a negative radius, which
! basically means that Creation or annihilation does nothing.
              RuserSky(i) = -1.0d0
         endif
      enddo

      endif

      if (NtargetSkyAdd /= 0) then

      Do i = 1,NtargetSkyAdd

         write(*,*) 'toto'
         write(*,*) Ttargetthreshold(i),kT/k_B,NtargetSkyAdd
         if((Ttargetthreshold(i).gt. kT/k_B).and.(i.le.NtargetSkyAdd)) then
            Call targetcreation(XtargetSky(i),YtargetSky(i),RtargetSky(i),.False.,0.0d0,1.0d0)
! Once created, we do not want the program to create a star-skyrmion
! at each temperature step. Thus we put a negative radius, which
! basically means that Creation or annihilation does nothing.
              RtargetSky(i) = -1.0d0
         endif
      enddo

      endif
      end subroutine kornevien
