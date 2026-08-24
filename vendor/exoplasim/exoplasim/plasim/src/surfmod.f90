      module surfmod
      use pumamod

      character(len=80) :: version = '12.06.2014 by Edi'

      integer, parameter :: nsurdim  = 100           ! max # of variables
      integer, parameter :: nsurunit =  35           ! read unit
      integer, parameter :: nlocunit =  36           ! read unit
      integer, parameter :: ncodunit =  37           ! write unit
      integer            :: nsurnum  =   0           ! total # of records
      integer            :: nfreefo  =   1           ! 1: free format
      character (len=16) :: ysfile   = 'surface.sra' ! surface file name
      character (len=16) :: ysurnam(nsurdim) = ' '   ! Name
      integer            :: nsurcod(nsurdim) =  -1   ! Code

!     namelist parameter

      integer :: noromax = NTRU

      real :: doro(NHOR) = 0.0     ! orography


!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(doro,nfreefo,noromax,nsurcod,nsurnum,version,ysfile,ysurnam)

      end module surfmod


!     ===================
!     SUBROUTINE SURFCODE
!     ===================

      subroutine surfcode(kcode,yn)
      use surfmod
      integer :: kcode
      character (len=*) :: yn

!     Build a table of codes and corresponding names
!     After that variables may be read from <surface.sra> or <surface.txt>
!     referencing their name

      nsurnum = nsurnum + 1
      if (nsurnum >= nsurdim) then
          call mpabort('NSURDIM in surfmod too low')
      endif
      nsurcod(nsurnum) = kcode ! Code
      ysurnam(nsurnum) = yn    ! Name
      return
      end subroutine surfcode


!     ==================
!     SUBROUTINE CHECKIO
!     ==================

      subroutine checkio(kstat,ym)
      integer :: kstat               ! status of last i/o
      character (len=* ) :: ym       ! calling routine
      character (len=60) :: ymessage ! complete abort message

      if (kstat /= 0) then
         write(ymessage,'(2A)') 'I/O error in ',ym
         call mpabort(ymessage)
      endif
      return
      end


!     =========================
!     SUBROUTINE GET_SURF_ARRAY
!     =========================

      subroutine get_surf_array(yn,pa,kdim,klot,kread)
      use surfmod

      character (len=*) :: yn  ! Array name
      real :: pa(kdim,klot)    ! array to receive values
      integer :: icode         ! code to read
      integer :: kdim          ! dim of array
      integer :: klot          ! number of arrays
      integer :: kread         ! flag to indicate success
      integer :: ilot          ! arrays read
      integer :: io            ! iostat
      integer :: ih(8)         ! header
      character (len=18) :: yf ! filename
      character (len=16) :: yc ! formatted array name
      logical :: lex           ! exists

      kread = 0 ! Initialize as "not read"
      icode = 0
      ilot  = 0
      do j = 1 , nsurnum
         if (yn == trim(ysurnam(j))) then
            icode = nsurcod(j)
            yc = ysurnam(j)
            exit
         endif
      enddo

      if (icode == 0) then
         write(nud,'(" *** unknown array [",A,"] in get_surf_array")') yn
         return
      endif

      call code_surf_file(icode,yf)
      inquire(file=yf,exist=lex)
      if (lex) then
         open (ncodunit,file=yf,form='formatted')
         if (klot == 1) then ! single array (no annual cycle)
            read (ncodunit,*,iostat=io) ih(:)
            call checkio(io,'get_surf_array 1')
            call check_surf_header(ih,yf)
            read (ncodunit,*,iostat=io) pa(:,1)
            call checkio(io,'get_surf_array 2')
            ilot = 1
         else                ! annual cycle (12 or 14 months)
            do jlot = 1 , klot
               read (ncodunit,*,iostat=io) ih(:)
               if (io /= 0) exit ! end-of-file
               call check_surf_header(ih,yf)
               read (ncodunit,*,iostat=io) pa(:,jlot)
               call checkio(io,'get_surf_array 3')
               ilot = jlot ! arrays read so far
            enddo
            if (ilot < klot) then ! read less arrays than expected
               if (ilot == 12 .and. klot == 14) then ! do cyclic expansion
                  do jlot = 12,1,-1
                     pa(:,jlot+1) = pa(:,jlot) ! Shift year to indices 2-13
                  enddo
                  pa(:,14) = pa(:, 2) ! Copy January
                  pa(:, 1) = pa(:,13) ! Copy December
               elseif (ilot == 1 .and. klot == 14) then ! copy mean to all months
                  do jlot = 2,14
                     pa(:,jlot) = pa(:,1)
                  enddo
               else
                  write(nud,*) 'Error reading file ',trim(yf)
                  write(nud,*) 'Expected ',klot,' arrays'
                  write(nud,*) 'Found    ',ilot,' arrays'
                  call mpabort('Wrong surface file')
               endif
            endif
         endif ! (ilot < klot) 
         close(ncodunit)
         kread = 1
         write(nud,'(" * Read ",A10,"    <",A,"> *")') yc,yf
         if (ilot < klot) then
          write(nud, &
          '(" * Expanded ",A10,"    ",I3," to",I3," months *")') &
          yc,ilot,klot
         endif
      else
         write(nud,'(" * Init ",A10," [code =",I4,"] internally *")') yc,icode
      endif
      return
      end subroutine get_surf_array


!     ============================
!     SUBROUTINE CHECK_SURF_HEADER
!     ============================

      subroutine check_surf_header(kh,yf)
      use surfmod
      integer :: kh(8)              ! the eight-word SRA header just read
      character (len=*) :: yf       ! the file it came from

!     THE HEADER CARRIES THE GRID, AND IT USED TO BE READ AND DISCARDED.
!     Word 5 is NLON and word 6 is NLAT, written by every producer of these
!     files. get_surf_array read them into ih(:) and never looked, so a field
!     on another rung's grid was taken by the list-directed read that follows:
!     an oversized file silently truncated to the running model's NUGP, an
!     undersized one an I/O error with no grid named. Nothing downstream can
!     recover the grid afterwards, because the array it lands in is dimensioned
!     by the COMPILE-TIME NLAT.
!
!     A named refusal rather than a warning, because the alternative is a run
!     that integrates a rearranged planet and reports nothing. The project side
!     already refuses this way -- read_sra compares the whole header -- and this
!     is the model saying it for itself. world-fuh.

      if (kh(5) /= NLON .or. kh(6) /= NLAT) then
         write(nud,*) 'Surface file ',trim(yf),' is on the wrong grid'
         write(nud,*) 'Header says NLON =',kh(5),' NLAT =',kh(6)
         write(nud,*) 'This model is    NLON =',NLON,' NLAT =',NLAT
         call mpabort('Surface file on the wrong grid')
      endif

      return
      end subroutine check_surf_header


!     =========================
!     SUBROUTINE CODE_SURF_FILE
!     =========================

      subroutine code_surf_file(kcode,yfn)
      use pumamod
      integer :: kcode
      character (len=18) :: yfn

      write(yfn,'("N",I3.3,"_surf_",I4.4,".sra")') NLAT,kcode
      return
      end

      
!     ======================
!     SUBROUTINE SURFACE_INI
!     ======================

      subroutine surface_ini
      use surfmod

      integer :: ih(8)         ! current header
      integer :: il(8)         ! last    header
      integer :: icmon
      character (len=18) :: yf ! file name
      character (len=20) :: yformat = "(8E12.6)"
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(yformat)
      real    :: zy(NUGP,0:13) ! code with annual cycle
      logical :: lm(0:13)      ! month read flag
      logical :: lex           ! file exists

!     Define code - arrayname relationship
 
      call surfcode( 129,'doro'    )
      call surfcode( 172,'dls'     )

!     miscmod arrays

      call surfcode( 130,'dtnudge' )
      call surfcode(  22,'dfnudge' )

!     oceanmod arrays

      call surfcode( 169,'yclsst'  )
      call surfcode( 172,'yls'     )  ! same as dls
      call surfcode( 903,'yfsst'   )

!     icemod arrays

      call surfcode( 169,'xclsst'  )  ! same as zclsst
      call surfcode( 172,'xls'     )  ! same as dls
      call surfcode( 210,'xclicec' )
      call surfcode( 211,'xcliced' )
      call surfcode( 709,'xflxice' )

!     simba arrays

      call surfcode( 200,'dlai'    ) ! leaf area index
      call surfcode( 304,'dcveg'   ) ! carbon in biomass
      call surfcode( 305,'dcsoil'  ) ! carbon in soil
      call surfcode(1604,'dgrow'   ) ! biomass growth
      call surfcode(1605,'dagg'    ) ! above ground growth
      call surfcode(1606,'dsc'     ) ! stomatal conductance
      call surfcode(1607,'dmr'     ) ! max. roughness

!     landmod arrays

      call surfcode( 173,'dz0clim' )
      call surfcode(1730,'dz0climo')
      call surfcode( 174,'dalbcl'  )   ! background albedo
!     CODES 175 AND 176 MEAN DIFFERENT THINGS IN AND OUT, and nothing renames
!     them because .sra files already staged and postprocessed output already
!     written both depend on the numbering they have.
!
!     IN: the staged BACKGROUND albedo of the bare surface, per band, which is
!     what dalbclim1 and dalbclim2 are initialised from and what landstep falls
!     back to where there is no snow.
!
!     OUT: outmod.f90 writes dalb -- the albedo of the surface as it actually
!     is, snow, forest, glacier and sea ice included -- as code 175, and
!     pyburn.py names it `alb`/`surface_albedo`. Code 176 is not written out at
!     all. So an output 175 read as an input 175 is a background albedo that
!     has had a cryosphere blended into it.
      call surfcode( 175,'dalbcl1' )   ! background albedo (<0.75 um)
      call surfcode( 176,'dalbcl2' )   ! background albedo (>0.75 um)
      call surfcode(1740,'dalbcls' )   ! albedo for bare soil
      call surfcode(1741,'dalbclv' )   ! albedo for vegetation
      call surfcode( 232,'dglac'   )   
      call surfcode( 212,'dforest' )
      call surfcode( 229,'dwmax'   )
      call surfcode( 209,'dtclsoil')
      call surfcode( 169,'dtcl'    )   ! same as xclsst
      call surfcode( 140,'dwcl'    )

!     radmod arrays

      call surfcode( 237,'dqo3cl'  )   ! climatological ozone
!     Prescribed aerosol column optical depth, band 1. DUST-11, extended to N
!     species by CLIM-39. Species s is code 1810+s. 1811 keeps the name it has
!     so that .sra files already on disk stay readable; the rest are numbered.
      call surfcode(1811,'ddustcol')   ! prescribed species 1, band-1 column AOD
      call surfcode(1812,'ddustcol2')  ! prescribed species 2
      call surfcode(1813,'ddustcol3')  ! prescribed species 3
      call surfcode(1814,'ddustcol4')  ! prescribed species 4

!     aeromod arrays: the in-model dust emission's source map. DUST-3.
!
!     Three and no fewer, and the split is decided by where each quantity
!     enters the flux rather than by tidiness. Two of the offline chain's
!     terms are linear prefactors and can be multiplied together outside the
!     model; the third enters inside the nonlinearity and cannot.
!
!     A composite folding the drag partition into dsrcw would be WRONG: the
!     partition multiplies u* before it is squared, compared against a
!     threshold and raised to the fragmentation exponent, so it does not
!     commute out to a prefactor.
!
!     All three are pure functions of terrain, lithology, the lake solution
!     and the soil, written by aeolian/scripts/build_dust_source_fields.py.
      call surfcode(1801,'dsrcw'  )   ! erodible fraction x clipped clay fraction
      call surfcode(1802,'ddrage' )   ! MB95 drag partition of the erodible bed
      call surfcode(1803,'dwpr'   )   ! Fecan residual soil moisture w', percent

!     Scan start_data for codes and store sequence number
!     Write surf-code files if not existent

      icmon = 0
      ih(:) = 0
      il(:) = 0
      lm(:) = .false.

!     Look first for 'surface.sra' written in free format
!     otherwise  for 'surface.txt' written in (8E12.6)

      open(nsurunit,file=ysfile,form='formatted',status='old',iostat=io)
      if (io /= 0) then ! Compatibility mode for old format surface file
         close(nsurunit)
         ysfile  = 'surface.txt'
         nfreefo = 0
         open(nsurunit,file=ysfile,form='formatted',status='old',iostat=io)
         if (io /= 0) return ! Neither <surface.sra> nor <surface.txt>
      endif
      do
         read (nsurunit,*,IOSTAT=io) ih(:)
         if (io /= 0) ih(:) = 0

!        write surface code file

         if (ih(1) /= il(1) .and. il(1) > 0) then
            call code_surf_file(il(1),yf)
            inquire(file=yf,exist=lex)
            if (.not. lex) then
               if (icmon == 12 .and. lm(12) .and. .not. lm(0)) then
                  zy(:,0) = zy(:,12)
                  lm(0) = .true.
               endif
               if (icmon == 12 .and. lm(1) .and. .not. lm(13)) then
                  zy(:,13) = zy(:,1)
                  lm(13) = .true.
               endif
               open(ncodunit,file=yf,form='formatted')
               do jmon = 0 , 13
                  if (lm(jmon)) then
                     write (ncodunit,'(8i10)') il(1:2),jmon*100,il(4:8)
                     write (ncodunit,'(4e16.6)') zy(:,jmon)
                  endif
               enddo ! jmon
               close(ncodunit)
               write(nud,'(" * Created <",A,"> from <",A,">  *")') & 
                     yf,trim(ysfile)
            endif ! lex
            icmon = 0
            lm(:) = .false.
         endif ! ih(1)

         if (io /= 0) exit ! end-of-file
         imon = ih(3)
         if (imon >= 100) imon = mod(imon/100,100)
         if (imon >= 0 .and. imon <= 13) then
            if (nfreefo == 1) then 
               read (nsurunit,*,IOSTAT=io) zy(:,imon)
            else
               read (nsurunit,yformat,IOSTAT=io) zy(:,imon)
            endif
            call checkio(io,'surface_ini')
            icmon = icmon + 1
            lm(imon) = .true.
            il(:) = ih(:)
         else
            write(nud,'(A)') "*** Illegal month in header ***"
            write(nud,'(8I10)') ih(:)
         endif
      enddo
      close(nsurunit)

      return
      end subroutine surface_ini


!     ==================
!     SUBROUTINE SURFINI
!     ==================

      subroutine surfini
      use surfmod
!
!     initialize surface parameter
!
      namelist/surfmod_nl/nspinit,noromax

      if (mypid == NROOT) then
       open(11,file=surfmod_namelist)
       read(11,surfmod_nl)
       close(11)
       write(nud,'(/,"***********************************************")')
       write(nud,'("* SURFMOD ",a35," *")') trim(version)
       write(nud,'("***********************************************")')
       if (naqua /= 0) then
       write(nud,'("* AQUA planet mode - ignoring land data       *")')
       endif
       write(nud,'("* Namelist SURFMOD_NL from <surfmod_namelist> *")')
       write(nud,'("***********************************************")')
       write(nud,surfmod_nl)
      endif

      
      call glacierprep

!     Aqua planet settings

      if (nrestart == 0 .and. naqua /= 0) then
         n_sea_points = NUGP   ! all gridpoints are water
         dls(:)  = 0.0         ! land/sea mask ro water
         doro(:) = 0.0         ! gridpoint orography
         so(:)   = 0.0         ! spectral  orography
         sp(:)   = 0.0         ! spectral  pressure
         spm(:)  = 0.0         ! spectral  pressure scattered
      endif
      
!     Desert planet settings

      if (nrestart == 0 .and. ndesert /= 0) then
         n_sea_points = 0      ! No gridpoints are water
         dls(:)  = 1.0         ! land/sea mask ro water
         doro(:) = 0.0         ! gridpoint orography (flat)
         so(:)   = 0.0         ! spectral  orography
         sp(:)   = 0.0         ! spectral  pressure
         spm(:)  = 0.0         ! spectral  pressure scattered
      endif

      if (nrestart == 0 .and. naqua == 0 .and. ndesert == 0) then ! need to read start data
         call mpsurfgp('doro',doro,NHOR,1)
         call mpsurfgp('dls' ,dls ,NHOR,1)
         if (npro == 1) then ! print only in single core runs
            write(nud,'(/,"Topography read from surface file")')
            write(nud,'("Maximum: ",f10.2," [m]")') maxval(doro) / ga
            write(nud,'("Minimum: ",f10.2," [m]")') minval(doro) / ga
            write(nud,'("Mean:    ",f10.2," [m]")') ugpmean(doro) / ga
         endif
         doro(:) = doro(:) * oroscale  ! Scale orography

         if (nglspec .eq. 0) then
         
!        Compute spectral orography

         call gp2fc(doro,NLON,NLPP)
         call fc2sp(doro,so)
         call mpsum(so,1)
         if (npro == 1) then ! print only in single core runs
            call sp2fc(so,doro)
            call fc2gp(doro,nlon,nlpp)
            write(nud,'(/,"Topography after spectral fitting")')
            write(nud,'("Maximum: ",f10.2," [m]")') maxval(doro) / ga
            write(nud,'("Minimum: ",f10.2," [m]")') minval(doro) / ga
            write(nud,'("Mean:    ",f10.2," [m]")') ugpmean(doro) / ga
         endif
         if (mypid == NROOT) then
            so(:) = so(:) / (cv*cv)
            if (noromax < NTRU) then
             jr=-1
             do jm=0,NTRU
              do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if(jn > noromax) then
                so(jr)=0.
                so(ji)=0.
               endif
              enddo
             enddo
            endif ! (noromax < NTRU)

!           Initialize surface pressure

           if (nspinit > 0) then
              sp(:) = -so(:)*cv*cv / (gascon * tgr)
           endif
        endif ! (mypid == NROOT)
        call mpscsp(sp,spm,1)

        endif
        
!*       adjust land sea mask  (workaround)

         dls(:)=AMAX1(dls(:),0.)
         dls(:)=AMIN1(dls(:),1.)

      endif ! (nrestart == 0)

                            call landini  ! land module
                            call glacierini(noromax) !glacier module
      if (nveg > 0)         call vegini   ! vegetation module
      if (n_sea_points > 0) call seaini   ! sea module

      return
      end subroutine surfini

!     ===================
!     SUBROUTINE SURFSTEP
!     ===================

      subroutine surfstep
      use surfmod

      if (naqua == 0)       call landstep
                            call glacierstep
      if (n_sea_points > 0) call seastep

      return
      end subroutine surfstep

!     ===================
!     SUBROUTINE SURFSTOP
!     ===================

      subroutine surfstop
      use surfmod

                            call landstop
                            call glacierstop
      if (nveg > 0)         call vegstop
      if (n_sea_points > 0) call seastop

      return
      end subroutine surfstop
