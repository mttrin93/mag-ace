!
! various utils
!
!----------------------------------------------------------------------------------------------------------
      
subroutine aceinit
   
    ! subroutine sets up required arrays for ace and reads parameter file
    !   
    use core
    use rw
    use global
    use functionparameters
    use modneigh
    use modtime
    implicit none
    logical doread, readheaderonly, last, converged, firstread, dostop, dofit
    integer nint
    integer neleread
    character(2) ele1, ele2
    character(2) eleread(1:nelements)
    integer nelemap(1:nelements)
    integer, allocatable :: nelemapplan(:,:), nelereadplan(:)
    character(2), allocatable:: elereadplan(:,:)
    integer i, j, k, l, m, n
    integer nele1, nele2, nread, nreadmax

    ! us c0 for setting doneinit -> may be moved somewhere else later
    if (.not.allocated(c0)) then
        doneinit = .false.
    endif

    if (.not.doneinit) then
        nreadmax = 1000
      
        ! check array 
        if (.not.allocated(element)) then
            print *,'Array elements needs to be available for initialization.'
            stop
        endif
        ! check if elements are ordered alphabetically
        do n = 1, nelements - 1
            ele1 = element(n)
            ele2 = element(n+1)
            if (ele1 == ele2 ) then
                print *,'Found two identical elements in element list. Stopping.'
                stop
            endif
            if ( ele1 > ele2 ) then
                print *,'Element list needs to be ordered alphabetically. Stopping.'
                stop
            endif         
        enddo

        tunpack = 0.d0
        tprepA = 0.d0
        tprepB = 0.d0
        tforce = 0.d0
        tradbase = 0.d0
        tradfunc = 0.d0
        tradpack = 0.d0
        tgetsh = 0.d0
        tpackA = 0.d0
        tcopyA = 0.d0

        xcountatom = 0.d0
        
        if (dorun) then
            ! no output for dorun
            debuglevel = 0
            verblevel = 0
            doforce = .true.
            dotorque = .true.
            doquick = .false.
            dotrain = .false.
            dotest = .false.
        endif

        allocate( lambda(1:2,1:nelements,1:nelements) )
        allocate( cut(1:2,1:nelements,1:nelements) )
        allocate( dcut(1:2,1:nelements,1:nelements) )
        allocate( ecut(1:nelements))
        allocate( decut(1:nelements))
      
        ! collect header information from parameter file and read available parameters
        allocate( elereadplan(1:nelements,1:nreadmax) )
        allocate( nelemapplan(1:nelements,1:nreadmax) )
        allocate( nelereadplan(1:nreadmax))
        nread = 0
        do n = 1, nelements
            nread = nread + 1
            nelereadplan(nread) = 1
            elereadplan(1,nread) = element(n)
            nelemapplan(1,nread) = n
        enddo
        if (nelements >= 2 ) then
            do n = 1, nelements-1
                do m = n+1, nelements
                    nread = nread + 1
                    nelereadplan(nread) = 2
                    elereadplan(1,nread) = element(n)
                    elereadplan(2,nread) = element(m)               
                    nelemapplan(1,nread) = n
                    nelemapplan(2,nread) = m
                enddo
            enddo
        endif
        if (nelements >= 3 ) then
            do n = 1, nelements-2
                do m = n+1, nelements-1
                    do k = m+1,nelements
                        nread = nread + 1
                        nelereadplan(nread) = 3
                        elereadplan(1,nread) = element(n)
                        elereadplan(2,nread) = element(m)
                        elereadplan(3,nread) = element(k)                  
                        nelemapplan(1,nread) = n
                        nelemapplan(2,nread) = m
                        nelemapplan(3,nread) = k                  
                    enddo
                enddo
            enddo
        endif
        if (nelements >= 4 ) then
            do n = 1, nelements-3
                do m = n+1, nelements-2
                    do k = m+1,nelements-1
                        do j = k+1,nelements
                            nread = nread + 1
                            nelereadplan(nread) = 4
                            elereadplan(1,nread) = element(n)
                            elereadplan(2,nread) = element(m)
                            elereadplan(3,nread) = element(k)
                            elereadplan(4,nread) = element(j)
                            nelemapplan(1,nread) = n
                            nelemapplan(2,nread) = m
                            nelemapplan(3,nread) = k
                            nelemapplan(4,nread) = j
                        enddo
                    enddo
                enddo
            enddo
        endif
        if (nelements >= 5 ) then
            do n = 1, nelements-4
                do m = n+1, nelements-3
                    do k = m+1,nelements-2
                        do j = k+1,nelements-1
                            do l = j+1,nelements
                                nread = nread + 1
                                nelereadplan(nread) = 5
                                elereadplan(1,nread) = element(n)
                                elereadplan(2,nread) = element(m)
                                elereadplan(3,nread) = element(k)
                                elereadplan(4,nread) = element(j)
                                elereadplan(5,nread) = element(l)                     
                                nelemapplan(1,nread) = n
                                nelemapplan(2,nread) = m
                                nelemapplan(3,nread) = k
                                nelemapplan(4,nread) = j
                                nelemapplan(5,nread) = l                        
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        endif
        if (nelements >= 6 ) then
            print *,'Please update aceini to read more than 5 elements.'
            stop
        endif

        if (nread > nreadmax) then
            print *,'Code should have crashed ...'
            print *,'Please increase nreadmax in aceini to more than: ', nread
            print *,'Then recompile.'
            stop
        endif

        ! now read all required files A.XC.in, B.XC.in, C.XC.in,..., AB.XC.in, AC.XC.in, ...., ABC.XC.in, ...    

        firstread = .true.
        do n = 1, nread
            neleread = nelereadplan(n)
            eleread(1:neleread) = elereadplan(1:neleread,n)
            nelemap(1:neleread)  = nelemapplan(1:neleread,n)
                
            doread = .true.
            last = .false.
            converged = .false.
            nint = 0
            readheaderonly = .true.
            call rwpara(doread,readheaderonly,last,converged,nint,neleread,eleread,nelemap)
        
            ! set basis parameters from first read
            if (firstread) then
                firstread = .false.
                nradbase_atomic = nradbaset_atomic
                nradbase_magnetic = nradbaset_magnetic
                nradial2_magnetic = nradial2t_magnetic
                l2max_magnetic = l2maxt_magnetic
                nradial3_atomic = nradial3t_atomic
                nradial3_magnetic = nradial3t_magnetic
                l3max_magnetic = l3maxt_magnetic
                l3max_atomic = l3maxt_atomic
                nradial4_atomic = nradial4t_atomic
                nradial4_magnetic = nradial4t_magnetic
                l4max_atomic = l4maxt_atomic
                l4max_magnetic = l4maxt_magnetic
                nradial5_atomic = nradial5t_atomic
                nradial5_magnetic = nradial5t_magnetic
                l5max_atomic = l5maxt_atomic
                l5max_magnetic = l5maxt_magnetic         
                nradial6_atomic = nradial6t_atomic
                nradial6_magnetic = nradial6t_magnetic         
                l6max_atomic = l6maxt_atomic
                l6max_magnetic = l6maxt_magnetic
                ndensity = ndensityt
                npot = npott
                nfemb = nfembt
                radtype_atomic = radtypet_atomic
                radtype_magnetic = radtypet_magnetic
                nradbase = max(nradbase_atomic,nradbase_magnetic)
                nradial = max(nradial3_atomic,nradial2_magnetic)
                lmax = max(l3max_atomic,l2max_magnetic)    
                nradial4 = max(nradial4_atomic,nradial4_magnetic)
                l4max = max(l4max_atomic,l4max_magnetic)
                nradial5 = max(nradial5_atomic,nradial5_magnetic)
                l5max = max(l5max_atomic,l5max_magnetic)
                nradial6 = max(nradial6_atomic,nradial6_magnetic)
                l6max = max(l6max_atomic,l6max_magnetic)
            
                ! check for limits
                if ((l3max_atomic<l4max_atomic).or.(l3max_atomic<l5max_atomic) &
                                                & .or.(l3max_atomic<l6max_atomic)) then
                    print *,' l3max_atomic needs to be largest.'
                    stop
                endif
                if ((l2max_magnetic<l3max_magnetic).or.(l2max_magnetic<l4max_magnetic) &
                        & .or.(l2max_magnetic<l5max_magnetic).or.(l2max_magnetic<l6max_magnetic)) then
                    print *,' l2max_magnetic needs to be largest.'
                    stop
                endif
                if ((lmax<l4max).or.(lmax<l5max).or.(lmax<l6max).or.(lmax<l7max)) then
                    print *,' lmax needs to be largest.'
                endif
                if ((nradial3_atomic<nradial4_atomic).or.(nradial3_atomic<nradial5_atomic) &
                                                & .or.(nradial3_atomic<nradial6_atomic)) then
                    print *,' nradial_atomic needs to be largest.'
                    stop
                endif
                if ((nradial2_magnetic<nradial3_magnetic).or.(nradial2_magnetic<nradial4_magnetic).or. &
                        & (nradial2_magnetic<nradial5_magnetic).or.(nradial2_magnetic<nradial6_magnetic)) then
                    print *,' nradial_magnetic needs to be largest.'
                    stop
                endif         
                if ((nradial<nradial4).or.(nradial<nradial5).or.(nradial<nradial6).or.(nradial<nradial7)) then
                    print *,' nradial needs to be largest.'
                    stop
                endif
                if (nradial3_atomic>nradbase_atomic) then
                    print *,' nradial_atomic < nradbase_atomic required'
                    stop
                endif
                if (nradial2_magnetic>nradbase_magnetic) then
                    print *,' nradial_magnetic < nradbase_magnetic required'
                    stop
                endif
                if (nradial>nradbase) then
                    print *,' nradial < nradbase required'
                    stop
                endif

                ! now is the time to set up all required tables
                if (verblevel > 0 ) then
                    write(*,'(A)') 'Preparing tables'
                endif
                ! set up tables required for dealing with spherical harmonics and related
                call setuptables
            
                if (.not. allocated(c0)) then
                    allocate(crad_atomic(1:nradbase_atomic,1:nradial3_atomic, &
                                    &   0:l3max_atomic,1:nelements,1:nelements))
                    crad_atomic = 0.d0
                    do nele1 = 1, nelements
                        do nele2 = 1, nelements
                            do m = 1, nradial3_atomic
                                crad_atomic(m,m,0:l3max_atomic,nele1,nele2) = 1.d0
                            enddo
                        enddo
                    enddo            
                    ! crad for magnetic radial functions
                    allocate(crad_magnetic(1:nradbase_magnetic,1:nradial2_magnetic,  &
                                    &   0:l2max_magnetic,1:nelements,1:nelements))
                    crad_magnetic = 0.d0
                    do nele1 = 1, nelements
                        do nele2 = 1, nelements
                            do m = 1, nradial2_magnetic
                                crad_magnetic(m,m,0:l2max_magnetic,nele1,nele2) = 1.d0
                            enddo
                        enddo
                    enddo
                    ! allocate logicals for crad
                    allocate(fitcrad_atomic(1:nradbase_atomic,1:nradial3_atomic, &
                                                &   0:l3max_atomic,1:nelements,1:nelements))
                    allocate(fitcrad_magnetic(1:nradbase_magnetic,1:nradial2_magnetic,  &
                                                &   0:l2max_magnetic,1:nelements,1:nelements))
                    fitcrad_atomic = .false.
                    fitcrad_magnetic = .false.

                    ! initialize parameters for core repulsion
                    allocate(corepara(1:2,1:nelements,1:nelements))
                    ! initialize with random value for prefactor and exponent
                    corepara = 3.d0
                    allocate(fitcorepara(1:2,1:nelements,1:nelements))
                    fitcorepara = .true.
                    
                    allocate(c0(1:nvar1,1:ndensity,1:nelements))
                    c0 = 0.d0
                    allocate(fitc0(1:nvar1,1:ndensity,1:nelements))
                    fitc0 = .false.
                    allocate(c2(1:nvar2,1:ndensity,1:nelements))
                    c2 = 0.d0
                    allocate(fitc2(1:nvar2,1:ndensity,1:nelements))
                    fitc2 = .false.
                    allocate(c3(1:nvar3,1:ndensity,1:nelements))
                    c3 = 0.d0
                    allocate(fitc3(1:nvar3,1:ndensity,1:nelements))
                    fitc3 = .false.
                    allocate(c4(1:nvar4,1:ndensity,1:nelements))
                    c4 = 0.d0
                    allocate(fitc4(1:nvar4,1:ndensity,1:nelements))
                    fitc4 = .false.
                    allocate(c5(1:nvar5,1:ndensity,1:nelements))
                    c5 = 0.d0
                    allocate(fitc5(1:nvar5,1:ndensity,1:nelements))
                    fitc5 = .false.
                    allocate(c6(1:nvar6,1:ndensity,1:nelements))
                    c6 = 0.d0
                    allocate(fitc6(1:nvar6,1:ndensity,1:nelements))
                    fitc6 = .false.
                endif
            
                if (.not. allocated(p1)) then
                    nvarp1 = 2
                    allocate(p1(1:nvarp1,1:ndensity,1:nelements))
                    p1 = 1.d0
                    allocate(fitp1(1:nvarp1,1:ndensity,1:nelements))
                    fitp1 = .false.
                endif        
            
            else
                ! all other required files need to have the same basis settings (except for cut, dcut and lambda)
                dostop = .false.
                if (nradbase_atomic /= nradbaset_atomic ) dostop = .true.
                if (nradbase_magnetic /= nradbaset_magnetic ) dostop = .true.
                if (nradial2_magnetic /= nradial2t_magnetic ) dostop = .true.
                if (l2max_magnetic /= l2maxt_magnetic ) dostop = .true.
                if (nradial3_atomic /= nradial3t_atomic ) dostop = .true.
                if (nradial3_magnetic /= nradial3t_magnetic ) dostop = .true.
                if (l3max_atomic /= l3maxt_atomic ) dostop = .true.
                if (l3max_magnetic /= l3maxt_magnetic ) dostop = .true.
                if (nradial4_atomic /= nradial4t_atomic ) dostop = .true.
                if (nradial4_magnetic /= nradial4t_magnetic ) dostop = .true.
                if (l4max_atomic /= l4maxt_atomic ) dostop = .true.
                if (l4max_magnetic /= l4maxt_magnetic ) dostop = .true.         
                if (nradial5_atomic /= nradial5t_atomic ) dostop = .true.
                if (nradial5_magnetic /= nradial5t_magnetic ) dostop = .true.
                if (l5max_atomic /= l5maxt_atomic ) dostop = .true.
                if (l5max_magnetic /= l5maxt_magnetic ) dostop = .true.
                if (nradial6_atomic /= nradial6t_atomic ) dostop = .true.
                if (nradial6_magnetic /= nradial6t_magnetic ) dostop = .true.
                if (l6max_atomic /= l6maxt_atomic ) dostop = .true.
                if (l6max_magnetic /= l6maxt_magnetic ) dostop = .true.         
                if (ndensity /= ndensityt ) dostop = .true.
                if (npot /= npott ) dostop = .true.
                if (nfemb /= nfembt ) dostop = .true.
                if (radtype_atomic /= radtypet_atomic ) dostop = .true.
                if (radtype_magnetic /= radtypet_magnetic ) dostop = .true.
                if (dostop) then
                    print *,'Inconsistent basis requirements.'
                    stop
                endif
            endif
            if (n == nread) then
                ! make sure that only the current model parameters are fitted, not previous ones
                ! i.e. for a fit of ABC, do not touch fits for A, B, C, AB, AC, BC
                
                ! overwrite all previous settings
                fitc0 = .false.
                fitc2 = .false.
                fitc3 = .false.
                fitc4 = .false.
                fitc5 = .false.
                fitc6 = .false.
                fitcrad_atomic = .false.
                fitcrad_magnetic = .false.
                ! next need to free fit parameters for the current model
                ! this is done by reading current model input, i.e. ABC.XC.in or in getnparatot
            endif

            ! here we have read the header once and all tables are in place, so we can continue reading the complete parameter file 
            doread = .true.
            last = .false.
            converged = .false.
            nint = 0
            readheaderonly = .false.
            call rwpara(doread,readheaderonly,last,converged,nint,neleread,eleread,nelemap)
        enddo
   
        if (.not.dotrainrad) then
            fitcrad_atomic = .false.
            fitcrad_magnetic = .false.
            if (dotrainrad) then
                print *,'Something wrong with too many logical flags.'
                stop
            endif
        endif

        ! check cutoff consistency   
        do nele1 = 1, nelements
            do nele2 = 1, nelements
                if ( dabs(cut(1,nele1,nele2) - cut(1,nele2,nele1)) > 1.d-10 ) then
                    print *,'Asymmetric atomic cutoffs. Stopping.'
                    stop
                endif
                if ( dabs(cut(2,nele1,nele2) - cut(2,nele2,nele1)) > 1.d-10 ) then
                    print *,'Asymmetric magnetic cutoffs. Stopping.'
                    stop
                endif
                if (cut(1,nele1,nele2) > cutoff) then
                    print *,'Global cutoff too short.'
                    print *,'global cutoff:', cutoff
                    print *,'cutoff: ', element(nele1),'-', element(nele2),':', cut(1,nele1,nele2)
                    stop
                endif
                if (cut(2,nele1,nele2) > cutoff_magnetic) then
                    print *,'Global cutoff too short.'
                    print *,'global cutoff:', cutoff_magnetic
                    print *,'cutoff: ', element(nele1),'-', element(nele2),':', cut(2,nele1,nele2)
                    stop
                endif
            enddo
        enddo
    
#ifdef USELOOKUP
        call setuplookuptables
#endif
    
        ! this completes initialization      
        doneinit = .true.

    endif
   
end subroutine aceinit
    
!----------------------------------------------------------------------------------------------------------

subroutine rwpara(doread,readheaderonly,last,converged,nint,neleread,eleread,nelemap)
  
    use core
    use rw
    use global
    use functionparameters
    use tables
    use modneigh
    implicit none
    integer n, nd, np, nr, nr1, nr2, nr3, nr4, nr5, nr6, na1, na, na2, na3
    integer nr0_p, nr1_p, nr2_p, nr3_p, nr4_p, nr5_p, nr6_p, n0_p, n1_p, n2_p, n3_p, n4_p, n5_p
    integer l0_p, l1_p, l2_p, l3_p, l4_p, l5_p, l6_p, lint01_p, lint23_p, lint45_p, lint0123_p
    integer l, l1, l2, l3, l4, l5, l6, lint, lint1, lint2, lint3, lint12, lint34, lint56
    integer m, nv, nclu, norder, nint, nele, nele1, nele2
    integer year, month, day, hour, minutes, seconds
    integer neleread
    integer nelemap(1:nelements)
    logical doread, last, converged, readheaderonly, fileisthere, check
    logical fitclu, fitpara, fitclu_p
    logical match, got1, got2, got3, got4, got5, got6
    character(40) str1, str2, str3, str4, xcread
    character(len=40) compound
    character(len=40) file
    character(12) paratype
    double precision xclu, xpara, xclu_p, smearingread
    character(20) str
    character(2) eleread(1:nelements), eletest(1:nelements), ele, ele1, ele2
    integer nelet
    character(2) elet(100)
    integer date_time(8)
    character(10) b(3)
    character(12) test

    nelet = 100

    ! check if elements are ordered alphabetically
    do n = 1, neleread - 1
        ele1 = eleread(n)
        ele2 = eleread(n+1)
        if (ele1 == ele2 ) then
            print *,'Found two identical elements in element list. Stopping.'
            stop
        endif
        if ( ele1 > ele2 ) then
            print *,'Element list needs to be ordered alphabetically. Stopping.'
            stop
        endif
    enddo
  
  
    if (neleread == 1) then
        compound = trim(eleread(1))
    elseif (neleread == 2) then
        compound = trim(eleread(1))//trim(eleread(2))
    elseif (neleread == 3) then
        compound = trim(eleread(1))//trim(eleread(2))//trim(eleread(3))
    elseif (neleread == 4) then
        compound = trim(eleread(1))//trim(eleread(2))//trim(eleread(3))//trim(eleread(4))
    elseif (neleread == 5) then
        compound = trim(eleread(1))//trim(eleread(2))//trim(eleread(3))//trim(eleread(4))//trim(eleread(5))
    elseif (neleread >= 6 ) then
        print *,'Please update rwpara to read more than 5 elements.'
        stop
    endif
  
    if (doread) then
        file = trim(compound)//'.'//trim(XC)//'.in'
        inquire(FILE=file,EXIST=fileisthere)
        if (.not.fileisthere) then
            print *,'File not found: ', trim(file)
            stop
        endif
        open(UNIT=1,FILE=file,STATUS='OLD')
    elseif (last) then
        file = trim(compound)//'.'//trim(XC)//'.out'
        open(UNIT=1,FILE=file,STATUS='UNKNOWN')
    else
        if ( nint > 0 ) then
            write (str, *) nint
            str = adjustl(str)
            file = trim(compound)//'.'//trim(XC)//'.'//trim(str)
            open(UNIT=1,FILE=file,STATUS='UNKNOWN')
        elseif ( nint == 0 ) then
            file = trim(compound)//'.'//trim(XC)//'.int'        
            open(UNIT=1,FILE=file,STATUS='UNKNOWN')
        else
            file = trim(compound)//'.'//trim(XC)//'.cur'        
            open(UNIT=1,FILE=file,STATUS='UNKNOWN')
        endif
    endif
  
    if(doread) then
        read(1,*)
    else
        if (converged) then
            write(1,*) ' Parameters (converged)'
        elseif (last) then
            write(1,*) ' Parameters (last)'
        else
            write(1,*) ' Parameters (intermediate)'
        endif
    endif
  
    if (doread) then
        read(1,*) str, year, month, day, hour, minutes, seconds
        if ((verblevel > 0 ).and.(readheaderonly)) then
            write(*,'(A,A)') 'Reading parameter file for: ', compound
            write(*,'(A,i4.4,A,i2.2,A,i2.2,A,i2.2,A,i2.2,A,i2.2)') 'Created y',year,'m',month,'d',day,' ',hour &
                    &,':',minutes,':',seconds
        endif
    else
        ! get current time
        call date_and_time(b(1), b(2), b(3), date_time)
        year    = date_time(1)
        month   = date_time(2)
        day     = date_time(3)
        hour    = date_time(5)
        minutes = date_time(6)
        seconds = date_time(7)
        write(1,*) 'created', year, month, day, hour, minutes, seconds
    endif
    if (doread) then    
        read(1,*) str, eletest(1:neleread)
        do n = 1, neleread
            if ( eletest(n) /= eleread(n) ) then
            print *,'Inconsistent element information in file: ', trim(compound)
            stop
            endif
        enddo
        read(1,*) str, XCread
        if (trim(XCread) /= trim(XC)) then
            print *,'Inconsistent XC functional in file: ', trim(compound)
            stop
        endif
        read(1,*) str, smearingread
        if (dabs(smearing-smearingread) > 1.d-6) then
            print *,'Inconsistent gauss-smearing in file: ', trim(compound)
            stop
        endif     
    else
        write(1,*) 'Elements ', eleread(1:neleread)
        write(1,*) 'XC ', xc
        write(1,*) 'gauss-smearing ', smearing
    endif
    if (doread) then
        if (neleread > 2) then
            read(1,*)
            read(1,*)
        elseif (neleread == 1) then
            read(1,*) str, lambda(1,nelemap(1),nelemap(1)), lambda(2,nelemap(1),nelemap(1))    
            read(1,*) str1, str2, str3, str4, cut(1,nelemap(1),nelemap(1)),  &  
            & cut(2,nelemap(1),nelemap(1)), dcut(1,nelemap(1),nelemap(1)),  &          
            & dcut(2,nelemap(1),nelemap(1)), ecut(nelemap(1)), decut(nelemap(1))
        elseif  (neleread == 2) then
            read(1,*) str1, str2, lambda(1,nelemap(1),nelemap(2)), lambda(2,nelemap(1),nelemap(2)),  &              
                &  lambda(1,nelemap(2),nelemap(1)), lambda(2,nelemap(2),nelemap(1))
            read(1,*) str1, str2, cut(1,nelemap(1),nelemap(2)), cut(2,nelemap(1),nelemap(2)),  &
                &  dcut(1,nelemap(1),nelemap(2)), dcut(2,nelemap(1),nelemap(2))
            cut(1,nelemap(2),nelemap(1)) =  cut(1,nelemap(1),nelemap(2))
            cut(2,nelemap(2),nelemap(1)) =  cut(2,nelemap(1),nelemap(2))
            dcut(1,nelemap(2),nelemap(1)) = dcut(1,nelemap(1),nelemap(2))
            dcut(2,nelemap(2),nelemap(1)) = dcut(2,nelemap(1),nelemap(2))
        else
            write(*,*) 'Should not be here. Stopping.'
            stop
        endif
    else
        if (neleread > 2) then
            write(1,*) ' for unary and binary only '
            write(1,*) ' for unary and binary only '
        elseif (neleread == 1) then
            write(1,*) 'lambda ', lambda(1,nelemap(1),nelemap(1)), lambda(2,nelemap(1),nelemap(1))
            write(1,*) 'rcut dcut ecut decut ',cut(1,nelemap(1),nelemap(1)),  & 
                    &  cut(2,nelemap(1),nelemap(1)), dcut(1,nelemap(1),nelemap(1)),   & 
                    &  dcut(2,nelemap(1),nelemap(1)), ecut(nelemap(1)), decut(nelemap(1))
        elseif  (neleread == 2) then
            write(1,*) 'lambdAB lambdaBA ',lambda(1,nelemap(1),nelemap(2)), lambda(2,nelemap(1),nelemap(2)),  &              
                    &  lambda(1,nelemap(2),nelemap(1)), lambda(2,nelemap(2),nelemap(1))
            write(1,*) 'cut dcut', cut(1,nelemap(1),nelemap(2)), cut(2,nelemap(1),nelemap(2)),  &
                    &  dcut(1,nelemap(1),nelemap(2)), dcut(2,nelemap(1),nelemap(2))
        else
            write(*,*) 'Should not be here. Stopping.'
            stop
        endif
    endif
    ! read further basis settings
    if (doread) then
        read(1,*) str1, str2, nradbaset_atomic, nradbaset_magnetic, nradial2t_magnetic, l2maxt_magnetic
        read(1,*) str1, str2, nradial3t_atomic, nradial3t_magnetic, l3maxt_atomic, l3maxt_magnetic
        read(1,*) str1, str2, nradial4t_atomic, nradial4t_magnetic, l4maxt_atomic, l4maxt_magnetic
        read(1,*) str1, str2, nradial5t_atomic, nradial5t_magnetic, l5maxt_atomic, l5maxt_magnetic
        read(1,*) str1, str2, nradial6t_atomic, nradial6t_magnetic, l6maxt_atomic, l6maxt_magnetic
        read(1,*) str1, ndensityt
        read(1,*) str1, npott, nfembt, radtypet_atomic, radtypet_magnetic     
    else
        write(1,*) 'nradbase l2max',  nradbase_atomic, nradbase_magnetic, nradial2_magnetic, l2max_magnetic
        write(1,*) 'nradial3 l3max ', nradial3_atomic, nradial3_magnetic, l3max_atomic, l3max_magnetic
        write(1,*) 'nradial4 l4max ', nradial4_atomic, nradial4_magnetic, l4max_atomic, l4max_magnetic
        write(1,*) 'nradial5 l5max ', nradial5_atomic, nradial5_magnetic, l5max_atomic, l5max_magnetic
        write(1,*) 'nradial6 l6max ', nradial6_atomic, nradial6_magnetic, l6max_atomic, l6max_magnetic 
        write(1,*) 'ndensity ', ndensity
        write(1,*) 'npot-nfemb-radtype ', npot, nfemb, radtype_atomic, radtype_magnetic
    endif
  
    if ((doread).and.(readheaderonly)) then
        goto 10000
    endif
  
  if (doread) then

     ! read all parameters - read only what we need, so that even simpler models can be read from more complex data
     do
        m = 0
        read(1,*,iostat=m) paratype, norder, nd
        if ( m == -1 ) then
           ! last line
           goto 1000
        endif

        if (trim(paratype) == 'clu' ) then
           ! read cluster contributions
           ! read order 1 contributions
           
           if (norder == 1) then
              read(1,*,iostat=m) elet(1), nr, xclu, fitclu
              if ( m == -1 ) then
                 ! last line
                 goto 1000
              endif
              if (neleread == 1) then
                 if ( nd <= ndensity ) then
                    do nele = 1, nelements
                        do n = 1, nvar1
                            match = .true.
                            ele = element(nele)
                            if (trim(ele) /= trim(elet(1))) match = .false.
                            if (nr /= b1index(n)) match = .false.
                            if (match) then
                               c0(n,nd,nele) = xclu
                               fitc0(n,nd,nele) = fitclu
                            endif
                        enddo
                    enddo
                 endif
              endif
           endif
           
           ! read order 2 contributions
           if (norder == 2) then
              read(1,*,iostat=m) elet(1), elet(2), nr0_p, nr1_p, nr1, l0_p, xclu, fitclu
              if ( m == -1 ) then
                 ! last line
                 goto 1000
              endif
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (.not.check) then
                 print *,'Read incorrect elements for interaction of order: ', norder
                 print *,'Read: ',elet(1:norder)
                 print *,'Expected number of elements: ', neleread
                 print *,'Expected: ', eleread(1:neleread)
                 stop
              endif
              if ( nd <= ndensity ) then
                 do nele = 1, nelements
                    do n = 1, nvar2
                       match = .true.
                       ele = element(nele)
                       if (trim(ele) /= trim(elet(1))) match = .false.
                       ele = element(b2occ(n))
                       if (trim(ele) /= trim(elet(2))) match = .false.
                       if (nr0_p /= b2index(1,n)) match = .false.
                       if (nr1_p /= b2index(2,n)) match = .false.
                       if (nr1 /= b2index(3,n)) match = .false.
                       if (l0_p /= b2index(4,n)) match = .false.
                       if ( match ) then
                          c2(n,nd,nele) = xclu
                          fitc2(n,nd,nele) = fitclu
                       endif
                    enddo
                 enddo
              endif
           endif
           
           ! read order 3 contributions
           if (norder == 3) then
              read(1,*,iostat=m) elet(1), elet(2), elet(3), nr1, nr2, nr0_p, nr1_p, & 
              &     nr2_p, l1, l0_p, l1_p, l2_p, xclu, fitclu
              if ( m == -1 ) then
                 ! last line
                 goto 1000
              endif
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (.not.check) then
                 print *,'Read incorrect elements for interaction of order: ', norder
                 print *,'Read: ',elet(1:norder)
                 print *,'Expected number of elements: ', neleread
                 print *,'Expected: ', eleread(1:neleread)
                 stop
              endif
              if ( nd <= ndensity ) then
                 do nele = 1, nelements
                    do n = 1, nvar3
                       match = .true.
                       ele = element(nele)
                       if (trim(ele) /= trim(elet(1))) match = .false.
                       do nele2 = 1, 2
                          ele = element(b3occ(nele2,n))
                          if (trim(ele) /= trim(elet(nele2+1))) match = .false.
                       enddo
                       if (nr1 /= b3index(1,n)) match = .false.
                       if (nr2 /= b3index(2,n)) match = .false.
                       if (nr0_p /= b3index(3,n)) match = .false.
                       if (nr1_p /= b3index(4,n)) match = .false.
                       if (nr2_p /= b3index(5,n)) match = .false.
                       if (l1    /= b3index(6,n)) match = .false. 
                       if (l0_p  /= b3index(7,n)) match = .false.
                       if (l1_p  /= b3index(8,n)) match = .false.
                       if (l2_p  /= b3index(9,n)) match = .false.
                       
                       if ( match ) then
                          c3(n,nd,nele) = xclu
                          fitc3(n,nd,nele) = fitclu
                       endif
                    enddo
                 enddo
              endif
           endif
           
           ! read order 4 contributions
           if (norder == 4) then
              read(1,*,iostat=m) elet(1), elet(2), elet(3), elet(4), nr1, nr2, nr3, nr0_p, &
              &   nr1_p, nr2_p, nr3_p, l1, l2, l3, l0_p, l1_p, l2_p, l3_p, lint01_p, xclu, fitclu
              if ( m == -1 ) then
                 ! last line
                 goto 1000
              endif
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (.not.check) then
                 print *,'Read incorrect elements for interaction of order: ', norder
                 print *,'Read: ',elet(1:norder)
                 print *,'Expected number of elements: ', neleread
                 print *,'Expected: ', eleread(1:neleread)
                 stop
              endif
              if ( nd <= ndensity ) then
                 do nele = 1, nelements
                    do n = 1, nvar4
                       match = .true.
                       ele = element(nele)
                       if (trim(ele) /= trim(elet(1))) match = .false.
                       do nele2 = 1, 3
                          ele = element(b4occ(nele2,n))
                          if (trim(ele) /= trim(elet(nele2+1))) match = .false.
                       enddo
                       if (nr1 /= b4index(1,n)) match = .false.
                       if (nr2 /= b4index(2,n)) match = .false.
                       if (nr3 /= b4index(3,n)) match = .false. 
                       if (nr0_p /= b4index(4,n)) match = .false.
                       if (nr1_p /= b4index(5,n)) match = .false.
                       if (nr2_p /= b4index(6,n)) match = .false.
                       if (nr3_p /= b4index(7,n)) match = .false.
                       if (l1  /= b4index(8,n)) match = .false.
                       if (l2  /= b4index(9,n)) match = .false.
                       if (l3  /= b4index(10,n)) match = .false. 
                       if (l0_p  /= b4index(11,n)) match = .false.
                       if (l1_p  /= b4index(12,n)) match = .false.
                       if (l2_p  /= b4index(13,n)) match = .false.
                       if (l3_p  /= b4index(14,n)) match = .false.
                       if (lint01_p /= b4index(15,n)) match = .false.
                       if ( match ) then
                          c4(n,nd,nele) = xclu
                          fitc4(n,nd,nele) = fitclu
                       endif
                    enddo
                 enddo
              endif
           endif
           
           ! read order 5 contributions
           if (norder == 5) then
              read(1,*,iostat=m) elet(1), elet(2), elet(3), elet(4), elet(5), nr1, nr2, nr3,    &
              &  nr4, nr0_p, nr1_p, nr2_p, nr3_p, nr4_p, l1, l2, l3, l4, l0_p, l1_p, l2_p,      &
              &  l3_p, l4_p, lint12, lint01_p, lint23_p, xclu, fitclu
              if ( m == -1 ) then
                 ! last line
                 goto 1000
              endif
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (.not.check) then
                 print *,'Read incorrect elements for interaction of order: ', norder
                 print *,'Read: ',elet(1:norder)
                 print *,'Expected number of elements: ', neleread
                 print *,'Expected: ', eleread(1:neleread)
                 stop
              endif
              if ( nd <= ndensity ) then
                 do nele = 1, nelements
                    do n = 1, nvar5
                       match = .true.
                       ele = element(nele)
                       if (trim(ele) /= trim(elet(1))) match = .false.
                       do nele2 = 1, 4
                          ele = element(b5occ(nele2,n))
                          if (trim(ele) /= trim(elet(nele2+1))) match = .false.
                       enddo
                       if (nr1 /= b5index(1,n)) match = .false.
                       if (nr2 /= b5index(2,n)) match = .false.
                       if (nr3 /= b5index(3,n)) match = .false.
                       if (nr4 /= b5index(4,n)) match = .false.
                       if (nr0_p /= b5index(5,n)) match = .false. 
                       if (nr1_p /= b5index(6,n)) match = .false. 
                       if (nr2_p /= b5index(7,n)) match = .false. 
                       if (nr3_p /= b5index(8,n)) match = .false. 
                       if (nr4_p /= b5index(9,n)) match = .false. 
                       if (l1  /= b5index(10,n)) match = .false.
                       if (l2  /= b5index(11,n)) match = .false.
                       if (l3  /= b5index(12,n)) match = .false.
                       if (l4  /= b5index(13,n)) match = .false.
                       if (l0_p  /= b5index(14,n)) match = .false.
                       if (l1_p  /= b5index(15,n)) match = .false.
                       if (l2_p  /= b5index(16,n)) match = .false.
                       if (l3_p  /= b5index(17,n)) match = .false.
                       if (l4_p  /= b5index(18,n)) match = .false.
                       if (lint12  /= b5index(19,n)) match = .false.
                       if (lint01_p  /= b5index(20,n)) match = .false.
                       if (lint23_p  /= b5index(21,n)) match = .false.
                       if ( match ) then
                          c5(n,nd,nele) = xclu
                          fitc5(n,nd,nele) = fitclu
                       endif
                    enddo
                 enddo
              endif
           endif
           
           ! read order 6 contributions
           if (norder == 6) then
              read(1,*,iostat=m) elet(1), elet(2), elet(3), elet(4), elet(5), elet(6), &
              &  nr1, nr2, nr3, nr4, nr5, nr0_p, nr1_p, nr2_p, nr3_p, nr4_p, nr5_p, l1, l2,  &
              &  l3, l4, l5, l0_p, l1_p, l2_p, l3_p, l4_p, l5_p, lint12, lint34,       &
              &  lint01_p, lint23_p, lint45_p, xclu, fitclu
              if ( m == -1 ) then
                 ! last line
                 goto 1000
              endif
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (.not.check) then
                 print *,'Read incorrect elements for interaction of order: ', norder
                 print *,'Read: ',elet(1:norder)
                 print *,'Expected number of elements: ', neleread
                 print *,'Expected: ', eleread(1:neleread)
                 stop
              endif
              if ( nd <= ndensity ) then
                 do nele = 1, nelements
                    do n = 1, nvar6
                       match = .true.
                       ele = element(nele)
                       if (trim(ele) /= trim(elet(1))) match = .false.
                       do nele2 = 1, 5
                          ele = element(b6occ(nele2,n))
                          if (trim(ele) /= trim(elet(nele2+1))) match = .false.
                       enddo
                       if (nr1 /= b6index(1,n)) match = .false.
                       if (nr2 /= b6index(2,n)) match = .false.
                       if (nr3 /= b6index(3,n)) match = .false.
                       if (nr4 /= b6index(4,n)) match = .false.
                       if (nr5 /= b6index(5,n)) match = .false.
                       if (nr0_p /= b6index(6,n)) match = .false.
                       if (nr1_p /= b6index(7,n)) match = .false.
                       if (nr2_p /= b6index(8,n)) match = .false.
                       if (nr3_p /= b6index(9,n)) match = .false.
                       if (nr4_p /= b6index(10,n)) match = .false.
                       if (nr5_p /= b6index(11,n)) match = .false.
                       if (l1  /= b6index(12,n)) match = .false.
                       if (l2  /= b6index(13,n)) match = .false.
                       if (l3  /= b6index(14,n)) match = .false.
                       if (l4  /= b6index(15,n)) match = .false.
                       if (l5  /= b6index(16,n)) match = .false.
                       if (l0_p  /= b6index(17,n)) match = .false.
                       if (l1_p  /= b6index(18,n)) match = .false.
                       if (l2_p  /= b6index(19,n)) match = .false.
                       if (l3_p  /= b6index(20,n)) match = .false.
                       if (l4_p  /= b6index(21,n)) match = .false.
                       if (l5_p  /= b6index(22,n)) match = .false.
                       if (lint12  /= b6index(23,n)) match = .false.
                       if (lint34  /= b6index(24,n)) match = .false.
                       if (lint01_p  /= b6index(25,n)) match = .false.
                       if (lint23_p  /= b6index(26,n)) match = .false.
                       if (lint45_p  /= b6index(27,n)) match = .false.
                       if ( match ) then
                          c6(n,nd,nele) = xclu
                          fitc6(n,nd,nele) = fitclu
                       endif
                    enddo
                 enddo
              endif
           endif
           
        elseif (trim(paratype) == 'par') then
           if (norder < 1) then
              print *,'norder > 1 required for par.'
              stop
           elseif (norder > nvarp1) then
              print *,'norder < nvarp1 required for par.'
              stop
           endif
           ! read parameters in functional
           read(1,*,iostat=m) elet(1), xpara, fitpara
           if ( m == -1 ) then
              ! last line
              goto 1000
           endif
           if (neleread == 1) then
              check = .true.
              if (elet(1) /= eleread(1)) check = .false.
              if (check) then
                 do nele1 = 1, nelements
                    if (element(nele1) == elet(1)) then
                       if (( nd <= ndensity ).and.( norder <= nvarp1 )) then
                          p1(norder,nd,nele1) = xpara
                          fitp1(norder,nd,nele1) = fitpara
                       endif
                    endif
                 enddo
              endif
           endif

           
           
        elseif (trim(paratype) == 'rad_atomic') then
           !           ! rad parameters in radial function
           read(1,*,iostat=m) elet(1), elet(2), nr1, nr2, l1, xclu, fitclu
           if ( m == -1 ) then
              ! last line
              goto 1000
           endif
           if (neleread == 1) then
              check = .true.
              if (elet(1) /= eleread(1)) check = .false.
              if (elet(2) /= eleread(1)) check = .false.
           elseif (neleread == 2) then
              check = .false.
              if ((elet(1) == eleread(1)).and.(elet(2) == eleread(2))) check = .true.
              if ((elet(2) == eleread(1)).and.(elet(1) == eleread(2))) check = .true.
           elseif (neleread > 2) then
              print *,'Found radial parameter in multi-body parameter file. Stopping.'
              stop
           endif
           if (.not.check) then
              print *,'Error when reading radial parameters.'
              stop
           endif
           do nele1 = 1, nelements
              do nele2 = 1, nelements
                 if ((element(nele1) == elet(1)).and.(element(nele2) == elet(2))) then
                    if ((nr1>=1).and.(nr1<=nradbase_atomic).and.(nr2>=1)                  &
            &        .and.(nr2<=nradial3_atomic).and.(l1>=0).and.(l1<=l3max_atomic)) then
                       crad_atomic(nr1,nr2,l1,nele1,nele2) = xclu
                       fitcrad_atomic(nr1,nr2,l1,nele1,nele2) = fitclu
                       if (.not.dotrainrad) then
                          fitcrad_atomic(nr1,nr2,l1,nele1,nele2) = .false.
                       endif
                    endif
                 endif
              enddo
           enddo        
           
        elseif (trim(paratype) == 'rad_magnetic') then
           !           ! rad parameters in radial function
           read(1,*,iostat=m) elet(1), elet(2), nr0_p, nr1_p, l0_p, xclu_p, fitclu_p
           if ( m == -1 ) then
              ! last line
              goto 1000
           endif
           if (neleread == 1) then
              check = .true.
              if (elet(1) /= eleread(1)) check = .false.
              if (elet(2) /= eleread(1)) check = .false.
           elseif (neleread == 2) then
              check = .false.
              if ((elet(1) == eleread(1)).and.(elet(2) == eleread(2))) check = .true.
              if ((elet(2) == eleread(1)).and.(elet(1) == eleread(2))) check = .true.
           elseif (neleread > 2) then
              print *,'Found radial parameter in multi-body parameter file. Stopping.'
              stop
           endif
           if (.not.check) then
              print *,'Error when reading radial parameters.'
              stop
           endif
           do nele1 = 1, nelements
              do nele2 = 1, nelements
                  if ((element(nele1) == elet(1)).and.(element(nele2) == elet(2))) then
                    if ((nr0_p>=1).and.(nr0_p<=nradbase_magnetic).and.(nr1_p>=1).and.   &
                &        (nr1_p<=nradial2_magnetic).and.(l0_p>=0).and.(l0_p<=l2max_magnetic)) then
                        crad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele2) = xclu_p
                        fitcrad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele2) = fitclu_p
                        if (.not.dotrainrad) then
                            fitcrad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele2) = .false.
                        endif
                    endif
                  endif
               enddo
            enddo
        elseif (trim(paratype) == 'core') then
           if (norder < 1) then
              print *,'norder > 0 required for core.'
              stop
           endif
           if (norder > 2) then
              print *,'norder < 3 required for core.'
              stop
           endif
           !           ! rad parameters for hard core repulsion
           read(1,*,iostat=m) elet(1), elet(2), xpara, fitpara
           if ( m == -1 ) then
              ! last line
              goto 1000
           endif
           if (neleread == 1) then
              check = .true.
              if (elet(1) /= eleread(1)) check = .false.
              if (elet(2) /= eleread(1)) check = .false.
           elseif (neleread == 2) then
              check = .false.
              if ((elet(1) == eleread(1)).and.(elet(2) == eleread(2))) check = .true.
              if ((elet(2) == eleread(1)).and.(elet(1) == eleread(2))) check = .true.
           elseif (neleread > 2) then
              print *,'Found hard core repulsion parameters in multi-body parameter file. Stopping.'
              stop
           endif
           if (.not.check) then
              print *,'Error when reading radial parameters.'
              stop
           endif
           do nele1 = 1, nelements
              do nele2 = 1, nelements
                 if ((element(nele1) == elet(1)).and.(element(nele2) == elet(2))) then
                    corepara(norder,nele1,nele2) = xpara
                    fitcorepara(norder,nele1,nele2) = fitpara
                 endif
              enddo
           enddo
        else
        endif
        
    enddo
     
1000 continue
     
  else
     ! write parameters
     if (nelements == 1) then
        nele = 1
        do nd = 1, ndensity
           do norder = 1, nvarp1
              write(1,*) 'par ', norder, nd
              write(1,*) element(nele), p1(norder,nd,nele), fitp1(norder,nd,nele)
           enddo
        enddo
     endif

     ! write parameters for core repulsion (two parameters for each pair of atoms)
     nd = 0
     np = 0
     if (neleread == 1) then
        do nele1 = 1, nelements
           if (eleread(1) == element(nele1)) then
              norder = 1
              write(1,*) 'core ', norder, nd
              write(1,*) eleread(1),' ', eleread(1), dabs(corepara(norder,nele1,nele1)), fitcorepara(norder,nele1,nele1)
              norder = 2
              write(1,*) 'core ', norder, nd
              write(1,*) eleread(1),' ', eleread(1), dabs(corepara(norder,nele1,nele1)), fitcorepara(norder,nele1,nele1)              
           endif
        enddo
     endif
     if (neleread == 2) then
        do nele1 = 1, nelements
           do nele2 = 1, nelements          
              if ((eleread(1) == element(nele1)).and.(eleread(2) == element(nele2))) then
                 norder = 1
                 write(1,*) 'core ', norder, nd
                 write(1,*) eleread(1),' ', eleread(2), dabs(corepara(norder,nele1,nele2)), fitcorepara(norder,nele1,nele2)
                 norder = 2
                 write(1,*) 'core ', norder, nd
                 write(1,*) eleread(1),' ', eleread(2), dabs(corepara(norder,nele1,nele2)), fitcorepara(norder,nele1,nele2)  
                 norder = 1
                 write(1,*) 'core ', norder, nd
                 write(1,*) eleread(2),' ', eleread(1), dabs(corepara(norder,nele2,nele1)), fitcorepara(norder,nele2,nele1)
                 norder = 2
                 write(1,*) 'core ', norder, nd
                 write(1,*) eleread(2),' ', eleread(1), dabs(corepara(norder,nele2,nele1)), fitcorepara(norder,nele2,nele1) 
              endif
           enddo
        enddo
     endif
     
     ! write parameters in radial function
     norder = 0
     nd = 0
     np = 0
     if (neleread == 1) then
        do nele1 = 1, nelements
           if (eleread(1) == element(nele1)) then
              do nr1 = 1, nradbase_atomic
                 do nr2 = 1, nradial3_atomic
                    do l1 = 0, l3max_atomic
                       write(1,*) 'rad_atomic ', norder, nd
                       write(1,*) eleread(1),' ', eleread(1), nr1, nr2, l1,   &
                            &    crad_atomic(nr1,nr2,l1,nele1,nele1),         &
                            &    fitcrad_atomic(nr1,nr2,l1,nele1,nele1)
                    enddo
                 enddo
              enddo
           endif
        enddo
        do nele1 = 1, nelements
           if (eleread(1) == element(nele1)) then
              do nr0_p = 1, nradbase_magnetic
                 do nr1_p = 1, nradial2_magnetic
                    do l0_p = 0, l2max_magnetic
                       write(1,*) 'rad_magnetic ', norder, nd
                       write(1,*) eleread(1),' ', eleread(1), nr0_p, nr1_p, l0_p, &
                            &    crad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele1),     &
                            &    fitcrad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele1)
                    enddo
                 enddo
              enddo
           endif
        enddo
     endif
     if (neleread == 2) then
        do nele1 = 1, nelements
           do nele2 = 1, nelements          
              if ((eleread(1) == element(nele1)).and.(eleread(2) == element(nele2))) then
                 do nr1 = 1, nradbase_atomic
                    do nr2 = 1, nradial3_atomic
                       do l1 = 0, l3max_atomic
                          write(1,*) 'rad_atomic ', norder, nd
                          write(1,*) eleread(1), eleread(2), nr1, nr2, l1, &
                               &  crad_atomic(nr1,nr2,l1,nele1,nele2),     &
                               &  fitcrad_atomic(nr1,nr2,l1,nele1,nele2)
                       enddo
                    enddo
                 enddo
                 do nr1 = 1, nradbase_atomic
                    do nr2 = 1, nradial3_atomic
                       do l1 = 0, l3max_atomic
                          write(1,*) 'rad_atomic ', norder, nd
                          write(1,*) eleread(2), ' ',eleread(1), nr1, nr2, l1, &
                               &   crad_atomic(nr1,nr2,l1,nele2,nele1),        &
                               &   fitcrad_atomic(nr1,nr2,l1,nele2,nele1)
                       enddo
                    enddo
                 enddo
              endif
           enddo
        enddo
        do nele1 = 1, nelements
           do nele2 = 1, nelements          
              if ((eleread(1) == element(nele1)).and.(eleread(2) == element(nele2))) then
                 do nr0_p = 1, nradbase_magnetic
                    do nr1_p = 1, nradial2_magnetic
                       do l0_p = 0, l2max_magnetic
                          write(1,*) 'rad_magnetic ', norder, nd
                          write(1,*) eleread(1), eleread(2), nr0_p, nr1_p, l0_p, &
                               &  crad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele2),   &
                               &  fitcrad_magnetic(nr0_p,nr1_p,l0_p,nele1,nele2)
                       enddo
                    enddo
                 enddo
                 do nr0_p = 1, nradbase_magnetic
                    do nr1_p = 1, nradial2_magnetic
                       do l0_p = 0, l2max_magnetic
                          write(1,*) 'rad_magnetic ', norder, nd
                          write(1,*) eleread(2), ' ',eleread(1), nr0_p, nr1_p, l0_p, &
                               &  crad_magnetic(nr0_p,nr1_p,l0_p,nele2,nele1),       &
                               &  fitcrad_magnetic(nr0_p,nr1_p,l0_p,nele2,nele1)
                       enddo
                    enddo
                 enddo
              endif
           enddo
        enddo
     endif
     
     ! write clusters

     ! write order 1 contributions
     if (nelements == 1 ) then
        norder = 1
        nele = 1
        !nclu = 1
        elet(1) = element(nele)
        do nd = 1, ndensity
            do nclu = 1, nvar1
                write(1,*) 'clu ', norder, nd
                write(1,*) elet(1), b1index(nclu), c0(nclu,nd,nele), fitc0(nclu,nd,nele)
            enddo
        enddo
     endif
     
     ! write order 2 contributions
     norder = 2
     do nd = 1, ndensity
        do nele = 1, nelements
           elet(1) = element(nele)
           do nclu = 1, nvar2
              elet(2) = element(b2occ(nvar2))
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (check) then
                 write(1,*) 'clu ', norder, nd
                 write(1,*) elet(1), ' ',elet(2), b2index(1,nclu), b2index(2,nclu), b2index(3,nclu), b2index(4,nclu),  &
                 &    c2(nclu,nd,nele), fitc2(nclu,nd,nele)
              endif
           enddo
        enddo
     enddo
     ! write order 3 contributions
     norder = 3
     do nd = 1, ndensity
        do nele = 1, nelements
           elet(1) = element(nele)
           do nclu = 1, nvar3
              elet(2:3) = element(b3occ(1:2,nvar3))
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (check) then
                 write(1,*) 'clu ', norder, nd
                 write(1,*) elet(1),' ',elet(2),' ',elet(3),               &
                      & b3index(1,nclu), b3index(2,nclu), b3index(3,nclu), &
                      & b3index(4,nclu), b3index(5,nclu), b3index(6,nclu), &
                      & b3index(7,nclu), b3index(8,nclu), b3index(9,nclu), &
                      & c3(nclu,nd,nele), fitc3(nclu,nd,nele)
              endif
           enddo
        enddo
     enddo
     ! write order 4 contributions
     norder = 4
     do nd = 1, ndensity
        do nele = 1, nelements
           elet(1) = element(nele)
           do nclu = 1, nvar4
              elet(2:4) = element(b4occ(1:3,nvar4))
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (check) then
                 write(1,*) 'clu ', norder, nd
                 write(1,*) elet(1),' ',elet(2),' ',elet(3),' ',elet(4), &
                      & b4index(1,nclu), b4index(2,nclu), b4index(3,nclu), b4index(4,nclu), &
                      & b4index(5,nclu), b4index(6,nclu), b4index(7,nclu), b4index(8,nclu), &
                      & b4index(9,nclu), b4index(10,nclu), b4index(11,nclu), b4index(12,nclu), &
                      & b4index(13,nclu), b4index(14,nclu), b4index(15,nclu),                  &
                      & c4(nclu,nd,nele), fitc4(nclu,nd,nele)
              endif
           enddo
        enddo
     enddo
     ! write order 5 contributions
     norder = 5
     do nd = 1, ndensity
        do nele = 1, nelements
           elet(1) = element(nele)
           do nclu = 1, nvar5
              elet(2:5) = element(b5occ(1:4,nvar5))
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (check) then
                 write(1,*) 'clu ', norder, nd
                 write(1,*)  elet(1),' ',elet(2),' ',elet(3),' ',elet(4),' ',elet(5), &
                      & b5index(1,nclu), b5index(2,nclu), b5index(3,nclu), b5index(4,nclu), &
                      & b5index(5,nclu), b5index(6,nclu), b5index(7,nclu), b5index(8,nclu), &
                      & b5index(9,nclu), b5index(10,nclu), b5index(11,nclu), b5index(12,nclu), &
                      & b5index(13,nclu), b5index(14,nclu), b5index(15,nclu), b5index(16,nclu), &
                      & b5index(17,nclu), b5index(18,nclu), b5index(19,nclu), b5index(20,nclu), &
                      & b5index(21,nclu), c5(nclu,nd,nele), fitc5(nclu,nd,nele)
              endif
           enddo
        enddo
     enddo
     ! write order 6 contributions
     norder = 6
     do nd = 1, ndensity
        do nele = 1, nelements
           elet(1) = element(nele)
           do nclu = 1, nvar6
              elet(2:6) = element(b6occ(1:5,nvar6))
              call elecheck(neleread,norder,nelements,nelet,elet,eleread,check)
              if (check) then
                 write(1,*) 'clu ', norder, nd
                 write(1,*) elet(1),' ',elet(2),' ',elet(3),' ',elet(4),' ',elet(5),' ',elet(6), &
                      & b6index(1,nclu), b6index(2,nclu), b6index(3,nclu), b6index(4,nclu), &
                      & b6index(5,nclu), b6index(6,nclu), b6index(7,nclu), b6index(8,nclu), &
                      & b6index(9,nclu), b6index(10,nclu),b6index(11,nclu),b6index(12,nclu), & 
                      & b6index(13,nclu), b6index(14,nclu),b6index(15,nclu),b6index(16,nclu), &
                      & b6index(17,nclu), b6index(18,nclu),b6index(19,nclu),b6index(20,nclu), &
                      & b6index(21,nclu), b6index(22,nclu),b6index(23,nclu),b6index(24,nclu), &
                      & b6index(25,nclu), b6index(26,nclu),b6index(27,nclu),                  &
                      & c6(nclu,nd,nele), fitc6(nclu,nd,nele)
              endif
           enddo
        enddo
     enddo
  endif

10000 continue
  
  close(1)
  
  if (verblevel > 0) then
     if (doread) then
        write(*,'(A)') 'Finished reading parameter file.'
     endif
  endif
  
     
end subroutine rwpara

!----------------------------------------------------------------------------------------------------------

subroutine elecheck(neleread,norder,nelements,nelet,elet,eleread,check)

    integer neleread, norder, nelements, nelet
    character(2) elet(1:nelet), eleread(1:nelements)
    logical check
    integer n, m, nele
    logical got(1:neleread)

    check = .true.

    got(1:neleread) = .false.
    nele = 0
    do n = 1, norder
        do m = 1, neleread
            if ( trim(elet(n)) == trim(eleread(m)) ) then
                got(m) = .true.
                nele = nele + 1
            endif
        enddo
    enddo

    if ( nele /= norder ) check = .false.
    do m = 1, neleread
        if (.not.got(m)) check = .false.
    enddo

end subroutine elecheck

!----------------------------------------------------------------------------------------------------------

subroutine allocateforprepAndBandforce
    use global
    use functionparameters
    use modprepA
    use modneigh
    implicit none
    integer ij

    if (n2atom > natomsmax) then
        print *,'Incorrect array dimension specifications.'
        stop
    endif
  
    if (allocated(abasepack)) then
        if (n2atom > size(abasepack,6) ) then
            deallocate(abasepack)
            deallocate(abase2pack)
            deallocate(abase2pack_magnetic)
            deallocate(abase2hcpack)    
        endif
    endif
    if (.not. allocated(abasepack)) then
        allocate(abasepack(1:nradial3_atomic,1:lmaxsq_atomic,  &
                &   1:nradial3_magnetic,1:lmaxsq_magnetic,1:nelements,1:natomsmax))
        allocate(abase2pack(1:nradbase_atomic,1:nradial2_magnetic, &
                &   1:lmaxsq_magnetic,1:nelements,1:natomsmax))
        allocate(abase2pack_magnetic(1:nradial2_magnetic,1:lmaxsq_magnetic, 1:natomsmax))
        allocate(abase2hcpack(1:natomsmax))
    endif
  
    if ((allocated(frpack_atomic)).and.(allocated(frpack_magnetic))) then
        ! last bond
        ij = n2stop(n2atom)
        if (ij > size(frpack_atomic,3) ) then
            deallocate(gradshpack_atomic)
            deallocate(gradshpack_magnetic)
            deallocate(gradshpack_magnetic_0)
            deallocate(gradfrpack_atomic)
            deallocate(gradfrpack_magnetic)
            deallocate(gradfrpack_magnetic_0)
            deallocate(gradfr2pack_atomic)
            deallocate(gradfr2hcpack)
            deallocate(shpack_atomic)
            deallocate(shpack_magnetic)
            deallocate(shpack_magnetic_0)
            deallocate(frpack_atomic)
            deallocate(frpack_magnetic)
            deallocate(grpack_magnetic_0)
            deallocate(gradgrpack_magnetic_0)
            deallocate(frpack_magnetic_0)
            deallocate(fr2pack_atomic)
            deallocate(fr2hcpack)
            if (ij > nei2countmax) then
                print *,'Incorrect array dimension specifications.'
                stop
            endif
        endif
    endif
    if (.not. allocated(frpack_atomic)) then    
        allocate(gradshpack_atomic(1:3,1:lmaxsq_atomic,1:nei2countmax))
        allocate(gradshpack_magnetic(1:3,1:lmaxsq_magnetic,1:nei2countmax))
        allocate(gradshpack_magnetic_0(1:3,1:lmaxsq_magnetic,1:natomsmax))
        allocate(gradfrpack_atomic(1:nradial3_atomic,0:l3max_atomic,1:nei2countmax))
        allocate(gradfrpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,1:nei2countmax))
        allocate(gradfrpack_magnetic_0(1:nradial2_magnetic,0:l2max_magnetic,1:natomsmax))
        allocate(gradfr2pack_atomic(1:nradbase_atomic,1:nei2countmax))
        allocate(gradfr2hcpack(1:nei2countmax))
        allocate(shpack_atomic(1:lmaxsq_atomic,1:nei2countmax))
        allocate(shpack_magnetic(1:lmaxsq_magnetic,1:nei2countmax))
        allocate(shpack_magnetic_0(1:lmaxsq_magnetic,1:natomsmax))
        allocate(frpack_atomic(1:nradial3_atomic,0:l3max_atomic,1:nei2countmax))
        allocate(frpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,1:nei2countmax))
        allocate(grpack_magnetic_0(1:nradbase_magnetic,1:natomsmax))
        allocate(gradgrpack_magnetic_0(1:nradbase_magnetic,1:natomsmax))
        allocate(frpack_magnetic_0(1:nradial2_magnetic,0:l2max_magnetic,1:natomsmax))
        allocate(fr2pack_atomic(1:nradbase_atomic,1:nei2countmax))
        allocate(fr2hcpack(1:nei2countmax))
        
    endif

    if((allocated(B1pack)).and.(n2atom > size(B1pack,2))) then
        deallocate(B1pack)
    endif
  
    if(.not.allocated(B1pack)) then
        allocate(B1pack(1:nvar1,1:natomsmax))
    endif
  
    if (nradial2_magnetic > 0) then
        if((allocated(B2pack)).and.(n2atom > size(B2pack,2))) then
            deallocate(B2pack)
            deallocate(B2hcpack)
            deallocate(db2pack)
            deallocate(db2pack_0)
        endif    
        if(.not.allocated(B2pack)) then
            allocate(B2pack(1:nvar2,1:natomsmax))
            allocate(B2hcpack(1:natomsmax))
            allocate(db2pack(-l2max_magnetic:l2max_magnetic,1:nvar2,1:natomsmax))
            allocate(db2pack_0(-l2max_magnetic:l2max_magnetic,1:nvar2,1:natomsmax))
        endif
    endif
  
    if ((nradial3_atomic > 0).or.(nradial3_magnetic > 0)) then
        if((allocated(B3pack)).and.(n2atom > size(B3pack,2))) then
            deallocate(B3pack)
            deallocate(db3pack)
            deallocate(db3pack_0)
        endif
        if (.not.allocated(B3pack)) then
            allocate(B3pack(1:nvar3,1:natomsmax))
            allocate(db3pack(1:2,-l3max_atomic:l3max_atomic,-l3max_magnetic:l3max_magnetic,1:nvar3,1:natomsmax))
            allocate(db3pack_0(-l3max_magnetic:l3max_magnetic,1:nvar3,1:natomsmax))
        endif
    endif
  
    if ((nradial4_atomic > 0).or.(nradial4_magnetic > 0)) then    
        if((allocated(B4pack)).and.(n2atom > size(B4pack,2))) then
            deallocate(B4pack)
            deallocate(db4pack)
            deallocate(db4pack_0)
        endif
        if(.not.allocated(B4pack)) then
            allocate(B4pack(1:nvar4,1:natomsmax))
            allocate(db4pack(1:3,-l4max_atomic:l4max_atomic,-l4max_magnetic:l4max_magnetic,1:nvar4,1:natomsmax))
            allocate(db4pack_0(-l4max_magnetic:l4max_magnetic,1:nvar4,1:natomsmax))
        endif
    endif
    if ((nradial5_atomic > 0).or.(nradial5_magnetic > 0)) then    
        if((allocated(B5pack)).and.(n2atom > size(B5pack,2))) then
            deallocate(B5pack)
            deallocate(db5pack)
            deallocate(db5pack_0)
        endif
        if(.not.allocated(B5pack)) then
            allocate(B5pack(1:nvar5,1:natomsmax))
            allocate(db5pack(1:4,-l5max_atomic:l5max_atomic,-l5max_magnetic:l5max_magnetic,1:nvar5,1:natomsmax))
            allocate(db5pack_0(-l5max_magnetic:l5max_magnetic,1:nvar5,1:natomsmax))
        endif
    endif
    if ((nradial6_atomic > 0).or.(nradial6_magnetic > 0)) then    
        if((allocated(B6pack)).and.(n2atom > size(B6pack,2))) then
            deallocate(B6pack)
            deallocate(db6pack)
            deallocate(db6pack_0)
        endif
        if(.not.allocated(B6pack)) then
            allocate(B6pack(1:nvar6,1:natomsmax))
            allocate(db6pack(1:5,-l6max_atomic:l6max_atomic,-l6max_magnetic:l6max_magnetic,1:nvar6,1:natomsmax))
            allocate(db6pack_0(-l6max_magnetic:l6max_magnetic,1:nvar6,1:natomsmax))
        endif
    endif
  
    if ((allocated(eicalc)).and.(n2atom>size(eicalc,1))) then
        deallocate(fcalcvec)
    endif
    if (.not. allocated(eicalc)) then
        allocate( eicalc(1:natomsmax))
    endif
  
    if (doforce) then
        if ((allocated(fcalcvec)).and.(n2atom>size(fcalcvec,2))) then
            deallocate(fcalcvec)
        endif
        if (.not. allocated(fcalcvec)) then
            allocate( fcalcvec(1:3,1:natomsmax))
        endif     
        if ((allocated(fijcalc)).and.(nei2countmax>size(fijcalc,2))) then
            deallocate(fijcalc)
        endif
        if (.not. allocated(fijcalc)) then
            allocate( fijcalc(1:3,1:nei2countmax))
        endif
        
    endif
  
    if (dotorque) then 
        if ((allocated(tcalcvec)).and.(n2atom>size(tcalcvec,2))) then
            deallocate(tcalcvec)
        endif
        if (.not. allocated(tcalcvec)) then
            allocate( tcalcvec(1:3,1:natomsmax))
        endif 
        
        if ((allocated(tijcalc0)).and.(n2atom>size(tijcalc0,2))) then
            deallocate(tijcalc0)
        endif
        if (.not. allocated(tijcalc0)) then
            allocate( tijcalc0(1:3,1:natomsmax))
        endif
        
        if ((allocated(tijcalc)).and.(nei2countmax>size(tijcalc,2))) then
            deallocate(tijcalc)
        endif
        if (.not. allocated(tijcalc)) then
            allocate( tijcalc(1:3,1:nei2countmax))
        endif
            
    endif
  
end subroutine allocateforprepAndBandforce

!----------------------------------------------------------------------------------------------------------

subroutine setuptables

    use core
    use global
    use functionparameters
    use tables
    use modneigh
#ifdef MPIPARALLEL
    use mpi
    use mpiglobal
#endif
    implicit none
    double precision, parameter :: pi = 3.14159265358979323846264d0
    integer n, m, lmm, lpm, k, nr, nrad
    integer na, na1, na2, na3
    integer l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, ltest
    integer j1, j2, j3, j4, j5
    logical index, second
    integer(kind=8) nfac
    double precision xfac
    logical sym12, sym23, sym34, sym45
    logical sym123, sym234, sym345
    logical sym123a45, sym12a345, sym12a45, sym12a34, sym23a45
    logical sym1234, sym2345, sym12345
    integer h(1:5)
    integer i1, i2, i3, i4, i5
    integer m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11
    integer maddmax_atomic, maddmax_magnetic, maddmax
    integer mlint1, mlint2, mlint3, mlint1234
    integer nbase, nele, nele1, nele2, nele3, nele4, nele5, nele6, nele7
    integer nc1, nc2, nc3, nc4, nc5, nc6, nc7
    integer kmax_atomic, kmax_magnetic
    integer k_atomic1, k_magnetic1, k_atomic2, k_magnetic2, k_atomic3, k_magnetic3
    integer k_atomic4, k_magnetic4, k_atomic5, k_magnetic5, k_atomic6, k_magnetic6
    integer L, l1r, l2r, l3r, l4r, l5r, l6r, l7r, l8r, l9r, l10r, l11r
    integer nr1, nr2, nr3, nr4, nr5, nr6, nr7, nr8, nr9, nr10, nr11, nr12
    integer nr1r, nr2r, nr3r, nr4r, nr5r, nr6r, nr7r, nr8r, nr9r, nr10r, nr11r, nr12r
    integer nele1r, nele2r, nele3r, nele4r, nele5r, nele6r
    integer lint, lint1, lint2, lint3, lint4, lint5, lintr, lint1r, lint2r, lint3r, lint4r, lint5r
    integer mint1, mint2, mint3, mint4, mint5, mint6, mint1234, mint0123
    integer nvar5tmp, ndeg, nuni, lsum1, lsum2, numlintr, nint
    integer n2v(1:2), n3v(1:3)
    integer nc1_atomic, nc2_atomic, nc3_atomic, nc4_atomic, nc5_atomic, nc6_atomic, nc7_atomic
    integer nc1_magnetic, nc2_magnetic, nc3_magnetic, nc4_magnetic, nc5_magnetic, nc6_magnetic, nc7_magnetic

    integer, allocatable :: buni(:)
    integer, allocatable :: b5occtmp(:,:), b5indextmp(:,:), b5ninttmp(:)
    logical readsuccess

    integer nmax
    integer nrank
    integer numlint, numlint_atomic, numlint_magnetic, nveclint
    integer, allocatable :: vecn(:), vecn_atomic(:), vecn_magnetic(:), vecnr(:), vecnr_atomic(:), vecnr_magnetic(:)
    integer, allocatable :: vecl(:), vecl_atomic(:), vecl_magnetic(:), arrlint(:,:), arrlint_atomic(:,:), arrlint_magnetic(:,:)
    logical isordered_atomic, isordered_magnetic, isbasis_atomic, isbasis_magnetic, isordered, isbasis

    character(8) atomic, magnetic

    atomic = 'atomic'
    magnetic = 'magnetic'

    ! precomputes tables

    ! angular bounds for atomic spherical harmonics
    lmaxsq_atomic = (l3max_atomic+1)**2
    k = 0
    do l = 0, l3max_atomic
        do m = 0, l
            k = k + 1
        enddo
    enddo
    lmaxhalf_atomic= k
  
    ! angular bounds for magnetic spherical harmonics
    lmaxsq_magnetic = (l2max_magnetic+1)**2
    k = 0
    do l = 0, l2max_magnetic
        do m = 0, l
            k = k + 1
        enddo
    enddo
    lmaxhalf_magnetic= k

    ! global angular bounds
    lmaxsq = max(lmaxsq_atomic,lmaxsq_magnetic)
    lmaxhalf = max(lmaxhalf_atomic,lmaxhalf_magnetic)
  
    !-------------------------------------------------------------------------
  
    !precompute required factorials
    allocate(fac(0:8*lmax+1))

    xfac = 1.d0
    fac(0) = xfac
    do l = 0, 8*lmax
        xfac = dble(l+1)*xfac
        fac(l+1) = xfac
        if ( fac(l+1) < 0.d0 ) then
            print *,'Stopping. Double precision overflow for fac calculation.'
            stop
        endif
    enddo  

    !shindex for atomic part 
    allocate(shindex_atomic(0:l3max_atomic,-l3max_atomic:l3max_atomic))
    shindex_atomic = 0
    allocate(shback_atomic(1:lmaxsq_atomic))
    shback_atomic = 0
  
    k = 1
    l = 0
    m = 0
    shindex_atomic(l,m) = k
    shback_atomic(k) = l
    ! pack m >= 0 first
    if (l3max_atomic >= 1) then
        do l = 1, l3max_atomic
            do m = 0, l
                k = k + 1
                shindex_atomic(l,m) = k
                shback_atomic(k) = l
                !        shbackz(k) = m
            enddo
        enddo
        do l = 1, l3max_atomic
            do m = -l, -1
                k = k + 1
                shindex_atomic(l,m) = k
                shback_atomic(k) = l
                !        shbackz(k) = m
            enddo
        enddo
    endif

    !shindex for magnetic part
    allocate(shindex_magnetic(0:l2max_magnetic,-l2max_magnetic:l2max_magnetic))
    shindex_magnetic = 0
    allocate(shback_magnetic(1:lmaxsq_magnetic))
    shback_magnetic = 0
  
    k = 1
    l = 0
    m = 0
    shindex_magnetic(l,m) = k
    shback_magnetic(k) = l
    ! pack m >= 0 first
    if (l2max_magnetic >= 1) then
        do l = 1, l2max_magnetic
            do m = 0, l
                k = k + 1
                shindex_magnetic(l,m) = k
                shback_magnetic(k) = l
                !        shbackz(k) = m
            enddo
        enddo
        do l = 1, l2max_magnetic
            do m = -l, -1
                k = k + 1
                shindex_magnetic(l,m) = k
                shback_magnetic(k) = l
                !        shbackz(k) = m
            enddo
        enddo
    endif
  
    ! set up table of Clebsch-Gordan coefficients and Wigner3j symbols
    call setupcg
            
    ! number of variables for radial functions
    if (verblevel > 0 ) then
        write(*,*) 'number of variables for atomic radial functions:', &
        &  nradbase_atomic*nradial3_atomic*(l3max_atomic+1)*nelements*nelements
    endif
    if (verblevel > 0 ) then
        write(*,*) 'number of variables for magnetic radial functions:', &
        &  nradbase_magnetic*nradial2_magnetic*(l2max_magnetic+1)*nelements*nelements
    endif

    !-------------------------------------------------------------
    ! read bace.in if available
    !-------------------------------------------------------------

    ! !!reading is done by all MPI processes simultaneously
    call readacebase(readsuccess)
    
    ! the purely magnetic one-body contribution is ignored
    
    !-------------------------------------------------------------
    ! number of variables for c1
    !-------------------------------------------------------------
    if (.not.readsuccess) then
  
        second = .false.
1050 continue
        if(second) then
            allocate( b1occ(1:nvar1) )
            allocate( b1index(1:nvar1) )
        endif
  
        nvar1 = 0
        do nele = 1, nelements
            do nr = 1, nradbase_magnetic       
                            
                nvar1 = nvar1 + 1
                if (second) then
                    b1occ(nvar1) = nele    ! mu_i
                    b1index(nvar1) = nr    ! n0'
                endif
                        
            enddo
        enddo
        if (.not.second) then
            second = .true.
            goto 1050
        endif 
  
    endif
    !  if (verblevel > 0 ) then
    write(*,*) 'number of variables for 0-body clusters:', nvar1
    !  endif  
  
    !-------------------------------------------------------------
    ! number of variables for c2
    !-------------------------------------------------------------

    second = .false.
2100 continue
    if (second) then
        allocate( b2lindex(1:nvar2l) )
    endif
    
    nvar2l = 0
    do l = 0, l2max_magnetic
        nvar2l = nvar2l + 1
        if (second) then
            b2lindex(nvar2l) = l       ! l0' (=l1')
        endif
    enddo
    if(.not.second) then
        second = .true.
        goto 2100
    endif

    if (.not.readsuccess) then
    
        second = .false.
2000 continue
        if (second) then
            allocate( b2occ(1:nvar2) )
            allocate( b2index(1:4,1:nvar2) )
        endif
  
        nrank = 2
        nmax = 50
        if (.not.second) then
            allocate( vecn(1:nrank) )
            allocate( vecnr(1:nrank) )
            allocate( vecl(1:nrank) )
            allocate( arrlint(1:nmax,1:nrank) )
        endif
  
        nvar2 = 0
        do nele1 = 1, nelements
            do nr1 = 1, nradial2_magnetic       
                do nr2 = 1, nradial2_magnetic
                    do nr3 = 1, nradbase_atomic
                        do k = 1, nvar2l
                        
                            l = b2lindex(k)
                            
                            vecn(1) = nele1
                            vecn(2) = nele1
                            vecnr(1) = nr1
                            vecnr(2) = nr2
                            vecl(1) = l
                            vecl(2) = l

                            call bcreate(magnetic,nrank,nmax,vecn,vecnr,vecl,numlint,arrlint,isordered,isbasis)
                        
                            if ((isordered).and.(isbasis)) then
                                index = .true.
                            else
                                index = .false.
                            endif                    
                        
                            if (index) then
                                nvar2 = nvar2 + 1
                                if (second) then
                                    b2occ(nvar2) = nele1      ! mu_1
                                    b2index(1,nvar2) = nr1    ! n0'
                                    b2index(2,nvar2) = nr2    ! n1'
                                    b2index(3,nvar2) = nr3    ! n1
                                    b2index(4,nvar2) = l      ! l0' (l1=0)
                                endif
                            endif
                    
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if (.not.second) then
            second = .true.
            goto 2000
        endif 
  
        if (second) then
            deallocate( vecn )
            deallocate( vecnr )
            deallocate( vecl )
            deallocate( arrlint )
        endif
  
    endif

    !  if (verblevel > 0 ) then
    write(*,*) 'number of variables for 1-body clusters:', nvar2
    !  endif

    allocate( b2map(1:nvar2) )
    do n = 1, nvar2
        b2map(n) = b2index(4,n) + 1
    enddo
  
    ! number of variables for a given angular momentum
    maddmax = 0
    allocate(madd2(nvar2l))
    do k = 1, nvar2l
        l1 = b2lindex(k)  
        madd2(k) = 0
        do m1 = -l1, l1
            madd2(k) = madd2(k) + 1
        enddo
        if (madd2(k) > maddmax) then
            maddmax = madd2(k)
        endif
    enddo
  
    allocate( b2lm(1:nvar2l,1:maddmax) )
    do k = 1, nvar2l
        l1 = b2lindex(k)
        m = 0
        do m1 = -l1, l1
            m = m + 1
            b2lm(k,m) = m1
        enddo
    enddo
  
    !---------------------------------------------------------
    ! number of variables for c3
    !---------------------------------------------------------
    ! first atomic, then atomic indices are stored in the 
    ! global arrays
    ! when it is possible the global arrays are splitted in 
    ! magnetic and atomic component
 
    ! angular part
    second = .false.
3000 continue
    if (second) then
        allocate( b3lindex(1:4,1:nvar3l) )
    endif
  
    nvar3l = 0 
    do l1 = 0, l3max_atomic               ! l1 (=l2)
        do l2 = 0, l3max_magnetic          ! l0'
            do l3 = 0, l3max_magnetic       ! l1'
                do l4 = 0, l3max_magnetic    ! l2' (=L01')
                    index = .false.
                    if ( ( l2 + l3 >= l4) .and. ( abs(l2-l3) <= l4 ) .and. (mod(l2+l3+l4,2) == 0 )) then
                        index = .true.
                    endif
                    if (index) then
                        nvar3l = nvar3l + 1
                        if (second) then
                            b3lindex(1,nvar3l) = l1   ! l1  (=l2)
                            b3lindex(2,nvar3l) = l2   ! l0'
                            b3lindex(3,nvar3l) = l3   ! l1'
                            b3lindex(4,nvar3l) = l4   ! l2' (=L01')
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
    
    if(.not.second) then
        second = .true.
        goto 3000
    endif
  
    ! angular part for atomic part only
    second = .false.
3050 continue
    if (second) then
        allocate( b3lindex_atomic(1:nvar3l_atomic) )
    endif  
    nvar3l_atomic = 0 
    do l1 = 0, l3max_atomic
        nvar3l_atomic = nvar3l_atomic + 1
        if (second) then
            b3lindex_atomic(nvar3l_atomic) = l1   ! l1  (=l2)
        endif
    enddo  
    if(.not.second) then
        second = .true.
        goto 3050
    endif
  
    ! angular part for magnetic part only  
    second = .false.
3055 continue
    if (second) then
        allocate( b3lindex_magnetic(1:3,1:nvar3l_magnetic) )
    endif  
    nvar3l_magnetic = 0 
    do l2 = 0, l3max_magnetic
        do l3 = 0, l3max_magnetic
            do l4 = 0, l3max_magnetic
                index = .false.
                if ( ( l2 + l3 >= l4) .and. ( abs(l2-l3) <= l4 ) .and. (mod(l2+l3+l4,2) == 0 )) then
                    index = .true.
                endif
                if (index) then
                    nvar3l_magnetic = nvar3l_magnetic + 1
                    if (second) then
                        b3lindex_magnetic(1,nvar3l_magnetic) = l2   ! l0'
                        b3lindex_magnetic(2,nvar3l_magnetic) = l3   ! l1'
                        b3lindex_magnetic(3,nvar3l_magnetic) = l4   ! l2' (=L01')
                    endif
                endif
            enddo
        enddo
    enddo  
    if(.not.second) then
        second = .true.
        goto 3055
    endif

    if (.not.readsuccess) then
  
        ! including radial part  
        second = .false.
3100 continue
        if (second) then
            allocate( b3occ_tmp(1:2,1:nvar3_tmp) )
            allocate( b3index_tmp(1:9,1:nvar3_tmp) )
        endif
  
        nvar3_tmp = 0
        ! nele1, nele2: neighboring atoms
        do nele1 = 1, nelements                             ! mu_1
            do nele2 = 1, nelements                          ! mu_2
                do nr1 = 1, nradial3_atomic                   ! n1
                    do nr2 = 1, nradial3_atomic                ! n2
                        do nr3 = 1, nradial3_magnetic           ! n0'
                            do nr4 = 1, nradial3_magnetic        ! n1'
                                do nr5 = 1, nradial3_magnetic     ! n2'
                                    do l1 = 0, l3max_atomic
                                        do l2 = 0, l3max_magnetic
                                            do l3 = 0, l3max_magnetic
                                                do l4 = 0, l3max_magnetic

!                        do k = 1, nvar3l
!                           l1 = b3lindex(1,k)        ! l1  (=l2)
!                           l2 = b3lindex(2,k)        ! l0'
!                           l3 = b3lindex(3,k)        ! l1'
!                           l4 = b3lindex(4,k)        ! l2' (=L01')
                          
                          k_atomic1 = nradial3_atomic*l1 + nr1
                          k_atomic2 = nradial3_atomic*l1 + nr2
                          
                          k_magnetic1 = nradial3_magnetic*l3 + nr4
                          k_magnetic2 = nradial3_magnetic*l4 + nr5
                          
                          kmax_atomic = nradial3_atomic*(l3max_atomic + 1)
                          kmax_magnetic = nradial3_magnetic*(l3max_magnetic + 1)
                          
                          nc1 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele1+  &
                          &  (kmax_magnetic+1)*k_atomic1 + k_magnetic1
                          nc2 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele2+  &
                          &  (kmax_magnetic+1)*k_atomic2 + k_magnetic2
                                                                                              
                          if ( (nc1 >= nc2) ) then
                            index = .true.
                          else
                            index = .false.
                          endif
                          
                          if (index) then
                             if (k_atomic1 .lt. k_atomic2) then
                                print*, 'Error! the atomic basis is not ordered lexicographically'
                                stop
                             endif
                          endif
                                
                          if (index) then
                              nvar3_tmp = nvar3_tmp + 1
                              if (second) then
                                    b3occ_tmp(1,nvar3_tmp) = nele1   ! mu_1
                                    b3occ_tmp(2,nvar3_tmp) = nele2   ! mu_2
                                    b3index_tmp(1,nvar3_tmp) = nr1   ! n1
                                    b3index_tmp(2,nvar3_tmp) = nr2   ! n2
                                    b3index_tmp(3,nvar3_tmp) = nr3   ! n0'
                                    b3index_tmp(4,nvar3_tmp) = nr4   ! n1'
                                    b3index_tmp(5,nvar3_tmp) = nr5   ! n2'
                                    b3index_tmp(6,nvar3_tmp) = l1    ! l1  (=l2)
                                    b3index_tmp(7,nvar3_tmp) = l2    ! l0'
                                    b3index_tmp(8,nvar3_tmp) = l3    ! l1'
                                    b3index_tmp(9,nvar3_tmp) = l4    ! l2' (=L01')
                              endif
                          endif
                          
                                                enddo
                                            enddo
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if(.not.second) then
            second = .true.
            goto 3100
        endif
  

        second = .false.
        3200 continue
        if (second) then
            allocate( b3occ(1:2,1:nvar3) )
            allocate( b3index(1:9,1:nvar3) )
        endif
    
        nrank = 2
        nmax = 50
        if (.not.second) then
            allocate( vecn_atomic(1:nrank) )
            allocate( vecn_magnetic(1:nrank+1) )
            allocate( vecnr_atomic(1:nrank) )
            allocate( vecnr_magnetic(1:nrank+1) )
            allocate( vecl_atomic(1:nrank) )
            allocate( vecl_magnetic(1:nrank+1) )
            allocate( arrlint_atomic(1:nmax,1:nrank) )
            allocate( arrlint_magnetic(1:nmax,1:nrank+1) )
        endif
  
        nvar3 = 0  
        do nvar = 1, nvar3_tmp
        
            vecn_atomic(1) = b3occ_tmp(1,nvar)               ! mu_1
            vecn_atomic(2) = b3occ_tmp(2,nvar)               ! mu_2 
            
            vecnr_atomic(1) = b3index_tmp(1,nvar)     ! n1
            vecnr_atomic(2) = b3index_tmp(2,nvar)     ! n2
            vecl_atomic(1) = b3index_tmp(6,nvar)      ! l1  
            vecl_atomic(2) = b3index_tmp(6,nvar)      ! l2

            call bcreate(atomic,nrank,nmax,vecn_atomic,vecnr_atomic,vecl_atomic,numlint_atomic,  &
            &   arrlint_atomic,isordered_atomic,isbasis_atomic)
                
            vecn_magnetic(1) = vecn_atomic(1)               ! mu_0
            vecn_magnetic(2) = vecn_atomic(1)               ! mu_1          
            vecn_magnetic(3) = vecn_atomic(2)               ! mu_2    
                
            vecnr_magnetic(1) = b3index_tmp(3,nvar)   ! n0'
            vecnr_magnetic(2) = b3index_tmp(4,nvar)   ! n1'
            vecnr_magnetic(3) = b3index_tmp(5,nvar)   ! n2'
            vecl_magnetic(1) = b3index_tmp(7,nvar)    ! l0'
            vecl_magnetic(2) = b3index_tmp(8,nvar)    ! l1'
            vecl_magnetic(3) = b3index_tmp(9,nvar)    ! l2' (=L01') 

            call bcreate(magnetic,nrank+1,nmax,vecn_magnetic,vecnr_magnetic,vecl_magnetic,numlint_magnetic,  &   
                                &  arrlint_magnetic,isordered_magnetic,isbasis_magnetic) 
            
            if ((isordered_atomic).and.(isbasis_atomic).and.(isbasis_magnetic)) then
                index = .true.
            else
                index = .false.
            endif
                
            if (index) then
            nvar3 = nvar3 + 1
                if (second) then
                    b3occ(1,nvar3) = vecn_atomic(1)         ! mu_1
                    b3occ(2,nvar3) = vecn_atomic(2)         ! mu_2
                    b3index(1,nvar3) = vecnr_atomic(1)      ! n1
                    b3index(2,nvar3) = vecnr_atomic(2)      ! n2
                    b3index(3,nvar3) = vecnr_magnetic(1)    ! n0'
                    b3index(4,nvar3) = vecnr_magnetic(2)    ! n1'
                    b3index(5,nvar3) = vecnr_magnetic(3)    ! n2'
                    b3index(6,nvar3) = vecl_atomic(1)       ! l1  (=l2)
                    b3index(7,nvar3) = vecl_magnetic(1)     ! l0'
                    b3index(8,nvar3) = vecl_magnetic(2)     ! l1'
                    b3index(9,nvar3) = vecl_magnetic(3)     ! l2' (=L01')
                endif
            endif     
            
        enddo
        if(.not.second) then
            second = .true.
            goto 3200
        endif

        if (second) then
            deallocate( vecn_atomic )
            deallocate( vecn_magnetic )
            deallocate( vecnr_atomic )
            deallocate( vecnr_magnetic )
            deallocate( vecl_atomic )
            deallocate( vecl_magnetic )
            deallocate( arrlint_atomic )
            deallocate( arrlint_magnetic )
        endif
  

    endif 
  
!  if (verblevel > 0 ) then
    write(*,*) 'number of variables for 3-body clusters:', nvar3
!  endif
  
    allocate( b3map_atomic(1:nvar3) )
    do n = 1, nvar3
        l1r = b3index(6,n)
        do k = 1, nvar3l_atomic
            l1 = b3lindex_atomic(k)
            if (l1 == l1r) then
                b3map_atomic(n) = k
            endif
        enddo
    enddo
  
    allocate( b3map_magnetic(1:nvar3) )
    do n = 1, nvar3
        l2r = b3index(7,n)
        l3r = b3index(8,n)
        l4r = b3index(9,n)
        do k = 1, nvar3l_magnetic
            l2 = b3lindex_magnetic(1,k)
            l3 = b3lindex_magnetic(2,k)
            l4 = b3lindex_magnetic(3,k)
            if ((l2 ==l2r).and.(l3==l3r).and.(l4==l4r)) then
                b3map_magnetic(n) = k
            endif
        enddo
    enddo
    
    ! number of ms for the atomic part
    maddmax_atomic = 0
    allocate(madd3_atomic(nvar3l_atomic))
    do k = 1, nvar3l_atomic
        l1 = b3lindex_atomic(k)  ! l1  (=l2)     
        madd3_atomic(k) = 0
        do m1 = -l1, l1
            madd3_atomic(k) = madd3_atomic(k) + 1
        enddo
        if (madd3_atomic(k) > maddmax_atomic) then
            maddmax_atomic = madd3_atomic(k)
        endif
    enddo
  
    ! number of ms for the magnetic part
    maddmax_magnetic = 0
    allocate(madd3_magnetic(nvar3l_magnetic))
    do k = 1, nvar3l_magnetic
        l2 = b3lindex_magnetic(1,k)  ! l0'
        l3 = b3lindex_magnetic(2,k)  ! l1'
        l4 = b3lindex_magnetic(3,k)  ! l2' (=L01')     
        madd3_magnetic(k) = 0
        do m2 = -l2, l2
            do m3 = -l3, l3
                m4 = - (m2 + m3)
                if ( abs(m4) <= l4 ) then
                    madd3_magnetic(k) = madd3_magnetic(k) + 1
                endif
            enddo
        enddo
        if (madd3_magnetic(k) > maddmax_magnetic) then
            maddmax_magnetic = madd3_magnetic(k)
        endif
    enddo
    
    ! list of allowed values of ms for atomic part
    allocate( b3lm_atomic(1:nvar3l_atomic,1:maddmax_atomic) )
    do k = 1, nvar3l_atomic
        l1 = b3lindex_atomic(k)  ! l1  (=l2)     
        m = 0
        do m1 = -l1, l1
            m = m + 1
            b3lm_atomic(k,m) = m1   ! m1 (=-m2)
        enddo
        if ( m /= madd3_atomic(k) ) then
            print *,'Something wrong here.'
            stop
        endif
    enddo 
    
    ! list of allowed values of ms for magnetic part
    allocate( b3lm_magnetic(1:nvar3l_magnetic,1:maddmax_magnetic,1:3) )
    do k = 1, nvar3l_magnetic
        l2 = b3lindex_magnetic(1,k)  ! l0'
        l3 = b3lindex_magnetic(2,k)  ! l1'
        l4 = b3lindex_magnetic(3,k)  ! l2' (=L01')     
        m = 0
        do m2 = -l2, l2
            do m3 = -l3, l3
                m4 = -(m2 + m3)
                if ( abs(m4) <= l4 ) then
                    m = m + 1
                    b3lm_magnetic(k,m,1) = m2   ! m0'
                    b3lm_magnetic(k,m,2) = m3   ! m1'
                    b3lm_magnetic(k,m,3) = m4   ! m2' (=-(m0'+m1'))
                endif
            enddo
        enddo
        if ( m /= madd3_magnetic(k) ) then
            print *,'Something wrong here.'
            stop
        endif
    enddo 

    !---------------------------------------------------------
    ! number of variables for c4
    !---------------------------------------------------------
    ! angular part for atomic part only
    second = .false.
4050 continue
    if (second) then
        allocate( b4lindex_atomic(1:3,1:nvar4l_atomic) )
    endif
    nvar4l_atomic = 0 
    do l1 = 0, l4max_atomic                              ! l1
        do l2 = 0, l4max_atomic                           ! l2
            do l3 = 0, l4max_atomic                        ! l3
                index = .false.
                ! condition on atomic indices
                if ( ( l1 + l2 >= l3) .and. ( abs(l1-l2) <= l3 ) .and. (mod(l1+l2+l3,2) == 0 )) then
                    index = .true.
                endif
                
                if (index) then
                nvar4l_atomic = nvar4l_atomic + 1
                    if (second) then
                        b4lindex_atomic(1,nvar4l_atomic) = l1
                        b4lindex_atomic(2,nvar4l_atomic) = l2
                        b4lindex_atomic(3,nvar4l_atomic) = l3
                    endif
                endif                          
            enddo
        enddo
    enddo
    if(.not.second) then
        second = .true.                     
        goto 4050
    endif
  

    if (.not.readsuccess) then
  
        ! find basis functions  
        second = .false.
4100 continue
        if (second) then
            allocate( b4occ_tmp(1:3,1:nvar4_tmp) )
            allocate( b4index_tmp(1:14,1:nvar4_tmp) )
        endif
  
        nvar4_tmp = 0
        do nele1 = 1, nelements                                                        ! mu_1
            do nele2 = 1, nelements                                                     ! mu_2
                do nele3 = 1, nelements                                                  ! mu_3
                    do nr1 = 1, nradial4_atomic                                        ! n1
                        do nr2 = 1, nradial4_atomic                                     ! n2
                            do nr3 = 1, nradial4_atomic                                  ! n3
                                do nr4 = 1, nradial4_magnetic                             ! n0'
                                    do nr5 = 1, nradial4_magnetic                          ! n1'
                                        do nr6 = 1, nradial4_magnetic                       ! n2'
                                            do nr7 = 1, nradial4_magnetic                    ! n3'
                                                do l1 = 0, l4max_atomic
                                                    do l2 = 0, l4max_atomic
                                                        do l3 = 0, l4max_atomic
                                                            do l4 = 0, l4max_magnetic
                                                                do l5 = 0, l4max_magnetic
                                                                    do l6 = 0, l4max_magnetic
                                                                        do l7 = 0, l4max_magnetic

                                    
                                    k_atomic1 = nradial4_atomic*l1 + nr1
                                    k_atomic2 = nradial4_atomic*l2 + nr2
                                    k_atomic3 = nradial4_atomic*l3 + nr3
                                    
                                    k_magnetic1 = nradial4_magnetic*l5 + nr5 
                                    k_magnetic2 = nradial4_magnetic*l6 + nr6
                                    k_magnetic3 = nradial4_magnetic*l7 + nr7
                                    
                                    kmax_atomic = nradial4_atomic*(l4max_atomic + 1)
                                    kmax_magnetic = nradial4_magnetic*(l4max_magnetic + 1)
                                    
                                    nc1 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele1 +  &
                                    &  (kmax_magnetic+1)*k_atomic1 + k_magnetic1
                                    nc2 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele2 +  &
                                    &  (kmax_magnetic+1)*k_atomic2 + k_magnetic2 
                                    nc3 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele3 +  &
                                    &  (kmax_magnetic+1)*k_atomic3 + k_magnetic3 
                                    

                                    ! take only ordered sequence
                                    if (( nc1 >= nc2 ).and.( nc2 >= nc3 )) then
                                       index = .true.
                                    else
                                       index = .false.
                                    endif
                                    
                                    if (index) then
                                        if ((k_atomic1 .lt. k_atomic2).or.(k_atomic2 .lt. k_atomic3)) then
                                            print*, 'Error! the atomic basis is not ordered lexicographically'
                                            stop
                                        endif
                                    endif
                                                            
                                    ! now everything should be o.k. for coupling                                      
                                    if (index) then
                                        nvar4_tmp = nvar4_tmp + 1
                                        if (second) then
                                            b4occ_tmp(1,nvar4_tmp) = nele1
                                            b4occ_tmp(2,nvar4_tmp) = nele2
                                            b4occ_tmp(3,nvar4_tmp) = nele3
                                            b4index_tmp(1,nvar4_tmp) = nr1   ! n1
                                            b4index_tmp(2,nvar4_tmp) = nr2   ! n2
                                            b4index_tmp(3,nvar4_tmp) = nr3   ! n3
                                            b4index_tmp(4,nvar4_tmp) = nr4   ! n0' 
                                            b4index_tmp(5,nvar4_tmp) = nr5   ! n1'
                                            b4index_tmp(6,nvar4_tmp) = nr6   ! n2'
                                            b4index_tmp(7,nvar4_tmp) = nr7   ! n3'
                                            b4index_tmp(8,nvar4_tmp) = l1    ! l1
                                            b4index_tmp(9,nvar4_tmp) = l2    ! l2
                                            b4index_tmp(10,nvar4_tmp) = l3   ! l3
                                            b4index_tmp(11,nvar4_tmp) = l4   ! l0'
                                            b4index_tmp(12,nvar4_tmp) = l5   ! l1'
                                            b4index_tmp(13,nvar4_tmp) = l6   ! l2'
                                            b4index_tmp(14,nvar4_tmp) = l7   ! l3'
                                        endif
                                    endif 
                                                                        enddo
                                                                    enddo
                                                                enddo
                                                            enddo
                                                        enddo
                                                    enddo
                                                enddo
                                            enddo
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo  
        if(.not.second) then
            second = .true.
            goto 4100
        endif
  
  
        second = .false.
4200 continue
        if (second) then
            allocate( b4occ(1:3,1:nvar4) )
            allocate( b4index(1:15,1:nvar4) )
            allocate( b4nint(1:nvar4) )
            b4nint = 0
        endif
  
        nrank = 3
        nmax = 50
        if (.not.second) then
            allocate( vecn_atomic(1:nrank) )
            allocate( vecn_magnetic(1:nrank+1) )
            allocate( vecnr_atomic(1:nrank) )
            allocate( vecnr_magnetic(1:nrank+1) )
            allocate( vecl_atomic(1:nrank) )
            allocate( vecl_magnetic(1:nrank+1) )
            allocate( arrlint_atomic(1:nmax,1:nrank) )
            allocate( arrlint_magnetic(1:nmax,1:nrank+1) )
        endif
  
        nvar4 = 0
        do nvar = 1, nvar4_tmp
            
            vecn_atomic(1) = b4occ_tmp(1,nvar)               ! mu_1
            vecn_atomic(2) = b4occ_tmp(2,nvar)               ! mu_2 
            vecn_atomic(3) = b4occ_tmp(3,nvar)               ! mu_3

            vecnr_atomic(1) = b4index_tmp(1,nvar)     ! n1
            vecnr_atomic(2) = b4index_tmp(2,nvar)     ! n2
            vecnr_atomic(3) = b4index_tmp(3,nvar)     ! n3
            vecl_atomic(1) = b4index_tmp(8,nvar)      ! l1  
            vecl_atomic(2) = b4index_tmp(9,nvar)      ! l2
            vecl_atomic(3) = b4index_tmp(10,nvar)     ! l3

            call bcreate(atomic,nrank,nmax,vecn_atomic,vecnr_atomic,vecl_atomic,numlint_atomic,  &
                            &   arrlint_atomic,isordered_atomic,isbasis_atomic)
            
            if(.not.isordered_atomic) then
                print*, 'Something wrong with rank = 3.'
                stop
            endif
            
            vecn_magnetic(1) = vecn_atomic(1)               ! mu_0
            vecn_magnetic(2) = vecn_atomic(1)               ! mu_1 
            vecn_magnetic(3) = vecn_atomic(2)               ! mu_2
            vecn_magnetic(4) = vecn_atomic(3)               ! mu_3
                
            vecnr_magnetic(1) = b4index_tmp(4,nvar)   ! n0'
            vecnr_magnetic(2) = b4index_tmp(5,nvar)   ! n1'
            vecnr_magnetic(3) = b4index_tmp(6,nvar)   ! n2'
            vecnr_magnetic(4) = b4index_tmp(7,nvar)   ! n3'
            vecl_magnetic(1) = b4index_tmp(11,nvar)    ! l0'
            vecl_magnetic(2) = b4index_tmp(12,nvar)    ! l1'
            vecl_magnetic(3) = b4index_tmp(13,nvar)    ! l2' 
            vecl_magnetic(4) = b4index_tmp(14,nvar)    ! l3'

            call bcreate(magnetic,nrank+1,nmax,vecn_magnetic,vecnr_magnetic,vecl_magnetic,numlint_magnetic,  &   
                            &  arrlint_magnetic,isordered_magnetic,isbasis_magnetic) 
            
            if( isordered_atomic.and.isbasis_atomic.and.isbasis_magnetic ) then
                
                do nint = 1, numlint_magnetic
                    nvar4 = nvar4 + 1
                    if (second) then
                        b4occ(1:3,nvar4) = vecn_atomic(1:3)
                        b4index(1:3,nvar4) = vecnr_atomic(1:3)
                        b4index(4:7,nvar4) = vecnr_magnetic(1:4)
                        b4index(8:10,nvar4) = vecl_atomic(1:3)
                        b4index(11:14,nvar4) = vecl_magnetic(1:4)
                        b4index(15,nvar4) = arrlint_magnetic(nint,1)
                        if (nint == 1) then
                            b4nint(nvar4) = numlint_magnetic
                        else
                            b4nint(nvar4) = 0
                        endif
                    endif
                enddo
                
            endif
            
        enddo
        if(.not.second) then
            second = .true.
            goto 4200
        endif

        if (second) then
            deallocate( vecn_atomic )
            deallocate( vecn_magnetic )
            deallocate( vecnr_atomic )
            deallocate( vecnr_magnetic )
            deallocate( vecl_atomic )
            deallocate( vecl_magnetic )
            deallocate( arrlint_atomic )
            deallocate( arrlint_magnetic )
        endif
  
    endif
  
  
    ! check for nint 1
    do n = 1, nvar4
        do m = 1, b4nint(nvar4) - 1
            if (b4nint(n+m) /= 0 ) then
                print *,'Something wrong with order of internal couplings in b4nint.'
                stop
            endif
        enddo
    enddo
  
    !  if (verblevel > 0 ) then
    write(*,*) 'number of variables for 4-body clusters:', nvar4
    !  endif
      
    allocate( b4map_atomic(1:nvar4) )
    b4map_atomic = 0  
    do n = 1, nvar4
        l1r = b4index(8,n) 
        l2r = b4index(9,n)
        l3r = b4index(10,n)
        do k = 1, nvar4l_atomic
            l1 = b4lindex_atomic(1,k) 
            l2 = b4lindex_atomic(2,k)
            l3 = b4lindex_atomic(3,k)
            if ((l1 == l1r).and.(l2 == l2r).and.(l3 == l3r)) then
                b4map_atomic(n) = k            
            endif
        enddo
    enddo
  
    ! b5map groups together all b5index vectors that have the same angular part
    allocate( b4map_magnetic(1:nvar4) )
    b4map_magnetic = 0

    do n = 1, nvar4
        l1 = b4index(11,n) 
        l2 = b4index(12,n)
        l3 = b4index(13,n)
        l4 = b4index(14,n)
        lint = b4index(15,n)
        do m = 1, nvar4
            l1r = b4index(11,m) 
            l2r = b4index(12,m)
            l3r = b4index(13,m)
            l4r = b4index(14,m)
            lintr = b4index(15,m)
            if ((l1r == l1).and.(l2r == l2).and.(l3r == l3).and.(l4r == l4).and.(lintr == lint)) then
                if ( b4map_magnetic(n) >= 0) then
                    ! take: increase index at n by 1
                    b4map_magnetic(n) = b4map_magnetic(n) + 1
                    ! burn m: has been found already
                    if (m > n) then
                        b4map_magnetic(m) = -n
                    endif
                    ! as both n and m start couting from 1, this may not happen
                    if (m < n) then
                        print *,'Should not be here.'
                        stop
                    endif              
                endif
            endif
        enddo
    enddo
    m = 0
    l = 0
    do n = 1, nvar4
        ! check if every n has been collected
        if (b4map_magnetic(n) == 0) then
            print *,'This may not happen.'
            stop
        endif 
        if (b4map_magnetic(n) > 0) then
            m = m + 1
            l = l + b4map_magnetic(n)
        endif   
    enddo
    ! number of different angular function vectors l1, l2, l3, l4, lint
    nvar4l_magnetic = m
    if ( nvar4 /= l) then
        print *,'Error in checksum for nvar4.'
        stop
    endif

!  print *,'nvar5l', nvar5l

    ! fill b5lindex
    allocate( b4lindex_magnetic(1:5,1:nvar4l_magnetic))
    m = 0
    do n = 1, nvar4
        if (b4map_magnetic(n) > 0) then
            m = m + 1
            b4lindex_magnetic(1:5,m) = b4index(11:15,n)
    !        print *, m, b5lindex(1:5,m)
        endif
    enddo
    ! finally update b5map
    m = 0
    do n = 1, nvar4
        if (b4map_magnetic(n) > 0) then
            m = m + 1
            b4map_magnetic(n) = m
        endif
    enddo
    do n = 1, nvar4
        if (b4map_magnetic(n) < 0) then
            m = abs(b4map_magnetic(n))
            b4map_magnetic(n) = b4map_magnetic(m)
        endif
    enddo  
    ! validate
    do n = 1, nvar4
        m = b4map_magnetic(n)
        k = sum(abs(b4lindex_magnetic(1:5,m)-b4index(11:15,n)))
        if (k > 0 ) then
            print *,'Something wrong with b4lindex_magnetic.'
            stop
        endif
    enddo  
  
  
    ! madd4 for atomic part only
    maddmax_atomic = 0
    allocate(madd4_atomic(nvar4l_atomic))
    do k = 1, nvar4l_atomic
        l1 = b4lindex_atomic(1,k)   ! l1
        l2 = b4lindex_atomic(2,k)   ! l2
        l3 = b4lindex_atomic(3,k)   ! l3
        
        !  finire qui con con i possibili valori di m e m'
        madd4_atomic(k) = 0
        do m1 = -l1, l1                      ! m1
            do m2 = -l2, l2                   ! m2
                m3 = -(m1 + m2)
                if ( abs(m3) <=l3 ) then
                    madd4_atomic(k) = madd4_atomic(k) + 1
                endif
            enddo
        enddo
        if (madd4_atomic(k) > maddmax_atomic) then
            maddmax_atomic = madd4_atomic(k)
        endif
    enddo
  
    ! madd4 for magnetic part only
    maddmax_magnetic = 0
    allocate(madd4_magnetic(nvar4l_magnetic))
    do k = 1, nvar4l_magnetic
        l4 = b4lindex_magnetic(1,k)   ! l0'
        l5 = b4lindex_magnetic(2,k)   ! l1'
        l6 = b4lindex_magnetic(3,k)   ! l2'
        l7 = b4lindex_magnetic(4,k)   ! l3'
        L  = b4lindex_magnetic(5,k)   ! L01'
        
        !  finire qui con con i possibili valori di m e m'
        madd4_magnetic(k) = 0
        do M = -L, L                            ! M01' = (m0'+m1') and M' = -m2'-m3'
            do m4 = -l4, l4                ! m0'
                do m5 = -l5, l5             ! m1'
                    do m6 = -l6, l6          ! m2'
                        do m7 = -l7, l7       ! m3'
                            if ( (M - (m4+m5) == 0).and.(M + (m6+m7) == 0) ) then
                                madd4_magnetic(k) = madd4_magnetic(k) + 1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if (madd4_magnetic(k) > maddmax_magnetic) then
            maddmax_magnetic = madd4_magnetic(k)
        endif
    enddo


    ! b4lm for atomic part only
    allocate( b4lm_atomic(1:nvar4l_atomic,1:maddmax_atomic,1:3) )
    do k = 1, nvar4l_atomic
        l1 = b4lindex_atomic(1,k)   ! l1
        l2 = b4lindex_atomic(2,k)   ! l2
        l3 = b4lindex_atomic(3,k)   ! l3
        n = 0
        do m1 = -l1, l1                      ! m1
            do m2 = -l2, l2                   ! m2
                m3 = -(m1 + m2)
                if ( abs(m3) <=l3 ) then
                    n = n + 1
                    b4lm_atomic(k,n,1) = m1
                    b4lm_atomic(k,n,2) = m2
                    b4lm_atomic(k,n,3) = m3
                endif
            enddo
        enddo
        if ( n /= madd4_atomic(k) ) then
            print *,'Something wrong here.'
            stop
        endif
    enddo
  
    ! b4lm for magnetic part only
    allocate( b4lm_magnetic(1:nvar4l_magnetic,1:maddmax_magnetic,1:4) )
    do k = 1, nvar4l_magnetic
        l4 = b4lindex_magnetic(1,k)   ! l0'
        l5 = b4lindex_magnetic(2,k)   ! l1'
        l6 = b4lindex_magnetic(3,k)   ! l2'
        l7 = b4lindex_magnetic(4,k)   ! l3'
        L  = b4lindex_magnetic(5,k)   ! L01'

        n = 0
        do M = -L, L                      ! M01' = (m0'+m1') and M' = -m2'-m3'
            do m4 = -l4, l4                ! m0'
                do m5 = -l5, l5             ! m1'
                    do m6 = -l6, l6          ! m2'
                        do m7 = -l7, l7       ! m3'
                            if ( (M - (m4+m5) == 0).and.(M + (m6+m7) == 0) ) then
                                n = n + 1
                                b4lm_magnetic(k,n,1) = m4
                                b4lm_magnetic(k,n,2) = m5
                                b4lm_magnetic(k,n,3) = m6
                                b4lm_magnetic(k,n,4) = m7
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( n /= madd4_magnetic(k) ) then
            print *,'Something wrong here.'
            stop
        endif
    enddo
  
    !---------------------------------------------------------
    ! number of variables for c5
    !---------------------------------------------------------
    
    if (.not.readsuccess) then
  
        ! find basis functions  
        second = .false.
5100 continue
        if (second) then
            allocate( b5occ_tmp(1:4,1:nvar5_tmp) )
            allocate( b5index_tmp(1:18,1:nvar5_tmp) )
        endif
  
        ! coupling scheme (b1*b2)_lint(b3*b4)_lint

        nvar5_tmp = 0
        do nele1 = 1, nelements                                                        ! mu_1
            do nele2 = 1, nelements                                                     ! mu_2
                do nele3 = 1, nelements                                                  ! mu_3
                    do nele4 = 1, nelements                                               ! mu_4
                        do nr1 = 1, nradial5_atomic                                     ! n1
                            do nr2 = 1, nradial5_atomic                                  ! n2
                                do nr3 = 1, nradial5_atomic                               ! n3
                                    do nr4 = 1, nradial5_atomic                            ! n4
                                        do nr5 = 1, nradial5_magnetic                       ! n0'
                                            do nr6 = 1, nradial5_magnetic                    ! n1'
                                                do nr7 = 1, nradial5_magnetic                 ! n2'
                                                    do nr8 = 1, nradial5_magnetic              ! n3'
                                                        do nr9 = 1, nradial5_magnetic           ! n4'                                     
                                                            do l1 = 0, l5max_atomic                                        ! l1
                                                                do l2 = 0, l5max_atomic                                    ! l2
                                                                    do l3 = 0, l5max_atomic                                ! l3
                                                                        do l4 = 0, l5max_atomic                            ! l4
                                                                            do l5 = 0, l5max_magnetic                      ! l0'
                                                                                do l6 = 0, l5max_magnetic                  ! l1'
                                                                                    do l7 = 0, l5max_magnetic              ! l2'
                                                                                        do l8 = 0, l5max_magnetic          ! l3'
                                                                                            do l9 = 0, l5max_magnetic      ! l4'
                                               
                                                
                                                k_atomic1 = nradial5_atomic*l1 + nr1
                                                k_atomic2 = nradial5_atomic*l2 + nr2
                                                k_atomic3 = nradial5_atomic*l3 + nr3
                                                k_atomic4 = nradial5_atomic*l4 + nr4
                                                
                                                k_magnetic1 = nradial5_magnetic*l6 + nr6
                                                k_magnetic2 = nradial5_magnetic*l7 + nr7
                                                k_magnetic3 = nradial5_magnetic*l8 + nr8
                                                k_magnetic4 = nradial5_magnetic*l9 + nr9
                                                
                                                kmax_atomic = nradial5_atomic*(l5max_atomic + 1)
                                                kmax_magnetic = nradial5_magnetic*(l5max_magnetic + 1)
                                            
                                                nc1 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele1 +  &
                                                &  (kmax_magnetic+1)*k_atomic1 + k_magnetic1
                                                nc2 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele2 +  &
                                                &  (kmax_magnetic+1)*k_atomic2 + k_magnetic2 
                                                nc3 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele3 +  &
                                                &  (kmax_magnetic+1)*k_atomic3 + k_magnetic3 
                                                nc4 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele4 +  &
                                                &  (kmax_magnetic+1)*k_atomic4 + k_magnetic4 
                                                
                                                index = .true.
                                                
                                                ! take only ordered sequence
                                                if (( nc1 < nc2 ).or.( nc2 < nc3 ).or.( nc3 < nc4 )) then
                                                    index = .false.
                                                endif
                                            
                                                if (index) then
                                                    if ((k_atomic1 .lt. k_atomic2).or.(k_atomic2 .lt. k_atomic3).or.  & 
                                                    &    (k_atomic3 .lt. k_atomic4)) then
                                                        print*, 'Error! the atomic basis is not ordered lexicographically'
                                                        stop
                                                    endif
                                                endif                                            
                                                
                                                ! now everything should be o.k. for coupling                                      
                                                if (index) then
                                                    nvar5_tmp = nvar5_tmp + 1
                                                    if (second) then
                                                        b5occ_tmp(1,nvar5_tmp) = nele1   ! mu_1
                                                        b5occ_tmp(2,nvar5_tmp) = nele2   ! mu_2
                                                        b5occ_tmp(3,nvar5_tmp) = nele3   ! mu_3
                                                        b5occ_tmp(4,nvar5_tmp) = nele4   ! mu_4
                                                        b5index_tmp(1,nvar5_tmp) = nr1   ! n1
                                                        b5index_tmp(2,nvar5_tmp) = nr2   ! n2
                                                        b5index_tmp(3,nvar5_tmp) = nr3   ! n3
                                                        b5index_tmp(4,nvar5_tmp) = nr4   ! n4
                                                        b5index_tmp(5,nvar5_tmp) = nr5   ! n0'
                                                        b5index_tmp(6,nvar5_tmp) = nr6   ! n1'
                                                        b5index_tmp(7,nvar5_tmp) = nr7   ! n2'
                                                        b5index_tmp(8,nvar5_tmp) = nr8   ! n3'
                                                        b5index_tmp(9,nvar5_tmp) = nr9   ! n4'
                                                        b5index_tmp(10,nvar5_tmp) = l1   ! l1
                                                        b5index_tmp(11,nvar5_tmp) = l2   ! l2
                                                        b5index_tmp(12,nvar5_tmp) = l3   ! l3
                                                        b5index_tmp(13,nvar5_tmp) = l4   ! l4
                                                        b5index_tmp(14,nvar5_tmp) = l5   ! l0'
                                                        b5index_tmp(15,nvar5_tmp) = l6   ! l1'
                                                        b5index_tmp(16,nvar5_tmp) = l7   ! l2'
                                                        b5index_tmp(17,nvar5_tmp) = l8   ! l3'
                                                        b5index_tmp(18,nvar5_tmp) = l9   ! l4'                                                    
                                                    endif
                                                endif   
                                                
                                                                                            enddo
                                                                                        enddo
                                                                                    enddo
                                                                                enddo
                                                                            enddo
                                                                        enddo
                                                                    enddo
                                                                enddo
                                                            enddo
                                                        enddo
                                                    enddo
                                                enddo
                                            enddo
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo 
        if(.not.second) then
            second = .true.
            goto 5100
        endif
  
        second = .false.
5200 continue
        if (second) then
            allocate( b5occ(1:4,1:nvar5) )
            allocate( b5index(1:21,1:nvar5) )
            allocate( b5nint(1:nvar5) )
            b5nint = 0
        endif
  
        nrank = 4
        nmax = 100
        if (.not.second) then
            allocate( vecn_atomic(1:nrank) )
            allocate( vecn_magnetic(1:nrank+1) )
            allocate( vecnr_atomic(1:nrank) )
            allocate( vecnr_magnetic(1:nrank+1) )
            allocate( vecl_atomic(1:nrank) )
            allocate( vecl_magnetic(1:nrank+1) )
            allocate( arrlint_atomic(1:nmax,1:nrank) )
            allocate( arrlint_magnetic(1:nmax,1:nrank+1) )
        endif
  
        nvar5 = 0
        do nvar = 1, nvar5_tmp
            
            vecn_atomic(1) = b5occ_tmp(1,nvar)        ! mu_1
            vecn_atomic(2) = b5occ_tmp(2,nvar)        ! mu_2 
            vecn_atomic(3) = b5occ_tmp(3,nvar)        ! mu_3
            vecn_atomic(4) = b5occ_tmp(4,nvar)        ! mu_4

            vecnr_atomic(1) = b5index_tmp(1,nvar)     ! n1
            vecnr_atomic(2) = b5index_tmp(2,nvar)     ! n2
            vecnr_atomic(3) = b5index_tmp(3,nvar)     ! n3
            vecnr_atomic(4) = b5index_tmp(4,nvar)     ! n4
            vecl_atomic(1) = b5index_tmp(10,nvar)     ! l1  
            vecl_atomic(2) = b5index_tmp(11,nvar)     ! l2
            vecl_atomic(3) = b5index_tmp(12,nvar)     ! l3
            vecl_atomic(4) = b5index_tmp(13,nvar)     ! l4

            call bcreate(atomic,nrank,nmax,vecn_atomic,vecnr_atomic,vecl_atomic,numlint_atomic,  &
                            &   arrlint_atomic,isordered_atomic,isbasis_atomic)
            
            if(.not.isordered_atomic) then
                print*, 'Something wrong with rank = 4.'
                stop
            endif
            
            vecn_magnetic(1) = vecn_atomic(1)               ! mu_0
            vecn_magnetic(2) = vecn_atomic(1)               ! mu_1 
            vecn_magnetic(3) = vecn_atomic(2)               ! mu_2
            vecn_magnetic(4) = vecn_atomic(3)               ! mu_3
            vecn_magnetic(5) = vecn_atomic(4)               ! mu_4
                
            vecnr_magnetic(1) = b5index_tmp(5,nvar)    ! n0'
            vecnr_magnetic(2) = b5index_tmp(6,nvar)    ! n1'
            vecnr_magnetic(3) = b5index_tmp(7,nvar)    ! n2'
            vecnr_magnetic(4) = b5index_tmp(8,nvar)    ! n3'
            vecnr_magnetic(5) = b5index_tmp(9,nvar)    ! n4'
            vecl_magnetic(1) = b5index_tmp(14,nvar)    ! l0'
            vecl_magnetic(2) = b5index_tmp(15,nvar)    ! l1'
            vecl_magnetic(3) = b5index_tmp(16,nvar)    ! l2' 
            vecl_magnetic(4) = b5index_tmp(17,nvar)    ! l3'
            vecl_magnetic(5) = b5index_tmp(18,nvar)    ! l4'

            call bcreate(magnetic,nrank+1,nmax,vecn_magnetic,vecnr_magnetic,vecl_magnetic,numlint_magnetic,  &   
                            &  arrlint_magnetic,isordered_magnetic,isbasis_magnetic) 
            
            if( isbasis_atomic.and.isbasis_magnetic ) then
                        
                numlint = max(numlint_atomic,numlint_magnetic)
                do nint = 1, numlint
                    nvar5 = nvar5 + 1
                    if (second) then
                        b5occ(1:4,nvar5) = vecn_atomic(1:4)
                        if ( nint.le.numlint_atomic ) then
                            b5index(1:4,nvar5) = vecnr_atomic(1:4)
                            b5index(10:13,nvar5) = vecl_atomic(1:4)
                            b5index(19,nvar5) = arrlint_atomic(nint,1)
                        else
                            b5index(1:4,nvar5) = 0
                            b5index(10:13,nvar5) = 0
                            b5index(19,nvar5) = 0
                        endif
                        if ( nint.le.numlint_magnetic ) then
                            b5index(5:9,nvar5) = vecnr_magnetic(1:5)
                            b5index(14:18,nvar5) = vecl_magnetic(1:5)
                            b5index(20:21,nvar5) = arrlint_magnetic(nint,1:2)
                        else
                            b5index(5:9,nvar5) = 0
                            b5index(14:18,nvar5) = 0
                            b5index(20:21,nvar5) = 0
                        endif
                        if (nint == 1) then
                            b5nint(nvar5) = numlint
                        else
                            b5nint(nvar5) = 0
                        endif
                    endif
                enddo
                
            endif
            
        enddo
        if(.not.second) then
            second = .true.
            goto 5200
        endif

        if (second) then
            deallocate( vecn_atomic )
            deallocate( vecn_magnetic )
            deallocate( vecnr_atomic )
            deallocate( vecnr_magnetic )
            deallocate( vecl_atomic )
            deallocate( vecl_magnetic )
            deallocate( arrlint_atomic )
            deallocate( arrlint_magnetic )
        endif
  
        ! check for nint 1
        do n = 1, nvar5
            do m = 1, b5nint(nvar5) - 1
                if (b5nint(n+m) /= 0 ) then
                    print *,'Something wrong with order of internal couplings in b5nint.'
                    stop
                endif
            enddo
        enddo
  
        ! ordering internal couplings to takes smallest lint2+lint3 first (only magnetic lint)
        do n = 1, nvar5
            do m = 1, b5nint(n) - 2
                do k = m+1, b5nint(n) - 1
                    lsum1 = b5index(20,n+m) + b5index(21,n+m)
                    lsum2 = b5index(20,n+k) + b5index(21,n+k)
                    if (lsum1 > lsum2) then
                        n2v(1:2) = b5index(20:21,n+m)
                        b5index(20:21,n+m) = b5index(20:21,n+k)
                        b5index(20:21,n+k) = n2v(1:2)
                    endif
                enddo
            enddo
        enddo
  
    endif

    !  if (verblevel > 0 ) then
    write(*,*) 'number of variables for 5-body clusters:', nvar5
    !  endif
      
    ! b5map for atomic part only
    allocate( b5map_atomic(1:nvar5) )
    b5map_atomic = 0
  
    do n = 1, nvar5
        l1 = b5index(10,n) 
        l2 = b5index(11,n)
        l3 = b5index(12,n)
        l4 = b5index(13,n)
        lint = b5index(19,n)
        do m = 1, nvar5
            l1r = b5index(10,m) 
            l2r = b5index(11,m)
            l3r = b5index(12,m)
            l4r = b5index(13,m)
            lintr = b5index(19,m)
            if ((l1r == l1).and.(l2r == l2).and.(l3r == l3).and.(l4r == l4).and.(lintr == lint)) then
                if ( b5map_atomic(n) >= 0) then
                    ! take: increase index at n by 1
                    b5map_atomic(n) = b5map_atomic(n) + 1
                    ! burn m: has been found already
                    if (m > n) then
                        b5map_atomic(m) = -n
                    endif
                    ! as both n and m start couting from 1, this may not happen
                    if (m < n) then
                        print *,'Should not be here.'
                        stop
                    endif              
                endif
            endif
        enddo
    enddo
    m = 0
    l = 0
    do n = 1, nvar5
        ! check if every n has been collected
        if (b5map_atomic(n) == 0) then
            print *,'This may not happen.'
            stop
        endif 
        if (b5map_atomic(n) > 0) then
            m = m + 1
            l = l + b5map_atomic(n)
        endif   
    enddo
    ! number of different angular function vectors l1, l2, l3, l4, lint
    nvar5l_atomic = m
    if ( nvar5 /= l) then
        print *,'Error in checksum for nvar5.'
        stop
    endif

    ! fill b5lindex
    allocate( b5lindex_atomic(1:5,1:nvar5l_atomic))
    m = 0
    do n = 1, nvar5
        if (b5map_atomic(n) > 0) then
            m = m + 1
            b5lindex_atomic(1:4,m) = b5index(10:13,n)
            b5lindex_atomic(5,m) = b5index(19,n)
    !        print *, m, b5lindex(1:5,m)
        endif
    enddo
    ! finally update b5map
    m = 0
    do n = 1, nvar5
        if (b5map_atomic(n) > 0) then
            m = m + 1
            b5map_atomic(n) = m
        endif
    enddo
    do n = 1, nvar5
        if (b5map_atomic(n) < 0) then
            m = abs(b5map_atomic(n))
            b5map_atomic(n) = b5map_atomic(m)
        endif
    enddo  
    ! validate
    do n = 1, nvar5
        m = b5map_atomic(n)
        k = sum(abs(b5lindex_atomic(1:4,m)-b5index(10:13,n))) + &
            &   abs(b5lindex_atomic(5,m)-b5index(19,n))
        if (k > 0 ) then
            print *,'Something wrong with b5lindex_atomic.'
            stop
        endif
    enddo
  

    ! b5map for magnetic part only
    allocate( b5map_magnetic(1:nvar5) )
    b5map_magnetic = 0
    
    do n = 1, nvar5
        l1 = b5index(14,n) 
        l2 = b5index(15,n)
        l3 = b5index(16,n)
        l4 = b5index(17,n)
        l5 = b5index(18,n)     
        lint1 = b5index(20,n)
        lint2 = b5index(21,n)     
        do m = 1, nvar5
            l1r = b5index(14,m) 
            l2r = b5index(15,m)
            l3r = b5index(16,m)
            l4r = b5index(17,m)
            l5r = b5index(18,m)     
            lint1r = b5index(20,m)
            lint2r = b5index(21,m)
            if ((l1r == l1).and.(l2r == l2).and.(l3r == l3).and.(l4r == l4).and.(l5r == l5).and. &
                                & (lint1r == lint1).and.(lint2r == lint2)) then
                if ( b5map_magnetic(n) >= 0) then
                    ! take: increase index at n by 1
                    b5map_magnetic(n) = b5map_magnetic(n) + 1
                    ! burn m: has been found already
                    if (m > n) then
                        b5map_magnetic(m) = -n
                    endif
                    ! as both n and m start counting from 1, this may not happen
                    if (m < n) then
                        print *,'Should not be here.'
                        stop
                    endif              
                endif
            endif
        enddo
    enddo
    m = 0
    l = 0
    do n = 1, nvar5
        ! check if every n has been collected
        if (b5map_magnetic(n) == 0) then
            print *,'This may not happen.'
            stop
        endif 
        if (b5map_magnetic(n) > 0) then
            m = m + 1
            l = l + b5map_magnetic(n)
        endif   
    enddo
    ! number of different angular function vectors l1, l2, l3, l4, l5, lint1, lint2
    nvar5l_magnetic = m
    if ( nvar5 /= l) then
        print *,'Error in checksum for nvar5.'
        stop
    endif

    ! fill b5lindex_magnetic
    allocate( b5lindex_magnetic(1:7,1:nvar5l_magnetic))
    m = 0
    do n = 1, nvar5
        if (b5map_magnetic(n) > 0) then
            m = m + 1
            b5lindex_magnetic(1:5,m) = b5index(14:18,n)
            b5lindex_magnetic(6:7,m) = b5index(20:21,n)
        endif
    enddo
    ! finally update b6map
    m = 0
    do n = 1, nvar5
        if (b5map_magnetic(n) > 0) then
            m = m + 1
            b5map_magnetic(n) = m
        endif
    enddo
    do n = 1, nvar5
        if (b5map_magnetic(n) < 0) then
            m = abs(b5map_magnetic(n))
            b5map_magnetic(n) = b5map_magnetic(m)
        endif
    enddo  
    ! validate
    do n = 1, nvar5
        m = b5map_magnetic(n)
        k = sum(abs(b5lindex_magnetic(1:5,m)-b5index(14:18,n))) + sum(abs(b5lindex_magnetic(6:7,m)-b5index(20:21,n))) 
        if (k > 0 ) then
            print *,'Something wrong with b5lindex_magnetic.'
            stop
        endif
    enddo
  
  
    ! madd5 for atomic part only
    maddmax_atomic= 0
    allocate(madd5_atomic(nvar5l_atomic))
    do k = 1, nvar5l_atomic
        l1 = b5lindex_atomic(1,k)       ! l1
        l2 = b5lindex_atomic(2,k)       ! l2
        l3 = b5lindex_atomic(3,k)       ! l3
        l4 = b5lindex_atomic(4,k)       ! l4
        lint1 = b5lindex_atomic(5,k)   ! L12  (= L34) 
        
        !  finire qui con con i possibili valori di m e m'
        madd5_atomic(k) = 0
        do mint1 = -lint1, lint1                   ! M12 = -(m3+m4) and M = m1+m2
            do m1 = -l1, l1                            ! m1
                do m2 = -l2, l2                         ! m2
                    do m3 = -l3, l3                      ! m3
                        do m4 = -l4, l4                   ! m4
                            if ( (mint1 - (m1+m2) == 0).and.(mint1 + (m3+m4) == 0) ) then
                                madd5_atomic(k) = madd5_atomic(k) + 1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if (madd5_atomic(k) > maddmax_atomic) then
            maddmax_atomic = madd5_atomic(k)
        endif
    enddo
  
    ! madd5 for magnetic part only
    maddmax_magnetic = 0
    allocate(madd5_magnetic(nvar5l_magnetic))
    do k = 1, nvar5l_magnetic
        l5 = b5lindex_magnetic(1,k)       ! l0'
        l6 = b5lindex_magnetic(2,k)       ! l1'
        l7 = b5lindex_magnetic(3,k)       ! l2'
        l8 = b5lindex_magnetic(4,k)       ! l3'
        l9 = b5lindex_magnetic(5,k)       ! l4'  (= L0123')
        lint2 = b5lindex_magnetic(6,k)   ! L01'
        lint3 = b5lindex_magnetic(7,k)   ! L23'
        
        !  finire qui con con i possibili valori di m e m'
        madd5_magnetic(k) = 0
        do mint2 = -lint2, lint2                 ! M01'
            do mint3 = -lint3, lint3              ! M23'
                do m5 = -l5, l5                    ! m0'
                    do m6 = -l6, l6                 ! m1'
                        do m7 = -l7, l7              ! m2'
                            do m8 = -l8, l8           ! m3'
                                do m9 = -l9, l9        ! m4'
                                        mint4 = -m9
                                        if ( (mint2 - (m5+m6) == 0).and.(mint3 - (m7+m8) == 0)         &
                                            &   .and.(mint4 - (mint2+mint3) == 0) ) then
                                            madd5_magnetic(k) = madd5_magnetic(k) + 1
                                        endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if (madd5_magnetic(k) > maddmax_magnetic) then
            maddmax_magnetic = madd5_magnetic(k)
        endif
    enddo
    
    ! b5lm for atomic part only
    allocate( b5lm_atomic(1:nvar5l_atomic,1:maddmax_atomic,1:4) )
    do k = 1, nvar5l_atomic
        l1 = b5lindex_atomic(1,k)       ! l1
        l2 = b5lindex_atomic(2,k)       ! l2
        l3 = b5lindex_atomic(3,k)       ! l3
        l4 = b5lindex_atomic(4,k)       ! l4
        lint1 = b5lindex_atomic(5,k)   ! L12  (= L34) 
        
        n = 0
        do mint1 = -lint1, lint1                ! M12 = -(m3+m4) and M = m1+m2
            do m1 = -l1, l1                      ! m1
                do m2 = -l2, l2                   ! m2
                    do m3 = -l3, l3                ! m3
                        do m4 = -l4, l4             ! m4
                            if ( (mint1 - (m1+m2) == 0).and.(mint1 + (m3+m4) == 0) ) then
                                n = n + 1
                                b5lm_atomic(k,n,1) = m1
                                b5lm_atomic(k,n,2) = m2
                                b5lm_atomic(k,n,3) = m3
                                b5lm_atomic(k,n,4) = m4
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( n/= madd5_atomic(k) ) then
            print*, 'Something wrong here'
            stop
        endif
    enddo
  
    ! b5lm for magnetic part only
    allocate( b5lm_magnetic(1:nvar5l_magnetic,1:maddmax_magnetic,1:5) )
    do k = 1, nvar5l_magnetic
        l5 = b5lindex_magnetic(1,k)       ! l0'
        l6 = b5lindex_magnetic(2,k)       ! l1'
        l7 = b5lindex_magnetic(3,k)       ! l2'
        l8 = b5lindex_magnetic(4,k)       ! l3'
        l9 = b5lindex_magnetic(5,k)       ! l4'  (= L0123')
        lint2 = b5lindex_magnetic(6,k)    ! L01'
        lint3 = b5lindex_magnetic(7,k)    ! L23'
        
        n = 0
        do mint2 = -lint2, lint2                ! M01'
            do mint3 = -lint3, lint3              ! M23'
                do m5 = -l5, l5                ! m0'
                    do m6 = -l6, l6             ! m1'
                        do m7 = -l7, l7          ! m2'
                            do m8 = -l8, l8       ! m3'
                                do m9 = -l9, l9    ! m4'
                                    mint4 = -m9
                                    if ( (mint2 - (m5+m6) == 0).and.(mint3 - (m7+m8) == 0)         &
                                        &   .and.(mint4 - (mint2+mint3) == 0) ) then
                                        n = n + 1
                                        b5lm_magnetic(k,n,1) = m5
                                        b5lm_magnetic(k,n,2) = m6
                                        b5lm_magnetic(k,n,3) = m7
                                        b5lm_magnetic(k,n,4) = m8
                                        b5lm_magnetic(k,n,5) = m9
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( n/= madd5_magnetic(k) ) then
            print*, 'Something wrong here'
            stop
        endif
    enddo

    !---------------------------------------------------------
    ! number of variables for c6
    !---------------------------------------------------------

    if ((2*l6max_atomic > l3max_atomic) .or. (2*l6max_magnetic> l2max_magnetic)) then
        print*, '6-body couplings not worked for lmax > 6'
        print*, 'stopping'
        stop
    endif

     
    if (.not.readsuccess) then
  
    ! find basis functions  
        second = .false.
6100 continue
        if (second) then
            allocate( b6occ_tmp(1:5,1:nvar6_tmp) )
            allocate( b6index_tmp(1:22,1:nvar6_tmp) )
        endif
    
        ! coupling scheme ( (b1*b2)_lint1 (b3*b4)_lint2 )_lint3 b5
        !
        nvar6_tmp = 0
do nele1 = 1, nelements                                                        ! mu_1
    do nele2 = 1, nelements                                                     ! mu_2
        do nele3 = 1, nelements                                                  ! mu_3
            do nele4 = 1, nelements                                               ! mu_4
                do nele5 = 1, nelements                                            ! mu_5
                    do nr1 = 1, nradial6_atomic                                     ! n1
                        do nr2 = 1, nradial6_atomic                                  ! n2
                            do nr3 = 1, nradial6_atomic                               ! n3
                                do nr4 = 1, nradial6_atomic                            ! n4
                                    do nr5 = 1, nradial6_atomic                         ! n5
                                        do nr6 = 1, nradial6_magnetic                       ! n0'
                                            do nr7 = 1, nradial6_magnetic                    ! n1'
                                                do nr8 = 1, nradial6_magnetic                 ! n2'
                                                    do nr9 = 1, nradial6_magnetic              ! n3'
                                                        do nr10 = 1, nradial6_magnetic          ! n4'
                                                            do nr11 = 1, nradial6_magnetic       ! n5'
                                                                do l1 = 0, l6max_atomic      
                                                                    do l2 = 0, l6max_atomic
                                                                        do l3 = 0, l6max_atomic
                                                                            do l4 = 0, l6max_atomic
                                                                                do l5 = 0, l6max_atomic
                                                                                    do l6 = 0, l6max_magnetic
                                                                                        do l7 = 0, l6max_magnetic
                                                                                            do l8 = 0, l6max_magnetic
                                                                                                do l9 = 0, l6max_magnetic
                                                                                                    do l10 = 0, l6max_magnetic
                                                                                                        do l11 = 0, l6max_magnetic

                                                                                
                                    k_atomic1 = nradial6_atomic*l1 + nr1
                                    k_atomic2 = nradial6_atomic*l2 + nr2
                                    k_atomic3 = nradial6_atomic*l3 + nr3
                                    k_atomic4 = nradial6_atomic*l4 + nr4
                                    k_atomic5 = nradial6_atomic*l5 + nr5
                                    
                                    k_magnetic1 = nradial6_magnetic*l7 + nr7
                                    k_magnetic2 = nradial6_magnetic*l8 + nr8
                                    k_magnetic3 = nradial6_magnetic*l9 + nr9
                                    k_magnetic4 = nradial6_magnetic*l10 + nr10
                                    k_magnetic5 = nradial6_magnetic*l11 + nr11
                                    
                                    kmax_atomic = nradial5_atomic*(l5max_atomic + 1)
                                    kmax_magnetic = nradial5_magnetic*(l5max_magnetic + 1)
                                
                                    nc1 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele1 +  &
                                    &  (kmax_magnetic+1)*k_atomic1 + k_magnetic1
                                    nc2 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele2 +  &
                                    &  (kmax_magnetic+1)*k_atomic2 + k_magnetic2 
                                    nc3 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele3 +  &
                                    &  (kmax_magnetic+1)*k_atomic3 + k_magnetic3 
                                    nc4 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele4 +  &
                                    &  (kmax_magnetic+1)*k_atomic4 + k_magnetic4     
                                    nc5 = ((kmax_magnetic+1)*kmax_atomic + kmax_magnetic)*nele5 +  &
                                    &  (kmax_magnetic+1)*k_atomic5 + k_magnetic5  
                                    
                                    index = .true.

                                    ! take only ordered sequence
                                    if (( nc1 < nc2 ).or.( nc2 < nc3 ).or.( nc3 < nc4 ).or.( nc4 < nc5 )) then
                                        index = .false.
                                    endif

                                    if (index) then
                                        if ((k_atomic1 .lt. k_atomic2).or.(k_atomic2 .lt. k_atomic3).or.  &
                                        &   (k_atomic3 .lt. k_atomic4).or.(k_atomic4 .lt. k_atomic5)) then
                                            print*, 'Error! the atomic basis is not ordered lexicographically'
                                            stop
                                        endif
                                    endif                             

                                    ! now everything should be o.k. for coupling                                      
                                    if (index) then
                                        nvar6_tmp = nvar6_tmp + 1
                                        if (second) then
                                            b6occ_tmp(1,nvar6_tmp) = nele1   ! mu_1
                                            b6occ_tmp(2,nvar6_tmp) = nele2   ! mu_2
                                            b6occ_tmp(3,nvar6_tmp) = nele3   ! mu_3
                                            b6occ_tmp(4,nvar6_tmp) = nele4   ! mu_4
                                            b6occ_tmp(5,nvar6_tmp) = nele5   ! mu_5 
                                            b6index_tmp(1,nvar6_tmp) = nr1   ! n1
                                            b6index_tmp(2,nvar6_tmp) = nr2   ! n2
                                            b6index_tmp(3,nvar6_tmp) = nr3   ! n3
                                            b6index_tmp(4,nvar6_tmp) = nr4   ! n4
                                            b6index_tmp(5,nvar6_tmp) = nr5   ! n5
                                            b6index_tmp(6,nvar6_tmp) = nr6    ! n0'
                                            b6index_tmp(7,nvar6_tmp) = nr7    ! n1'
                                            b6index_tmp(8,nvar6_tmp) = nr8    ! n2'
                                            b6index_tmp(9,nvar6_tmp) = nr9    ! n3'
                                            b6index_tmp(10,nvar6_tmp) = nr10   ! n4'
                                            b6index_tmp(11,nvar6_tmp) = nr11  ! n5'
                                            b6index_tmp(12,nvar6_tmp) = l1   ! l1
                                            b6index_tmp(13,nvar6_tmp) = l2   ! l2
                                            b6index_tmp(14,nvar6_tmp) = l3   ! l3
                                            b6index_tmp(15,nvar6_tmp) = l4   ! l4
                                            b6index_tmp(16,nvar6_tmp) = l5   ! l5
                                            b6index_tmp(17,nvar6_tmp) = l6   ! l0'
                                            b6index_tmp(18,nvar6_tmp) = l7   ! l1'
                                            b6index_tmp(19,nvar6_tmp) = l8   ! l2'
                                            b6index_tmp(20,nvar6_tmp) = l9   ! l3'
                                            b6index_tmp(21,nvar6_tmp) = l10  ! l4'
                                            b6index_tmp(22,nvar6_tmp) = l11  ! l5'
                                        endif
                                    endif                     
                                                          
                                                                                
                                                                                                        enddo
                                                                                                    enddo
                                                                                                enddo
                                                                                            enddo
                                                                                        enddo
                                                                                    enddo
                                                                                enddo
                                                                            enddo
                                                                        enddo
                                                                    enddo
                                                                enddo
                                                            enddo
                                                        enddo
                                                    enddo
                                                enddo
                                            enddo
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo
        if(.not.second) then
            second = .true.
            goto 6100
        endif
  
        second = .false.
6200 continue
        if (second) then
            allocate( b6occ(1:5,1:nvar6) )
            allocate( b6index(1:27,1:nvar6) )
            allocate( b6nint(1:nvar6) )
            b6nint = 0
        endif
  
        nrank = 5
        nmax = 200
        if (.not.second) then
            allocate( vecn_atomic(1:nrank) )
            allocate( vecn_magnetic(1:nrank+1) )
            allocate( vecnr_atomic(1:nrank) )
            allocate( vecnr_magnetic(1:nrank+1) )
            allocate( vecl_atomic(1:nrank) )
            allocate( vecl_magnetic(1:nrank+1) )
            allocate( arrlint_atomic(1:nmax,1:nrank) )
            allocate( arrlint_magnetic(1:nmax,1:nrank+1) )
        endif
  
        nvar6 = 0
        do nvar = 1, nvar6_tmp
            
            vecn_atomic(1) = b6occ_tmp(1,nvar)        ! mu_1
            vecn_atomic(2) = b6occ_tmp(2,nvar)        ! mu_2 
            vecn_atomic(3) = b6occ_tmp(3,nvar)        ! mu_3
            vecn_atomic(4) = b6occ_tmp(4,nvar)        ! mu_4
            vecn_atomic(5) = b6occ_tmp(5,nvar)        ! mu_5

            vecnr_atomic(1) = b6index_tmp(1,nvar)     ! n1
            vecnr_atomic(2) = b6index_tmp(2,nvar)     ! n2
            vecnr_atomic(3) = b6index_tmp(3,nvar)     ! n3
            vecnr_atomic(4) = b6index_tmp(4,nvar)     ! n4
            vecnr_atomic(5) = b6index_tmp(5,nvar)     ! n5
            vecl_atomic(1) = b6index_tmp(12,nvar)     ! l1  
            vecl_atomic(2) = b6index_tmp(13,nvar)     ! l2
            vecl_atomic(3) = b6index_tmp(14,nvar)     ! l3
            vecl_atomic(4) = b6index_tmp(15,nvar)     ! l4
            vecl_atomic(5) = b6index_tmp(16,nvar)     ! l5

            call bcreate(atomic,nrank,nmax,vecn_atomic,vecnr_atomic,vecl_atomic,numlint_atomic,  &
            &   arrlint_atomic,isordered_atomic,isbasis_atomic)
            
            if(.not.isordered_atomic) then
                print*, 'Something wrong with rank = 5.'
                stop
            endif
            
            vecn_magnetic(1) = vecn_atomic(1)               ! mu_0
            vecn_magnetic(2) = vecn_atomic(1)               ! mu_1 
            vecn_magnetic(3) = vecn_atomic(2)               ! mu_2
            vecn_magnetic(4) = vecn_atomic(3)               ! mu_3
            vecn_magnetic(5) = vecn_atomic(4)               ! mu_4
            vecn_magnetic(6) = vecn_atomic(5)               ! mu_5
                
            vecnr_magnetic(1) = b6index_tmp(6,nvar)     ! n0'
            vecnr_magnetic(2) = b6index_tmp(7,nvar)     ! n1'
            vecnr_magnetic(3) = b6index_tmp(8,nvar)     ! n2'
            vecnr_magnetic(4) = b6index_tmp(9,nvar)     ! n3'
            vecnr_magnetic(5) = b6index_tmp(10,nvar)    ! n4'
            vecnr_magnetic(6) = b6index_tmp(11,nvar)    ! n5'
            vecl_magnetic(1) = b6index_tmp(17,nvar)    ! l0'
            vecl_magnetic(2) = b6index_tmp(18,nvar)    ! l1'
            vecl_magnetic(3) = b6index_tmp(19,nvar)    ! l2' 
            vecl_magnetic(4) = b6index_tmp(20,nvar)    ! l3'
            vecl_magnetic(5) = b6index_tmp(21,nvar)    ! l4'
            vecl_magnetic(6) = b6index_tmp(22,nvar)    ! l5'

            call bcreate(magnetic,nrank+1,nmax,vecn_magnetic,vecnr_magnetic,vecl_magnetic,numlint_magnetic,  &   
            &  arrlint_magnetic,isordered_magnetic,isbasis_magnetic) 
            
            if( isbasis_atomic.and.isbasis_magnetic ) then
                        
                numlint = max(numlint_atomic,numlint_magnetic)
                do nint = 1, numlint
                    nvar6 = nvar6 + 1
                    if (second) then
                        b6occ(1:5,nvar6) = vecn_atomic(1:5)
                        if ( nint.le.numlint_atomic ) then
                            b6index(1:5,nvar6) = vecnr_atomic(1:5)
                            b6index(12:16,nvar6) = vecl_atomic(1:5)
                            b6index(23:24,nvar6) = arrlint_atomic(nint,1:2)
                        else
                            b6index(1:5,nvar6) = 0
                            b6index(12:16,nvar6) = 0
                            b6index(23:24,nvar6) = 0
                        endif
                        if ( nint.le.numlint_magnetic ) then
                            b6index(6:11,nvar6) = vecnr_magnetic(1:6)
                            b6index(17:22,nvar6) = vecl_magnetic(1:6)
                            b6index(25:27,nvar6) = arrlint_magnetic(nint,1:3)
                        else
                            b6index(6:11,nvar6) = 0
                            b6index(17:22,nvar6) = 0
                            b6index(25:27,nvar6) = 0
                        endif
                        if (nint == 1) then
                            b6nint(nvar6) = numlint
                        else
                            b6nint(nvar6) = 0
                        endif
                    endif
                enddo
                
            endif
            
        enddo
        if(.not.second) then
            second = .true.
            goto 6200
        endif

        if (second) then
            deallocate( vecn_atomic )
            deallocate( vecn_magnetic )
            deallocate( vecnr_atomic )
            deallocate( vecnr_magnetic )
            deallocate( vecl_atomic )
            deallocate( vecl_magnetic )
            deallocate( arrlint_atomic )
            deallocate( arrlint_magnetic )
        endif  

        ! check for nint 1
        do n = 1, nvar6
            do m = 1, b6nint(nvar6) - 1
                if (b6nint(n+m) /= 0 ) then
                    print *,'Something wrong with order of internal couplings in b6nint.'
                    stop
                endif
            enddo
        enddo
  
        ! ordering internal couplings to take smallest lint1+lint2 first (only atomic lints)
        do n = 1, nvar6
            do m = 1, b6nint(n) - 2
                do k = m+1, b6nint(n) - 1
                    lsum1 = b6index(23,n+m) + b6index(24,n+m)
                    lsum2 = b6index(23,n+k) + b6index(24,n+k)
                    if (lsum1 > lsum2) then
                        n2v(1:2) = b6index(23:24,n+m)
                        b6index(23:24,n+m) = b6index(23:24,n+k)
                        b6index(23:24,n+k) = n2v(1:2)
                    endif
                enddo
            enddo
        enddo
        ! ordering internal couplings to take smallest lint3 + lint4 + lint5 first (only magnetic lints)
        ! capire cosa sarebbe sta roba che si scambiano di posizione i lints????????
        do n = 1, nvar6
            do m = 1, b6nint(n) - 2
                do k = m+1, b6nint(n) - 1
                    lsum1 =  b6index(25,n+m) + b6index(26,n+m) + b6index(27,n+m)
                    lsum2 =  b6index(25,n+k) + b6index(26,n+k) + b6index(27,n+k)
                    if (lsum1 > lsum2) then
                        n3v(1:3) =  b6index(25:27,n+m)
                        b6index(25:27,n+m) = b6index(25:27,n+k)
                        b6index(25:27,n+k) = n3v(1:3)
                    endif
                enddo
            enddo
        enddo

    endif
    
    !  if (verblevel > 0 ) then
    write(*,*) 'number of variables for 6-body clusters:', nvar6
    !  endif
  
    ! b6map for magnetic part only
    allocate( b6map_atomic(1:nvar6) )
    b6map_atomic = 0
 
    do n = 1, nvar6
        l1 = b6index(12,n) 
        l2 = b6index(13,n)
        l3 = b6index(14,n)
        l4 = b6index(15,n)
        l5 = b6index(16,n)     
        lint1 = b6index(23,n)
        lint2 = b6index(24,n)     
        do m = 1, nvar6
            l1r = b6index(12,m) 
            l2r = b6index(13,m)
            l3r = b6index(14,m)
            l4r = b6index(15,m)
            l5r = b6index(16,m)     
            lint1r = b6index(23,m)
            lint2r = b6index(24,m)
            if ((l1r == l1).and.(l2r == l2).and.(l3r == l3).and.(l4r == l4).and.(l5r == l5).and. &
                & (lint1r == lint1).and.(lint2r == lint2)) then
                if ( b6map_atomic(n) >= 0) then
                    ! take: increase index at n by 1
                    b6map_atomic(n) = b6map_atomic(n) + 1
                    ! burn m: has been found already
                    if (m > n) then
                        b6map_atomic(m) = -n
                    endif
                    ! as both n and m start counting from 1, this may not happen
                    if (m < n) then
                        print *,'Should not be here.'
                        stop
                    endif              
                endif
            endif
        enddo
    enddo
    m = 0
    l = 0
    do n = 1, nvar6
        ! check if every n has been collected
        if (b6map_atomic(n) == 0) then
            print *,'This may not happen.'
            stop
        endif 
        if (b6map_atomic(n) > 0) then
            m = m + 1
            l = l + b6map_atomic(n)
        endif   
    enddo
    ! number of different angular function vectors l1, l2, l3, l4, l5, lint1, lint2
    nvar6l_atomic = m
    if ( nvar6 /= l) then
        print *,'Error in checksum for nvar6.'
        stop
    endif

    ! fill b6lindex
    allocate( b6lindex_atomic(1:7,1:nvar6l_atomic))
    m = 0
    do n = 1, nvar6
        if (b6map_atomic(n) > 0) then
            m = m + 1
            b6lindex_atomic(1:5,m) = b6index(12:16,n)
            b6lindex_atomic(6:7,m) = b6index(23:24,n)
        endif
    enddo
    ! finally update b6map
    m = 0
    do n = 1, nvar6
        if (b6map_atomic(n) > 0) then
            m = m + 1
            b6map_atomic(n) = m
        endif
    enddo
    do n = 1, nvar6
        if (b6map_atomic(n) < 0) then
            m = abs(b6map_atomic(n))
            b6map_atomic(n) = b6map_atomic(m)
        endif
    enddo  
    ! validate
    do n = 1, nvar6
        m = b6map_atomic(n)
        k = sum(abs(b6lindex_atomic(1:5,m)-b6index(12:16,n))) + &
        &  sum(abs(b6lindex_atomic(6:7,m)-b6index(23:24,n)))
        if (k > 0 ) then
            print *,'Something wrong with b6lindex.'
            stop
        endif
    enddo
  
  
    ! b6map for magnetic part
    allocate( b6map_magnetic(1:nvar6) )
    b6map_magnetic = 0
    
    do n = 1, nvar6
        l1 = b6index(17,n) 
        l2 = b6index(18,n)
        l3 = b6index(19,n)
        l4 = b6index(20,n)
        l5 = b6index(21,n)
        l6 = b6index(22,n)
        lint1 = b6index(25,n)
        lint2 = b6index(26,n)     
        lint3 = b6index(27,n)
        do m = 1, nvar6
            l1r = b6index(17,m) 
            l2r = b6index(18,m)
            l3r = b6index(19,m)
            l4r = b6index(20,m)
            l5r = b6index(21,m)
            l6r = b6index(22,m)
            lint1r = b6index(25,m)
            lint2r = b6index(26,m)
            lint3r = b6index(27,m)
            if ((l1r == l1).and.(l2r == l2).and.(l3r == l3).and.(l4r == l4).and.(l5r == l5).and.(l6r == l6).and. &
                & (lint1r == lint1).and.(lint2r == lint2).and.(lint3r == lint3)) then
                if ( b6map_magnetic(n) >= 0) then
                    ! take: increase index at n by 1
                    b6map_magnetic(n) = b6map_magnetic(n) + 1
                    ! burn m: has been found already
                    if (m > n) then
                        b6map_magnetic(m) = -n
                    endif
                    ! as both n and m start counting from 1, this may not happen
                    if (m < n) then
                        print *,'Should not be here.'
                        stop
                    endif              
                endif
            endif
        enddo
    enddo
    m = 0
    l = 0
    do n = 1, nvar6
        ! check if every n has been collected
        if (b6map_magnetic(n) == 0) then
            print *,'This may not happen.'
            stop
        endif 
        if (b6map_magnetic(n) > 0) then
            m = m + 1
            l = l + b6map_magnetic(n)
        endif   
    enddo
    ! number of different angular function vectors l1, l2, l3, l4, l5, l6, lint1, lint2, lint3
    nvar6l_magnetic = m
    if ( nvar6 /= l) then
        print *,'Error in checksum for nvar6.'
        stop
    endif


    ! fill b6lindex
    allocate( b6lindex_magnetic(1:9,1:nvar6l_magnetic))
    m = 0
    do n = 1, nvar6
        if (b6map_magnetic(n) > 0) then
            m = m + 1
            b6lindex_magnetic(1:6,m) = b6index(17:22,n)
            b6lindex_magnetic(7:9,m) = b6index(25:27,n)
        endif
    enddo
    ! finally update b7map
    m = 0
    do n = 1, nvar6
        if (b6map_magnetic(n) > 0) then
            m = m + 1
            b6map_magnetic(n) = m
        endif
    enddo
    do n = 1, nvar6
        if (b6map_magnetic(n) < 0) then
            m = abs(b6map_magnetic(n))
            b6map_magnetic(n) = b6map_magnetic(m)
        endif
    enddo  
    ! validate
    do n = 1, nvar6
        m = b6map_magnetic(n)
        k = sum(abs(b6lindex_magnetic(1:6,m)-b6index(17:22,n))) + &
        & sum(abs(b6lindex_magnetic(7:9,m)-b6index(25:27,n))) 
        if (k > 0 ) then
            print *,'Something wrong with b6lindex_magnetic.'
            stop
        endif
    enddo
  
    !collect non-zero summation contributions for m0', m1', m2', m3', m1, m2, m3, M'
    maddmax_atomic = 0
    allocate(madd6_atomic(nvar6l_atomic))
    do k = 1, nvar6l_atomic
        l1 = b6lindex_atomic(1,k)       ! l1
        l2 = b6lindex_atomic(2,k)       ! l2
        l3 = b6lindex_atomic(3,k)       ! l3
        l4 = b6lindex_atomic(4,k)       ! l4
        l5 = b6lindex_atomic(5,k)       ! l5
        lint1 = b6lindex_atomic(6,k)   ! L12   
        lint2 = b6lindex_atomic(7,k)   ! L34
        
        madd6_atomic(k) = 0
        do mint1 = -lint1, lint1                    ! M12 
            do mint2 = -lint2, lint2                 ! M34
                do m1 = -l1, l1                            ! m1
                    do m2 = -l2, l2                         ! m2
                        do m3 = -l3, l3                      ! m3
                            do m4 = -l4, l4                   ! m4
                                do m5 = -l5, l5                ! m5

                                    mint1234 = -m5     ! M1234  (=-m5) 
                                    
                                    if ( (mint1 - (m1+m2) == 0).and.(mint2 - (m3+m4) == 0)                 &
                                        &   .and.(mint1234 - (mint1+mint2) == 0) ) then
                                        madd6_atomic(k) = madd6_atomic(k) + 1
                                    endif 

                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if (madd6_atomic(k) > maddmax_atomic) then
            maddmax_atomic = madd6_atomic(k)
        endif
    enddo
  
    maddmax_magnetic = 0
    allocate(madd6_magnetic(nvar6l_magnetic))
    do k = 1, nvar6l_magnetic
        l6 = b6lindex_magnetic(1,k)       ! l0'
        l7 = b6lindex_magnetic(2,k)       ! l1'
        l8 = b6lindex_magnetic(3,k)       ! l2'
        l9 = b6lindex_magnetic(4,k)       ! l3'
        l10 = b6lindex_magnetic(5,k)     ! l4'
        l11 = b6lindex_magnetic(6,k)     ! l5' 
        lint3 = b6lindex_magnetic(7,k)   ! L01'
        lint4 = b6lindex_magnetic(8,k)   ! L23'
        lint5 = b6lindex_magnetic(9,k)   ! L45'
        
    !      madd6_magnetic(k) = 0
        n = 0
        do mint3 = -lint3, lint3              ! M01'
            do mint4 = -lint4, lint4           ! M23'
                do mint5 = -lint5, lint5        ! M45'
                    do m6 = -l6, l6                     ! m0'
                        do m7 = -l7, l7                  ! m1'
                            do m8 = -l8, l8              ! m2'
                                do m9 = -l9, l9           ! m3'
                                    do m10 = -l10, l10     ! m4'
                                        do m11 = -l11, l11  ! m5' 

                                            mint0123 = -mint5  ! M0123' (=-M45')

                                            if ( (mint3 - (m6+m7) == 0).and.(mint4 - (m8+m9) == 0)       &
                                                &   .and.(mint5 - (m10+m11) == 0)                        &
                                                &   .and.(mint0123 - (mint3+mint4) == 0) ) then
    !                                             madd6_magnetic(k) = madd6_magnetic(k) + 1
                                                n = n + 1
                                            endif 

                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    !      if (madd6_magnetic(k) > maddmax_magnetic) then
        if ( n > maddmax_magnetic ) then 
    !         maddmax_magnetic = madd6_magnetic(k)
            maddmax_magnetic = n
        endif
        madd6_magnetic(k) = n
    enddo
    
    allocate( b6lm_atomic(1:nvar6l_atomic,1:maddmax_atomic,1:5) )
    do k = 1, nvar6l_atomic
        l1 = b6lindex_atomic(1,k)       ! l1
        l2 = b6lindex_atomic(2,k)       ! l2
        l3 = b6lindex_atomic(3,k)       ! l3
        l4 = b6lindex_atomic(4,k)       ! l4
        l5 = b6lindex_atomic(5,k)       ! l5
        lint1 = b6lindex_atomic(6,k)   ! L12   
        lint2 = b6lindex_atomic(7,k)   ! L34
        
        n = 0
        do mint1 = -lint1, lint1                    ! M12 
            do mint2 = -lint2, lint2                 ! M34
                do m1 = -l1, l1                            ! m1
                    do m2 = -l2, l2                         ! m2
                        do m3 = -l3, l3                      ! m3
                            do m4 = -l4, l4                   ! m4
                                do m5 = -l5, l5                ! m5

                                    mint1234 = -m5     ! M1234  (=-m5) 

                                    if ( (mint1 - (m1+m2) == 0).and.(mint2 - (m3+m4) == 0)                 &
                                                    &   .and.(mint1234 - (mint1+mint2) == 0) ) then
                                        n = n + 1
                                        b6lm_atomic(k,n,1) = m1    
                                        b6lm_atomic(k,n,2) = m2
                                        b6lm_atomic(k,n,3) = m3
                                        b6lm_atomic(k,n,4) = m4
                                        b6lm_atomic(k,n,5) = m5
                                    endif

                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( n/= madd6_atomic(k) ) then
            print*, 'Something wrong here'
            stop
        endif
    enddo
  
    allocate( b6lm_magnetic(1:nvar6l_magnetic,1:maddmax_magnetic,1:6) )
    do k = 1, nvar6l_magnetic
        l6 = b6lindex_magnetic(1,k)       ! l0'
        l7 = b6lindex_magnetic(2,k)       ! l1'
        l8 = b6lindex_magnetic(3,k)       ! l2'
        l9 = b6lindex_magnetic(4,k)       ! l3'
        l10 = b6lindex_magnetic(5,k)     ! l4'
        l11 = b6lindex_magnetic(6,k)     ! l5' 
        lint3 = b6lindex_magnetic(7,k)   ! L01'
        lint4 = b6lindex_magnetic(8,k)   ! L23'
        lint5 = b6lindex_magnetic(9,k)   ! L45'
        
        n = 0
        do mint3 = -lint3, lint3              ! M01'
            do mint4 = -lint4, lint4           ! M23'
                do mint5 = -lint5, lint5        ! M45'
                    do m6 = -l6, l6                     ! m0'
                        do m7 = -l7, l7                  ! m1'
                            do m8 = -l8, l8              ! m2'
                                do m9 = -l9, l9           ! m3'
                                    do m10 = -l10, l10     ! m4'
                                        do m11 = -l11, l11  ! m5' 

                                            mint0123 = -mint5  ! M0123' (=-M45')

                                            if ( (mint3 - (m6+m7) == 0).and.(mint4 - (m8+m9) == 0)                 &
                                                    &   .and.(mint5 - (m10+m11) == 0)                                  &
                                                    &   .and.(mint0123 - (mint3+mint4) == 0)                           &
                                                    &   .and.(mint0123 + (m10+m11) == 0)  ) then
                                                n = n + 1
                                                b6lm_magnetic(k,n,1) = m6
                                                b6lm_magnetic(k,n,2) = m7
                                                b6lm_magnetic(k,n,3) = m8
                                                b6lm_magnetic(k,n,4) = m9
                                                b6lm_magnetic(k,n,5) = m10
                                                b6lm_magnetic(k,n,6) = m11
                                            endif

                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if ( n/= madd6_magnetic(k) ) then
            print*, 'Something wrong here'
            stop
        endif
    enddo

    !---------------------------------------------------------------------------------------
   
    ! set up table of Clebsch-Gordan coefficients and Wigner3j symbols
  
    ! set up tables for higher couplings
    call setuphighercouplings
  
end subroutine setuptables

!----------------------------------------------------------------------------------------------------------

subroutine setupcg

    use tables
    use functionparameters
    use core
    implicit none

    integer J, j1, j2
    integer M, m1, m2
    integer k, n, l1, l2, l3
    integer :: l0_p, l1_p, l2_p, l3_p, m0_p, m1_p, m2_p, m3_p
    integer :: L01_p, l23_p, L, size_cg    

    integer j1pj2mJ, j1mj2pJ, j2mj1pJ, j1pj2pJp1, JpM, JmM, j1pm1, j1mm1, j2pm2, j2mm2
    integer j1pj2mJmk, j1mm1mk, j2pm2mk, Jmj2pm1pk, Jmj1mm2pk

    double precision prefac, xsum

    ! for checking
    double precision xw3j
    character(20) :: str1, str2, str3, str4
    integer m3, j3

    allocate( cgcoeff(0:2*lmax,0:2*lmax,0:4*lmax,-2*lmax:2*lmax,-2*lmax:2*lmax,-4*lmax:4*lmax) )

    cgcoeff = 0.d0

    do J = 0, 4*lmax
        do j1 = 0, 2*lmax
            do j2 = 0, 2*lmax
                
                ! basic rules for angular momentum addition   | j1 - j2 | <= J <= j1 + j2
                if ( ( j1+j2 >= J ) .and. ( abs(j1-j2) <= J ) ) then
                
                    do m1 = -j1, j1
                        do m2 = -j2, j2
                        
                            !only contributions if M = m1 + m2 and |M| <= J
                            M = m1 + m2
                        
                            if ( abs(M) <= J) then
                                
                                j1pj2mJ   =  j1 + j2 - J
                                j1mj2pJ   =  j1 - j2 + J
                                j2mj1pJ   = -j1 + j2 + J
                                j1pj2pJp1 =  j1 + j2 + J + 1
                                JpM       =  J + M
                                JmM       =  J - M
                                j1pm1     =  j1 + m1
                                j1mm1     =  j1 - m1
                                j2pm2     =  j2 + m2
                                j2mm2     =  j2 - m2
                                
                                ! \sqrt{\frac{(2J+1)(j_1+j_2-J)!(j_1-j_2+J)!(-j_1+j_2+J)!}{(j_1+j_2+J+1)!}}
                                prefac = ( dble(2*J+1)*fac(j1pj2mJ)*fac(j1mj2pJ)*fac(j2mj1pJ)/fac(j1pj2pJp1) ) &
                                    ! \sqrt{(J+M)!(J-M)!(j_1+m_1)!(j_1-m_1)!(j_2+m_2)!(j_2-m_2)!}      
                                    &  *( fac(JpM)*fac(JmM)*fac(j1pm1)*fac(j1mm1)*fac(j2pm2)*fac(j2mm2) ) 

                                if ( prefac < 0.d0 ) then
                                    print *,'Error in calculation of Clebsch-Gordan coefficients. Stopping.'
                                    stop
                                endif
                                
                                prefac = sqrt(prefac)
                                
                                xsum = 0.d0
                                do k = 0, 4*lmax
                                    
                                    j1pj2mJmk =  j1 + j2 - J - k
                                    j1mm1mk   =  j1 - m1 - k
                                    j2pm2mk   =  j2 + m2 - k
                                    Jmj2pm1pk =  J - j2 + m1 + k
                                    Jmj1mm2pk =  J - j1 - m2 + k
                                    
                                    if ( (j1pj2mJmk >= 0) .and. &
                                            & (j1mm1mk >= 0) .and. &
                                            & (j2pm2mk >= 0) .and. &
                                            & (Jmj2pm1pk >= 0) .and. &
                                            & (Jmj1mm2pk >= 0) ) then
                                
                                        ! &\sum_z \frac{(-1)^z}{z!(j_1+j_2-J-z)!(j_1-m_1-z)!(j_2+m_2-z)!(J-j_2+m_1+z)!(J-j_1-m_2+z)!} \,,
                                        xsum = xsum + &
                                            & dble( (-1)**k ) &
                                            & /( fac(k)*fac(j1pj2mJmk)*fac(j1mm1mk)*fac(j2pm2mk)*fac(Jmj2pm1pk)*fac(Jmj1mm2pk) )
                                        
                                    endif ! relevant contributions to sum
                                enddo ! sum k
                            
                                cgcoeff(j1,j2,J,m1,m2,M) = prefac*xsum
                                !w3j(j1,j2,J,m1,m2,-M) = dble( (-1)**(j1 - j2 + M) )/sqrt(dble( 2*J + 1))*cgcoeff(j1,j2,J,m1,m2,M)

                            endif ! M boundaries
                        
                        enddo ! m2
                    enddo ! m1
                
                endif  ! J boundaries
                
            enddo ! j2
        enddo ! j1
    enddo ! J

    
end subroutine setupcg

!----------------------------------------------------------------------------------------------------------

subroutine setuphighercouplings
      
    use functionparameters
    use tables
    implicit none
    integer low(2), up(2)
    integer n, l1, l2, l3, l4, l5, l6, l7, L, Lc1
    integer m0_p, m1_p, m2_p, m3_p, m4_p, m5_p
    integer L01_p, l23_p
    integer m1, m2, m3, m4, m5, m6, M, Mc1, k, M34
    integer lint1, lint2, lint3, lint1234
    integer mint1, mint2, mint3, mint1234
    integer l0_p, l1_p, l2_p, l3_p, l4_p, l5_p, lint12, lint34, lint01_p, lint23_p, lint45_p
    integer M12, M01_p, M23_p, M45_p, size_cg

    !redefinition of the three-body contribution with Clebsch-Gordan
    !-------------------------------------------------------------------
    ! one body-term (purely magnetic), this is not present in the normal
    ! ace, because there is no coupling at first-order term 
    !-------------------------------------------------------------------

    allocate( cg01_magnetic(1:nvar2l,-lmax:lmax) )
    cg01_magnetic = 0.d0
    
    do k = 1, nvar2l
        l0_p = b2lindex(k)       ! l0' (=l1')
        do n = 1, madd2(k)
            m0_p = b2lm(k,n)
            cg01_magnetic(k,m0_p) = cg01_magnetic(k,m0_p) + cgcoeff(l0_p,l0_p,0,m0_p,-m0_p,0)
        enddo
    enddo

    !--------------------------------------------------------------------
    !  two body term (the atomic part is coupled to zero)
    !--------------------------------------------------------------------

    allocate( cg02_atomic(1:nvar3l_atomic,-lmax:lmax) )
    allocate( cg02_magnetic(1:nvar3l_magnetic,-lmax:lmax,-lmax:lmax,-lmax:lmax) )
    cg02_atomic = 0.d0
    cg02_magnetic = 0.d0

    do k = 1, nvar3l_atomic
        l1 = b3lindex_atomic(k)    ! l1  (=l2)
        do n = 1, madd3_atomic(k)
            m1 = b3lm_atomic(k,n)     ! m1 (=-m2)            
                        
            cg02_atomic(k,m1) = cg02_atomic(k,m1) + cgcoeff(l1,l1,0,m1,-m1,0)

            if ( (dabs( cg02_atomic(k,m1) ) > 1.d-12).and.(m1.gt.l1)) then
                print *,'Something funny with cg03_atomic.'
                stop
            endif
        enddo
    enddo
        
    do k = 1, nvar3l_magnetic
        l0_p = b3lindex_magnetic(1,k)  ! l0'
        l1_p = b3lindex_magnetic(2,k)  ! l1'
        l2_p = b3lindex_magnetic(3,k)  ! l2' (=L01')
        do n = 1, madd3_magnetic(k)
            m0_p = b3lm_magnetic(k,n,1)   ! m0'                  
            m1_p = b3lm_magnetic(k,n,2)   ! m1'                  
            m2_p = b3lm_magnetic(k,n,3)   ! m2' (=-(m0'+m1'))    
                        
            cg02_magnetic(k,m0_p,m1_p,m2_p) = cg02_magnetic(k,m0_p,m1_p,m2_p) + &
                        &  cgcoeff(l0_p,l1_p,l2_p,m0_p,m1_p,-m2_p)* cgcoeff(l2_p,l2_p,0,m2_p,-m2_p,0)
                    
            if ( (dabs( cg02_magnetic(k,m0_p,m1_p,m2_p) ) > 1.d-12).and. (m0_p + m1_p + m2_p /= 0)) then
                print *,'Something funny with cg03.'
                stop
            endif
            if ( (dabs( cgcoeff(l0_p,l1_p,l2_p,m0_p,m1_p,m2_p) ) > 1.d-12) .and. &
                    &       (m0_p + m1_p - m2_p /= 0)) then
                print *,'Something funny with cgcoeff.'
                stop
            endif 
        enddo
    enddo

    !--------------------------------------------------------------------
    !  three body term 
    !--------------------------------------------------------------------
        
    allocate( cg03_atomic(1:nvar4l_atomic,-lmax:lmax,-lmax:lmax,-lmax:lmax) )
    allocate( cg03_magnetic(1:nvar4l_magnetic,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax) )
    cg03_atomic = 0.d0
    cg03_magnetic = 0.d0

    do k = 1, nvar4l_atomic
        l1 = b4lindex_atomic(1,k)   ! l1
        l2 = b4lindex_atomic(2,k)   ! l2
        l3 = b4lindex_atomic(3,k)   ! l3
        do n = 1, madd4_atomic(k)
            m1 = b4lm_atomic(k,n,1)     ! m1                  
            m2 = b4lm_atomic(k,n,2)     ! m2                  
            m3 = b4lm_atomic(k,n,3)     ! m3      
                
            cg03_atomic(k,m1,m2,m3) = cg03_atomic(k,m1,m2,m3) + cgcoeff(l1,l2,l3,m1,m2,-m3)*cgcoeff(l3,l3,0,m3,-m3,0)
        
            if ( (dabs( cg03_atomic(k,m1,m2,m3) ) > 1.d-12).and.(m1 + m2 + m3 /= 0) ) then
                print *,'Something funny with cg03_atomic.'
                stop
            endif
            if ( (dabs( cgcoeff(l1,l2,l3,m1,m2,m3) ) > 1.d-12) .and. (m1 + m2 - m3 /= 0)) then
                print *,'Something funny with cgcoeff.'
                stop
            endif                
        enddo
    enddo
        
    do k = 1, nvar4l_magnetic
        l0_p = b4lindex_magnetic(1,k)   ! l0'
        l1_p = b4lindex_magnetic(2,k)   ! l1'
        l2_p = b4lindex_magnetic(3,k)   ! l2'
        l3_p = b4lindex_magnetic(4,k)   ! l3'
        L  = b4lindex_magnetic(5,k)   ! L01'
        do n = 1, madd4_magnetic(k)
            m0_p = b4lm_magnetic(k,n,1)   ! m0'       
            m1_p = b4lm_magnetic(k,n,2)   ! m1'     
            m2_p = b4lm_magnetic(k,n,3)   ! m2'     
            m3_p = b4lm_magnetic(k,n,4)   ! m3'     

            M = m0_p + m1_p
                
            cg03_magnetic(k,m0_p,m1_p,m2_p,m3_p) = cg03_magnetic(k,m0_p,m1_p,m2_p,m3_p) +  & 
                    &    cgcoeff(l0_p,l1_p,L,m0_p,m1_p,M)*cgcoeff(l2_p,l3_p,L,m2_p,m3_p,-M)*cgcoeff(L,L,0,M,-M,0)

            if ( (dabs( cg03_magnetic(k,m0_p,m1_p,m2_p,m3_p) ) > 1.d-12) .and. (m0_p + m1_p + m2_p + m3_p /= 0) ) then
                print *,'Something funny with cg03_magnetic.'
                stop
            endif                
        enddo
    enddo
      
    allocate( cg04_atomic(1:nvar5l_atomic,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax) )
    allocate( cg04_magnetic(1:nvar5l_magnetic,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax) )
    cg04_atomic = 0.d0
    cg04_magnetic = 0.d0

    !redefinition of the last coupling with the Clebsch-Gordan coupled to zero
    ! initialize the required couplings
    do k = 1, nvar5l_atomic
        l1 = b5lindex_atomic(1,k)       ! l1
        l2 = b5lindex_atomic(2,k)       ! l2
        l3 = b5lindex_atomic(3,k)       ! l3
        l4 = b5lindex_atomic(4,k)       ! l4
        lint12 = b5lindex_atomic(5,k)     ! L12  (= L34) 
        do n = 1, madd5_atomic(k)
            m1 = b5lm_atomic(k,n,1)    ! m1
            m2 = b5lm_atomic(k,n,2)    ! m2
            m3 = b5lm_atomic(k,n,3)    ! m3
            m4 = b5lm_atomic(k,n,4)    ! m4
        
            M12 = m1 + m2
                
            if (( abs(m1+m2) > lint12).or.( abs(m3+m4) > lint12) ) then
                print *,'something wrong'
                stop
            endif
            if ( m1 + m2 + m3 + m4 /= 0 ) then
                print *,'something more wrong'
                stop
            endif
            
            cg04_atomic(k,m1,m2,m3,m4) = cg04_atomic(k,m1,m2,m3,m4) +  &
                    &    cgcoeff(l1,l2,lint12,m1,m2,M12)*cgcoeff(l3,l4,lint12,m3,m4,-M12)*cgcoeff(lint12,lint12,0,M12,-M12,0)
        
        enddo
    enddo
        
    do k = 1, nvar5l_magnetic
        l0_p = b5lindex_magnetic(1,k)        ! l0'
        l1_p = b5lindex_magnetic(2,k)        ! l1'
        l2_p = b5lindex_magnetic(3,k)        ! l2'
        l3_p = b5lindex_magnetic(4,k)        ! l3'
        l4_p = b5lindex_magnetic(5,k)        ! l4'  (= L0123')
        lint01_p = b5lindex_magnetic(6,k)    ! L01'
        lint23_p = b5lindex_magnetic(7,k)    ! L23'
        do n = 1, madd5_magnetic(k)
            m0_p = b5lm_magnetic(k,n,1)  ! m0'
            m1_p = b5lm_magnetic(k,n,2)  ! m1'
            m2_p = b5lm_magnetic(k,n,3)  ! m2'
            m3_p = b5lm_magnetic(k,n,4)  ! m3'
            m4_p = b5lm_magnetic(k,n,5)  ! m4'
            
            M01_p = m0_p + m1_p
            M23_p = m2_p + m3_p  
                
            if (( abs(m0_p+m1_p) > lint01_p).or.( abs(m2_p+m3_p) > lint23_p) ) then
                print *,'something wrong'
                stop
            endif
            if ( m0_p + m1_p + m2_p + m3_p + m4_p /= 0 ) then
                print *,'something more wrong'
                stop
            endif
            
            cg04_magnetic(k,m0_p,m1_p,m2_p,m3_p,m4_p) = cg04_magnetic(k,m0_p,m1_p,m2_p,m3_p,m4_p) +                    &
            &    cgcoeff(l0_p,l1_p,lint01_p,m0_p,m1_p,M01_p)*cgcoeff(l2_p,l3_p,lint23_p,m2_p,m3_p,M23_p)*              &    
            &    cgcoeff(lint01_p,lint23_p,l4_p,M01_p,M23_p,-m4_p)*cgcoeff(l4_p,l4_p,0,-m4_p,m4_p,0)
            
        enddo
    enddo

    ! when l6max_atomic is equal to 2 and l6max_atomic to zero segmentation fault
    allocate( cg05_atomic(1:nvar6l_atomic,-l6max_atomic:l6max_atomic,-l6max_atomic:l6max_atomic, &
    &  -l6max_atomic:l6max_atomic,-l6max_atomic:l6max_atomic,-l6max_atomic:l6max_atomic) )

    ! when l6max_atomic is equal to 4 and l6max_magnetic to zero segmentation fault 
    allocate( cg05_magnetic(1:nvar6l_magnetic,-l6max_magnetic:l6max_magnetic,-l6max_magnetic:    &
    &  l6max_magnetic,-l6max_magnetic:l6max_magnetic,-l6max_magnetic:l6max_magnetic,            &
    &  -l6max_magnetic:l6max_magnetic,-l6max_magnetic:l6max_magnetic) )
    
    cg05_atomic = 0.d0
    cg05_magnetic = 0.d0
    !redefinition of the last coupling with the Clebsch-Gordan coupled to zero
    ! initialize the required couplings
    ! generalized coupling for atomic part
    do k = 1, nvar6l_atomic
        l1 = b6lindex_atomic(1,k)       ! l1
        l2 = b6lindex_atomic(2,k)       ! l2
        l3 = b6lindex_atomic(3,k)       ! l3
        l4 = b6lindex_atomic(4,k)       ! l4
        l5 = b6lindex_atomic(5,k)       ! l5
        lint12 = b6lindex_atomic(6,k)    ! L12   
        lint34 = b6lindex_atomic(7,k)    ! L34
        do n = 1, madd6_atomic(k)
            m1 = b6lm_atomic(k,n,1)    ! m1
            m2 = b6lm_atomic(k,n,2)    ! m2
            m3 = b6lm_atomic(k,n,3)    ! m3
            m4 = b6lm_atomic(k,n,4)    ! m4
            m5 = b6lm_atomic(k,n,5)    ! m5
            
            M12 = m1 + m2
            M34 = m3 + m4
                !print*, m1, m2, M12, m3, m4, M34, m5,         
            if (( abs(m1+m2) > lint12).or.( abs(m3+m4) > lint34)) then
                print *,'something wrong'
                stop
            endif
            if ( m1 + m2 + m3 + m4 + m5 /= 0 ) then
                print *,'something more wrong'
                stop
            endif
            
            cg05_atomic(k,m1,m2,m3,m4,m5) = cg05_atomic(k,m1,m2,m3,m4,m5) +  &
            &    cgcoeff(l1,l2,lint12,m1,m2,M12)*cgcoeff(l3,l4,lint34,m3,m4,M34)*cgcoeff(lint12,lint34,l5,M12,M34,-m5)*   &
            &    cgcoeff(l5,l5,0,-m5,m5,0)

        enddo
    enddo
    
    do k = 1, nvar6l_magnetic
        l0_p = b6lindex_magnetic(1,k)       ! l0'
        l1_p = b6lindex_magnetic(2,k)       ! l1'
        l2_p = b6lindex_magnetic(3,k)       ! l2'
        l3_p = b6lindex_magnetic(4,k)       ! l3'
        l4_p = b6lindex_magnetic(5,k)      ! l4'
        l5_p = b6lindex_magnetic(6,k)      ! l5' 
        lint01_p = b6lindex_magnetic(7,k)   ! L01'
        lint23_p = b6lindex_magnetic(8,k)   ! L23'
        lint45_p = b6lindex_magnetic(9,k)   ! L45'
        do n = 1, madd6_magnetic(k)
            m0_p = b6lm_magnetic(k,n,1)  ! m0'
            m1_p = b6lm_magnetic(k,n,2)  ! m1'
            m2_p = b6lm_magnetic(k,n,3)  ! m2'
            m3_p = b6lm_magnetic(k,n,4)  ! m3'
            m4_p = b6lm_magnetic(k,n,5)  ! m4'
            m5_p = b6lm_magnetic(k,n,6)  ! m5'

            M01_p = m0_p + m1_p
            M23_p = m2_p + m3_p 
            M45_p = m4_p + m5_p
                
            if (( abs(m0_p+m1_p) > lint01_p).or.( abs(m2_p+m3_p) > lint23_p).or.  &
                        &   ( abs(m4_p+m5_p) > lint45_p) ) then
                print *,'something wrong'
                stop
            endif
            if ( m0_p + m1_p + m2_p + m3_p + m4_p + m5_p /= 0 ) then
                print *,'something more wrong'
                stop
            endif
            
            cg05_magnetic(k,m0_p,m1_p,m2_p,m3_p,m4_p,m5_p) = cg05_magnetic(k,m0_p,m1_p,m2_p,m3_p,m4_p,m5_p) +          &
            &    cgcoeff(l0_p,l1_p,lint01_p,m0_p,m1_p,M01_p)*cgcoeff(l2_p,l3_p,lint23_p,m2_p,m3_p,M23_p)*              &    
            &    cgcoeff(l4_p,l5_p,lint45_p,m4_p,m5_p,M45_p)*cgcoeff(lint01_p,lint23_p,lint45_p,M01_p,M23_p,-M45_p)*   &    
            &    cgcoeff(lint45_p,lint45_p,0,-M45_p,M45_p,0)
            
        enddo
    enddo
      
end subroutine setuphighercouplings

!----------------------------------------------------------------------------------------------------------

subroutine readacebase(readsuccess)

    use modneigh
    use functionparameters
    use tables
    implicit none
    logical readsuccess

    logical exist
    integer nv, nread1, nread2, nread3, nread4

    inquire(FILE="bace.in", EXIST=exist)

    if (.not. exist) then
        readsuccess = .false.
    else
        open(UNIT=201,FILE='bace.in',STATUS='OLD')

        readsuccess = .true.
        read(201,*)  nread1
        if (nread1 /= nelements) then
            readsuccess = .false.
        endif
        read(201,*)  nread1, nread2, nread3, nread4
        if ((nread1 /= nradbase_atomic).or.(nread2 /= nradbase_magnetic).or.  &
                    &   (nread3 /= nradial2_magnetic).or.(nread4 /= l2max_magnetic)) then
            readsuccess = .false.
        endif
        read(201,*)  nread1, nread2, nread3, nread4
        if ((nread1 /= nradial3_atomic).or.(nread2 /= nradial3_magnetic).or. &
                    &   (nread3 /= l3max_atomic).or.(nread4 /= l3max_magnetic)) then
            readsuccess = .false.
        endif
        read(201,*)  nread1, nread2, nread3, nread4
        if ((nread1 /= nradial4_atomic).or.(nread2 /= nradial4_magnetic).or. &
                    &   (nread3 /= l4max_atomic).or.(nread4 /= l4max_magnetic)) then
            readsuccess = .false.
        endif
        read(201,*)  nread1, nread2, nread3, nread4
        if ((nread1 /= nradial5_atomic).or.(nread2 /= nradial5_magnetic).or. &
                    &   (nread3 /= l5max_atomic).or.(nread4 /= l5max_magnetic)) then
            readsuccess = .false.
        endif
        read(201,*)  nread1, nread2, nread3, nread4
        if ((nread1 /= nradial6_atomic).or.(nread2 /= nradial6_magnetic).or. &
                    &   (nread3 /= l6max_atomic).or.(nread4 /= l6max_magnetic)) then
            readsuccess = .false.
        endif
        
        if (readsuccess) then
            read(201,*)  nread1, nvar1
            read(201,*)  nread1, nvar2
            read(201,*)  nread1, nvar3
            read(201,*)  nread1, nvar4
            read(201,*)  nread1, nvar5
            read(201,*)  nread1, nvar6
        !         read(201,*)  nread1, nvar7
            
            allocate( b1occ(1:nvar1) )
            allocate( b1index(1:nvar1) )
            
            allocate( b2occ(1:nvar2) )
            allocate( b2index(1:4,nvar2) )
            
            allocate( b3occ(1:2,1:nvar3) )
            allocate( b3index(1:9,1:nvar3) )
            
            allocate( b4nint(1:nvar4) )
            allocate( b4occ(1:3,1:nvar4) )
            allocate( b4index(1:15,1:nvar4) )
            
            allocate(b5nint(1:nvar5))
            allocate(b5occ(1:4,1:nvar5))
            allocate(b5index(1:21,1:nvar5))
            
            allocate(b6nint(1:nvar6))
            allocate(b6occ(1:5,1:nvar6))
            allocate(b6index(1:27,1:nvar6))
        
        
            do nv = 1, nvar1
                read(201,*) b1occ(nv)
                read(201,*) b1index(nv)
            enddo
            do nv = 1, nvar2
                read(201,*) b2occ(nv)
                read(201,*) b2index(1:4,nv)
            enddo
            do nv = 1, nvar3
                read(201,*) b3occ(1:2,nv)
                read(201,*) b3index(1:9,nv)
            enddo
            do nv = 1, nvar4
                read(201,*) b4nint(nv)
                read(201,*) b4occ(1:3,nv)
                read(201,*) b4index(1:15,nv)
            enddo
            do nv = 1, nvar5
                read(201,*) b5nint(nv)
                read(201,*) b5occ(1:4,nv)
                read(201,*) b5index(1:21,nv)
            enddo
            do nv = 1, nvar6
                read(201,*) b6nint(nv)
                read(201,*) b6occ(1:5,nv)
                read(201,*) b6index(1:27,nv)
            enddo
        endif
        
        close(201)
        if (readsuccess) then
            print *,'Read basis from file bace.in'
        else
            print *,'Basis functions in file bace.in not compatible. Re-generating bace.in'
        endif
    endif

end subroutine readacebase

!----------------------------------------------------------------------------------------------------------

subroutine bcreate(dof,nrank,nmax,vecn,vecnr,vecl,numlint,arrlint,isordered,isbasis)

    use functionparameters
    use tables
    use modneigh
    use modrandAB
    implicit none
    integer nrank
    integer nmax
    integer numlint

    logical isordered, isbasis, index, reducebasis

    integer vecn(1:nrank)
    integer vecnr(1:nrank)
    integer vecl(1:nrank)
    integer arrlint(1:nmax,1:nrank)

    integer veck(1:nrank)

    integer nrmvec(1:nmax)

    double precision, allocatable :: amatref(:,:)
    double precision, allocatable :: ct4(:,:,:,:)
    double precision, allocatable :: ct5(:,:,:,:,:)
    double precision, allocatable :: ct6(:,:,:,:,:,:)

    double complex xiy


    integer nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8
    integer nr1, nr2, nr3, nr4, nr5, nr6, nr7, nr8
    integer nele1, nele2, nele3, nele4, nele5, nele6, nele7, nele8
    integer k, k1,k2,k3,k4,k5,k6,k7
    integer l, l1,l2,l3,l4,l5,l6,l7
    integer m, m1,m2,m3,m4,m5,m6,m7
    integer lint, lint1, lint2, lint3, lint4, lint12, lint34, lint56, lint1234
    integer mint, mint1, mint2, mint3, mint4, mint12, mint34, mint56, mint1234
    integer lintmax
    integer n, nlint, num

    integer nsample, ntest, ndotest

    character(8) dof
     
    ntest = 1000
    nsample = 20

    call randAB(dof,ntest,nsample)

    do n = 1, nrank
        if ((vecn(n) < 1) .or. (vecn(n) > nelements)) then
            print *,'Element vector outside allowed values.'
            stop
        endif
        if ((vecnr(n) < 1) .or. (vecnr(n) > nradial)) then
            print *,'Radial vector outside allowed values.'
            stop
        endif
        if ((vecl(n) < 0) .or. (vecl(n) > lmax)) then
            print *,'Angular momentum vector outside allowed values.'
            stop
        endif
        !     veck(n) =  ( (nelements + 1 - vecn(n) )*(nradial+1) + vecnr(n) )*(lmax+1) + vecl(n)
        ! new order: element first, then angualar momenta, then radial functions
        veck(n) =  ( (nelements + 1 - vecn(n) )*(lmax+1) + vecl(n) )*(nradial+1) + vecnr(n)
    enddo


    select case (nrank)
        
    case(2) ! rank 2, 3-body
        
        numlint = 1
        arrlint = 0
        
        l1 = vecl(1)
        l2 = vecl(2)
        nc1 = veck(1)
        nc2 = veck(2)
        
        isordered = .false.
        isbasis = .false.
        
        if ( l1 >= l2 ) then
            isbasis = .true.
        endif
        
        if ( nc1 >= nc2 ) then
            isordered = .true.
        endif

     
    case(3) ! rank 3, 4-body
     
        nlint = 0
        numlint = 1

        l1 = vecl(1)
        l2 = vecl(2)
        l3 = vecl(3)
        nc1 = veck(1)
        nc2 = veck(2)
        nc3 = veck(3)

        isbasis = .false.
        isordered = .false.

        if ( ( l1 + l2 >= l3) .and. ( abs(l1-l2) <= l3 ) .and. (mod(l1+l2+l3,2) == 0 )) then
            isbasis = .true.
        endif

        if ( (nc1 >= nc2) .and.( nc2 >= nc3 ) ) then
            isordered = .true.
        endif


    case(4)     ! rank 4, 5-body
     
     !     coupling ( (n1 n2)lint (n3 n4)lint )0
     
        numlint = 0
        arrlint = 0
        
        l1 = vecl(1)
        l2 = vecl(2)
        l3 = vecl(3)
        l4 = vecl(4)
        nr1 = vecnr(1)
        nr2 = vecnr(2)
        nr3 = vecnr(3)
        nr4 = vecnr(4)
        nele1 = vecn(1)
        nele2 = vecn(2)
        nele3 = vecn(3)
        nele4 = vecn(4)
        nc1 = veck(1)
        nc2 = veck(2)
        nc3 = veck(3)
        nc4 = veck(4)
        
        lintmax = 2*maxval(vecl(1:4))
        
        isbasis = .false.
        
        isordered = .true.
        ! take only ordered sequence
        if (( nc1 < nc2 ).or.( nc2 < nc3 ).or.( nc3 < nc4 )) then
            isordered = .false.
        endif
        
        ! loop over possible couplings
        do lint = 0, lintmax
            index = .true.
            ! make sure that coupling is possible at all
            if ((lint > l1+l2) .or. (lint > l3+l4)) then
                index = .false.
            endif
            if ((lint < abs(l1-l2)) .or. (lint < abs(l3-l4))) then
                index = .false.
            endif
            ! remove pairwise couplings that are not possible by inversion
            if ((nc1 == nc2) .and. (mod(lint,2) /= 0)) then
                index = .false.
            endif
            if ((nc3 == nc4) .and. (mod(lint,2) /= 0)) then
                index = .false.
            endif
            if (index) then
                numlint = numlint + 1
                if (numlint > nmax) then
                    print *,'Increase nmax for rank = 4 /nbody = 5 and run again.'
                    stop
                endif
                arrlint(numlint,1) = lint
            endif
        enddo
        if (numlint > 0) then
            isbasis = .true.
        endif
        
        ! remove couplings that are not possible by inversion
        if ( mod(l1+l2+l3+l4,2) /= 0) then
        isbasis = .false.
        endif

        reducebasis = .false.
        if (numlint > 1) then
            if (isbasis) then
                ! need to reduce basis for two or more identical basis functions
                do k1 = 2, 4
                    do k2 = 1, k1 -1
                        if (veck(k1) == veck(k2)) then
                        reducebasis = .true.
                        endif
                    enddo
                enddo
            endif
        endif
        
        if (reducebasis) then
        
            ! this requires CGs to be in place up to 2*lmax
                
            allocate( ct4(-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax))
            ndotest = 4*numlint
            if (ndotest > ntest) then
                print *,'Increase ntest in bcreate.'
                stop
            endif
            allocate(amatref(1:ndotest,1:nmax))
            ! setup coupling coefficients
            ct4 = 0.d0
            amatref = 0.d0
            
            do num = 1, numlint

                ! build ct4 for num
                L = arrlint(num,1)
                do m1 = -l1, l1
                    do m2 = -l2, l2
                        do m3 = -l3, l3
                            do m4 = -l4, l4
                                if (m1 + m2 + m3 + m4 == 0) then
                                    M = m1 + m2
                                    if ( abs(M) <= L ) then
                                        ct4(m1,m2,m3,m4) =  ct4(m1,m2,m3,m4) + &
                                        & cgcoeff(L,L,0,M,-M,0)*cgcoeff(l1,l2,L,m1,m2,M)*cgcoeff(l3,l4,L,m3,m4,-M)
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                enddo

                ! build testvector for num

                do n = 1, ndotest
                    do m1 = -l1, l1
                        do m2 = -l2, l2
                            do m3 = -l3, l3
                                do m4 = -l4, l4
                                    if (m1 + m2 + m3 + m4 == 0 ) then
                                
                                        if (trim(dof) == 'atomic') then
                                            k1 = shindex_atomic(l1,m1)
                                            k2 = shindex_atomic(l2,m2)
                                            k3 = shindex_atomic(l3,m3)
                                            k4 = shindex_atomic(l4,m4)
                                        elseif (trim(dof) == 'magnetic') then
                                            k1 = shindex_magnetic(l1,m1)
                                            k2 = shindex_magnetic(l2,m2)
                                            k3 = shindex_magnetic(l3,m3)
                                            k4 = shindex_magnetic(l4,m4)
                                        endif                            
                                        
                                            xiy = ct4(m1,m2,m3,m4)*randABvec(n,nele1,nr1,k1)*randABvec(n,nele2,nr2,k2) &
                                                & *randABvec(n,nele3,nr3,k3)*randABvec(n,nele4,nr4,k4) 
                                            
                                            amatref(n,num) = amatref(n,num) + dreal(xiy) 
                                    endif
                                
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        
    !         call redbasis(nrank, ndotest, nmax, numlint, amatref, arrlint)
        
            deallocate(ct4)
            deallocate(amatref)
        endif

    case(5)     ! rank 5, 6-body
     
        !     coupling ( ( (n1 n2)lint12 (n3 n4)lint34 )lint1234 n5) 0 with lint1234 = l5

        numlint = 0
        arrlint = 0

        l1 = vecl(1)
        l2 = vecl(2)
        l3 = vecl(3)
        l4 = vecl(4)
        l5 = vecl(5)
        nr1 = vecnr(1)
        nr2 = vecnr(2)
        nr3 = vecnr(3)
        nr4 = vecnr(4)
        nr5 = vecnr(5)
        nele1 = vecn(1)
        nele2 = vecn(2)
        nele3 = vecn(3)
        nele4 = vecn(4)
        nele5 = vecn(5)
        nc1 = veck(1)
        nc2 = veck(2)
        nc3 = veck(3)
        nc4 = veck(4)
        nc5 = veck(5)
     
        lintmax = 2*maxval(vecl(1:5))
        
        isbasis = .false.
        
        isordered = .true.
        ! take only ordered sequence
        if (( nc1 < nc2 ).or.( nc2 < nc3 ).or.( nc3 < nc4 ).or.( nc4 < nc5 )) then
            isordered = .false.
        endif
     
        ! loop over possible couplings
        do lint12 = 0, lintmax
            do lint34 = 0, lintmax
                index = .true.
                ! make sure that coupling is possible at all
                if ((lint12 > l1+l2) .or. (lint34 > l3+l4)) then
                    index = .false.
                endif
                if ((lint12 < abs(l1-l2)) .or. (lint34 < abs(l3-l4))) then
                    index = .false.
                endif
                if ( (l5 > lint12+lint34) .or. (l5 < abs(lint12-lint34)) ) then
                    index = .false.
                endif
                ! remove pairwise couplings that are not possible by inversion
                if ((nc1 == nc2) .and. (mod(lint12,2) /= 0)) then
                    index = .false.
                endif
                if ((nc3 == nc4) .and. (mod(lint34,2) /= 0)) then
                    index = .false.
                endif
                if ((nc1 == nc3) .and.(nc2 == nc4) .and. (lint12 == lint34) .and. &
                        & (mod(l5,2) /= 0)) then
                    index = .false.
                endif
                if (index) then
                    numlint = numlint + 1
                    if (numlint > nmax) then
                        print *,'Increase nmax for rank = 5 /nbody = 6 and run again.'
                        stop
                    endif
                    arrlint(numlint,1) = lint12
                    arrlint(numlint,2) = lint34
                endif
            enddo
        enddo
        if (numlint > 0) then
            isbasis = .true.
        endif

        ! remove couplings that are not possible by inversion
        if ( mod(l1+l2+l3+l4+l5,2) /= 0) then
            isbasis = .false.
        endif

        reducebasis = .false.
        if (numlint > 1) then
            if (isbasis) then
            ! need to reduce basis for two or more identical basis functions
                do k1 = 2, 5
                    do k2 = 1, k1 - 1
                        if (veck(k1) == veck(k2)) then
                            reducebasis = .true.
                        endif
                    enddo
                enddo
            endif
        endif
        
        if (reducebasis) then
            
            ! this requires CGs to be in place up to 2*lmax
            
            allocate( ct5(-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax))
            ndotest = 4*numlint
            if (ndotest > ntest) then
                print *,'Increase ntest in bcreate.'
                stop
            endif
            allocate(amatref(1:ndotest,1:nmax))
            ! setup coupling coefficients
            ct5 = 0.d0
            amatref = 0.d0
            
            do num = 1, numlint
            
                ! build ct5 for num
                lint12 = arrlint(num,1)
                lint34 = arrlint(num,2)
                lint1234 = l5
                do m1 = -l1, l1
                    do m2 = -l2, l2
                        do m3 = -l3, l3
                            do m4 = -l4, l4
                                do m5 = -l5, l5
                                    if (m1 + m2 + m3 + m4 + m5 == 0) then
                                        mint12 = m1 + m2
                                        mint34 = m3 + m4
                                        mint1234 = -m5
                                        if ((abs(mint12) <= lint12).and.(abs(mint34) <= lint34)) then
                                            ct5(m1,m2,m3,m4,m5) =  ct5(m1,m2,m3,m4,m5) + &
                                    & cgcoeff(lint1234,lint1234,0,mint1234,-mint1234,0)*cgcoeff(l1,l2,lint12,m1,m2,mint12)* &
                                    & cgcoeff(l3,l4,lint34,m3,m4,mint34)*cgcoeff(lint12,lint34,lint1234,mint12,mint34,mint1234)
                                        endif
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo

                ! build testvector for num

                do n = 1, ndotest
                    do m1 = -l1, l1
                        do m2 = -l2, l2
                            do m3 = -l3, l3
                                do m4 = -l4, l4
                                    do m5 = -l5, l5
                                        if (m1 + m2 + m3 + m4 + m5 == 0 ) then
                                                                
                                            if (trim(dof) == 'atomic') then
                                                k1 = shindex_atomic(l1,m1)
                                                k2 = shindex_atomic(l2,m2)
                                                k3 = shindex_atomic(l3,m3)
                                                k4 = shindex_atomic(l4,m4)
                                                k5 = shindex_atomic(l5,m5)
                                            elseif (trim(dof) == 'magnetic') then
                                                k1 = shindex_magnetic(l1,m1)
                                                k2 = shindex_magnetic(l2,m2)
                                                k3 = shindex_magnetic(l3,m3)
                                                k4 = shindex_magnetic(l4,m4)
                                                k5 = shindex_magnetic(l5,m5)
                                            endif   
                                                                    
                                            xiy = ct5(m1,m2,m3,m4,m5)*randABvec(n,nele1,nr1,k1)*randABvec(n,nele2,nr2,k2) &
                                                & *randABvec(n,nele3,nr3,k3)*randABvec(n,nele4,nr4,k4)*randABvec(n,nele5,nr5,k5)
                                            
                                            amatref(n,num) = amatref(n,num) + dreal(xiy) 
                                        endif
                                    
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            
    !         call redbasis(nrank, ndotest, nmax, numlint, amatref, arrlint)
            
            deallocate(ct5)
            deallocate(amatref)
        endif

    case(6)     ! rank 6, 7-body
     
        ! coupling scheme ( ( (n1 n2) lint12 (n3 n4) lint34 ) lint1234  (n5 n6) lint56 ) 0
        ! requires lint56 = lint1234
        
        numlint = 0
        arrlint = 0
        
        l1 = vecl(1)
        l2 = vecl(2)
        l3 = vecl(3)
        l4 = vecl(4)
        l5 = vecl(5)
        l6 = vecl(6)
        nr1 = vecnr(1)
        nr2 = vecnr(2)
        nr3 = vecnr(3)
        nr4 = vecnr(4)
        nr5 = vecnr(5)
        nr6 = vecnr(6)
        nele1 = vecn(1)
        nele2 = vecn(2)
        nele3 = vecn(3)
        nele4 = vecn(4)
        nele5 = vecn(5)
        nele6 = vecn(6)
        nc1 = veck(1)
        nc2 = veck(2)
        nc3 = veck(3)
        nc4 = veck(4)
        nc5 = veck(5)
        nc6 = veck(6)
        
        lintmax = 2*maxval(vecl(1:6))
        
        isbasis = .false.
        
        isordered = .true.
        ! take only ordered sequence
        if (( nc1 < nc2 ).or.( nc2 < nc3 ).or.( nc3 < nc4 ).or.( nc4 < nc5 ).or.( nc5 < nc6 )) then
            isordered = .false.
        endif
     
        ! loop over possible couplings
        do lint12 = 0, lintmax
            do lint34 = 0, lintmax
                do lint56 = 0, lintmax
                    index = .true.
                    ! make sure that coupling is possible at all
                    if ((lint12 > l1+l2) .or. (lint34 > l3+l4) .or. (lint56 > l5+l6)) then
                        index = .false.
                    endif
                    if ((lint12 < abs(l1-l2)) .or. (lint34 < abs(l3-l4)) .or. (lint56 < abs(l5-l6))) then
                        index = .false.
                    endif
                    if ( (lint56 > lint12+lint34) .or. (lint56 < abs(lint12-lint34)) ) then
                        index = .false.
                    endif
                    ! remove pairwise couplings that are not possible by inversion
                    if ((nc1 == nc2) .and. (mod(lint12,2) /= 0)) then
                        index = .false.
                    endif
                    if ((nc3 == nc4) .and. (mod(lint34,2) /= 0)) then
                        index = .false.
                    endif
                    if ((nc5 == nc6) .and. (mod(lint56,2) /= 0)) then
                        index = .false.
                    endif
                    if ((nc1 == nc3) .and.(nc2 == nc4) .and. (lint12 == lint34) .and. &
                        & (mod(lint56,2) /= 0)) then
                        index = .false.
                    endif
                    if (index) then
                        numlint = numlint + 1
                        if (numlint > nmax) then
                            print *,'Increase nmax for rank = 6 /nbody = 7 and run again.'
                            stop
                        endif
                        arrlint(numlint,1) = lint12
                        arrlint(numlint,2) = lint34
                        arrlint(numlint,3) = lint56              
                    endif
                enddo
            enddo
        enddo
        if (numlint > 0) then
            isbasis = .true.
        endif
        
        ! remove couplings that are not possible by inversion
        if ( mod(l1+l2+l3+l4+l5+l6,2) /= 0) then
            isbasis = .false.
        endif
     
        reducebasis = .false.
        if (numlint > 1) then
            if (isbasis) then
            ! need to reduce basis for two or more identical basis functions
                do k1 = 2, 6
                    do k2 = 1, k1 - 1
                        if (veck(k1) == veck(k2)) then
                            reducebasis = .true.
                        endif
                    enddo
                enddo
            endif
        endif
     
        if (reducebasis) then
        
            ! this requires CGs to be in place up to 2*lmax
            
            allocate( ct6(-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax,-lmax:lmax))
            ndotest = 4*numlint
            if (ndotest > ntest) then
                print *,'Increase ntest in bcreate.'
                stop
            endif
            allocate(amatref(1:ndotest,1:nmax))
            ! setup coupling coefficients
            ct6 = 0.d0
            amatref = 0.d0
        
            do num = 1, numlint
            
                ! build ct6 for num
                lint12 = arrlint(num,1)
                lint34 = arrlint(num,2)
                lint56 = arrlint(num,3)
                do m1 = -l1, l1
                    do m2 = -l2, l2
                        do m3 = -l3, l3
                            do m4 = -l4, l4
                                do m5 = -l5, l5
                                    do m6 = -l6, l6                         
                                        if (m1 + m2 + m3 + m4 + m5 + m6 == 0) then
                                            mint12 = m1 + m2
                                            mint34 = m3 + m4
                                            mint56 = m5 + m6
                                            if ((abs(mint12) <= lint12).and.(abs(mint34) <= lint34) &
                                                        & .and.(abs(mint56) <= lint56)) then
                                            ct6(m1,m2,m3,m4,m5,m6) =  ct6(m1,m2,m3,m4,m5,m6) + &
                                                    & cgcoeff(lint56,lint56,0,-mint56,mint56,0)* &
                                                    & cgcoeff(l1,l2,lint12,m1,m2,mint12)* &
                                                    & cgcoeff(l3,l4,lint34,m3,m4,mint34)* &
                                                    & cgcoeff(lint12,lint34,lint56,mint12,mint34,-mint56)* &
                                                    & cgcoeff(l5,l6,lint56,m5,m6,mint56)
                                        
                                            endif
                                        endif
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
           
                ! build testvector for num
                
                do n = 1, ndotest
                    do m1 = -l1, l1
                        do m2 = -l2, l2
                            do m3 = -l3, l3
                                do m4 = -l4, l4
                                    do m5 = -l5, l5
                                        do m6 = -l6, l6
                                            if (m1 + m2 + m3 + m4 + m5 + m6 == 0 ) then
                                        
                                                if (trim(dof) == 'atomic') then
                                                    k1 = shindex_atomic(l1,m1)
                                                    k2 = shindex_atomic(l2,m2)
                                                    k3 = shindex_atomic(l3,m3)
                                                    k4 = shindex_atomic(l4,m4)
                                                    k5 = shindex_atomic(l5,m5)
                                                    k6 = shindex_atomic(l6,m6)
                                                elseif (trim(dof) == 'magnetic') then
                                                    k1 = shindex_magnetic(l1,m1)
                                                    k2 = shindex_magnetic(l2,m2)
                                                    k3 = shindex_magnetic(l3,m3)
                                                    k4 = shindex_magnetic(l4,m4)
                                                    k5 = shindex_magnetic(l5,m5)
                                                    k6 = shindex_magnetic(l6,m6)
                                                endif  

                                                xiy = ct6(m1,m2,m3,m4,m5,m6)*randABvec(n,nele1,nr1,k1)* &
                                                        & randABvec(n,nele2,nr2,k2)*randABvec(n,nele3,nr3,k3)* &
                                                        & randABvec(n,nele4,nr4,k4)*randABvec(n,nele5,nr5,k5)* &
                                                        & randABvec(n,nele6,nr6,k6)
                                                
                                                amatref(n,num) = amatref(n,num) + dreal(xiy) 
                                            endif
                                        
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
                
!         call redbasis(nrank, ndotest, nmax, numlint, amatref, arrlint)
        
            deallocate(ct6)
            deallocate(amatref)
        endif
     
    case default
        print *,'Should not be here: rank = ',nrank,' not implemented.'
        stop
   
    end select

end subroutine bcreate

!----------------------------------------------------------------------------------------------------------

subroutine randAB(dof,ntest,nsample)

    use functionparameters
    use tables
    use modneigh
    use modrandAB
    implicit none
    integer ntest, nsample
    double precision ebond(1:3)
    double precision r, x, x1, x2, x3
    integer lsq, nr, nele
    integer l, n, k , k1, m, ns, nsr, nrp
    double complex shvec(1:(lmax+1)**2)
    double complex dshvec(1:3,1:(lmax+1)**2)
    double precision cheb1(0:nradial)
    double precision dcheb1(0:nradial)
    character(8) dof
  
    lsq =  (lmax+1)**2
    
    ! loop only if randSHvec has not been set up before
    if ( allocated(randABvec) ) then
        if ( size( randABvec,2) /= nelements ) then
            deallocate(randABvec)
        elseif ( size( randABvec,3) /= nradial ) then
            deallocate(randABvec)
        elseif ( size( randABvec,4) /= lsq ) then
            deallocate(randABvec)
        endif
    endif
  
    if (.not.allocated(randABvec)) then
        allocate(randABvec(1:ntest,1:nelements,1:nradial,1:lsq))
        randABvec = (0.d0,0.d0)

        if (nsample <= 1) then
            print *,'!!!! Need nsample > 1 to include permutation invariance.'
            print *,'!!!! Setting nsample = 20.'
            nsample = 20
        endif
     
        do n = 1, ntest
            do nele = 1, nelements
                call random_number(x1)
                nsr = floor(x1*dble(nsample)) + 1
                do ns = 1, nsr
                
                    call random_number(x1)
                    call random_number(x2)
                    call random_number(x3)
                    x1 = 2.d0*x1 - 1.d0
                    x2 = 2.d0*x2 - 1.d0
                    x3 = 2.d0*x3 - 1.d0
                    
                    r = sqrt(x1**2 + x2**2 + x3**2)
                    if (r < 1.d-9) then
                        x1 = 1.d0
                        r = sqrt(x1**2 + x2**2 + x3**2)
                    endif
                
                    ebond(1) = x1/r
                    ebond(2) = x2/r
                    ebond(3) = x3/r
                    
                    call random_number(x1)
                    r = 2.d0*x1 - 1.d0
                    if ( dabs(r) > 1.d0 ) then
                        print *,'Something wrong in randAB.'
                        stop
                    endif
                            
                    ! set up Chebyshev polynomials
                    !               nrp = nradial + 1
                    nrp = nradial
                    call calcCheb(nrp,r,cheb1,dcheb1)

                    ! set up spherical harmonics shvec(1:lsq)
                    !call getshhalf(ebond,lmax,lsq,shvec,dshvec)
                    call compute_ylm_cart(dof,ebond(1),ebond(2),ebond(3),lmax,lsq,shvec,dshvec)
                    do l = 1, lmax
                        do m = 1, l
                            x = dble((-1)**(m))
                            
                            if (trim(dof) == 'atomic') then
                                k = shindex_atomic(l,-m)
                                k1 = shindex_atomic(l,m)
                            elseif (trim(dof) == 'magnetic') then
                                k = shindex_magnetic(l,-m)
                                k1 = shindex_magnetic(l,m)
                            endif 
                    
                            shvec(k) = x*conjg( shvec(k1) )
                        enddo
                    enddo
                
                    ! put sh into vector
                    do nr = 1, nradial
                        do l = 1, lsq
                            randABvec(n,nele,nr,l) = randABvec(n,nele,nr,l) + cheb1(nr)*shvec(l)
                        enddo
                    enddo
                
                enddo
            enddo
        enddo

   endif

end subroutine randAB

!----------------------------------------------------------------------------------------------------------
