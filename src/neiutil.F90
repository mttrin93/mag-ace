!
! neigborlist utils
!------------------------------

 subroutine packstructures
!
! subroutine builds neighbours and packs all structures and their neighbourlist into larger list
!
    use global
    use neilistglobal        
    implicit none
    integer, allocatable :: nr2(:)
    double precision, allocatable :: dr2(:), er2(:,:)
    double precision, allocatable :: dr2_mom(:), er2_mom(:,:)

    integer i,m, n
    integer natompack, npack, npackdelta, natomboxpack
    integer nei2pack
    integer nei2packalloc
    integer natomboxpackalloc
    integer nei2count
    integer nei2atomallmax

    double precision printcount, dprintcount
    double precision distminmin, distminmax, vol

    character(60) nameminmin, nameminmax

#ifndef MPIPARALLEL
    write(*,'(A)') 'Start packing structures'
#endif

    rcut = cutoff

    do3 = .false.
    do4 = .false.

    if ((do3).or.(do4)) then
        print *,'3-atom and 4-atom neighbor list not activated.'
        stop
    endif

    natomsmaxcheck = 0
    natomboxmax = 0
    nei2countmax = 0
    nei2atomallmax = 0
    nei2count = -12345

    ! deallocte everything if it is already there
    if (allocated(strucaddress)) deallocate(strucaddress)
    if (allocated(dmin)) deallocate(dmin)
    if (allocated(dminmax)) deallocate(dminmax)
    if (allocated(atvol)) deallocate(atvol)
    if (allocated(natomboxpackstart)) deallocate(natomboxpackstart)
    if (allocated(natomboxpackstop)) deallocate(natomboxpackstop)
    if (allocated(originpack)) deallocate(originpack)
    if (allocated(n2startpack)) deallocate(n2startpack)
    if (allocated(n2stoppack)) deallocate(n2stoppack)
    if (allocated(n2listpack)) deallocate(n2listpack)
    if (allocated(d2listpack)) deallocate(d2listpack)
    if (allocated(d2listpack_mom)) deallocate(d2listpack_mom)
    if (allocated(e2listpack)) deallocate(e2listpack)
    if (allocated(e2listpack_mom)) deallocate(e2listpack_mom)

    if (allocated(d2listpack_momi)) deallocate(d2listpack_momi)
    if (allocated(e2listpack_momi)) deallocate(e2listpack_momi)


    ! index for structures
    allocate(strucaddress(1:nstruc,1:4))
    do n = 1, nstruc
        m = natoms(n)
        if (m > natomsmaxcheck) natomsmaxcheck = m
    enddo

    ! total number of atoms in all cells: ntotatoms

    npack = 100000
    npackdelta = 10000

    natomboxpackalloc = npack
    nei2packalloc = npack

    allocate(dmin(1:nstruc))
    allocate(dminmax(1:nstruc))
    allocate(atvol(1:nstruc))
    atvol = 0.d0

    allocate(natomboxpackstart(1:nstruc))
    allocate(natomboxpackstop(1:nstruc))
    allocate(originpack(1:natomboxpackalloc))

    allocate(n2startpack(1:ntotatoms))
    allocate(n2stoppack(1:ntotatoms))

    allocate(n2listpack(nei2packalloc))
    allocate(d2listpack(nei2packalloc))
    allocate(e2listpack(3,nei2packalloc))
    allocate(d2listpack_mom(nei2packalloc))
    allocate(e2listpack_mom(3,nei2packalloc))

    allocate(d2listpack_momi(1:ntotatoms))
    allocate(e2listpack_momi(1:3,1:ntotatoms))


    natompack = 0
    natomboxpack = 0
    nei2pack = 0

    printcount = 0.d0
    dprintcount = 0.1d0

    distminmin = 1.d10
    distminmax = 0.d0

    do n = 1, nstruc
    if (allocated(posin) ) deallocate(posin)
        natom = natoms(n)
        allocate(posin(1:3,1:natom))   
    if (allocated(magin) ) deallocate(magin)
        allocate(magin(1:3,1:natom))
    do i = 1, 3
        posin(i,1:natom) = atompos(i,1:natom,n)
        magin(i,1:natom) = magmom(i,1:natom,n)
    enddo
    periodic = periodiccell(n)
    cell(1:3,1:3) = unitcell(1:3,1:3,n)

    if (periodic) then
        call cellvol(cell,vol)
        atvol(n) = vol/float(natom)
    endif
        
    call neilist(do3,do4)        

    ! pack origin
    if (natomboxpack + natombox > natomboxpackalloc) then
        allocate(nr2(1:natomboxpackalloc))
        nr2(1:natomboxpack) = originpack(1:natomboxpack)
        deallocate(originpack)
        allocate(originpack(1:natomboxpack+natombox+npackdelta))
        originpack(1:natomboxpack) = nr2(1:natomboxpack)
        deallocate(nr2)
    endif
    natomboxpackstart(n) = natomboxpack + 1
    natomboxpackstop(n) = natomboxpack + natombox
    originpack(natomboxpack+1:natomboxpack+natombox) = origin(1:natombox)
    
    
    d2listpack_momi(natompack+1:natompack+natom) = dist_momi(1:natom)
    e2listpack_momi(1:3,natompack+1:natompack+natom) = eist_momi(1:3,1:natom)
    
    
    ! pack atomic pairs
    strucaddress(n,1) = natompack + 1
    n2startpack(natompack+1:natompack+natom) = nrangestart(1:natom)
    n2stoppack(natompack+1:natompack+natom) = nrangestop(1:natom)
    nei2count = nrangestop(natom)
    if (nei2pack + nei2count > nei2packalloc ) then
        allocate(nr2(1:nei2packalloc))
        allocate(dr2(1:nei2packalloc))
        allocate(dr2_mom(1:nei2packalloc))
        allocate(er2(1:3,1:nei2packalloc))
        allocate(er2_mom(1:3,1:nei2packalloc))
        nr2(1:nei2packalloc) = n2listpack(1:nei2packalloc)
        dr2(1:nei2packalloc) = d2listpack(1:nei2packalloc)
        dr2_mom(1:nei2packalloc) = d2listpack_mom(1:nei2packalloc)
        er2(1:3,1:nei2packalloc) = e2listpack(1:3,1:nei2packalloc)
        er2_mom(1:3,1:nei2packalloc) = e2listpack_mom(1:3,1:nei2packalloc)
        deallocate(n2listpack)
        deallocate(d2listpack)
        deallocate(d2listpack_mom)
        deallocate(e2listpack)
        deallocate(e2listpack_mom)
        allocate(n2listpack(1:nei2packalloc+nei2count+npackdelta))
        allocate(d2listpack(1:nei2packalloc+nei2count+npackdelta))
        allocate(d2listpack_mom(1:nei2packalloc+nei2count+npackdelta))
        allocate(e2listpack(1:3,1:nei2packalloc+nei2count+npackdelta))
        allocate(e2listpack_mom(1:3,1:nei2packalloc+nei2count+npackdelta))
        n2listpack(1:nei2packalloc) = nr2(1:nei2packalloc)
        d2listpack(1:nei2packalloc) = dr2(1:nei2packalloc)
        d2listpack_mom(1:nei2packalloc) = dr2_mom(1:nei2packalloc)
        e2listpack(1:3,1:nei2packalloc) = er2(1:3,1:nei2packalloc)
        e2listpack_mom(1:3,1:nei2packalloc) = er2_mom(1:3,1:nei2packalloc)
        deallocate(nr2)
        deallocate(dr2)
        deallocate(dr2_mom)
        deallocate(er2)
        deallocate(er2_mom)
    endif
    strucaddress(n,2) = nei2pack + 1
    if (nei2count > 0 ) then
        n2listpack(nei2pack+1:nei2pack+nei2count) = origin(nlist(1:nei2count))
        d2listpack(nei2pack+1:nei2pack+nei2count) = dist(1:nei2count)
        d2listpack_mom(nei2pack+1:nei2pack+nei2count) = dist_mom(1:nei2count)
        e2listpack(1:3,nei2pack+1:nei2pack+nei2count) = eist(1:3,1:nei2count)
        e2listpack_mom(1:3,nei2pack+1:nei2pack+nei2count) = eist_mom(1:3,1:nei2count)
    endif

    dmin(n) = distmin
    dminmax(n) = distminmax2
    if ( (natom > 1) .or. (periodiccell(n))) then
        if (distmin < distminmin) then
            distminmin = distmin
            nameminmin = name(n)
        endif
        if (distmin > distminmax) then
            distminmax = distmin
            nameminmax = name(n)
        endif
    endif
    
    natompack = natompack + natom
    if (natomboxpack + natombox > natomboxpackalloc) then
        natomboxpackalloc = natomboxpack + natombox + npackdelta
    endif
    natomboxpack = natomboxpack + natombox
    if (nei2atommax > nei2atomallmax) nei2atomallmax = nei2atommax
    if (natombox > natomboxmax) natomboxmax = natombox
    if (nei2count > nei2countmax) nei2countmax = nei2count
    if (nei2pack + nei2count > nei2packalloc ) then
        nei2packalloc = nei2packalloc + nei2count + npackdelta
    endif
    nei2pack = nei2pack + nei2count
    
    call deallocate_neilistglobal

#ifndef MPIPARALLEL   
    if ( printcount <= float(n)/float(nstruc) ) then
        write(*,'(A,f6.1,A,i7,A,i7,A)') 'Packed ', float(n)/float(nstruc)*100.d0,'% (',n,' of ',nstruc,')'
        printcount = printcount + dprintcount
    endif
#endif
    
    enddo

    nei2atommax =  nei2atomallmax

    if (natomsmax < natomsmaxcheck) then
        print *,'Inconsistent count for number of atoms.'
        stop
    endif

#ifndef MPIPARALLEL
    write(*,'(A)') 'Done packing structures'
    write(*,'(A)') ' ---------- Summary from packing structures ----------' 
    write(*,'(A,i7)') 'max number of atoms in any structure                      ', natomsmaxcheck
    write(*,'(A,i7)') 'max number of atoms in any structure + periodic boundary  ', natomboxmax
    write(*,'(A,i7)') 'max number of pairs for any atom                          ', nei2atommax
    write(*,'(A,f6.2)') 'average number of pairs for an atom                       ', dble(nei2pack)/dble(natompack)
    write(*,'(A,f12.5,A,A)') 'shortest distance in any structure                   ', distminmin,'  ',nameminmin
    write(*,'(A,f12.5,A,A)') 'longest shortest distance in any structure           ', distminmax,'  ',nameminmax
    if (abs(distminmax - rcut) < 1.d-9 ) then
        m = 0
        do n = 1, nstruc
            if ( abs(dmin(n) - rcut) < 1.d-9 ) then
                m = m + 1
            endif
        enddo
    write(*,'(A,i5,A,i5,A)') 'WARNING: ', m,' of ', nstruc,' structures with smallest distance outside cutoff'
    endif
#endif

end subroutine packstructures

!------------------------------------------------------------------------------------      
      subroutine unpackstructure(n)
!-----------------------------------------------------------------------------------

! gets n2list from memory
                
    use global       
    use neilistglobal
    use modneigh
    !use mpiglobal
    implicit none
    integer n,k

    integer pack1count, pack2count, boxstart, boxstop, nbox
    integer nei2count

    ! allocate arrays in order that they can hold any of the nstruc structures
    ! add + 1 everywhere as arguments may be zero
    if (.not. allocated(n2start)) allocate(n2start(natomsmax+1))
    if (.not. allocated(n2stop)) allocate(n2stop(natomsmax+1))

    if (.not. allocated(aty)) allocate(aty(natomsmax+1))
    if (.not. allocated(amom)) allocate(amom(3,natomsmax+1))
    if (.not. allocated(ol)) allocate(ol(natomboxmax+1))

    if (.not. allocated(aty2l)) allocate(aty2l(nei2countmax+1))
    if (.not. allocated(n2l)) allocate(n2l(nei2countmax+1))
    if (.not. allocated(d2l)) allocate(d2l(nei2countmax+1))
    if (.not. allocated(d2l_mom)) allocate(d2l_mom(nei2countmax+1))
    if (.not. allocated(e2l)) allocate(e2l(3,nei2countmax+1))
    if (.not. allocated(e2l_mom)) allocate(e2l_mom(3,nei2countmax+1))

    if (.not. allocated(d2l_momi)) allocate(d2l_momi(natomsmax+1))
    if (.not. allocated(e2l_momi)) allocate(e2l_momi(3,natomsmax+1))

    strucname = name(n)
    natom = natoms(n)
    n2atom = natom
    !!! change this to all atoms, including ghost atoms
    aty(1:natom) = occpos(1:natom,n)
    amom(1:3,1:natom) = magmom(1:3,1:natom,n)
    boxstart = natomboxpackstart(n)
    boxstop = natomboxpackstop(n)
    nbox = boxstop - boxstart + 1
    ol(1:nbox) = originpack(boxstart:boxstop)

    pack1count = strucaddress(n,1)

    n2start(1:natom) = n2startpack(pack1count:pack1count+natom-1)
    n2stop(1:natom) = n2stoppack(pack1count:pack1count+natom-1)
    nei2count = n2stop(natom)
    pack2count = strucaddress(n,2)

    d2l_momi(1:natom) = d2listpack_momi(pack1count:pack1count+natom-1)
    e2l_momi(1:3,1:natom) = e2listpack_momi(1:3,pack1count:pack1count+natom-1)

    if (nei2count > 0 ) then
        n2l(1:nei2count) = n2listpack(pack2count:pack2count+nei2count-1)  
        d2l(1:nei2count) = d2listpack(pack2count:pack2count+nei2count-1)
        d2l_mom(1:nei2count) = d2listpack_mom(pack2count:pack2count+nei2count-1)
        e2l(1:3,1:nei2count) = e2listpack(1:3,pack2count:pack2count+nei2count-1)
        e2l_mom(1:3,1:nei2count) = e2listpack_mom(1:3,pack2count:pack2count+nei2count-1)
        aty2l(1:nei2count) = aty(ol(n2l(1:nei2count)))
    endif

end subroutine unpackstructure
!-------------------------------------------------------------------------------
