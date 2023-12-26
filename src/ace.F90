
! ace interface for libcall 
!
!  2 parts
!  1. acelibini: initializes arrays, can be called as often as one likes, but once is sufficient
!  2. acelib: computes energies and forces

!-----------------------------------------------------------------------------

subroutine acelibini(nele,rcut,rcut_magnetic,smear,xcnumber,elenumbers)
  
    use iso_c_binding
    use global
    use modneigh
    use functionparameters
    implicit none
    integer(kind=c_int) nele  
    real(kind=c_double) rcut, rcut_magnetic
    real(kind=c_double) smear
    integer(kind=c_int) xcnumber
    integer(kind=c_int) elenumbers(nele)

    character(kind=c_char,len=60) jobname

    character(2) char2a, char2b
    character(2) elementlist(1:104)

    integer n, k , l, m, nr

    elementlist = (/'X ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na', &
        & 'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn',&
        & 'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr',&
        & 'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',&
        & 'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
        & 'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir',&
        & 'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',&
        & 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/)

    jobname = 'libcall '

    cutoff = rcut
    cutoff_magnetic = rcut_magnetic

    smearing = smear

    if (nele /= 1 ) then
        print *,'Multiple elements require fixing acelibini call first.'
        stop
    endif

    if (xcnumber == 1 ) then
        xc = "pbe"
    elseif (xcnumber == 2 ) then
        xc = "lda"
    else
        print *,'Unknown XC functional received by acelibini.'
        stop
    endif  

    nelements = nele
    if (.not.allocated(element)) then
        allocate(element(1:nelements))
    elseif (size(element,1) /= nelements) then
        deallocate(element)
        allocate(element(1:nelements))
    endif
    do n = 1, nelements
        element(n) = elementlist(elenumbers(n))
    enddo

    dorun = .true.

    call aceinit ! aceinit can be called multiple times

end subroutine acelibini

!---------------------------------------------------------------------------

subroutine acelib(natom,nall,nei,nlist,occ,occlist,nstart,nstop,eout,dlist,elist,flist,   &
        &    dlist_momi,elist_momi,dlist_mom,elist_mom,tlist,tlist0)
  
    use iso_c_binding
    use global
    use modneigh
    use functionparameters
    use modforces
    implicit none    
    integer(kind=c_int) natom
    integer(kind=c_int) nall
    integer(kind=c_int) nei    
    integer(kind=c_int) nlist(nei)
    integer(kind=c_int) occ(natom)
    integer(kind=c_int) occlist(nei)
    integer(kind=c_int) nstart(natom)
    integer(kind=c_int) nstop(natom)
    real(kind=c_double) eout(natom)
    real(kind=c_double) dlist(nei)
    real(kind=c_double) elist(3*nei)
    real(kind=c_double) flist(3*nei)    

    ! magnetic atom and neighbor-lists, list for torques
    real(kind=c_double) dlist_momi(natom)
    real(kind=c_double) elist_momi(3*natom)
    real(kind=c_double) dlist_mom(nei)
    real(kind=c_double) elist_mom(3*nei)
    real(kind=c_double) tlist(3*nei) 
    real(kind=c_double) tlist0(3*natom)
   
    ! problem with passing too many strings ...    
    character(kind=c_char,len=60) jobname

    integer nei2count
    integer atomk, ij, k, n, m, atomi

    ! copy everything to global variables 

    jobname = 'libcall'

    strucname = jobname

    n2atom = natom
    ! for array allocation
    natomsmax = nall
    !    natomsmax = nmax

    if (.not.allocated(aty)) then
        allocate(aty(1:natom))
    elseif (size(aty,1) /= natom) then
        deallocate(aty)
        allocate(aty(1:natom))
    endif
    aty = occ

    ! magnetic part: atomic list
    if (.not.allocated(d2l_momi)) then
        allocate(d2l_momi(1:natom))
    elseif (size(d2l_momi,1) /= natom) then
        deallocate(d2l_momi)
        allocate(d2l_momi(1:natom))
    endif          
    if (.not.allocated(e2l_momi)) then
        allocate(e2l_momi(1:3,1:natom))
    elseif (size(e2l_momi,2) /= natom) then
        deallocate(e2l_momi)
        allocate(e2l_momi(1:3,1:natom))
    endif
    d2l_momi = dlist_momi
    k = 0
    do n = 1, natom
        do m = 1, 3
            k = k + 1
            e2l_momi(m,n) = elist_momi(k)
        enddo
    enddo 

    nei2count = nei
    if (.not.allocated(n2l)) then
        allocate(n2l(1:nei2count))
    elseif (size(n2l,1) /= nei2count) then
        deallocate(n2l)
        allocate(n2l(1:nei2count))
    endif
    n2l = nlist
    if (.not.allocated(aty2l)) then
        allocate(aty2l(1:nei2count))
    elseif (size(aty2l,1) /= nei2count) then
        deallocate(aty2l)
        allocate(aty2l(1:nei2count))
    endif
    aty2l = occlist

    ! atomic part: need to populate distances and bond directions
    if (.not.allocated(d2l)) then
        allocate(d2l(1:nei2count))
    elseif (size(d2l,1) /= nei2count) then
        deallocate(d2l)
        allocate(d2l(1:nei2count))
    endif          
    if (.not.allocated(e2l)) then
        allocate(e2l(1:3,1:nei2count))
    elseif (size(e2l,2) /= nei2count) then
        deallocate(e2l)
        allocate(e2l(1:3,1:nei2count))
    endif
    d2l = dlist
    k = 0
    do n = 1, nei2count
        do m = 1, 3
            k = k + 1
            e2l(m,n) = elist(k)
        enddo
    enddo

    ! magnetic part: need to populate magnitudes and moments directions
    if (.not.allocated(d2l_mom)) then
        allocate(d2l_mom(1:nei2count))
    elseif (size(d2l_mom,1) /= nei2count) then
        deallocate(d2l_mom)
        allocate(d2l_mom(1:nei2count))
    endif          
    if (.not.allocated(e2l_mom)) then
        allocate(e2l_mom(1:3,1:nei2count))
    elseif (size(e2l_mom,2) /= nei2count) then
        deallocate(e2l_mom)
        allocate(e2l_mom(1:3,1:nei2count))
    endif
    d2l_mom = dlist_mom
    k = 0
    do n = 1, nei2count
        do m = 1, 3
            k = k + 1
            e2l_mom(m,n) = elist_mom(k)
        enddo
    enddo

    if (.not.allocated(n2start)) then
        allocate(n2start(1:n2atom))
    elseif (size(n2start,1) /= n2atom) then
        deallocate(n2start)
        allocate(n2start(1:n2atom))
    endif
        if (.not.allocated(n2stop)) then
        allocate(n2stop(1:n2atom))
    elseif (size(n2stop,1) /= n2atom) then
        deallocate(n2stop)
        allocate(n2stop(1:n2atom))
    endif   
    n2start = nstart
    n2stop = nstop

    if (.not.allocated(fcalcvec)) then
        allocate(fcalcvec(1:3,1:n2atom))
    elseif (size(fcalcvec,2) /= n2atom) then
        deallocate(fcalcvec)
        allocate(fcalcvec(1:3,1:n2atom))
    endif

    if (.not.allocated(tcalcvec)) then
        allocate(tcalcvec(1:3,1:n2atom))
    elseif (size(tcalcvec,2) /= n2atom) then
        deallocate(tcalcvec)
        allocate(tcalcvec(1:3,1:n2atom))
    endif

    if (.not.allocated(eicalc)) then
        allocate(eicalc(n2atom))
    elseif (size(eicalc,1) /= n2atom) then
        deallocate(eicalc)
        allocate(eicalc(1:n2atom))
    endif

    if ((allocated(fijcalc)).and.(nei2count/=size(fijcalc,2))) then
        deallocate(fijcalc)
    endif
    if (.not. allocated(fijcalc)) then
        allocate( fijcalc(1:3,1:nei2count))
    endif

    if ((allocated(tijcalc)).and.(nei2count/=size(tijcalc,2))) then
        deallocate(tijcalc)
    endif
    if (.not. allocated(tijcalc)) then
        allocate( tijcalc(1:3,1:nei2count))
    endif

    if ((allocated(tijcalc0)).and.(n2atom/=size(tijcalc0,2))) then
        deallocate(tijcalc0)
    endif
    if (.not. allocated(tijcalc0)) then
        allocate( tijcalc0(1:3,1:n2atom))
    endif

    ! array boundaries
    nei2countmax = nei2count

    !-------------------------------------------------------------
    ! here everything is in place for energy and force computation
    !-------------------------------------------------------------

    call allocateforprepAndBandforce ! check if arrays need updates

    call prepA         ! prepares atomic base
    call prepB         ! prepares ACE basis functions
    call forceweights  ! calculates forces and energy part I
    call forcecompute  ! calculates forces and energy part II
    call torquecompute ! calculates magnetic gradients part III

    eout(1:n2atom) = eicalc(1:n2atom)

    ! copy everything in 1d array for transfer to C
    k = 0
    do atomk = 1, n2atom
        do ij = n2start(atomk), n2stop(atomk)
            do n = 1, 3
                k = k + 1
                flist(k) = fijcalc(n,ij)
            enddo
        enddo
    enddo

    ! copy everything in 1d array for transfer to C
    k = 0
    do atomk = 1, n2atom
        do ij = n2start(atomk), n2stop(atomk)
            do n = 1, 3
                k = k + 1
                tlist(k) = tijcalc(n,ij)
            enddo
        enddo
    enddo

    k = 0
    do atomk = 1, n2atom
        do n = 1, 3
            k = k + 1
            tlist0(k) = tijcalc0(n,atomk)
        enddo
    enddo

    
end subroutine acelib

!---------------------------------------------------------

