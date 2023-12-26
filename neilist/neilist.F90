
!--------------------------------------------------------
module neilistglobal
    implicit none

    logical periodic
    integer natom, nalloc
    double precision rcut, r3cut, r4cut
    double precision, parameter :: eps = 1.d-9
    double precision, allocatable :: pos(:,:), posin(:,:)
    double precision, allocatable :: magin(:,:), mag(:,:)
    character(2), allocatable :: species(:)
    ! double precision a1(3), a2(3), a3(3)
    double precision cell(3,3)
    integer, allocatable :: nlist(:)
    double precision, allocatable :: dist(:)
    double precision, allocatable :: eist(:,:)
    integer, allocatable :: nrangestart(:), nrangestop(:)
    integer, allocatable :: n3list(:,:)
    double precision, allocatable :: d3ist(:,:)
    double precision, allocatable :: e3ist(:,:,:)
    double precision, allocatable :: dist_mom(:)
    double precision, allocatable :: eist_mom(:,:)

    double precision, allocatable :: dist_momi(:)
    double precision, allocatable :: eist_momi(:,:)

    integer, allocatable :: n3rangestart(:), n3rangestop(:)
    integer, allocatable :: n4list(:,:)
    double precision, allocatable :: d4ist(:,:)
    double precision, allocatable :: e4ist(:,:,:)
    integer, allocatable :: n4rangestart(:), n4rangestop(:)
    integer nei2atommax, nei3atommax, nei4atommax
    double precision distmin, distminmax2


    integer, allocatable :: nshift(:,:)
    integer, allocatable :: origin(:)

    integer hc1, hc2, hc3, lc1, lc2, lc3
    integer natombox
    integer nresize

    contains
    !--------------------------------------------------------------
        subroutine deallocate_neilistglobal
    !--------------------------------------------------------------
            if (allocated(pos)) deallocate(pos)
            if (allocated(posin)) deallocate(posin)
            if (allocated(mag)) deallocate(mag)
            if (allocated(magin)) deallocate(magin)
            if (allocated(species)) deallocate(species)
            if (allocated(nshift)) deallocate(nshift)          
            if (allocated(origin)) deallocate(origin)
            ! two body list          
            if (allocated(nlist)) deallocate(nlist)
            if (allocated(dist)) deallocate(dist)
            if (allocated(eist)) deallocate(eist)
            if (allocated(dist_mom)) deallocate(dist_mom)
            if (allocated(eist_mom)) deallocate(eist_mom)
            
            if (allocated(dist_momi)) deallocate(dist_momi)
            if (allocated(eist_momi)) deallocate(eist_momi)          
            
            if (allocated(nrangestart)) deallocate(nrangestart)          
            if (allocated(nrangestop)) deallocate(nrangestop)          
            ! three body list          
            if (allocated(n3list)) deallocate(n3list)          
            if (allocated(d3ist)) deallocate(d3ist)
            if (allocated(e3ist)) deallocate(e3ist)          
            if (allocated(n3rangestart)) deallocate(n3rangestart)          
            if (allocated(n3rangestop)) deallocate(n3rangestop)
            ! four body list          
            if (allocated(n4list)) deallocate(n4list)
            if (allocated(d4ist)) deallocate(d4ist)
            if (allocated(e4ist)) deallocate(e4ist)          
            if (allocated(n4rangestart)) deallocate(n4rangestart)          
            if (allocated(n4rangestop)) deallocate(n4rangestop) 
            
        end subroutine deallocate_neilistglobal
            
end module neilistglobal

!----------------------------------------------------------------------------
module grid
!---------------------------------------------------------------------------
        
    implicit none
    !     number of gridboxes
    integer ngridbox
    integer nglow1, nglow2, nglow3, ngup1, ngup2, ngup3
    !     range of atoms in gridboxes
    integer, allocatable :: natgbstart(:), natgbstop(:)
    !     natgb: i = natgb(m): atom m is in gridbox number i
    integer, allocatable :: natgb(:)
    !     gbnat: m = gbnat(i): position i is atom m
    integer, allocatable :: gbnat(:)
    !     short-range neighbour-list
    integer ngbnei
    integer, allocatable ::  ngbrangestart(:), ngbrangestop(:)
    integer, allocatable :: ngblist(:)
    !     long-range neighbour-list
    integer ngb2nei
    integer, allocatable ::  ngb2rangestart(:), ngb2rangestop(:)
    integer, allocatable :: ngb2list(:)
    !
    contains
    
        !--------------------------------------------------------------
        subroutine allocate_grid(natombox,ngridbox,nmax)
        !--------------------------------------------------------------
            implicit none
            integer, intent(in):: natombox,ngridbox,nmax
            !       atom <-> gridbox lists     
            allocate(natgbstart(1:ngridbox))
            allocate(natgbstop(1:ngridbox))
            allocate(natgb(1:natombox))
            allocate(gbnat(1:natombox))
            !       short-range 
            allocate(ngblist(1:nmax))
            allocate(ngbrangestart(1:ngridbox))
            allocate(ngbrangestop(1:ngridbox))
            ngblist       = 0
            ngbrangestart = 0
            ngbrangestop  = 0

        end subroutine allocate_grid
        
        !--------------------------------------------------------------
        subroutine deallocate_grid
        !--------------------------------------------------------------
            implicit none
            if (allocated(natgbstart))     deallocate(natgbstart)
            if (allocated(natgbstop))      deallocate(natgbstop)
            if (allocated(natgb))          deallocate(natgb)
            if (allocated(gbnat))          deallocate(gbnat)
            if (allocated(ngblist))        deallocate(ngblist)
            if (allocated(ngbrangestart))  deallocate(ngbrangestart)
            if (allocated(ngbrangestop))   deallocate(ngbrangestop)

        end subroutine deallocate_grid
        
!--------------------------------------------------------------
end module grid
!---------------------------------------------------------------

!--------------------------------------------------------
subroutine neilist(do3,do4)
!--------------------------------------------------------
! global variables in neilistglobal required as input
! in: natom: number of atoms in unit cell or cluster
! in: periodic: true for periodic cells, false for cluster  
! in: posin(1:natom,1:3) : atomic positions
! in: species(1:natom)
! in: cell(1:3,1:3) : unit cell   cell(1:3,1) contains the first unit cell vector, etc.
! in: rcut : cut-off for neighbour list
! in: magin(1:3,1:natom) : atomic magnetic moments

! out: natombox: number of atoms in cell plus (if required) periodic boundary
! out: pos(1:natombox,1:3) : positions of all atoms inclusive periodic boundary
! out: mag(1:natombox,1:3) : magnetic moments of all atoms inclusive periodic boundary
! out: nrangestart(1:natombox), nrangestop(1:natombox):
!    neighbour index for atom n ranges from i = nrangestart(n) to nrangestop(n)
! !!! for atoms natom+1 to natombox only the neighbours in the real space cluster are evaluated, i.e. not all neighbours in the periodic crystal !! 
! out: nlist(i): neighbours of atom n
! out: nshift(1:natombox,1:3) unit cell index for each atom
!  integer shifts nhisft(n,1:3) in terms of the unit cell vectors required
!  to reach from the unit cell at the origin to the unit cell of atom n
! out: origin(1:natombox)
!   k = origin(n) gives the number of atom k that is obtained when atom n is shifted
!   in the original unit cell
! out: nresize: number of array resize operations  

    implicit none
    logical do3, do4
    
    call buildcluster
    call buildgrid
    call refreshgrid
    call refreshneighbours
    if (do3) then
    call threepoint
    if (do4) then
        call fourpoint
    endif
    endif

    return
end subroutine neilist

!-----------------------------------------------------------------------------
subroutine buildcluster
!-----------------------------------------------------------------------------  

    use neilistglobal
    use grid  
    implicit none  
    integer n,k,k1,k2,k3,n1,n2,n3
    double precision atompos(1:3,1:natom),cellinv(1:3,1:3),magmom(1:3,1:natom)
    double precision v(3), v1(3), v2(3), v3(3),x, r1(3), r2(3), rout, dxhc1, dxlc1
    double precision dxhc2, dxlc2, dxhc3, dxlc3, dxhcc1, dxlcc1, dxhcc2, dxlcc2, dxhcc3, dxlcc3
    integer, allocatable :: origintmp(:), nshifttmp(:,:)
    double precision, allocatable :: postmp(:,:),magtmp(:,:)

    if (rcut < 1.d-6) then
    print *,'value for rcut to small.'
    stop
    endif
    if (rcut > 100.d0) then
    print *,'value for rcut to large.'
    stop
    endif


    if (periodic) then
    !     1. calculate size of cluster by adding a layer of thickness
    !        rcut parallel to the faces of the unit cell
    
        v1(1:3) = cell(1:3,1)
        v2(1:3) = cell(1:3,2)
        v3(1:3) = cell(1:3,3)
        
        !     x is the thickness of the cell orthogonal to the v1-v2 plane
        call north(v1,v2,v)
        x = sum(v3(1:3)*v(1:3))
        if (x < 1.d-3) then
            print *,'Extreme shape of unit cell. Stopping.'
            stop
        endif
        dxhc3 =  1.0d0 + rcut/x
        dxlc3 = -dxhc3 + 1.0d0    
        hc3 = ceiling(dxhc3)
        lc3 = floor(dxlc3)
    
        !     x is the thickness of the cell orthogonal to the v3-v1 plane
        call north(v3,v1,v)
        x = sum(v2(1:3)*v(1:3))
        if (x < 1.d-3) then
            print *,'Extreme shape of unit cell. Stopping.'
            stop
        endif
        dxhc2 =  1.0d0 + rcut/x
        dxlc2 = -dxhc2 + 1.0d0    
        hc2 = ceiling(dxhc2)
        lc2 = floor(dxlc2)
        
        !     x is the thickness of the cell orthogonal to the v2-v3 plane      
        call north(v2,v3,v)
        x = sum(v1(1:3)*v(1:3))
        if (x < 1.d-3) then
            print *,'Extreme shape of unit cell. Stopping.'
            stop
        endif
        dxhc1 =  1.0d0 + rcut/x
        dxlc1 = -dxhc1 + 1.0d0  
        hc1 = ceiling(dxhc1)
        lc1 = floor(dxlc1)
    
        !     2. evaluate size of orthorombic box that contains unit cell
        !        plus distance rcutmax around unit cell 
        
        do n = 1, 3
            rout = 0.d0
            do n1 = 0, 1
                do n2 = 0, 1
                    do n3 = 0, 1
                    x = rcut +  DBLE(n1)*cell(n,1) + DBLE(n2)*cell(n,2) + DBLE(n3)*cell(n,3)
                    if ( x .gt. rout ) then
                        rout = x
                    endif
                    enddo
                enddo
            enddo
            v1(n) = rout     
        enddo
    
        !     now in negative direction
        
        do n = 1, 3
            rout = 0.d0
            do n1 = 0, 1
                do n2 = 0, 1
                    do n3 = 0, 1
                    x = -rcut +  DBLE(n1)*cell(n,1) + DBLE(n2)*cell(n,2) + DBLE(n3)*cell(n,3)
                    if ( x .lt. rout) then
                        rout = x
                    endif
                    enddo
                enddo
            enddo
            v2(n) = rout     
        enddo
        
        dxhcc1 = v1(1)
        dxhcc2 = v1(2)     
        dxhcc3 = v1(3)
        
        dxlcc1 = v2(1)
        dxlcc2 = v2(2)     
        dxlcc3 = v2(3)
        
        !number of gridboxes in positive/negative direction
        !(gridboxes oriented along positive axis)
        
        ngup1 = floor(v1(1)/rcut + eps)
        ngup2 = floor(v1(2)/rcut + eps)     
        ngup3 = floor(v1(3)/rcut + eps)
        
        nglow1 = floor(v2(1)/rcut - eps)
        nglow2 = floor(v2(2)/rcut - eps)    
        nglow3 = floor(v2(3)/rcut - eps)
    
    else
   ! settings for non-periodic calculation
   ! evaluate orthorhmobic box that contains the cluster


        v1(1:3) = -1.d10
        v2(1:3) = 1.d10

        do n = 1, natom
            if (posin(1,n) > v1(1) ) v1(1) = posin(1,n)
            if (posin(2,n) > v1(2) ) v1(2) = posin(2,n)
            if (posin(3,n) > v1(3) ) v1(3) = posin(3,n)
            if (posin(1,n) < v2(1) ) v2(1) = posin(1,n)
            if (posin(2,n) < v2(2) ) v2(2) = posin(2,n)
            if (posin(3,n) < v2(3) ) v2(3) = posin(3,n)
        enddo

        dxhcc1 = v1(1)
        dxhcc2 = v1(2)     
        dxhcc3 = v1(3)

        dxlcc1 = v2(1)
        dxlcc2 = v2(2)     
        dxlcc3 = v2(3)

        dxhc1 = 0.d0
        dxlc1 = 0.d0
        dxhc2 = 0.d0
        dxlc2 = 0.d0
        dxhc3 = 0.d0
        dxlc3 = 0.d0
   
        ngup1 = floor(v1(1)/rcut + eps)
        ngup2 = floor(v1(2)/rcut + eps)     
        ngup3 = floor(v1(3)/rcut + eps)
        
        nglow1 = floor(v2(1)/rcut - eps)
        nglow2 = floor(v2(2)/rcut - eps)    
        nglow3 = floor(v2(3)/rcut - eps)
    endif   
   
    if (periodic) then
    !convert positions from cartesian into direct coordinates
        call inv3x3(cell,cellinv)
        do n = 1, natom
            v1(1:3) = posin(1:3,n)
            do k = 1, 3
                v(k) = sum( cellinv(k,1:3)*v1(1:3) )
            enddo
            !       put atoms inside box
            do k=1,3
                v(k)  = v(k) - DBLE( floor(v(k) + eps) )
            enddo
            atompos(1:3,n) = v(1:3)
        enddo
    endif

    ! guess natombox (upper limit)
    if (periodic) then
        natombox = (hc1-lc1)*(hc2-lc2)*(hc3-lc3)*natom
    else
        natombox = natom
    endif

    allocate( origintmp(natombox) )
    allocate( postmp(3,natombox) )
    allocate( magtmp(3,natombox) )
    allocate( nshifttmp(3,natombox))

    if (periodic) then
    ! build cluster including periodic images   
        k = 0
        do n = 1, natom
            k = k + 1
            !       store direct coordinates first, convert to cartesian later
            r1(1:3) = atompos(1:3,n)
            r2(1:3)    = cell(1:3,1)*r1(1) + cell(1:3,2)*r1(2) + cell(1:3,3)*r1(3)
            postmp(1:3,n) = r2(1:3)
            magtmp(1:3,n) = magin(1:3,n)
            origintmp(k) = n
            nshifttmp(1,k) = 0
            nshifttmp(2,k) = 0
            nshifttmp(3,k) = 0
        enddo
        do k1 = lc1,hc1
            do k2 = lc2,hc2
                do k3 = lc3,hc3
                    if ((k1.eq.0).and.(k2.eq.0).and.(k3.eq.0)) goto 3000
                    do n = 1, natom
                    r1(1) = atompos(1,n) + dble(k1)
                    r1(2) = atompos(2,n) + dble(k2)
                    r1(3) = atompos(3,n) + dble(k3)
                    if ((r1(1).gt.dxlc1) .and. (r1(1).le.dxhc1) .and. &
                            &            (r1(2).gt.dxlc2) .and. (r1(2).le.dxhc2) .and. &
                            &            (r1(3).gt.dxlc3) .and. (r1(3).le.dxhc3) ) then
                        
                        !     now we take this atoms for the cluster
                        !               transform to cartesian coordinates
                        r2(1:3) = cell(1:3,1)*r1(1) + cell(1:3,2)*r1(2) + cell(1:3,3)*r1(3)
                        if( (r2(1).gt.dxlcc1) .and. (r2(1).le.dxhcc1) .and. &
                            &              (r2(2).gt.dxlcc2) .and. (r2(2).le.dxhcc2) .and. &
                            &              (r2(3).gt.dxlcc3) .and. (r2(3).le.dxhcc3)) then
                            !                 here we know that atom is inside of the 1. and 2. box
                            !                 around the original unit cell.
                            
                            k = k + 1
                            if (k>natombox) then
                                print *,'Too many atoms. Stopping.'
                            endif
                            postmp(1:3,k) = r2(1:3)
                            magtmp(1:3,k) = magtmp(1:3,n)
                            origintmp(k) = n
                            nshifttmp(1,k) = k1
                            nshifttmp(2,k) = k2
                            nshifttmp(3,k) = k3
                        endif
                    endif
                    enddo
        3000        continue
                enddo
            enddo
        enddo
    else
! non-periodic: only use the atoms of the cluster   
        k = 0
        do n = 1, natom
            k = k + 1
            postmp(1:3,n) = posin(1:3,n)
            magtmp(1:3,n) = magin(1:3,n)
            origintmp(k) = n
            nshifttmp(1,k) = 0
            nshifttmp(2,k) = 0
            nshifttmp(3,k) = 0
        enddo
    endif

    natombox = k

    if (.not. periodic ) then
        if (natom /= natombox ) then
            print *,'Strange things happen. Stopping.'
            stop
        endif
    endif

    allocate( pos(3,natombox) )
    allocate( mag(3,natombox) )
    allocate( nshift(3,natombox) )
    allocate( origin(natombox))

    do k = 1, 3
        pos(k,1:natombox)    =  postmp(k,1:natombox)
        mag(k,1:natombox)    =  magtmp(k,1:natombox)
        nshift(k,1:natombox) =  nshifttmp(k,1:natombox)
    enddo
        origin(1:natombox) = origintmp(1:natombox)

    deallocate( postmp )
    deallocate( magtmp )
    deallocate( origintmp )
    deallocate( nshifttmp )

    return
end subroutine buildcluster

!----------------------------------------------------------------------
      subroutine buildgrid
!----------------------------------------------------------------------
!     natgb: specify the number of an atom in a gridbox, 
!     get the number that the atom has in the cluster
!     specify number of atom in cluster,
!     get the number of the gridbox in which the atom is located
!----------------------------------------------------------------------
use neilistglobal
use grid        
implicit none
integer n,h1,h2,h3,nmax,neicount,k1,k2,k3,ng1, ng2, ng3
integer, allocatable :: nreshape(:)
integer, parameter :: nsafe_bg = 1000000
integer, parameter :: nsize = 10

call deallocate_grid

ng1 = ngup1 - nglow1 + 1
ng2 = ngup2 - nglow2 + 1
ng3 = ngup3 - nglow3 + 1
ngridbox = ng1*ng2*ng3
nmax = nsize*ngridbox 

call allocate_grid(natombox,ngridbox,nmax)

!     I. neighbourlist for gridboxes with rcut
neicount = 0
n = 0
do h1 = nglow1, ngup1
   do h2 = nglow2, ngup2
      do h3 = nglow3, ngup3
         n = n + 1
         neicount = neicount + 1
         if (neicount .ge. nmax) then
            allocate( nreshape(nmax) )
            nreshape = ngblist
            deallocate(ngblist)
            allocate(ngblist(nmax+nsafe_bg))
            ngblist(1:nmax) = nreshape(1:nmax)
            deallocate(nreshape)
            nmax = nmax + nsafe_bg
            nresize = nresize + 1
         endif
         ngbrangestart(n) = neicount
         ngblist(neicount) = n
         
         !           neighbours are all values of hi +- 1 (except boundary boxes)
         do k1 = h1-1, h1+1
            !             avoid outside regions
            if ((k1.ge.nglow1).and.(k1.le.ngup1)) then
               do k2 = h2-1,h2+1
                  if ((k2.ge.nglow2).and.(k2.le.ngup2)) then
                     do k3 = h3-1,h3+1
                        if ((k3.ge.nglow3).and.(k3.le.ngup3)) then
                           if ((k1.ne.h1).or.(k2.ne.h2).or.(k3.ne.h3)) then
                              neicount = neicount + 1
                              if (neicount .ge. nmax) then
                                 allocate( nreshape(nmax) )
                                 nreshape = ngblist
                                 deallocate(ngblist)
                                 allocate(ngblist(nmax+nsafe_bg))
                                 ngblist(1:nmax) = nreshape(1:nmax)
                                 deallocate(nreshape)
                                 nmax = nmax + nsafe_bg
                                 nresize = nresize + 1
                              endif
                              ngblist(neicount) = ng3*( ng2*(k1-nglow1) + (k2-nglow2) ) + (k3-nglow3) + 1
                           endif
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo
         ngbrangestop(n) = neicount
      enddo
   enddo
enddo
ngbnei = neicount
allocate( nreshape(nmax) )
nreshape = ngblist
deallocate(ngblist)
allocate(ngblist(neicount))
ngblist(1:neicount) = nreshape(1:neicount)
deallocate(nreshape)
!     
return
end subroutine buildgrid

!-----------------------------------------------------------------------     
      subroutine refreshgrid
!-----------------------------------------------------------------------
!     refresh atoms in boxes, 0(N) scheme
!-----------------------------------------------------------------------
use neilistglobal
use grid        
implicit none
integer atomi, n, m, ng1, ng2, ng3
integer ng(3)
double precision invrcut

ng1 = -nglow1 + ngup1 + 1
ng2 = -nglow2 + ngup2 + 1
ng3 = -nglow3 + ngup3 + 1

!
!     atoms in gridboxes ( this is all O(N) )
natgbstart = 0
natgbstop = 0

!     first run: sort all atoms in gridboxes and
!     count the number of atoms in the gridboxes,
!     generate inverse list in the second run.
invrcut = 1.d0/rcut

!     find origin of gridbox for each atom
do atomi = 1, natombox
   do n=1,3 
      ng(n) = floor( pos(n,atomi)*invrcut )
   enddo
!       calculate number of gridbox
   n = ng3*(ng2*(ng(1)-nglow1) + (ng(2)-nglow2)) + (ng(3)-nglow3) + 1
!       store atom in gridbox
   gbnat(atomi) = n
!       temporarily store number of atoms in each gridbox in natgbstop
   natgbstop(n) =  natgbstop(n) + 1
enddo

!     generate startvector natgbstart
m = 1
do n = 1, ngridbox
   natgbstart(n) = m
   m = m + natgbstop(n)
enddo

!     second run, generate natgbstop and natgb 
natgbstop = natgbstart - 1
do atomi = 1, natombox
   n = gbnat(atomi)
   natgbstop(n) = natgbstop(n) + 1
   natgb(natgbstop(n)) = atomi
enddo

!DEVELOP
!     consistency check I
do n = 1, ngridbox-1
   m = natgbstart(n+1) - natgbstop(n)
   if (m .ne. 1) then
      print *,'refreshgrid: wrong number of atoms in gridbox. Stopping.'
      stop
   endif
enddo
!     consistency check II
m = sum( natgbstop(1:ngridbox) - natgbstart(1:ngridbox)) + ngridbox
if( m .ne. natombox ) then
   print *,'refreshgrid: wrong total number of atoms in gridboxes. Stopping.'
   stop
endif
return
end subroutine refreshgrid

!--------------------------------------------------------------------------      
      subroutine refreshneighbours
!-----------------------------------------------------------------------
!     builds:
!     - short ranged neighbourlist
!     - inverse neighbourlist
!     - originnlist
!     - long neighbourlist
!-----------------------------------------------------------------------
use neilistglobal
use grid        
implicit none
integer natomboxmax,neicount,n1,n2,l,l1,l2,k2,m2,nbondmax
integer ndelta
integer, allocatable :: nreshape(:)
double precision, allocatable :: dreshape(:), d2reshape(:,:)
double precision, allocatable :: dreshape_mom(:), d2reshape_mom(:,:)
double precision rcutsq, r1(3),r2(3),d, dminat, mag_mom, mag_momi, mom1(3),mom2(3)

if (allocated(nrangestart)) deallocate(nrangestart)
if (allocated(nrangestop))  deallocate(nrangestop)
if (allocated(nlist))       deallocate(nlist)
if (allocated(dist))       deallocate(dist)
if (allocated(eist))       deallocate(eist)
if (allocated(dist_mom))   deallocate(dist_mom)
if (allocated(eist_mom))   deallocate(eist_mom)

if (allocated(dist_momi))   deallocate(dist_momi)
if (allocated(eist_momi))   deallocate(eist_momi)
!     value of nbondmax is not really important
!     (must be larger than 0)
nbondmax= 100
natomboxmax = natom*nbondmax
allocate(nrangestart(1:natom))
allocate(nrangestop(1:natom))
allocate(nlist(1:natomboxmax))
allocate(dist(1:natomboxmax))
allocate(eist(1:3,1:natomboxmax))
allocate(dist_mom(1:natomboxmax))
allocate(eist_mom(1:3,1:natomboxmax))

allocate(dist_momi(1:natom))
allocate(eist_momi(1:3,1:natom))

rcutsq = rcut*rcut

neicount = 0
nei2atommax = 0
distmin = rcut
distminmax2 = 0.d0
!     build neighbourlist of atoms in unit cell
!     an atom is NOT its own neighbour
ndelta = 10000
do n1 = 1, natom
   dminat = 1.d12
   l1 = gbnat(n1)
   r1(1:3) = pos(1:3,n1)
   
   mom1(1:3) = mag(1:3,n1)
   mag_momi = sqrt(sum(mom1(1:3)*mom1(1:3)) ) 
   
   nrangestart(n1) = neicount + 1
!       inner loop: loop all atoms in neighbouring gridboxes
   do k2 = ngbrangestart(l1), ngbrangestop(l1)
      l2 = ngblist(k2)
      do m2 =  natgbstart(l2), natgbstop(l2)
         n2 = natgb(m2)
         if ( n1 .ne. n2) then
            r2(1:3) = pos(1:3,n2)
            r2(1:3) = r2(1:3) - r1(1:3)
            mom2(1:3) = mag(1:3,n2)
            call fastdist(r2,rcut,rcutsq,l)
            if ( l .eq. 0 ) then
               d = sqrt(sum(r2(1:3)*r2(1:3)) )
               mag_mom = sqrt(sum(mom2(1:3)*mom2(1:3)) )
               neicount = neicount + 1
               if ( neicount .ge. natomboxmax) then
                  allocate( nreshape(natomboxmax) )
                  allocate( dreshape(natomboxmax) )
                  allocate( dreshape_mom(natomboxmax) )
                  allocate( d2reshape(1:3,natomboxmax) )
                  allocate( d2reshape_mom(1:3,natomboxmax) )
                  nreshape(1:natomboxmax) = nlist(1:natomboxmax)
                  dreshape(1:natomboxmax) = dist(1:natomboxmax)
                  dreshape_mom(1:natomboxmax) = dist_mom(1:natomboxmax)
                  d2reshape(1:3,1:natomboxmax) = eist(1:3,1:natomboxmax)
                  d2reshape_mom(1:3,1:natomboxmax) = eist_mom(1:3,1:natomboxmax)
                  deallocate(nlist)
                  deallocate(dist)
                  deallocate(dist_mom)
                  deallocate(eist)
                  deallocate(eist_mom)
                  allocate(nlist(natomboxmax+ndelta))
                  allocate(dist(natomboxmax+ndelta))
                  allocate(dist_mom(natomboxmax+ndelta))
                  allocate(eist(3,natomboxmax+ndelta))
                  allocate(eist_mom(3,natomboxmax+ndelta))
                  nlist(1:natomboxmax) = nreshape(1:natomboxmax)
                  dist(1:natomboxmax) = dreshape(1:natomboxmax)
                  dist_mom(1:natomboxmax) = dreshape_mom(1:natomboxmax)
                  eist(1:3,1:natomboxmax) = d2reshape(1:3,1:natomboxmax)
                  eist_mom(1:3,1:natomboxmax) = d2reshape_mom(1:3,1:natomboxmax)
                  deallocate(nreshape)
                  deallocate(dreshape)
                  deallocate(dreshape_mom)
                  deallocate(d2reshape)
                  deallocate(d2reshape_mom)
                  natomboxmax = natomboxmax + ndelta     
                  nresize = nresize + 1
               endif
               if ( sum(abs(r2(1:3))) .lt. eps) then
                  print *,'Atoms with identical positions.'
                  print *, r1, r2
                  stop
               endif
               nlist(neicount) = n2
               dist(neicount) = d
               dist_mom(neicount) = mag_mom
               if (d < distmin) then
                  distmin = d
               endif
               if (d < dminat ) then
                  dminat = d
               endif
               eist(1:3,neicount) = r2(1:3)/d
               if (mag_mom.eq.0.d0) then
                  eist_mom(1:2,neicount) = 0.d0
                  eist_mom(3,neicount) = 1.d0
               else
                  eist_mom(1:3,neicount) = mom2(1:3)/mag_mom
               endif
            endif
         endif
      enddo
   enddo 
   
   dist_momi(n1) = mag_momi
   if (mag_momi.eq.0.d0) then
      eist_momi(1:2,n1) = 0.d0
      eist_momi(3,n1) = 1.d0
   else
      eist_momi(1:3,n1) = mom1(1:3)/mag_momi
   endif
    
   nrangestop(n1) = neicount
   if (nrangestop(n1) - nrangestart(n1) + 1 > nei2atommax ) then
      nei2atommax = nrangestop(n1) - nrangestart(n1) + 1
   endif
   if (dminat > distminmax2 ) then
      distminmax2 = dminat
   endif
enddo

call deallocate_grid

return
end subroutine refreshneighbours
!----------------------------------------------------------------------
      subroutine north(v1,v2,vn)
!--------------------------------------------------------------
!     calculate orthonormal vector
!--------------------------------------------------------------
implicit none
double precision, intent(in)::  v1(3), v2(3)
double precision, intent(out):: vn(3)
double precision x

vn(1) = v1(2)*v2(3) - v1(3)*v2(2)
vn(2) = v1(3)*v2(1) - v1(1)*v2(3)
vn(3) = v1(1)*v2(2) - v1(2)*v2(1)

x = sqrt(sum(vn(1:3)*vn(1:3)))

if(x.lt.1.d-9) then
   print *,'Very short vector. Stopping.'
   stop
endif

vn = vn/x

return
end subroutine north
     

!--------------------------------------------------------------
      subroutine fastdist(v,rcut,rcutsq,ncut)
!--------------------------------------------------------------
implicit none
double precision, intent(in)::  v(1:3), rcut, rcutsq
integer, intent(out):: ncut
double precision x

ncut = 1

if (abs(v(1)) .lt. rcut) then
   if (abs(v(2)) .lt. rcut) then
      if (abs(v(3)) .lt. rcut) then
     
         x = sum(v(1:3)*v(1:3))
         if ( x .lt. rcutsq) then
            ncut = 0
         endif
     
      endif
   endif
endif

return
end subroutine fastdist
!--------------------------------------------------------------
      subroutine inv3x3(a,ainv)
!--------------------------------------------------------------
implicit none
double precision,intent(in) :: a(1:3,1:3)
double precision,intent(out):: ainv(1:3,1:3)
double precision det
integer n,m
double precision x

ainv(1,1) =     a(3,3)*a(2,2)-a(3,2)*a(2,3)  
ainv(1,2) =  -( a(3,3)*a(1,2)-a(3,2)*a(1,3) )
ainv(1,3) =     a(2,3)*a(1,2)-a(2,2)*a(1,3)

ainv(2,1) =  -( a(3,3)*a(2,1)-a(3,1)*a(2,3) ) 
ainv(2,2) =     a(3,3)*a(1,1)-a(3,1)*a(1,3)
ainv(2,3) =  -( a(2,3)*a(1,1)-a(2,1)*a(1,3) )

ainv(3,1) =     a(3,2)*a(2,1)-a(3,1)*a(2,2)
ainv(3,2) =  -( a(3,2)*a(1,1)-a(3,1)*a(1,2) )
ainv(3,3) =     a(2,2)*a(1,1)-a(2,1)*a(1,2)

det = a(1,1)*(a(3,3)*a(2,2)-a(3,2)*a(2,3)) &
&     -a(2,1)*(a(3,3)*a(1,2)-a(3,2)*a(1,3)) &
&     +a(3,1)*(a(2,3)*a(1,2)-a(2,2)*a(1,3))

if ( abs(det).gt.1.d-9) then
   ainv = ainv/det
else
   print *,'Problem inv3x3: ZERO DETERMINANT.'
   stop
endif
    
!     test
     
do n = 1, 3
   do m = 1, 3
      x = sum( ainv(n,1:3)*a(1:3,m) )
      if ( n .eq. m) then
         if ( abs(x - 1.d0) .gt. 1.d-9) then
            print *,'Problem 1 in inv3x3.'
            stop
         endif
      else
         if ( abs(x) .gt. 1.d-9) then
            print *,'Problem 2 in inv3x3.'
            stop
         endif
      endif
   enddo
enddo
     
return
end subroutine inv3x3
!--------------------------------------------------------------
subroutine threepoint
!--------------------------------------------------------------
use neilistglobal      
implicit none
integer iat, kat, jat, nbondmax, ij, ik, l, natomboxmax, ndelta, neicount, nei2count
double precision d, rcutsq
double precision r1(3), r2(3)
integer, allocatable :: nreshape(:,:)
double precision, allocatable :: dreshape(:,:), d2reshape(:,:,:)

if (r3cut > rcut ) then
   print *,'r3cut must be smaller or equal to rcut.'
   stop
endif
rcutsq = r3cut*r3cut

nbondmax= 100
natomboxmax = natom*nbondmax
ndelta = 10000
if (allocated(n3rangestart)) deallocate(n3rangestart)
if (allocated(n3rangestop))  deallocate(n3rangestop)
if (allocated(n3list))       deallocate(n3list)
if (allocated(d3ist))       deallocate(d3ist)
if (allocated(e3ist))       deallocate(e3ist)

if (nrangestop(natom) > 0) then
   nei2count = nrangestop(natom)
else
   nei2count = 1
endif
allocate(n3rangestart(1:nei2count))
allocate(n3rangestop(1:nei2count))
allocate(n3list(1:natomboxmax,1:2))
allocate(d3ist(1:natomboxmax,1:3))
allocate(e3ist(1:natomboxmax,1:3,1:3))


! loop all neighbours
neicount = 0
nei3atommax = 0
do iat = 1, natom
   do ij = nrangestart(iat), nrangestop(iat)
        n3rangestart(ij) = neicount + 1 
        if ( dist(ij) <= r3cut ) then
           jat = nlist(ij)
           r1(1:3) = pos(1:3,jat)
           do ik = nrangestart(iat), nrangestop(iat)
              if ( dist(ik) <= r3cut ) then            
                 kat = nlist(ik)
                 !take triangles with jat /= kat         
                 if (jat /= kat) then
                    r2(1:3) = pos(1:3,kat) - r1(1:3)
                    call fastdist(r2,r3cut,rcutsq,l)
                    if ( l .eq. 0 ) then
                       d = sqrt( sum(r2(1:3)*r2(1:3)) )
                       neicount = neicount + 1
                       
                       if ( neicount .ge. natomboxmax) then
                          allocate( nreshape(natomboxmax,2) )
                          allocate( dreshape(natomboxmax,3) )
                           allocate( d2reshape(natomboxmax,3,3) )
                          nreshape(1:natomboxmax,1:2) = n3list(1:natomboxmax,1:2)                        
                          dreshape(1:natomboxmax,1:3) = d3ist(1:natomboxmax,1:3)
                          d2reshape(1:natomboxmax,1:3,1:3) = e3ist(1:natomboxmax,1:3,1:3)
                        deallocate(n3list)
                        deallocate(d3ist)
                        deallocate(e3ist)
                        allocate(n3list(natomboxmax+ndelta,2))     
                        allocate(d3ist(natomboxmax+ndelta,3))
                        allocate(e3ist(natomboxmax+ndelta,3,3))
                        n3list(1:natomboxmax,1:2) = nreshape(1:natomboxmax,1:2)                        
                        d3ist(1:natomboxmax,1:3) = dreshape(1:natomboxmax,1:3)
                        e3ist(1:natomboxmax,1:3,1:3) = d2reshape(1:natomboxmax,1:3,1:3) 
                        deallocate(nreshape)                        
                        deallocate(dreshape)
                        deallocate(d2reshape)
                        natomboxmax = natomboxmax + ndelta     
                        nresize = nresize + 1
                     endif
                     
                     n3list(neicount,1) = jat
                     n3list(neicount,2) = kat
                     if (ij == ik ) then
                        print *,'Something wrong here. Subroutine threepoint.'
                        stop
                     endif
                     d3ist(neicount,1) = dist(ij)
                     d3ist(neicount,2) = dist(ik)
                     d3ist(neicount,3) = d
                     e3ist(neicount,1,1:3) = eist(1:3,ij)  ! iat -> jat
                     e3ist(neicount,2,1:3) = eist(1:3,ik)  ! iat -> kat
                     e3ist(neicount,3,1:3) = r2(1:3)/d     ! jat -> kat                
                     
                  endif
               endif
            endif
         enddo
      endif
      n3rangestop(ij) = neicount
      if (n3rangestop(ij) - n3rangestart(ij) + 1 > nei3atommax ) then
         nei3atommax = n3rangestop(ij) - n3rangestart(ij) + 1
      endif
   enddo
enddo

end subroutine threepoint

!--------------------------------------------------------------
subroutine fourpoint
!--------------------------------------------------------------
use neilistglobal      
implicit none
integer iat, kat, jat, lat, nbondmax, ij, ijk, il, l1, l2, l3, natomboxmax, ndelta, neicount, nei3count, nei2count
double precision dil, djl, dkl, rcutsq
double precision r1(3), r2(3), r3(3)
double precision riat(3), rjat(3), rkat(3), rlat(3)
integer, allocatable :: nreshape(:,:)
double precision, allocatable :: dreshape(:,:), d2reshape(:,:,:)

if ( r4cut > r3cut ) then
   print *,'r4cut must be smaller or equal to r3cut.'
   stop
endif
rcutsq = r4cut*r4cut

nbondmax= 100
natomboxmax = natom*nbondmax
ndelta = 1000
if (allocated(n4rangestart)) deallocate(n4rangestart)
if (allocated(n4rangestop))  deallocate(n4rangestop)
if (allocated(n4list))       deallocate(n4list)
if (allocated(d4ist))       deallocate(d4ist)
if (allocated(e4ist))       deallocate(e4ist)
if (nrangestop(natom) > 0) then
   nei2count = nrangestop(natom)
else
   nei2count = 1
endif
if (n3rangestop(nei2count) > 0) then
   nei3count = n3rangestop(nei2count)
else
   nei3count = 1
endif

allocate(n4rangestart(1:nei3count))
allocate(n4rangestop(1:nei3count))
allocate(n4list(1:natomboxmax,1:3))
allocate(d4ist(1:natomboxmax,1:6))
allocate(e4ist(1:natomboxmax,1:6,1:3))

! loop all neighbours
neicount = 0
nei4atommax = 0
do iat = 1, natom
   riat(1:3) = pos(1:3,iat)
   do ij = nrangestart(iat), nrangestop(iat)
      do ijk = n3rangestart(ij), n3rangestop(ij)
         n4rangestart(ijk) = neicount + 1
         if ((d3ist(ijk,1) <= r4cut).and.(d3ist(ijk,2) <= r4cut).and.(d3ist(ijk,3) <= r4cut)) then
            jat = n3list(ijk,1)
            kat = n3list(ijk,2)
            rjat(1:3) = pos(1:3,jat)
            rkat(1:3) = pos(1:3,kat)
            
            do il = nrangestart(iat), nrangestop(iat)        
               lat = nlist(il)      
               if ( (iat /= lat) .and. (jat /= lat)  .and. (kat /= lat) ) then
                  rlat(1:3) = pos(1:3,lat)
                  r1(1:3) = rlat(1:3) - riat(1:3)
                  call fastdist(r1,r4cut,rcutsq,l1)
                  if ( l1 == 0 ) then
                     r2(1:3) = rlat(1:3) - rjat(1:3)
                     call fastdist(r2,r4cut,rcutsq,l2)
                     if ( l2 == 0 ) then
                        r3(1:3) = rlat(1:3) - rkat(1:3)
                        call fastdist(r3,r4cut,rcutsq,l3)
                        if ( l3 == 0 ) then
                           dil = sqrt(sum(r1(1:3)*r1(1:3)) )                       
                           djl = sqrt(sum(r2(1:3)*r2(1:3)) )
                           dkl = sqrt(sum(r3(1:3)*r3(1:3)) )
                           
                           neicount = neicount + 1
                           
                           if ( neicount >= natomboxmax) then
                              allocate( nreshape(natomboxmax,3) )
                              allocate( dreshape(natomboxmax,6) )
                              allocate( d2reshape(natomboxmax,6,3) )
                              nreshape(1:natomboxmax,1:3) = n4list(1:natomboxmax,1:3)
                              dreshape(1:natomboxmax,1:6) = d4ist(1:natomboxmax,1:6)
                              d2reshape(1:natomboxmax,1:6,1:3) = e4ist(1:natomboxmax,1:6,1:3)
                              deallocate(n4list)
                              deallocate(d4ist)
                              deallocate(e4ist)
                              allocate(n4list(natomboxmax+ndelta,3))
                              allocate(d4ist(natomboxmax+ndelta,6))
                              allocate(e4ist(natomboxmax+ndelta,6,3))
                              n4list(1:natomboxmax,1:3) = nreshape(1:natomboxmax,1:3)
                              d4ist(1:natomboxmax,1:6) = dreshape(1:natomboxmax,1:6)                  
                              e4ist(1:natomboxmax,1:6,1:3) = d2reshape(1:natomboxmax,1:6,1:3)
                              deallocate(nreshape)
                              deallocate(dreshape)
                              deallocate(d2reshape)
                              natomboxmax = natomboxmax + ndelta     
                              nresize = nresize + 1
                           endif
                           
                           n4list(neicount,1) = jat
                           n4list(neicount,2) = kat
                           n4list(neicount,3) = lat
                           
                           d4ist(neicount,1) = d3ist(ijk,1)       ! dij
                           d4ist(neicount,2) = d3ist(ijk,2)       ! dik
                           d4ist(neicount,3) = dil               ! dil
                           d4ist(neicount,4) = d3ist(ijk,3)       ! djk
                           d4ist(neicount,5) = djl               ! djl
                           d4ist(neicount,6) = dkl               ! dkl

                           e4ist(neicount,1,1:3) = e3ist(ijk,1,1:3)       ! i->j
                           e4ist(neicount,2,1:3) = e3ist(ijk,2,1:3)       ! i->k
                           e4ist(neicount,3,1:3) = r1(1:3)/dil               ! i->l
                           e4ist(neicount,4,1:3) = e3ist(ijk,3,1:3)       ! j->k
                           e4ist(neicount,5,1:3) = r2(1:3)/djl               ! j->l
                           e4ist(neicount,6,1:3) = r3(1:3)/dkl               ! k->l 
                           
                        endif
                     endif
                  endif
               endif
               
            enddo
         endif
         n4rangestop(ijk) = neicount
         if (n4rangestop(ijk) - n4rangestart(ijk) + 1 > nei4atommax ) then
            nei4atommax = n4rangestop(ijk) - n4rangestart(ijk) + 1
         endif
      enddo
   enddo
enddo

end subroutine fourpoint

!----------------------------------------------------------------------------

