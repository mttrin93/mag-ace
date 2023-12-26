
#ifdef MPIPARALLEL

module mpiglobal
    use mpi
    integer status(MPI_Status_size)
    integer ierr, procid, rootid, numprocs
    integer nstrucproc
    integer nlow, nup, nsize, tag, id
    integer send_data_tag
    logical exloop

end module mpiglobal

#endif 
!---------------------------------------------------------------------------

module core
    integer verblevel
    integer debuglevel
end module core

!---------------------------------------------------------------------------
      
module rw
    integer kmaxt
    integer nradbaset_atomic, nradbaset_magnetic, nradial2t_magnetic, l2maxt_magnetic
    integer nradial3t_atomic, l3maxt_atomic, nradial3t_magnetic, l3maxt_magnetic
    integer nradial4t_atomic, l4maxt_atomic, nradial4t_magnetic, l4maxt_magnetic
    integer nradial5t_atomic, l5maxt_atomic, nradial5t_magnetic, l5maxt_magnetic
    integer nradial6t_atomic, l6maxt_atomic, nradial6t_magnetic, l6maxt_magnetic
    integer nradial7t_atomic, l7maxt_atomic, nradial7t_magnetic, l7maxt_magnetic
    integer ndensityt, npott
    integer nfembt, radtypet, radtypet_atomic, radtypet_magnetic
end module rw

!-----------------------------------------------------------------------------

module solvermod
    integer npara
    double precision, allocatable :: paraall(:)
    integer, allocatable :: paramap(:)
end module solvermod

!---------------------------------------------------------------------------

module global
    implicit none

    integer nparatot, npararad
    integer npararadtot, nparaalltot

    character(20) task
    character(60) xc, solver

    logical do3, do4
    logical doforce, donumforce, iszerror, dotorque
    double precision maxforceerr, maxforceerrZ
    integer nmaxforceerr, nmaxforceerrZ
    logical dotest, dotrain, dotrainrad, dorun, doneinit, doquick, doprintbasis, dosetupbasis
    double precision forceweight
    double precision deltaenethreshold, deltaenethresholdlow, wlow
    logical takeperiodic

    integer ncount, mncount, pcount, natomsmaxtakeclu, natomsmaxtakebulk, nclimb
    integer nstrucfit, nstructest
    double precision fracclimb

    double precision jloss, jlossE, jlossF
    double precision xcountatom
    double precision sigmashake
    integer nshakemax

    ! energy and forces
    double precision ecalc
    ! cehck if this is required here
    double precision, allocatable ::  ecalcvec(:)
    double precision, allocatable ::  fcalcvec(:,:)
    double precision, allocatable ::  fijcalc(:,:)
    double precision, allocatable ::  tcalcvec(:,:)
    double precision, allocatable ::  tijcalc(:,:)
    double precision, allocatable ::  tijcalc0(:,:)
    double precision, allocatable ::  eicalc(:)
    double precision, allocatable ::  edfcut(:)
    double precision emincalc

    ! natom in neilistglobal -> change
    integer nstruc, nstructot, nfailed, natomsmax, natomsmaxcheck, ntotatoms, natomboxmax
    integer nei2countmax, nei3countmax, nei4countmax

    ! for MPI integration
    integer nstrucglobal, ntotatomsglobal
    integer, allocatable :: natomsglobal(:)

    integer, allocatable :: natoms(:)
    integer, allocatable :: strucaddress(:,:)

    integer, allocatable :: originpack(:)
    integer, allocatable :: ol(:)
    integer, allocatable :: natomboxpackstart(:)
    integer, allocatable :: natomboxpackstop(:)

    integer, allocatable :: n3start(:), n3stop(:)
    integer, allocatable :: n4start(:), n4stop(:) 
    integer, allocatable :: n2startpack(:), n2stoppack(:)
    integer, allocatable :: n3startpack(:), n3stoppack(:)
    integer, allocatable :: n4startpack(:), n4stoppack(:) 


    integer, allocatable :: n3l(:,:)
    integer, allocatable :: n4l(:,:)      
    integer, allocatable :: n2listpack(:)
    integer, allocatable :: n3listpack(:,:)
    integer, allocatable :: n4listpack(:,:)


    double precision, allocatable :: d3l(:,:)
    double precision, allocatable :: d4l(:,:)
    double precision, allocatable :: e3l(:,:,:)
    double precision, allocatable :: e4l(:,:,:)      
    double precision, allocatable :: d2listpack(:)
    double precision, allocatable :: d2listpack_mom(:)
    double precision, allocatable :: d2listpack_momi(:)
    double precision, allocatable :: d3listpack(:,:)
    double precision, allocatable :: d4listpack(:,:)
    double precision, allocatable :: e2listpack(:,:)
    double precision, allocatable :: e2listpack_mom(:,:)
    double precision, allocatable :: e2listpack_momi(:,:)
    double precision, allocatable :: e3listpack(:,:,:)
    double precision, allocatable :: e4listpack(:,:,:)
    double precision, allocatable :: dmin(:)
    double precision, allocatable :: dminmax(:)     
    double precision, allocatable :: dminglobal(:)
    double precision, allocatable :: dminmaxglobal(:)      
    double precision, allocatable :: atvol(:)

    double precision, allocatable :: fr2listpack(:,:)
    double precision, allocatable :: fr2sumpack(:,:)
    double precision, allocatable :: fr3sumpack(:,:,:,:)

    !      double precision, allocatable :: base2corrpack(:,:,:)
    double precision, allocatable :: b1pack(:,:)
    double precision, allocatable :: b2pack(:,:)
    double precision, allocatable :: b2hcpack(:)
    double precision, allocatable :: b3pack(:,:)
    double precision, allocatable :: b4pack(:,:)
    double precision, allocatable :: b5pack(:,:)
    double precision, allocatable :: b6pack(:,:)
    double precision, allocatable :: b7pack(:,:)      

    double precision, allocatable :: b1packall(:,:)
    double precision, allocatable :: b2packall(:,:)
    double precision, allocatable :: b2hcpackall(:)
    double precision, allocatable :: b3packall(:,:)
    double precision, allocatable :: b4packall(:,:)
    double precision, allocatable :: b5packall(:,:)
    double precision, allocatable :: b6packall(:,:)
    double precision, allocatable :: b7packall(:,:)      

    double complex, allocatable :: db2pack(:,:,:)
    double complex, allocatable :: db3pack(:,:,:,:,:)
    double complex, allocatable :: db3pack1(:,:,:,:)
    double complex, allocatable :: db3pack2(:,:,:,:)
    double complex, allocatable :: db4pack(:,:,:,:,:)     
    double complex, allocatable :: db5pack(:,:,:,:,:)
    double complex, allocatable :: db6pack(:,:,:,:,:)
    double complex, allocatable :: db7pack(:,:,:,:,:)  

    double complex, allocatable :: db2pack_0(:,:,:)
    double complex, allocatable :: db3pack_0(:,:,:)
    double complex, allocatable :: db4pack_0(:,:,:)     
    double complex, allocatable :: db5pack_0(:,:,:)
    double complex, allocatable :: db6pack_0(:,:,:)

    double precision cutoff, cutoff3, cutoff4, cutoff_magnetic

    double precision dcutoff
    double precision weightcut, weightexp, weightpre
    double precision deltawE, deltawF, deltawM
    double precision weightbulk, weight2, weight3, weight4
    double precision smearing

    ! linear and quadratic regularization weights
    double precision alphalin, alphasq

    ! for storing information of many structures
    double precision, allocatable :: atompos(:,:,:)
    integer, allocatable :: occpos(:,:)
    double precision, allocatable :: force(:,:,:)
    double precision, allocatable :: magmom(:,:,:)
    double precision, allocatable :: torque(:,:,:)
    double precision, allocatable :: unitcell(:,:,:)
    double precision, allocatable :: energy(:)
    !     double precision, allocatable :: r(:,:)  

    double precision, allocatable :: weightread(:)
    double precision, allocatable :: sqrtweightE(:), sqrtweightF(:,:)
    double precision, allocatable :: deglobal(:), dfglobal(:,:)
    double precision, allocatable :: eglobal(:), fglobal(:,:)
    logical, allocatable :: periodiccellglobal(:)
    double precision av, referr

    character(60) strucname
    character(60), allocatable :: name(:)
    character(60), allocatable :: nameglobal(:)
    character(60), allocatable :: runtype(:)      

    logical, allocatable :: groupstart(:)
    logical, allocatable :: forceavail(:)
    logical, allocatable :: magmomavail(:)
    logical, allocatable :: torqueavail(:)
    logical, allocatable :: periodiccell(:)

    double precision, allocatable :: eneuncorrected(:)
    double precision, allocatable :: enefreeene(:)
    double precision, allocatable :: eneTtozero(:)
    double precision, allocatable :: cputime(:)

    double precision oldrms(1:10), avrms

    double precision rho2att_0min, rho2att_0max
    double precision rho2rep_0min, rho2rep_0max

    double precision paraweightfunc(1:10)

    double precision wscoreRnl(1:3), totscoreRnl
    double precision wscoreRnl_atomic(1:3), wscoreRnl_magnetic(1:3) 
    double precision totscoreRnl_atomic, totscoreRnl_magnetic

      
end module global

!--------------------------------------------------------------

module functionparameters
        
    implicit none

    integer ndensity, npot, nexp, nfemb

    integer lmax, lmaxhalf, lmaxsq
    integer lmax_atomic, lmaxsq_atomic, lmaxhalf_atomic
    integer lmax_magnetic, lmaxsq_magnetic, lmaxhalf_magnetic
    integer nvar, nvar1, nvar2, nvar3, nvar4, nvar5, nvar6, nvar7
    integer nvar3_tmp, nvar4_tmp, nvar5_tmp, nvar6_tmp, nvar7_tmp
    integer nvar2l, nvar3l, nvar4l, nvar5l, nvar6l, nvar7l

    integer nvar3l_atomic, nvar3l_magnetic
    integer nvar4l_atomic, nvar4l_magnetic
    integer nvar5l_atomic, nvar5l_magnetic
    integer nvar6l_atomic, nvar6l_magnetic

    integer radtype, radtype_atomic, radtype_magnetic
    double precision radpara
    character(40) radialbasis
    double precision, allocatable :: lambda(:,:,:)

    double precision, allocatable :: c0(:,:,:)
    logical, allocatable :: fitc0(:,:,:)

    integer nradbase, nradial

    integer nradbase_atomic, nradbase_magnetic
    integer nradial2_magnetic, l2max_magnetic

    integer nradial3_atomic, l3max_atomic
    integer nradial3_magnetic, l3max_magnetic
    double precision, allocatable :: c2(:,:,:)
    logical, allocatable :: fitc2(:,:,:)
    double precision, allocatable :: c3(:,:,:)
    logical, allocatable :: fitc3(:,:,:)

    integer nradial4, l4max
    integer nradial4_atomic, l4max_atomic
    integer nradial4_magnetic, l4max_magnetic       
    double precision, allocatable :: c4(:,:,:)
    logical, allocatable :: fitc4(:,:,:)

    integer nradial5, l5max
    integer nradial5_atomic, l5max_atomic
    integer nradial5_magnetic, l5max_magnetic
    double precision, allocatable :: c5(:,:,:)
    logical, allocatable :: fitc5(:,:,:)

    integer nradial6, l6max
    integer nradial6_atomic, l6max_atomic
    integer nradial6_magnetic, l6max_magnetic
    double precision, allocatable :: c6(:,:,:)
    logical, allocatable :: fitc6(:,:,:)

    integer nradial7, l7max
    integer nradial7_atomic, l7max_atomic
    integer nradial7_magnetic, l7max_magnetic
    double precision, allocatable :: c7(:,:,:)
    logical, allocatable :: fitc7(:,:,:)

    integer nvarp1
    double precision, allocatable :: p1(:,:,:)
    logical, allocatable :: fitp1(:,:,:)

    double precision, allocatable :: crad(:,:,:,:,:)
    double precision, allocatable :: crad_atomic(:,:,:,:,:)
    double precision, allocatable :: crad_magnetic(:,:,:,:,:)
    logical, allocatable :: fitcrad(:,:,:,:,:)
    logical, allocatable :: fitcrad_atomic(:,:,:,:,:)
    logical, allocatable :: fitcrad_magnetic(:,:,:,:,:)

    double precision, allocatable :: corepara(:,:,:)
    logical, allocatable :: fitcorepara(:,:,:)

end module functionparameters

!-----------------------------------------------------------------------------------

module tables
                
    integer, allocatable:: shindex(:,:)
    integer, allocatable:: shindex_atomic(:,:)
    integer, allocatable:: shindex_magnetic(:,:)
    integer, allocatable:: shback(:)
    integer, allocatable:: shback_atomic(:)
    integer, allocatable:: shback_magnetic(:)

    integer, allocatable:: b1index(:)
    integer, allocatable:: b1occ(:)        

    integer, allocatable:: b2index(:,:)
    integer, allocatable:: b2lindex(:)
    integer, allocatable:: b2occ(:)
    integer, allocatable:: b2map(:)
    integer, allocatable:: madd2(:)
    integer, allocatable:: b2lm(:,:)

    integer, allocatable:: b3index(:,:)
    integer, allocatable:: b3index_tmp(:,:)
    integer, allocatable:: b3lindex(:,:)
    integer, allocatable:: b3occ(:,:)
    integer, allocatable:: b3occ_tmp(:,:)
    integer, allocatable:: b3map(:)
    integer, allocatable:: madd3(:)
    integer, allocatable:: b3lm(:,:,:)
    integer, allocatable:: b3lindex_atomic(:)
    integer, allocatable:: b3lindex_magnetic(:,:)
    integer, allocatable:: b3map_atomic(:)
    integer, allocatable:: b3map_magnetic(:)
    integer, allocatable:: b3lm_atomic(:,:)
    integer, allocatable:: b3lm_magnetic(:,:,:) 
    integer, allocatable:: madd3_atomic(:)
    integer, allocatable:: madd3_magnetic(:) 

    integer, allocatable:: b4index(:,:)
    integer, allocatable:: b4index_tmp(:,:)
    integer, allocatable:: b4lindex(:,:)
    integer, allocatable:: b4occ(:,:)
    integer, allocatable:: b4occ_tmp(:,:)
    integer, allocatable:: b4map(:)
    integer, allocatable:: madd4(:)
    integer, allocatable:: b4lm(:,:,:)
    integer, allocatable:: b4nint(:)
    integer, allocatable:: b4nint_tmp(:)
    integer, allocatable:: b4lindex_atomic(:,:)
    integer, allocatable:: b4lindex_magnetic(:,:)
    integer, allocatable:: b4map_atomic(:)
    integer, allocatable:: b4map_magnetic(:)
    integer, allocatable:: b4lm_atomic(:,:,:)
    integer, allocatable:: b4lm_magnetic(:,:,:) 
    integer, allocatable:: madd4_atomic(:)
    integer, allocatable:: madd4_magnetic(:)
        
    integer, allocatable:: b5index(:,:)
    integer, allocatable:: b5index_tmp(:,:)
    integer, allocatable:: b5lindex(:,:)
    integer, allocatable:: b5occ(:,:)
    integer, allocatable:: b5occ_tmp(:,:)
    integer, allocatable:: b5map(:)
    integer, allocatable:: madd5(:)
    integer, allocatable:: b5lm(:,:,:)
    integer, allocatable:: b5nint(:) 
    integer, allocatable:: b5nint_tmp(:)
    integer, allocatable:: b5lindex_atomic(:,:)
    integer, allocatable:: b5lindex_magnetic(:,:)
    integer, allocatable:: b5map_atomic(:)
    integer, allocatable:: b5map_magnetic(:)
    integer, allocatable:: b5lm_atomic(:,:,:)
    integer, allocatable:: b5lm_magnetic(:,:,:) 
    integer, allocatable:: madd5_atomic(:)
    integer, allocatable:: madd5_magnetic(:)

    integer, allocatable:: b6index(:,:)
    integer, allocatable:: b6index_tmp(:,:)
    integer, allocatable:: b6lindex(:,:)
    integer, allocatable:: b6occ(:,:)
    integer, allocatable:: b6occ_tmp(:,:)
    integer, allocatable:: b6map(:)
    integer, allocatable:: madd6(:)
    integer, allocatable:: b6lm(:,:,:)
    integer, allocatable:: b6nint(:)
    integer, allocatable:: b6nint_tmp(:)
    integer, allocatable:: b6lindex_atomic(:,:)
    integer, allocatable:: b6lindex_magnetic(:,:)
    integer, allocatable:: b6map_atomic(:)
    integer, allocatable:: b6map_magnetic(:)
    integer, allocatable:: b6lm_atomic(:,:,:)
    integer, allocatable:: b6lm_magnetic(:,:,:) 
    integer, allocatable:: madd6_atomic(:)
    integer, allocatable:: madd6_magnetic(:)

    integer, allocatable:: b7index(:,:)
    integer, allocatable:: b7index_tmp(:,:)
    integer, allocatable:: b7lindex(:,:)
    integer, allocatable:: b7occ(:,:)
    integer, allocatable:: b7occ_tmp(:,:)
    integer, allocatable:: b7map(:)
    integer, allocatable:: madd7(:)
    integer, allocatable:: b7lm(:,:,:)
    integer, allocatable:: b7nint(:) 
    integer, allocatable:: b7nint_tmp(:)

    double precision, allocatable :: fac(:)
    double precision, allocatable :: alm(:),blm(:),clm(:),cl(:),dl(:),el(:)
    double precision, allocatable :: plm(:),splm(:),dplm(:)

    double precision, allocatable :: cgcoeff(:,:,:,:,:,:)
    double precision, allocatable :: w3j(:,:,:,:,:,:)

    double precision, allocatable :: cg02(:,:)
    double precision, allocatable :: cg03(:,:,:,:)
    double precision, allocatable :: cg04(:,:,:,:,:)
    double precision, allocatable :: cg05(:,:,:,:,:,:)
    double precision, allocatable :: cg06(:,:,:,:,:,:,:)

    double precision, allocatable :: cg01_magnetic(:,:)
    double precision, allocatable :: cg02_atomic(:,:), cg02_magnetic(:,:,:,:)
    double precision, allocatable :: cg03_atomic(:,:,:,:), cg03_magnetic(:,:,:,:,:)
    double precision, allocatable :: cg04_atomic(:,:,:,:,:), cg04_magnetic(:,:,:,:,:,:)
    double precision, allocatable :: cg05_atomic(:,:,:,:,:,:), cg05_magnetic(:,:,:,:,:,:,:)
        
end module tables
      
!------------------------------------------------------------------------------------      

module modforces

    double precision, allocatable :: wf1(:,:)
    double complex, allocatable :: wf2(:,:,:,:,:)
    double complex, allocatable :: wf(:,:,:,:,:,:,:)
    double complex, allocatable :: wf_0(:,:,:,:)

end module modforces

!-------------------------------------------------------------------------------------

module lookuptables

    integer nlut
    double precision, allocatable :: lutfr(:,:,:,:,:)
    double precision, allocatable :: lutdfr(:,:,:,:,:)
    double precision, allocatable :: lutgr(:,:,:,:)
    double precision, allocatable :: lutdgr(:,:,:,:)
    double precision, allocatable :: luthcr(:,:,:)
    double precision, allocatable :: lutdhcr(:,:,:)

    double precision, allocatable :: lutfrs(:,:,:,:,:,:)        
    double precision, allocatable :: lutgrs(:,:,:,:,:)
    double precision, allocatable :: luthcrs(:,:,:,:)

    double precision, allocatable :: lutfrs_atomic(:,:,:,:,:,:)
    double precision, allocatable :: lutfrs_magnetic(:,:,:,:,:,:)
    double precision, allocatable :: lutgrs_atomic(:,:,:,:,:)
    double precision, allocatable :: lutgrs_magnetic(:,:,:,:,:)
    double precision, allocatable :: luthcrs_atomic(:,:,:,:)
    double precision, allocatable :: luthcrs_magnetic(:,:,:,:)

    double precision invrscalelookup
    double precision rscalelookup
    double precision invrscalelookup_atomic, invrscalelookup_magnetic  
    double precision rscalelookup_atomic, rscalelookup_magnetic

    double precision, allocatable :: scoreRnl(:)
    double precision, allocatable :: scoreRnl_atomic(:)
    double precision, allocatable :: scoreRnl_magnetic(:)

end module lookuptables

!----------------------------------------------------------------------------------------      

module modprepA

    double precision, allocatable :: gradfr2pack_atomic(:,:)
    double precision, allocatable :: fr2pack_atomic(:,:)
    double precision, allocatable :: gradfr2hcpack(:)
    double precision, allocatable :: fr2hcpack(:)      
    double precision, allocatable :: gradfrpack_atomic(:,:,:)
    double precision, allocatable :: gradfrpack_magnetic(:,:,:)
    double precision, allocatable :: gradfrpack_magnetic_0(:,:,:)
    double precision, allocatable :: frpack_atomic(:,:,:)
    double precision, allocatable :: frpack_magnetic(:,:,:)
    double precision, allocatable :: grpack_magnetic_0(:,:)
    double precision, allocatable :: gradgrpack_magnetic_0(:,:)
    double precision, allocatable :: frpack_magnetic_0(:,:,:)
    double complex, allocatable :: gradshpack_atomic(:,:,:)
    double complex, allocatable :: gradshpack_magnetic(:,:,:)
    double complex, allocatable :: gradshpack_magnetic_0(:,:,:)
    double complex, allocatable :: shpack_atomic(:,:)
    double complex, allocatable :: shpack_magnetic(:,:)
    double complex, allocatable :: shpack_magnetic_0(:,:)
    double precision, allocatable :: abase2hcpack(:)
    double complex, allocatable :: abasepack(:,:,:,:,:,:)
    double complex, allocatable :: abase2pack(:,:,:,:,:)
    double complex, allocatable :: abase2pack_magnetic(:,:,:)

end module modprepA

!-----------------------------------------------------------------------------------------
    
module modneigh

    integer nelements
    integer n2atom
    character(2), allocatable :: element(:)
    integer, allocatable :: n2start(:), n2stop(:)
    integer, allocatable :: n2l(:)
    integer, allocatable :: aty(:)
    integer, allocatable :: amom(:,:)
    integer, allocatable :: aty2l(:)
    double precision, allocatable :: d2l(:)
    double precision, allocatable :: d2l_momi(:)
    double precision, allocatable :: d2l_mom(:)
    double precision, allocatable :: e2l(:,:)
    double precision, allocatable :: e2l_momi(:,:)
    double precision, allocatable :: e2l_mom(:,:)
    double precision, allocatable :: cut(:,:,:)
    double precision, allocatable :: dcut(:,:,:)
    double precision, allocatable :: ecut(:)
    double precision, allocatable :: decut(:)      
      
end module modneigh

!--------------------------------------------------------------------------------------------

module modtime

    double precision tcpu, tcpuold
    double precision tunpack, tprepA, tprepB, tforce, ttotal, tstart
    double precision tradbase, tradfunc, tradpack, tgetsh, tpackA, tcopyA

end module modtime

!--------------------------------------------------------------------------------------------    

module modrandAB

    double complex, allocatable :: randABvec(:,:,:,:)

end module modrandAB

!------------------------------------------------------------------------
