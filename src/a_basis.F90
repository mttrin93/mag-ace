!
!  set up atomic base
!
!-----------------------------------------------------------
subroutine prepA

    use modprepA
    use modneigh
    use modtime
    use functionparameters, only: lmax, l2max_magnetic, l3max_atomic, l3max_magnetic, lmaxsq, &
    &  lmaxsq_atomic, lmaxsq_magnetic, lmaxhalf_atomic, lmaxhalf_magnetic, nradbase, nradbase_atomic, &  
    &  nradbase_magnetic, nradial, nradial3_atomic, nradial2_magnetic, nradial3_magnetic, crad, & 
    &  crad_atomic, crad_magnetic, lambda, corepara, radtype_atomic, radtype_magnetic 
    use tables, only: shindex_atomic,shindex_magnetic, shback_atomic, shback_magnetic

    implicit none
    double precision x, mag_mom, mi
    double complex x1, y1
    integer ii
    double precision ebond(1:3)
    double precision mom(1:3),momi(1:3)
    double complex shvec_atomic(1:lmaxsq_atomic), dshvec_atomic(1:3,1:lmaxsq_atomic)
    double complex shvec_magnetic(1:lmaxsq_magnetic), dshvec_magnetic(1:3,1:lmaxsq_magnetic)
    double complex shvec_magnetic_i(1:lmaxsq_magnetic), dshvec_magnetic_i(1:3,1:lmaxsq_magnetic)
    double complex abase(1:nradial3_atomic,1:lmaxsq_atomic, &
    &  1:nradial3_magnetic,1:lmaxsq_magnetic,1:nelements)       
    double complex abase2_magnetic(1:nradial2_magnetic,1:lmaxsq_magnetic)
    double complex abase2(1:nradbase_atomic,1:nradial2_magnetic,1:lmaxsq_magnetic,1:nelements)
    double precision abase2hc
    double precision fr_atomic(1:nradial3_atomic,0:l3max_atomic),    &
    &    dfr_atomic(1:nradial3_atomic,0:l3max_atomic)
    double precision fr_magnetic(1:nradial2_magnetic,0:l2max_magnetic), &
    &    dfr_magnetic(1:nradial2_magnetic,0:l2max_magnetic)
    double precision fr_magnetic_i(1:nradial2_magnetic,0:l2max_magnetic), &
    &    dfr_magnetic_i(1:nradial2_magnetic,0:l2max_magnetic)
    double precision gr_atomic(1:nradbase_atomic), dgr_atomic(1:nradbase_atomic)
    double precision gr_magnetic(1:nradbase_magnetic), dgr_magnetic(1:nradbase_magnetic)
    double precision gr_magnetic_i(1:nradbase_magnetic), dgr_magnetic_i(1:nradbase_magnetic)
    double precision cr, dcr
    
    integer len, n, nat, nei, ij, ijstart, ijstop, ji, k, l, l1, l2, k1, k2, k3, k4
    integer pack1count, pack2count, atomi, nr, nr1, nsh, nsh1, nr2
    
    integer nei2count, elei, elej
    integer i, m, m1, m2

    double precision t, told
    double precision eps, eps2, abs2
    double precision lam(2), cu(2), dc(2), lc, pc, lami(2), cui(2), dci(2)
    
    character(8) atomic, magnetic
    
    atomic = 'atomic'
    magnetic = 'magnetic'
    
    eps = 1.d-9
    eps2 = eps*eps
    
    do atomi = 1, n2atom
        elei = aty(atomi)
        abase = dcmplx(0.d0,0.d0)
        abase2 = dcmplx(0.d0,0.d0)
        abase2_magnetic = dcmplx(0.d0,0.d0)
        abase2hc = 0.d0
        ijstart = n2start(atomi)
        ijstop = n2stop(atomi)

        ! constuction of Aimu_in_0'l_0'm_0'
        mi = d2l_momi(atomi)
        momi(1:3) = e2l_momi(1:3,atomi)

#ifdef USELOOKUP
        call lookupRadspline_magnetic(mi,nradbase_magnetic,nradial2_magnetic,l2max_magnetic,elei,elei, &
                    &   fr_magnetic_i,dfr_magnetic_i,gr_magnetic_i,dgr_magnetic_i)
#else
        lami(1:2) = lambda(1:2,elei,elei)
        cui(1:2) = cut(1:2,elei,elei)
        dci(1:2) = dcut(1:2,elei,elei)   

        call radbase(radtype_magnetic,lami(2),nradbase_magnetic,cui(2),dci(2),mi,gr_magnetic_i,dgr_magnetic_i)
        call radfunc(nradbase_magnetic,nelements,elei,elei,nradial2_magnetic,  &
                    &  l2max_magnetic,crad_magnetic,gr_magnetic_i,dgr_magnetic_i,fr_magnetic_i,dfr_magnetic_i)
#endif       

        ! spherical harmonics for magnetic indices for atom i
        call compute_ylm_cart(magnetic,momi(1),momi(2),momi(3),l2max_magnetic,lmaxhalf_magnetic,shvec_magnetic_i,dshvec_magnetic_i)


        ! construction of abase2_magnetic: Aimu_in_0'l_0'm_0'
        do nr = 1, nradial2_magnetic
            do nsh = 1, lmaxhalf_magnetic               
                abase2_magnetic(nr,nsh) = fr_magnetic_i(nr,shback_magnetic(nsh))*shvec_magnetic_i(nsh)
            enddo
        enddo

        ! update missing half of abase2_magnetic
        do nr = 1, nradial2_magnetic
            do l = 1, l2max_magnetic
                do m = 1, l
                    k = shindex_magnetic(l,-m)
                    k1 = shindex_magnetic(l,m)
                    x = dble((-1)**(m))
                    abase2_magnetic(nr,k) = x*conjg( abase2_magnetic(nr,k1) )
                enddo
            enddo
        enddo
       
        grpack_magnetic_0(1:nradbase_magnetic,atomi)=gr_magnetic_i(1:nradbase_magnetic)
        gradgrpack_magnetic_0(1:nradbase_magnetic,atomi)=dgr_magnetic_i(1:nradbase_magnetic)
        frpack_magnetic_0(1:nradial2_magnetic,0:l2max_magnetic,atomi) = &
                                &  fr_magnetic_i(1:nradial2_magnetic,0:l2max_magnetic) 
        gradfrpack_magnetic_0(1:nradial2_magnetic,0:l2max_magnetic,atomi) = &
                                &  dfr_magnetic_i(1:nradial2_magnetic,0:l2max_magnetic)

        shpack_magnetic_0(1:lmaxhalf_magnetic,atomi) = shvec_magnetic_i(1:lmaxhalf_magnetic)
        gradshpack_magnetic_0(1:3,1:lmaxhalf_magnetic,atomi) = dshvec_magnetic_i(1:3,1:lmaxhalf_magnetic)

        do l = 1, l2max_magnetic
            do m = 1, l
                x = dble((-1)**(m))
                k = shindex_magnetic(l,-m)
                k1 = shindex_magnetic(l,m)
                shpack_magnetic_0(k,atomi) = x*conjg( shvec_magnetic_i(k1) )
                gradshpack_magnetic_0(1:3,k,atomi) = x*conjg( dshvec_magnetic_i(1:3,k1) )
            enddo
        enddo

        do ij = ijstart, ijstop

            x = d2l(ij)
            ebond(1:3) = e2l(1:3,ij)
            elej = aty2l(ij)
            mag_mom = d2l_mom(ij)               ! magnitude magnetic moment
            mom(1:3) = e2l_mom(1:3,ij)    ! unit vector that points in the magnetic moment direction 

            ! looking up radial functions and gradients          
            call CPU_TIME(t)
            told = t

#ifdef USELOOKUP
            if (nradbase_atomic > 0) then
                call lookupRadspline_atomic(x,nradbase_atomic,nradial3_atomic,l3max_atomic,elei,elej,  &
                            &   fr_atomic,dfr_atomic,gr_atomic,dgr_atomic,cr,dcr)
            endif

            call lookupRadspline_magnetic(mag_mom,nradbase_magnetic,nradial2_magnetic,  & 
                            &   l2max_magnetic,elei,elej,fr_magnetic,dfr_magnetic,gr_magnetic,dgr_magnetic)          
#else
            lam(1:2) = lambda(1:2,elei,elej)
            cu(1:2) = cut(1:2,elei,elej)
            dc(1:2) = dcut(1:2,elei,elej)

            ! call radbase and radfunc for the atomic part only
            call radbase(radtype_atomic,lam(1),nradbase_atomic,cu(1),dc(1),x,gr_atomic,dgr_atomic)
            call radfunc(nradbase_atomic,nelements,elei,elej,nradial3_atomic,l3max_atomic,  &
                            &  crad_atomic,gr_atomic,dgr_atomic,fr_atomic,dfr_atomic)
                    
            ! call radbase and radfunc for the magnetic part only
            call radbase(radtype_magnetic,lam(2),nradbase_magnetic,cu(2),dc(2),mag_mom,gr_magnetic,dgr_magnetic)
            call radfunc(nradbase_magnetic,nelements,elei,elej,nradial2_magnetic,  &
                            &  l2max_magnetic,crad_magnetic,gr_magnetic,dgr_magnetic,fr_magnetic,dfr_magnetic)    
                    
            ! hard core repulsion
            pc = corepara(1,elei,elej)
            lc = corepara(2,elei,elej)
            call radcore(x,pc,lc,cu(1),cr,dcr)
#endif
          
            call CPU_TIME(t)
            tradfunc = tradfunc + t - told
            told = t          

            ! set up spherical harmonics
            call CPU_TIME(t)
            told = t 

            ! spherical harmonics for atomic indices
            call compute_ylm_cart(atomic,ebond(1),ebond(2),ebond(3),l3max_atomic,lmaxhalf_atomic,shvec_atomic,dshvec_atomic)

            ! spherical harmonics for magnetic indices
            call compute_ylm_cart(magnetic,mom(1),mom(2),mom(3),l2max_magnetic,lmaxhalf_magnetic,shvec_magnetic,dshvec_magnetic)

            call CPU_TIME(t)
            tgetsh = tgetsh + t - told
            told = t

            ! store radial functions and their gradient
            call CPU_TIME(t)
            told = t
            ! vectors instead of loops

            ! construction of abase2: Aimu_1n_100n_1'l_1'-m0'
            do nr = 1, nradial2_magnetic
                do l = 0, l2max_magnetic
                    do m = 0, l
                        k = shindex_magnetic(l,m)
                        abase2(1:nradbase_atomic,nr,k,elej) = abase2(1:nradbase_atomic,nr,k,elej) + &
                                & gr_atomic(1:nradbase_atomic)*fr_magnetic(nr,shback_magnetic(k))*shvec_magnetic(k)
                    enddo
                enddo
            enddo
          
            abase2hc = abase2hc + cr
            fr2pack_atomic(1:nradbase_atomic,ij) = gr_atomic(1:nradbase_atomic)
            gradfr2pack_atomic(1:nradbase_atomic,ij) = dgr_atomic(1:nradbase_atomic)
            fr2hcpack(ij) = cr
            gradfr2hcpack(ij) = dcr

            frpack_atomic(1:nradial3_atomic,0:l3max_atomic,ij) = fr_atomic(1:nradial3_atomic,0:l3max_atomic)
            frpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,ij) = fr_magnetic(1:nradial2_magnetic,0:l2max_magnetic)

            gradfrpack_atomic(1:nradial3_atomic,0:l3max_atomic,ij) = dfr_atomic(1:nradial3_atomic,0:l3max_atomic)
            gradfrpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,ij) = dfr_magnetic(1:nradial2_magnetic,0:l2max_magnetic)

            ! store only half of the atomic and magnetic spherical harmonics
            shpack_atomic(1:lmaxhalf_atomic,ij) = shvec_atomic(1:lmaxhalf_atomic)
            shpack_magnetic(1:lmaxhalf_magnetic,ij) = shvec_magnetic(1:lmaxhalf_magnetic)
            gradshpack_atomic(1:3,1:lmaxhalf_atomic,ij) = dshvec_atomic(1:3,1:lmaxhalf_atomic)
            gradshpack_magnetic(1:3,1:lmaxhalf_magnetic,ij) = dshvec_magnetic(1:3,1:lmaxhalf_magnetic)


            do l = 1, l2max_magnetic
                do m = 1, l
                    x = dble((-1)**(m))
                    k = shindex_magnetic(l,-m)
                    k1 = shindex_magnetic(l,m)
                    shpack_magnetic(k,ij) = x*conjg( shvec_magnetic(k1) )
                    gradshpack_magnetic(1:3,k,ij) = x*conjg( dshvec_magnetic(1:3,k1) )
                enddo
            enddo

            do l = 1, l3max_atomic
                do m = 1, l
                    x = dble((-1)**(m))
                    k = shindex_atomic(l,-m)
                    k1 = shindex_atomic(l,m)
                    shpack_atomic(k,ij) = x*conjg( shvec_atomic(k1) )
                    gradshpack_atomic(1:3,k,ij) = x*conjg( dshvec_atomic(1:3,k1) )
                enddo
            enddo


            call CPU_TIME(t)
            tradpack = tradpack + t - told
            told = t

            ! compute atomic base Aimunlmn'l'm'
            do nr1 = 1, nradial3_atomic
                do nr2 = 1, nradial3_magnetic
                    do l1 = 0, l3max_atomic
                        do l2 = 0, l3max_magnetic
                            do m2 = -l2, l2

                                m1 = 0
                                k1 = shindex_atomic(l1,m1)
                                k2 = shindex_magnetic(l2,m2)
                                x1 = fr_atomic(nr1,shback_atomic(k1))*fr_magnetic(nr2,shback_magnetic(k2))*shvec_magnetic(k2)
                                y1 = shvec_atomic(k1)
                            
                                abase(nr1,k1,nr2,k2,elej) = abase(nr1,k1,nr2,k2,elej) + x1*y1
                            
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            do nr1 = 1, nradial3_atomic
                do nr2 = 1, nradial3_magnetic
                    do l1 = 1, l3max_atomic
                        do l2 = 0, l3max_magnetic
                            do m1 = 1, l1
                                do m2 = 0, l2
                                                            
                                    k1 = shindex_atomic(l1,m1)
                                    k2 = shindex_magnetic(l2,m2)
                                    k3 = shindex_atomic(l1,-m1)

                                    x = dble((-1)**(m1))
                                    x1 = fr_atomic(nr1,shback_atomic(k1))*fr_magnetic(nr2,shback_magnetic(k2))*shvec_magnetic(k2)
                                    y1 = shvec_atomic(k1)
                                    
                                    abase(nr1,k1,nr2,k2,elej) = abase(nr1,k1,nr2,k2,elej) + x1*y1
                                    abase(nr1,k3,nr2,k2,elej) = abase(nr1,k3,nr2,k2,elej) + x1*conjg(y1)*x
                                
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
          
            call CPU_TIME(t)
            tpackA = tpackA + t - told
            told = t          

        enddo


        ! add prefactors and update missing half of abase and copy to storage
        call CPU_TIME(t)
        told = t
        do elej = 1, nelements  
        
            ! update missing half of abase
            do nr1 = 1, nradial3_atomic
                do nr2 = 1, nradial3_magnetic
                    do l1 = 1, l3max_atomic
                        do l2 = 1, l3max_magnetic
                            do m1 = 1, l1
                                do m2 = -l2, -1
                                
                                    k1 = shindex_atomic(l1,m1)
                                    k2 = shindex_atomic(l1,-m1)
                                    k3 = shindex_magnetic(l2,m2)
                                    k4 = shindex_magnetic(l2,-m2)
                                    
                                    x = dble((-1)**(m1+m2))
                                    
                                    abase(nr1,k1,nr2,k3,elej) = x*conjg(abase(nr1,k2,nr2,k4,elej)) 
                                    abase(nr1,k2,nr2,k3,elej) = x*conjg(abase(nr1,k1,nr2,k4,elej))
                                
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            ! update missing half of abase2
            do nr = 1, nradial2_magnetic
                do l = 1, l2max_magnetic
                    do m = 1, l
                        
                        k = shindex_magnetic(l,-m)
                        k1 = shindex_magnetic(l,m)
                        x = dble((-1)**(m))
                        abase2(1:nradbase_atomic,nr,k,elej) = x*conjg(abase2(1:nradbase_atomic,nr,k1,elej))
                        
                    enddo
                enddo
            enddo
            
        enddo
 
        abasepack(1:nradial3_atomic,1:lmaxsq_atomic,1:nradial3_magnetic,1:lmaxsq_magnetic,1:nelements,atomi) = &
        &   abase(1:nradial3_atomic,1:lmaxsq_atomic,1:nradial3_magnetic,1:lmaxsq_magnetic,1:nelements)

        abase2pack(1:nradbase_atomic,1:nradial2_magnetic,1:lmaxsq_magnetic,1:nelements,atomi) =   &
        &   abase2(1:nradbase_atomic,1:nradial2_magnetic,1:lmaxsq_magnetic,1:nelements) 

        abase2pack_magnetic(1:nradial2_magnetic,1:lmaxsq_magnetic,atomi) =  &
        &   abase2_magnetic(1:nradial2_magnetic,1:lmaxsq_magnetic)

        abase2hcpack(atomi) = abase2hc
        call CPU_TIME(t)
        tcopyA = tcopyA + t - told
        told = t
    enddo
    
end subroutine prepA
