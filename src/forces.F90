! routines for force and energy computation

!-----------------------------------------------------------------------------------------
      
subroutine forceweights
  
    ! compute weights for forces
    use core
    use global
    use tables
    use functionparameters
    use modforces
    use modneigh
    implicit none
    integer atomk, atomi, nele
    double precision energyi
    double precision rho(1:ndensity), dEdrho(1:ndensity)
    double complex dw_0(0:l2max_magnetic,0:l2max_magnetic,1:nradial2_magnetic)
    double precision dw1(1:nradbase_magnetic)
    double complex dw2(1:nradbase_atomic,0:l2max_magnetic,0:l2max_magnetic,1:nradial2_magnetic)
    double complex  dw(-l3max_atomic:l3max_atomic,0:l3max_atomic,1:nradial3_atomic,0:l3max_magnetic, &
                    &   0:l3max_magnetic,1:nradial3_magnetic)
    double precision ceff1(1:nvar1), ceff2(1:nvar2), ceff3(1:nvar3), ceff4(1:nvar4), ceff5(1:nvar5), &
                    &   ceff6(1:nvar6), ceff7(1:nvar7)
    integer pack1count, pack2count
    integer n, i, ij, ijstart, ijstop, k, l, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11
    integer m, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, mat, nat, mati,k1
    integer nclu, nei, nr, nr1, nr2, nr3, nr4, nr5, nr6, nr7, nr8, nr9, nr10, nr11
    double precision x, pm, invr
    double precision t, told, dt
    double complex xiyvec(1:3), wei, xiy, wi, wk, pmwi, xiyvec2(1:3)
    double precision revec(1:3), imvec(1:3), rew, imw, w
    double precision edfi

    ! for debugging
    double precision eu, ed, eps, xv(1:3), xvref(1:3), Fk2(1:3)
    double precision diff, denom
    integer nc

    if (allocated(wf2)) then
        if(size(wf2,5)/=natomsmax) then
            deallocate(wf1)
            deallocate(wf2)
            deallocate(wf)
            deallocate(wf_0)
        endif
    endif

    if (.not.allocated(wf1)) then
        allocate(wf1(1:nradbase_magnetic,1:natomsmax))
    endif    

    if (.not.allocated(wf2)) then
        allocate(wf2(1:nradbase_atomic,0:l2max_magnetic,0:l2max_magnetic, &
                    &  1:nradial2_magnetic,1:natomsmax))
    endif

    if (.not.allocated(wf)) then
        allocate(wf(-l3max_atomic:l3max_atomic,0:l3max_atomic,1:nradial3_atomic,&
                    &  0:l3max_magnetic,0:l3max_magnetic,1:nradial3_magnetic,1:natomsmax))
    endif

    if (.not.allocated(wf_0)) then
        allocate(wf_0(0:l2max_magnetic,0:l2max_magnetic,1:nradial2_magnetic,1:natomsmax))
    endif

    if (allocated(edfcut)) then
        if(n2atom > size(edfcut,1)) then
            deallocate(edfcut)
        endif
    endif
    if (.not.allocated(edfcut)) then
        allocate(edfcut(1:n2atom))
    endif    


    eps = 1.d-5

    if (verblevel >= 10 ) then
        call CPU_TIME(t)
        told = t
    endif

    do atomi = 1, n2atom
        nele = aty(atomi)
        
        ! 1. get effective expansion coefficients          
        call calcenergyatom(atomi,rho,dEdrho,edfi,energyi)
        
        edfcut(atomi) = edfi
        eicalc(atomi) = energyi

        do nclu = 1, nvar1
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c0(nclu,k,nele)
            enddo
            ceff1(nclu) = x
        enddo       
        do nclu = 1, nvar2
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c2(nclu,k,nele)
            enddo
            ceff2(nclu) = x
        enddo
        do nclu = 1, nvar3
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c3(nclu,k,nele)
            enddo
            ceff3(nclu) = x
        enddo
        do nclu = 1, nvar4
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c4(nclu,k,nele)
            enddo
            ceff4(nclu) = x
        enddo
        do nclu = 1, nvar5
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c5(nclu,k,nele)
            enddo
            ceff5(nclu) = x
        enddo
        do nclu = 1, nvar6
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c6(nclu,k,nele)
            enddo
            ceff6(nclu) = x
        enddo
        do nclu = 1, nvar7
            x = 0.d0
            do k = 1, ndensity
                x =  x + dEdrho(k)*c7(nclu,k,nele)
            enddo
            ceff7(nclu) = x
        enddo
        
        ! 2. multiply with dB/dA

        dw = (0.d0,0.d0)
        dw_0 = (0.d0,0.d0)
        dw2 = (0.d0,0.d0)
        dw1 = 0.d0

        ! 1-body
        do nclu = 1, nvar1
            nr = b1index(nclu)
            dw1(nr) = dw1(nr) + ceff1(nclu)
        enddo

        
        ! 2-body
        do nclu = 1, nvar2
            nr1 = b2index(1,nclu)  ! n0'
            nr2 = b2index(2,nclu)  ! n1'
            nr3 = b2index(3,nclu)  ! n1
            l = b2index(4,nclu)    ! l0' (=l1')
            do m = 0, l
                dw2(nr3,m,l,nr2) = dw2(nr3,m,l,nr2) + ceff2(nclu)*db2pack(m,nclu,atomi)
            enddo

            do m = 0, l
                dw_0(m,l,nr1) = dw_0(m,l,nr1) + ceff2(nclu)*db2pack_0(m,nclu,atomi)  
            enddo
        enddo
        
        ! 3-body
        do nclu = 1, nvar3
            nr1 = b3index(1,nclu)   ! n1
            nr2 = b3index(2,nclu)   ! n2
            nr3 = b3index(3,nclu)   ! n0'
            nr4 = b3index(4,nclu)   ! n1'
            nr5 = b3index(5,nclu)   ! n2'
            l1 = b3index(6,nclu)    ! l1  
            l2 = b3index(7,nclu)    ! l0' 
            l4 = b3index(8,nclu)    ! l1'
            l5 = b3index(9,nclu)    ! l2'
            do m1 = -l1, l1
                do m4 = 0, l4
                    dw(m1,l1,nr1,m4,l4,nr4) = dw(m1,l1,nr1,m4,l4,nr4) + ceff3(nclu)*db3pack(1,m1,m4,nclu,atomi)
                enddo
            enddo
            do m1 = -l1, l1
                do m5 = 0, l5
                    dw(m1,l1,nr2,m5,l5,nr5) = dw(m1,l1,nr2,m5,l5,nr5) + ceff3(nclu)*db3pack(2,m1,m5,nclu,atomi)
                enddo
            enddo  
            do m2 = 0, l2
                dw_0(m2,l2,nr3) = dw_0(m2,l2,nr3) + ceff3(nclu)*db3pack_0(m2,nclu,atomi) 
            enddo
        enddo
        
        ! 4-body
        do nclu = 1, nvar4
            nr1 = b4index(1,nclu)  ! n1
            nr2 = b4index(2,nclu)  ! n2
            nr3 = b4index(3,nclu)  ! n3
            nr4 = b4index(4,nclu)  ! n0'
            nr5 = b4index(5,nclu)  ! n1'
            nr6 = b4index(6,nclu)  ! n2'
            nr7 = b4index(7,nclu)  ! n3'
            l1 = b4index(8,nclu)   ! l1
            l2 = b4index(9,nclu)   ! l2
            l3 = b4index(10,nclu)   ! l3
            l4 = b4index(11,nclu)   ! l0'
            l5 = b4index(12,nclu)   ! l1'
            l6 = b4index(13,nclu)   ! l2'
            l7 = b4index(14,nclu)   ! l3'
            do m1 = -l1, l1
                do m5 = 0, l5
                    dw(m1,l1,nr1,m5,l5,nr5) = dw(m1,l1,nr1,m5,l5,nr5) + ceff4(nclu)*db4pack(1,m1,m5,nclu,atomi)
                enddo
            enddo
            do m2 = -l2, l2
                do m6 = 0, l6
                    dw(m2,l2,nr2,m6,l6,nr6) = dw(m2,l2,nr2,m6,l6,nr6) + ceff4(nclu)*db4pack(2,m2,m6,nclu,atomi)
                enddo
            enddo
            do m3 = -l3, l3
                do m7 = 0, l7
                    dw(m3,l3,nr3,m7,l7,nr7) = dw(m3,l3,nr3,m7,l7,nr7) + ceff4(nclu)*db4pack(3,m3,m7,nclu,atomi)
                enddo
            enddo 
            do m4 = 0, l4
                dw_0(m4,l4,nr4) = dw_0(m4,l4,nr4) + ceff4(nclu)*db4pack_0(m4,nclu,atomi) 
            enddo
        enddo
        
        ! 5-body
        do nclu = 1, nvar5
            nr1 = b5index(1,nclu)  ! n1
            nr2 = b5index(2,nclu)  ! n2
            nr3 = b5index(3,nclu)  ! n3
            nr4 = b5index(4,nclu)  ! n4
            nr5 = b5index(5,nclu)  ! n0'
            nr6 = b5index(6,nclu)  ! n1'
            nr7 = b5index(7,nclu)  ! n2'
            nr8 = b5index(8,nclu)  ! n3'
            nr9 = b5index(9,nclu)  ! n4'
            l1 = b5index(10,nclu)  ! l1
            l2 = b5index(11,nclu)  ! l2
            l3 = b5index(12,nclu)  ! l3
            l4 = b5index(13,nclu)  ! l4
            l5 = b5index(14,nclu)  ! l0'
            l6 = b5index(15,nclu)  ! l1'
            l7 = b5index(16,nclu)  ! l2'
            l8 = b5index(17,nclu)  ! l3'
            l9 = b5index(18,nclu)  ! l4'
            do m1 = -l1, l1
                do m6 = 0, l6
                    dw(m1,l1,nr1,m6,l6,nr6) = dw(m1,l1,nr1,m6,l6,nr6) + ceff5(nclu)*db5pack(1,m1,m6,nclu,atomi)
                enddo
            enddo
            do m2 = -l2, l2
                do m7 = 0, l7
                    dw(m2,l2,nr2,m7,l7,nr7) = dw(m2,l2,nr2,m7,l7,nr7) + ceff5(nclu)*db5pack(2,m2,m7,nclu,atomi)
                enddo
            enddo
            do m3 = -l3, l3
                do m8 = 0, l8
                    dw(m3,l3,nr3,m8,l8,nr8) = dw(m3,l3,nr3,m8,l8,nr8) + ceff5(nclu)*db5pack(3,m3,m8,nclu,atomi)
                enddo
            enddo
            do m4 = -l4, l4
                do m9 = 0, l9
                    dw(m4,l4,nr4,m9,l9,nr9) = dw(m4,l4,nr4,m9,l9,nr9) + ceff5(nclu)*db5pack(4,m4,m9,nclu,atomi)
                enddo
            enddo
            do m5 = 0, l5
                dw_0(m5,l5,nr5) = dw_0(m5,l5,nr5) + ceff5(nclu)*db5pack_0(m5,nclu,atomi) 
            enddo
        enddo

        ! 6-body
        do nclu = 1, nvar6
            nr1 = b6index(1,nclu)  ! n1
            nr2 = b6index(2,nclu)  ! n2
            nr3 = b6index(3,nclu)  ! n3
            nr4 = b6index(4,nclu)  ! n4
            nr5 = b6index(5,nclu)  ! n5
            nr6 = b6index(6,nclu)  ! n0'
            nr7 = b6index(7,nclu)  ! n1'
            nr8 = b6index(8,nclu)  ! n2'
            nr9 = b6index(9,nclu)  ! n3'
            nr10 = b6index(10,nclu)  ! n4'
            nr11 = b6index(11,nclu)  ! n5'
            l1 = b6index(12,nclu)  ! l1
            l2 = b6index(13,nclu)  ! l2
            l3 = b6index(14,nclu)  ! l3
            l4 = b6index(15,nclu)  ! l4
            l5 = b6index(16,nclu)  ! l5
            l6 = b6index(17,nclu)  ! l0'
            l7 = b6index(18,nclu)  ! l1'
            l8 = b6index(19,nclu)  ! l2'
            l9 = b6index(20,nclu)  ! l3'
            l10 = b6index(21,nclu)  ! l4'
            l11 = b6index(22,nclu)  ! l5'
            do m1 = -l1, l1
                do m7 = 0, l7
                    dw(m1,l1,nr1,m7,l7,nr7) = dw(m1,l1,nr1,m7,l7,nr7) + ceff6(nclu)*db6pack(1,m1,m7,nclu,atomi)
                enddo
            enddo
            do m2 = -l2, l2
                do m8 = 0, l8
                    dw(m2,l2,nr2,m8,l8,nr8) = dw(m2,l2,nr2,m8,l8,nr8) + ceff6(nclu)*db6pack(2,m2,m8,nclu,atomi)
                enddo
            enddo
            do m3 = -l3, l3
                do m9 = 0, l9
                    dw(m3,l3,nr3,m9,l9,nr9) = dw(m3,l3,nr3,m9,l9,nr9) + ceff6(nclu)*db6pack(3,m3,m9,nclu,atomi)
                enddo
            enddo
            do m4 = -l4, l4
                do m10 = 0, l10
                    dw(m4,l4,nr4,m10,l10,nr10) = dw(m4,l4,nr4,m10,l10,nr10) + ceff6(nclu)*db6pack(4,m4,m10,nclu,atomi)
                enddo
            enddo
            do m5 = -l5, l5
                do m11 = 0, l11
                    dw(m5,l5,nr5,m11,l11,nr11) = dw(m5,l5,nr5,m11,l11,nr11) + ceff6(nclu)*db6pack(5,m5,m11,nclu,atomi)
                enddo
            enddo
            do m6 = 0, l6
                dw_0(m6,l6,nr6) = dw_0(m6,l6,nr6) + ceff6(nclu)*db6pack_0(m6,nclu,atomi) 
            enddo
        enddo


        ! 3. store weights
        wf1(1:nradbase_magnetic,atomi) = dw1(1:nradbase_magnetic)
        
        wf2(1:nradbase_atomic,0:l2max_magnetic,0:l2max_magnetic, &
            &  1:nradial2_magnetic,atomi) = dw2(1:nradbase_atomic,0:l2max_magnetic, &
            &  0:l2max_magnetic,1:nradial2_magnetic)
                
        ! we use the full range of the atomic m index, but only half for m'
        wf(-l3max_atomic:l3max_atomic,0:l3max_atomic,1:nradial3_atomic,   &
            &  0:l3max_magnetic,0:l3max_magnetic,1:nradial3_magnetic,atomi) =    &
            &  dw(-l3max_atomic:l3max_atomic,0:l3max_atomic,1:nradial3_atomic,   &
            &  0:l3max_magnetic,0:l3max_magnetic,1:nradial3_magnetic)

        wf_0(0:l2max_magnetic,0:l2max_magnetic,1:nradial2_magnetic,atomi) = &
            &    dw_0(0:l2max_magnetic,0:l2max_magnetic,1:nradial2_magnetic)

    enddo



    if (verblevel >= 100 ) then
        call CPU_TIME(t)
        dt = (t - told)/dble(n2atom)

        told = t
        write(*,'(A,E7.2)') '   force: weight:    CPUtime/atom= ', dt
    endif

end subroutine forceweights

!-----------------------------------------------------------------------------------------

subroutine forcecompute

    use core
    use global
    use tables
    use functionparameters
    use modforces
    use modneigh
    use modprepA
    use modtime
    implicit none
    integer atomk, atomi, nele

    double precision ebond(1:3), emom(1:3), Fk(1:3), Fi(1:3)
    double complex Fkc(1:3)
    double precision dfr2_atomic(1:nradbase_atomic), dfr2hc
    double precision dfr_atomic(1:nradial3_atomic,0:l3max_atomic)
    double precision fr_atomic(1:nradial3_atomic,0:l3max_atomic), fr_magnetic(1:nradial2_magnetic,0:l2max_magnetic)
    double complex dshvec_atomic(1:3,1:lmaxsq_atomic)
    double complex shvec_atomic(1:lmaxsq_atomic), shvec_magnetic(1:lmaxsq_magnetic)
    integer pack1count, pack2count
    integer n, i, ij, ijstart, ijstop, k, l, l1, l2, l3, m, m1, m2, m3, mat, nat, mati,k1
    integer nclu, nei, nr, nr1, nr2, nr3, index_debugl, index_debugn
    double precision x, mom, pm, invr, invr_m
    double precision t, told, dt
    double complex xiyvec(1:3), xiyvec1(1:3), xiyvec2(1:3), wei, xiy, wi, wk, pmwi
    double precision revec(1:3), imvec(1:3), rew, imw, w
    double precision revec1(1:3), imvec1(1:3), revec2(1:3), imvec2(1:3)
    double complex cx1, cvec2(1:3)
    double precision x2, v2(3)
    double precision edf

    ! for debugging
    double precision eu, ed, eps, xv(1:3), xvref(1:3), Fk2(1:3)
    double precision diff, denom
    integer nc
    double complex im

    im = (0.d0,1.d0)

    if (.not.allocated(wf2)) then
        print *,'Something wrong here. wf2 and wf need to be allocated.'
        stop
    endif

    if (verblevel >= 10 ) then
        call CPU_TIME(t)
        told = t
    endif

    ! loop for force computation

    do atomk = 1, n2atom

        edf = edfcut(atomk)
                
        do ij = n2start(atomk), n2stop(atomk)
            Fk(1:3) = 0.d0
            atomi = n2l(ij)      
            ebond(1:3) = e2l(1:3,ij)
            x = d2l(ij)
            invr = 1.d0/x
            
            ! get gradients of basis functions from memory
            
            ! grad g_{k}(r_{ki})
            dfr2_atomic(1:nradbase_atomic) = gradfr2pack_atomic(1:nradbase_atomic,ij)
            dfr2hc =  gradfr2hcpack(ij)

            ! grad R_{nl}(r_{ki})
            dfr_atomic(1:nradial3_atomic,0:l3max_atomic) = &
            &   gradfrpack_atomic(1:nradial3_atomic,0:l3max_atomic,ij)

            ! R_{nl}(r_{ki})
            fr_atomic(1:nradial3_atomic,0:l3max_atomic) =  &
            &   frpack_atomic(1:nradial3_atomic,0:l3max_atomic,ij)

            ! M_{n'l'}(m_{k})
            fr_magnetic(1:nradial2_magnetic,0:l2max_magnetic) =  &
            &   frpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,ij)
            

            ! grad Y_{lm}(\widehat_{r}_{ki})
            dshvec_atomic(1:3,1:lmaxsq_atomic) = gradshpack_atomic(1:3,1:lmaxsq_atomic,ij)
            
            ! Y_{lm}(\widehat_{r}_{ki})
            shvec_atomic(1:lmaxsq_atomic) = shpack_atomic(1:lmaxsq_atomic,ij)
            
            ! Y_{l'm'}(\widehat_{m}_{k})
            shvec_magnetic(1:lmaxsq_magnetic) = shpack_magnetic(1:lmaxsq_magnetic,ij)

            ! loop basis functions

            ! first pair part only
            do nr1 = 1, nradbase_atomic
                do nr2 = 1, nradial2_magnetic
                    do l = 0, l2max_magnetic

                        v2(1:3) = dfr2_atomic(nr1)*ebond(1:3)
                        m = 0
                        wei = wf2(nr1,m,l,nr2,atomk)
                        rew = dreal(wei)
                        imw = dimag(wei)
                        k = shindex_magnetic(l,m)
                        xiyvec(1:3) = fr_magnetic(nr2,l)*shvec_magnetic(k)*v2(1:3)
                        revec(1:3) = dreal(xiyvec(1:3))
                        imvec(1:3) = dimag(xiyvec(1:3))
                        Fk(1:3) = Fk(1:3) + rew*revec(1:3) - imw*imvec(1:3)


                        do m = 1, l
                        
                            wei = wf2(nr1,m,l,nr2,atomk)
                            rew = dreal(wei)
                            imw = dimag(wei)
                            k = shindex_magnetic(l,m)
                            xiyvec(1:3) = fr_magnetic(nr2,l)*shvec_magnetic(k)*v2(1:3)
                            revec(1:3) = dreal(xiyvec(1:3))
                            imvec(1:3) = dimag(xiyvec(1:3))
                            Fk(1:3) = Fk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))
                        
                        enddo

                    enddo
                enddo
            enddo
                    
            ! hard core part
            Fk(1:3) = Fk(1:3) + edf*dfr2hc*ebond(1:3)

            ! now multi-body interactions
            do nr1 = 1, nradial3_atomic
                do nr2 = 1, nradial3_magnetic
                    do l1 = 0, l3max_atomic
                        do l2 = 0, l3max_magnetic
                    
                            x2 = 0.d0
                            wei = 0.d0
                            v2(1:3) = 0.d0
                            xiyvec(1:3) = 0.d0
                            revec(1:3) = 0.d0
                            imvec(1:3) = 0.d0
                            
                            ! m = 0, m' = 0
                            m1 = 0
                            m2 = 0
                            x2 = invr*fr_atomic(nr1,l1) ! R
                            v2(1:3) = dfr_atomic(nr1,l1)*ebond(1:3) ! DR * r_hat
                            wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                            rew = dreal(wei)
                            imw = dimag(wei)
                            k = shindex_atomic(l1,m1)
                            k1 = shindex_magnetic(l2,m2)

                            xiyvec(1:3) = shvec_magnetic(k1)*fr_magnetic(nr2,l2)* &
                                    &  ( v2(1:3)*shvec_atomic(k) + x2*dshvec_atomic(1:3,k) )
                            
                            revec(1:3) = dreal(xiyvec(1:3))
                            imvec(1:3) = dimag(xiyvec(1:3))
                                    
                            Fk(1:3) = Fk(1:3) + rew*revec(1:3) - imw*imvec(1:3)



                            do m1 = 1, l1
                                do m2 = 1, l2
                                
                                    ! m1 > 0, m' > 0
                                    x2 = invr*fr_atomic(nr1,l1) ! R
                                    v2(1:3) = dfr_atomic(nr1,l1)*ebond(1:3) ! DR * r_hat
                                    wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                                    rew = dreal(wei)
                                    imw = dimag(wei)
                                    k = shindex_atomic(l1,m1)
                                    k1 = shindex_magnetic(l2,m2)
                                        
                                    xiyvec(1:3) = shvec_magnetic(k1)*fr_magnetic(nr2,l2)* &
                                                &  ( v2(1:3)*shvec_atomic(k) + x2*dshvec_atomic(1:3,k) )

                                    revec(1:3) = dreal(xiyvec(1:3))
                                    imvec(1:3) = dimag(xiyvec(1:3))
                                                
                                    Fk(1:3) = Fk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3)) 
                                    

                                    
                                    ! m < 0, m' > 0
                                    x2 = invr*fr_atomic(nr1,l1) ! R
                                    v2(1:3) = dfr_atomic(nr1,l1)*ebond(1:3) ! DR * r_hat
                                    wei = wf(-m1,l1,nr1,m2,l2,nr2,atomk)
                                    rew = dreal(wei)
                                    imw = dimag(wei)
                                    k = shindex_atomic(l1,-m1)
                                    k1 = shindex_magnetic(l2,m2)
                                    
                                    xiyvec(1:3) = shvec_magnetic(k1)*fr_magnetic(nr2,l2)* &
                                                &  ( v2(1:3)*shvec_atomic(k) + x2*dshvec_atomic(1:3,k) )
                                    revec(1:3) = dreal(xiyvec(1:3))
                                    imvec(1:3) = dimag(xiyvec(1:3))
                                                
                                    Fk(1:3) = Fk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))
                                
                                
                                enddo
                            enddo
                        
                            do m1 = 1, l1
                            
                                ! m1 > 0 , m' = 0   
                                m2 = 0
                                x2 = invr*fr_atomic(nr1,l1) ! R
                                v2(1:3) = dfr_atomic(nr1,l1)*ebond(1:3) ! DR * r_hat
                                wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                                rew = dreal(wei)
                                imw = dimag(wei)
                                k = shindex_atomic(l1,m1)
                                k1 = shindex_magnetic(l2,m2)
                                
                                xiyvec(1:3) = shvec_magnetic(k1)*fr_magnetic(nr2,l2)* &
                                        &  ( v2(1:3)*shvec_atomic(k) + x2*dshvec_atomic(1:3,k) )
                                revec(1:3) = dreal(xiyvec(1:3))
                                imvec(1:3) = dimag(xiyvec(1:3))
                                        
                                Fk(1:3) = Fk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))

                            
                            enddo
                                            
                            do m2 = 1, l2
                            
                                ! m1 = 0, m' > 0  
                                m1 = 0
                                x2 = invr*fr_atomic(nr1,l1) ! R
                                v2(1:3) = dfr_atomic(nr1,l1)*ebond(1:3) ! DR * r_hat
                                wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                                rew = dreal(wei)
                                imw = dimag(wei)
                                k = shindex_atomic(l1,m1)
                                k1 = shindex_magnetic(l2,m2)
                                
                                xiyvec(1:3) = shvec_magnetic(k1)*fr_magnetic(nr2,l2)* &
                                                &  ( v2(1:3)*shvec_atomic(k) + x2*dshvec_atomic(1:3,k) )
                                revec(1:3) = dreal(xiyvec(1:3))
                                imvec(1:3) = dimag(xiyvec(1:3))
                                                
                                Fk(1:3) = Fk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))
                                
                            
                            enddo
                        
                ! m > 0 part, -m and +m are handled together
                        enddo
                    enddo
                enddo
            enddo
            ! store forces for each bond          
            fijcalc(1:3,ij) = Fk(1:3)
        enddo ! do ij
    enddo  ! do n2atom
        

    if (verblevel >= 100 ) then
        call CPU_TIME(t)
        dt = (t - told)/dble(n2atom)
        told = t
        write(*,'(A,E7.2)') '   force: grad mult: CPUtime/atom= ', dt
    endif
    
end subroutine forcecompute
  
!-----------------------------------------------------------------------------------------

subroutine torquecompute
    use core
    use global
    use tables
    use functionparameters
    use modforces
    use modneigh
    use modprepA
    use modtime
    implicit none
    integer atomk, atomi, nele

    double precision ebond(1:3), emom(1:3), Tk(1:3), Tk0(1:3)
    double precision dfr2hc, fr2_atomic(1:nradbase_atomic)
    double precision dfr_magnetic(1:nradial2_magnetic,0:l2max_magnetic)
    double precision fr_atomic(1:nradial3_atomic,0:l3max_atomic), fr_magnetic(1:nradial2_magnetic,0:l2max_magnetic)
    double precision fr_magnetic_k(1:nradial2_magnetic,0:l2max_magnetic), dfr_magnetic_k(1:nradial2_magnetic,0:l2max_magnetic)
    double precision dgr_magnetic_k(1:nradbase_magnetic)
    double complex dshvec_magnetic(1:3,1:lmaxsq_magnetic)
    double complex shvec_atomic(1:lmaxsq_atomic), shvec_magnetic(1:lmaxsq_magnetic)
    double complex shvec_magnetic_k(1:lmaxsq_magnetic), dshvec_magnetic_k(1:3,1:lmaxsq_magnetic)
    integer pack1count, pack2count
    integer n, i, ij, ijstart, ijstop, k, l, l1, l2, l3, m, m1, m2, m3, mat, nat, mati,k1
    integer nclu, nei, nr, nr1, nr2, nr3, index_debugl, index_debugn
    double precision x, mom, pm, invr, invr_m
    double precision t, told, dt
    double complex xiyvec(1:3), wei, xiy, wi, wk, pmwi, xiyvec2(1:3)
    double precision revec(1:3), imvec(1:3), rew, imw, w
    double complex cx1, cvec2(1:3)
    double precision x2, v2(3)
    double precision edf


    ! for debugging
    double precision eu, ed, eps, xv(1:3), xvref(1:3), Fk2(1:3)
    double precision diff, denom
    integer nc

    if (.not.allocated(wf_0)) then
        print *,'Something wrong here. wf2 and wf need to be allocated.'
        stop
    endif

    if (verblevel >= 10 ) then
        call CPU_TIME(t)
        told = t
    endif

    ! loop for force computation

    do atomk = 1, n2atom
        
        mom = d2l_momi(atomk)
        emom(1:3) = e2l_momi(1:3,atomk)
        invr_m = 1.d0/mom
        
                    
        ! M_{n'l'}_{m_{k}}
        fr_magnetic_k(1:nradial2_magnetic,0:l2max_magnetic) = &
        &  frpack_magnetic_0(1:nradial2_magnetic,0:l2max_magnetic,atomk)
        
        ! grad M_{n'l'}_{m_{k}}
        dfr_magnetic_k(1:nradial2_magnetic,0:l2max_magnetic) = &
        &  gradfrpack_magnetic_0(1:nradial2_magnetic,0:l2max_magnetic,atomk)
        
        ! grad g_{n'}_{m_{k}}
        dgr_magnetic_k(1:nradbase_magnetic) = gradgrpack_magnetic_0(1:nradbase_magnetic,atomk)
        
        ! Y_{l'm'}(\widehat_{m}_{k})
        shvec_magnetic_k(1:lmaxhalf_magnetic) = shpack_magnetic_0(1:lmaxhalf_magnetic,atomk)
        
        ! grad Y_{l'm'}(\widehat_{m}_{k})
        dshvec_magnetic_k(1:3,1:lmaxhalf_magnetic) = gradshpack_magnetic_0(1:3,1:lmaxhalf_magnetic,atomk)
                
        Tk0(1:3) = 0.d0
        
        ! first order contribution
        do nr = 1, nradbase_magnetic
            Tk0(1:3) = Tk0(1:3) + wf1(nr,atomk)*dgr_magnetic_k(nr)*emom(1:3)
        enddo
        
        do nr = 1, nradial2_magnetic
            do l = 0, l2max_magnetic
            
                v2(1:3) = dfr_magnetic_k(nr,l)*emom(1:3)
                x2 = invr_m*fr_magnetic_k(nr,l) ! R
                m = 0
                wei = wf_0(m,l,nr,atomk)
                rew = dreal(wei)
                imw = dimag(wei)
                k = shindex_magnetic(l,m)
                xiyvec(1:3) = shvec_magnetic_k(k)*v2(1:3) + dshvec_magnetic_k(1:3,k)*x2
                revec(1:3) = dreal(xiyvec(1:3))
                imvec(1:3) = dimag(xiyvec(1:3))
                Tk0(1:3) = Tk0(1:3) + rew*revec(1:3) - imw*imvec(1:3)


                do m = 1, l
                    
                    wei = wf_0(m,l,nr,atomk)
                    rew = dreal(wei)
                    imw = dimag(wei)
                    k = shindex_magnetic(l,m)
                    xiyvec(1:3) = shvec_magnetic_k(k)*v2(1:3) + dshvec_magnetic_k(1:3,k)*x2
                    revec(1:3) = dreal(xiyvec(1:3))
                    imvec(1:3) = dimag(xiyvec(1:3))
                    Tk0(1:3) = Tk0(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3)) 
                
                enddo
            
            enddo
        enddo
        
        do ij = n2start(atomk), n2stop(atomk)
            Tk(1:3) = 0.d0
            
            atomi = n2l(ij)      
            ebond(1:3) = e2l(1:3,ij)
            x = d2l(ij)
            invr = 1.d0/x
            
            emom(1:3) = e2l_mom(1:3,ij)
            mom = d2l_mom(ij)
            invr_m = 1.0/mom
            
            
            ! get gradients of basis functions from memory
                    
            ! g_{k}(r_{ki})
            fr2_atomic(1:nradbase_atomic) = fr2pack_atomic(1:nradbase_atomic,ij)

            ! R_{nl}(r_{ki})
            fr_atomic(1:nradial3_atomic,0:l3max_atomic) =  &
            &   frpack_atomic(1:nradial3_atomic,0:l3max_atomic,ij)

            ! M_{n'l'}(m_{k})
            fr_magnetic(1:nradial2_magnetic,0:l2max_magnetic) =  &
            &   frpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,ij)
            
            ! grad M_{n'l'}(m_{k})
            dfr_magnetic(1:nradial2_magnetic,0:l2max_magnetic) =  &
            &   gradfrpack_magnetic(1:nradial2_magnetic,0:l2max_magnetic,ij)
            
            ! Y_{lm}(\widehat_{r}_{ki})
            shvec_atomic(1:lmaxsq_atomic) = shpack_atomic(1:lmaxsq_atomic,ij)
            
            ! Y_{l'm'}(\widehat_{m}_{k})
            shvec_magnetic(1:lmaxhalf_magnetic) = shpack_magnetic(1:lmaxhalf_magnetic,ij)
            
            ! grad Y_{l'm'}(\widehat_{m}_{k})
            dshvec_magnetic(1:3,1:lmaxhalf_magnetic) = gradshpack_magnetic(1:3,1:lmaxhalf_magnetic,ij)

            ! loop basis functions

            ! first pair part only
            do nr1 = 1, nradbase_atomic
                do nr2 = 1, nradial2_magnetic
                    do l = 0, l2max_magnetic
                        
                        ! m = 0
                        m = 0
                        x2 = invr_m*fr_magnetic(nr2,l) ! M/m
                        v2(1:3) = dfr_magnetic(nr2,l)*emom(1:3) ! DM * m_hat
                        wei = wf2(nr1,m,l,nr2,atomk)
                        rew = dreal(wei)
                        imw = dimag(wei)
                        k1 = shindex_magnetic(l,m)
                        
                        xiyvec(1:3) = fr2_atomic(nr1)* &
                            &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )

                        revec(1:3) = dreal(xiyvec(1:3))
                        imvec(1:3) = dimag(xiyvec(1:3))
                                
                        Tk(1:3) = Tk(1:3) + rew*revec(1:3) - imw*imvec(1:3)

                        do m = 1, l
                                
                            ! m > 0    
                            wei = wf2(nr1,m,l,nr2,atomk)
                            rew = dreal(wei)
                            imw = dimag(wei)
                            k1 = shindex_magnetic(l,m)
                            
                            xiyvec(1:3) = fr2_atomic(nr1)* &
                                &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )
                            revec(1:3) = dreal(xiyvec(1:3))
                            imvec(1:3) = dimag(xiyvec(1:3))
                                    
                            Tk(1:3) = Tk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))
                            
                        enddo
                        
                    enddo
                enddo
            enddo

            ! now multi-body interactions
            do nr1 = 1, nradial3_atomic
                do nr2 = 1, nradial3_magnetic
                    do l1 = 0, l3max_atomic
                        do l2 = 0, l3max_magnetic
                    
                    
                            x2 = 0.d0
                            wei = 0.d0
                            v2(1:3) = 0.d0
                            xiyvec(1:3) = 0.d0
                            revec(1:3) = 0.d0
                            imvec(1:3) = 0.d0                       
                        
                            ! m = 0, m' = 0
                            m1 = 0
                            m2 = 0
                            x2 = invr_m*fr_magnetic(nr2,l2) ! M/m
                            v2(1:3) = dfr_magnetic(nr2,l2)*emom(1:3) ! DM * m_hat
                            wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                            rew = dreal(wei)
                            imw = dimag(wei)
                            k = shindex_atomic(l1,m1)
                            k1 = shindex_magnetic(l2,m2)
                            
                            xiyvec(1:3) = shvec_atomic(k)*fr_atomic(nr1,l1)* &
                                &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )
                            revec(1:3) = dreal(xiyvec(1:3))
                            imvec(1:3) = dimag(xiyvec(1:3))
                                    
                            Tk(1:3) = Tk(1:3) + rew*revec(1:3) - imw*imvec(1:3)


                            do m1 = 1, l1
                                do m2 = 1, l2
                                
                                    ! m > 0 , m' > 0  
                                    x2 = invr_m*fr_magnetic(nr2,l2) ! M/m
                                    v2(1:3) = dfr_magnetic(nr2,l2)*emom(1:3) ! DM * m_hat
                                    wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                                    rew = dreal(wei)
                                    imw = dimag(wei)
                                    k = shindex_atomic(l1,m1)
                                    k1 = shindex_magnetic(l2,m2)
                                    
                                    xiyvec(1:3) = shvec_atomic(k)*fr_atomic(nr1,l1)* &
                                    &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )
                                    revec(1:3) = dreal(xiyvec(1:3))
                                    imvec(1:3) = dimag(xiyvec(1:3))
                                            
                                    Tk(1:3) = Tk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))


                            
                                    ! m < 0, m' > 0
                                    x2 = invr_m*fr_magnetic(nr2,l2) ! M/m
                                    v2(1:3) = dfr_magnetic(nr2,l2)*emom(1:3) ! DM * m_hat
                                    wei = wf(-m1,l1,nr1,m2,l2,nr2,atomk)
                                    rew = dreal(wei)
                                    imw = dimag(wei)
                                    k = shindex_atomic(l1,-m1)
                                    k1 = shindex_magnetic(l2,m2)
                                    
                                    xiyvec(1:3) = shvec_atomic(k)*fr_atomic(nr1,l1)* &
                                    &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )
                                    revec(1:3) = dreal(xiyvec(1:3))
                                    imvec(1:3) = dimag(xiyvec(1:3))
                                            
                                    Tk(1:3) = Tk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3)) 


                                enddo
                            enddo
                        
                            do m1 = 1, l1
                            
                                ! m1 > 0, m' = 0
                                m2 = 0
                                x2 = invr_m*fr_magnetic(nr2,l2) ! M/m
                                v2(1:3) = dfr_magnetic(nr2,l2)*emom(1:3) ! DM * m_hat
                                wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                                rew = dreal(wei)
                                imw = dimag(wei)
                                k = shindex_atomic(l1,m1)
                                k1 = shindex_magnetic(l2,m2)
                                    
                                xiyvec(1:3) = shvec_atomic(k)*fr_atomic(nr1,l1)* &
                                &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )
                                revec(1:3) = dreal(xiyvec(1:3))
                                imvec(1:3) = dimag(xiyvec(1:3))
                                            
                                Tk(1:3) = Tk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))


                            enddo
                        
                            do m2 = 1, l2
                        
                                ! m = 0, m' > 0
                                m1 = 0
                                x2 = invr_m*fr_magnetic(nr2,l2) ! M/m
                                v2(1:3) = dfr_magnetic(nr2,l2)*emom(1:3) ! DM * m_hat
                                wei = wf(m1,l1,nr1,m2,l2,nr2,atomk)
                                rew = dreal(wei)
                                imw = dimag(wei)
                                k = shindex_atomic(l1,m1)
                                k1 = shindex_magnetic(l2,m2)
                                
                                xiyvec(1:3) = shvec_atomic(k)*fr_atomic(nr1,l1)* &
                                &  ( v2(1:3)*shvec_magnetic(k1) + x2*dshvec_magnetic(1:3,k1) )
                                revec(1:3) = dreal(xiyvec(1:3))
                                imvec(1:3) = dimag(xiyvec(1:3))
                                        
                                Tk(1:3) = Tk(1:3) + 2.d0*(rew*revec(1:3) - imw*imvec(1:3))

                            
                            enddo
                                        
                ! m > 0 part, -m and +m are handled together
                        enddo
                    enddo
                enddo
            enddo
            ! store torques for each bond          
            tijcalc(1:3,ij) = Tk(1:3)
        enddo ! do ij
        tijcalc0(1:3,atomk) = Tk0(1:3)
    enddo  ! do n2atom
        

    if (verblevel >= 100 ) then
        call CPU_TIME(t)
        dt = (t - told)/dble(n2atom)
        told = t
        write(*,'(A,E7.2)') '   torque: grad mult: CPUtime/atom= ', dt
    endif
    
end subroutine torquecompute
  
!-----------------------------------------------------------------------------------------

subroutine expansionvariables(atomi,rho)

    use global
    use functionparameters
    use tables
    use modneigh
    implicit none
    integer atomi, ns
    double precision rho(1:ndensity)
    double precision b1(1:nvar1), b2(1:nvar2), b3(1:nvar3), b4(1:nvar4), b5(1:nvar5), b6(1:nvar6), b7(1:nvar7)
    integer nd, np, nclu, nele
    integer atomj, atomiinlist

    integer k, n, bb
    double precision x
    double precision pX0(1:ndensity)
    double precision pX2(1:ndensity)
    double precision pX3(1:ndensity)
    double precision pX4(1:ndensity)
    double precision pX5(1:ndensity)
    double precision pX6(1:ndensity)
    double precision pX7(1:ndensity)

    nele = aty(atomi)

    b1(1:nvar1) = b1pack(1:nvar1,atomi)
    do nd = 1, ndensity
        pX0(nd)  = sum( c0(1:nvar1,nd,nele)*b1(1:nvar1) )
    enddo

    rho(1:ndensity) = pX0(1:ndensity)

    if (nradial2_magnetic > 0) then 
        b2(1:nvar2) = b2pack(1:nvar2,atomi)
        do nd = 1, ndensity
            pX2(nd) = sum( c2(1:nvar2,nd,nele)*b2(1:nvar2) )
        enddo
        rho(1:ndensity) = rho(1:ndensity) + pX2(1:ndensity)
    endif

    if ((nradial3_atomic > 0).or.(nradial3_magnetic > 0)) then
        b3(1:nvar3) = b3pack(1:nvar3,atomi)
        do nd = 1, ndensity
            pX3(nd) = sum( c3(1:nvar3,nd,nele)*b3(1:nvar3) )
        enddo
        rho(1:ndensity) = rho(1:ndensity) +  pX3(1:ndensity)
    endif

    if ((nradial4_atomic > 0).or.(nradial4_magnetic > 0)) then
        b4(1:nvar4) = b4pack(1:nvar4,atomi)
        do nd = 1, ndensity
            pX4(nd) = sum( c4(1:nvar4,nd,nele)*b4(1:nvar4) )
        enddo
        rho(1:ndensity) = rho(1:ndensity) +  pX4(1:ndensity)
    endif

    if ((nradial5_atomic > 0).or.(nradial5_magnetic > 0)) then
        b5(1:nvar5) = b5pack(1:nvar5,atomi)
        do nd = 1, ndensity
            pX5(nd) = sum( c5(1:nvar5,nd,nele)*b5(1:nvar5) )
        enddo
        rho(1:ndensity) = rho(1:ndensity) +  pX5(1:ndensity)
    endif

    if ((nradial6_atomic > 0).or.(nradial6_magnetic > 0)) then
        b6(1:nvar6) = b6pack(1:nvar6,atomi)
        do nd = 1, ndensity
            pX6(nd) = sum( c6(1:nvar6,nd,nele)*b6(1:nvar6) )
        enddo
        rho(1:ndensity) = rho(1:ndensity) +  pX6(1:ndensity)
    endif

end subroutine expansionvariables

!-----------------------------------------------------------------------------------------

subroutine embedding(energyi,rho,won,xpon,dEdrho,ndensity,nfemb,doforce)

    implicit none
    integer ndensity
    double precision energyi
    double precision rho(1:ndensity), won(1:ndensity), xpon(1:ndensity), dEdrho(1:ndensity)
    integer np, nd
    double precision F1, dF1, F2, dF2, F3, dF3, F4, dF4, F5, dF5, F6, dF6
    double precision F, dF, wp,xp, rho1, rho2, rho3, rho4, rho5, rho6
    integer nfemb
    logical doforce

    select case(ndensity)
    
        case(1)
        wp = won(1)
        xp = xpon(1)
        rho1 = rho(1)
        call Fexp(rho1, wp, xp, F, dF,nfemb,doforce)
        energyi = F
        dEdrho(1) = dF
    
        case(2)
        energyi = 0.d0
        do nd = 1, ndensity
            wp = won(nd)
            xp = xpon(nd)
            call Fexp(rho(nd), wp, xp, F, dF,nfemb,doforce)
            energyi = energyi + F
            dEdrho(nd) = dF
        enddo
        
        case default
            print *,'Energy function not implemented.'
            stop
    end select
        
end subroutine embedding
     
!-----------------------------------------------------------------------------------------

subroutine calcenergyatom(atomi,rho,dEdrho,edfi,energyi)  

    use global
    use functionparameters
    use modneigh
    implicit none
    integer atomi, ns, nele
    double precision rho(1:ndensity), dEdrho(1:ndensity), xpon(1:ndensity),won(1:ndensity)
    double precision energyi, edfi
    double precision e, ecutin, decutin, fcut, dfcut

    call expansionvariables(atomi,rho)

    nele = aty(atomi)
    won(1:ndensity) = p1(1,1:ndensity,nele)
    xpon(1:ndensity) = p1(2,1:ndensity,nele)

    call embedding(energyi,rho,won,xpon,dEdrho,ndensity,nfemb,doforce)

    if (nradial2_magnetic> 0) then
        ! energy of core repulsion
        e = b2hcpack(atomi)
        ! inner cutoff
        ecutin = ecut(nele)
        decutin = decut(nele)

        call innercutoff(e,ecutin,decutin,fcut,dfcut)

        ! store gradient for cut
        edfi = dfcut*energyi + 1.d0
        ! update gradients accordingly
        dEdrho(1:ndensity) = fcut*dEdrho(1:ndensity)
        ! modify energy
        energyi = fcut*energyi + e
        
    endif

end subroutine calcenergyatom

!-----------------------------------------------------------------------------------------

subroutine Fexp(rho, wpre, mexp, F, dF,nfemb,doforce)
  
    implicit none

    double precision rho, wpre, mexp, F, dF
    integer nfemb
    logical doforce
    double precision g, dg, omg, a, da, dy, y1, dy1, y2, dy2, srho
    double precision, parameter :: w = 10.d0 
    double precision, parameter :: eps = 1.d-10


    select case (nfemb)

    case(0)
        F = rho
        if(doforce) then
            dF = 1.d0
        end if
    case (1)
        !     
        ! F = wpre*sign(x)*(  ( ( 1 - exp(-w*x**2) )*abs(x) )^m +  m*exp(-w*x**2)*abs(x) )
        !
        
        if ( abs(rho) > eps ) then
        
            y1 = w*rho**2
            ! to avoid underflow
            if (y1 > 30.d0) then
                g = 0.d0
            else
                g  = exp(-y1) !exp(-w*x**2)
            endif
            omg = 1.d0 - g      !( 1 - exp(-w*x**2) )
            a  = abs(rho) !*abs(x)
            y1  =  ( omg*a )**mexp ! ( ( 1 - exp(-w*x**2) )*abs(x) )^m
            y2  = mexp*g*a
            F = wpre*sign(1.d0,rho)*( y1 + y2 )
        
            if (doforce) then
                dg = -2.d0*w*rho*g
                da = sign(1.d0,rho)
                if (dabs(y1) < eps) then
                    dy = 0.d0
                else
                    dy  =  mexp*y1/( omg*a )
                endif
                dy1 =  dy*( -dg*a + omg*da )
                dy2 = mexp*(dg*a + g*da)
                dF = wpre*sign(1.d0,rho)*( dy1 + dy2 )
            else
            dF = 0.d0
            endif
        else
            F = wpre*mexp*rho
            if (doforce) then
                dF = wpre*mexp
            else
                dF = 0.d0
            endif
        
        endif
    case (2)
        !     
        ! F = wpre*( sign(1+x)*(  ( ( 1 - exp(-w*(1+x)**2) )*abs(1+x) )^m +  m*exp(-w*(1+x)**2)*abs(1+x) ) - 1 )
        srho = 1.d0 + rho
        if ( abs(srho) > eps ) then
        
            y1 = w*srho**2
            ! to avoid underflow
            if (y1 > 30.d0) then
                g = 0.d0
            else
                g  = exp(-y1)
            endif
            omg = 1.d0 - g     
            a  = abs(srho)
            y1  =  ( omg*a )**mexp
            y2  = mexp*g*a
            F = wpre*(sign(1.d0,srho)*( y1 + y2) - 1.d0)
        
            if (doforce) then
                dg = -2.d0*w*srho*g
                da = sign(1.d0,srho)
                if (dabs(y1) < eps) then
                    dy = 0.d0
                else
                    dy  =  mexp*y1/( omg*a )
                endif
                dy1 =  dy*( -dg*a + omg*da )
                dy2 = mexp*(dg*a + g*da)
                dF = wpre*sign(1.d0,srho)*( dy1 + dy2 )
            endif
        else
            F = wpre*mexp*(srho - 1.d0)
            if (doforce) then
                dF = wpre*mexp
            endif
        
        endif
    case(3)
        if (doforce) then
            print *,'Gradient not yet implemented.'
            stop
        else
            a = abs(rho)
            if ( a > eps) then
                F = wpre*mexp*sign(1.d0,rho)*a*log(a)
            else
                F = 0.d0
            endif
        endif
        
    case(4)
        if (doforce) then
            print *,'Gradient not yet implemented.'
            stop
        else
            a = abs(rho)
            if ( a > eps) then
                F = wpre*mexp*sign(1.d0,rho)*tanh(a)
            else
                F = 0.d0
            endif
        endif
        
        case default
            print *,'Should not be here. Embedding function not implemented:', nfemb
        stop
        
    end select

end subroutine Fexp

!-----------------------------------------------------------------------------------------

subroutine innercutoff(e,ecutin,decutin,fcut,dfcut)

    implicit none
    ! computes inner cut-off based on energy
    double precision e, ecutin, decutin, fcut, dfcut, elow
    double precision, parameter :: pi = 3.14159265358979323846264338327950288419d0

    elow = ecutin - decutin
    if ( e >= ecutin ) then
        fcut = 0.d0
        dfcut = 0.d0
    elseif ( e <= elow ) then
        fcut = 1.d0
        dfcut = 0.d0       
    else
        fcut = 0.5d0*(1.d0 + cos( pi*(e-elow)/decutin ))
        dfcut = -0.5d0*sin( pi*(e-elow)/decutin )*pi/decutin
    endif

end subroutine  

!-----------------------------------------------------------------------------------------
