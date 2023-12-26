! routines for computation of spherical harmnics

!-----------------------------------------------------------------

subroutine pre_compute(dof,lmax,lmaxhalf)
    use tables, only: shindex_atomic,shindex_magnetic,alm,blm,clm,cl,dl,el
    implicit none

    double precision a,b,c,lsq,ld,l1,l2,l3,l4,msq
    integer l,m,k,lmaxhalf,lmax
    character(8) dof

    if (.not. allocated(alm)) allocate(alm(1:lmaxhalf))
    if (.not. allocated(blm)) allocate(blm(1:lmaxhalf))
    if (.not. allocated(cl)) allocate(cl(1:(lmax+1)))
    if (.not. allocated(dl)) allocate(dl(1:(lmax+1)))

    do l = 1, lmax

        lsq = l * l
        ld = 2 * l
        l1 = (4 * lsq - 1)
        l2 = lsq - ld + 1
        
        do m = 0, l-2

            msq = m * m
            a = sqrt(l1 / (lsq - msq))
            b = -sqrt((l2 - msq) / (4 * l2 - 1))

            if (trim(dof) == 'atomic') then
                k = shindex_atomic(l,m)
            elseif (trim(dof) == 'magnetic') then
                k = shindex_magnetic(l,m)
            endif
            
            alm(k) = a
            blm(k) = b
    
        enddo
        
    enddo

    do l = 1, lmax

    cl(l) = -sqrt(dble(1.d0 + 0.5d0 / dble(l)))
    dl(l) = sqrt(dble(2.d0 * (l - 1) + 3.d0))
    
    enddo

end subroutine
    
!--------------------------------------------------------------------------

subroutine compute_barplm(dof,lmax,lmaxhalf,rz)
    use tables, only: shindex_atomic,shindex_magnetic,alm,blm,dl,cl,plm,dplm
    implicit none

    ! requires -1 <= rz <= 1 , NO CHECKING IS PERFORMED !!!!!!!!!
    ! prefactors include 1/sqrt(2) factor compared to reference
    double precision rz
    double precision temp,sq3o1,sq3o2,sq1o4pi,sq4pi,sq3
    double precision t0,t0d,t0s,t, td, ts,tt,invst,plmlm1lm1
    integer l,m,k,k1,k2,lmax,lmaxhalf
    double precision Y00
    character(8) dof
    
    Y00 = 1.d0
    sq1o4pi = 0.28209479177387814347d0
    sq4pi = 3.54490770181103176384d0
    sq3 = 1.73205080756887719318d0
    sq3o2 = 1.22474487139158894067d0

    if (.not. allocated(plm) ) allocate(plm(1:lmaxhalf))
    if (.not. allocated(dplm) ) allocate(dplm(1:lmaxhalf))

    ! l=0, m=0
    ! plm(0, 0) = Y00/sq1o4pi; //= sq1o4pi
    l = 0
    m = 0
    
    if (trim(dof) == 'atomic') then
        k = shindex_atomic(l,m)
    elseif (trim(dof) == 'magnetic') then
        k = shindex_magnetic(l,m)
    endif
        
    plm(k) = Y00 
    dplm(k) = 0.d0

    if (lmax > 0) then

        ! l=1, m=0
        l = 1
        m = 0
        
        if (trim(dof) == 'atomic') then
            k = shindex_atomic(l,m)
        elseif (trim(dof) == 'magnetic') then
            k = shindex_magnetic(l,m)
        endif
        
        plm(k) = Y00 * sq3 * rz
        dplm(k) = Y00 * sq3

        ! l=1, m=1
        l = 1
        m = 1
        
        if (trim(dof) == 'atomic') then
            k = shindex_atomic(l,m)
        elseif (trim(dof) == 'magnetic') then
            k = shindex_magnetic(l,m)
        endif
        
        plm(k) = -sq3o2 * Y00
        dplm(k) = 0.0d0

        ! loop l = 2, lmax
        do l = 2, lmax
            do m = 0, l - 2
            
                if (trim(dof) == 'atomic') then
                    k = shindex_atomic(l,m)
                    k1 = shindex_atomic(l-1,m)
                    k2 = shindex_atomic(l-2,m)
                elseif (trim(dof) == 'magnetic') then
                    k = shindex_magnetic(l,m)
                    k1 = shindex_magnetic(l-1,m)
                    k2 = shindex_magnetic(l-2,m)
                endif

                plm(k) = alm(k) * (rz * plm(k1) + blm(k) * plm(k2))
                dplm(k) = alm(k) * (plm(k1) + rz * dplm(k1) + blm(k) * dplm(k2))
            enddo
            
            if (trim(dof) == 'atomic') then
                k = shindex_atomic(l-1,l-1)
                k1 = shindex_atomic(l,l-1)
                k2 = shindex_atomic(l,l)
            elseif (trim(dof) == 'magnetic') then
                k = shindex_magnetic(l-1,l-1)
                k1 = shindex_magnetic(l,l-1)
                k2 = shindex_magnetic(l,l)
            endif
            
            t = dl(l) * plm(k)
            plm(k1) = t * rz
            dplm(k1) = t
            plm(k2) = cl(l) * plm(k)
            dplm(k2) = 0.0
        enddo
    endif
    
end subroutine

!--------------------------------------------------------------------------

subroutine compute_ylm_cart(dof,rx,ry,rz,lmax,lmaxhalf,shvec,dshvec)

    ! requires rx^2 + ry^2 + rz^2 = 1 , NO CHECKING IS PERFORMED !!!!!!!!!
    use tables, only: shindex_atomic,shindex_magnetic,alm,blm,clm,dl,el,plm,splm,dplm
    !use functionparameters, only: lmaxhalf
    implicit none

    double precision rx,ry,rz
    double precision costheta,sintheta,cosphi,sinphi
    double precision eps
    double complex phase,phasem, mphasem1,dyx, dyy, dyz, rdy
    double precision c, s, ctcp, ctsp, c1, c2
    double precision s1, s2, two_cos_phi
    integer k,k1,m,lmax,l,lmaxhalf
    double complex shvec(1:(lmax+1)**2)
    double complex dshvec(1:3,1:(lmax+1)**2)
    double complex im
    character(8) dof

    im = dcmplx(0.d0,1.d0)
    
    call pre_compute(dof,lmax,lmaxhalf)

    !compute barplm
    call compute_barplm(dof,lmax,lmaxhalf,rz)
    
    phase = dcmplx(rx,ry)   
    

    !m = 0
    m = 0
    do l = 0, lmax

        if (trim(dof) == 'atomic') then
            k = shindex_atomic(l,m)
        elseif (trim(dof) == 'magnetic') then
            k = shindex_magnetic(l,m)
        endif
        
        shvec(k) = plm(k) + im * 0.d0

        dyz = dplm(k) + im * 0.d0 
        rdy = dreal(dyz) * rz + im * 0.d0

        dshvec(1,k) = -dreal(rdy) * rx + im * 0.d0
        dshvec(2,k) = -dreal(rdy) * ry + im * 0.d0
        dshvec(3,k) = (dreal(dyz) - dreal(rdy) * rz) + im * 0.d0
    enddo
    
    !m = 0
    m = 1
    do l = 1, lmax

        if (trim(dof) == 'atomic') then
            k = shindex_atomic(l,m)
        elseif (trim(dof) == 'magnetic') then
            k = shindex_magnetic(l,m)
        endif
        
        shvec(k) = phase * plm(k);

        dyx = plm(k) + im * 0.d0
        dyy = 0.0 + im * plm(k)
        dyz = dreal(phase) * dplm(k) + im * dimag(phase) * dplm(k)

        rdy = (rx * dreal(dyx) + rz * dreal(dyz)) + im * (ry * dimag(dyy) + rz * dimag(dyz))

        dshvec(1,k) = (dreal(dyx) - dreal(rdy) * rx) - im * dimag(rdy) * rx
        dshvec(2,k) = -dreal(rdy) * ry + im * (dimag(dyy) - dimag(rdy) * ry)
        dshvec(3,k) = (dreal(dyz) - dreal(rdy) * rz) + im * (dimag(dyz) - dimag(rdy) * rz)
        
    enddo

    ! m > 1
    phasem = phase
    do m = 2, lmax

        mphasem1 = dreal(phasem) * dble(m) + im * dimag(phasem) * dble(m)
        phasem = phasem * phase

        do l = m, lmax

            if (trim(dof) == 'atomic') then
                k = shindex_atomic(l,m)
            elseif (trim(dof) == 'magnetic') then
                k = shindex_magnetic(l,m)
            endif
        
            shvec(k) = dreal(phasem) * plm(k) + im * dimag(phasem) * plm(k)

            dyx = mphasem1 * plm(k)
            dyy = -dimag(dyx) + im * dreal(dyx)
            dyz = phasem * dplm(k)

            rdy = (rx * dreal(dyx) + ry * dreal(dyy) + rz * dreal(dyz)) + &
              &  im * (rx * dimag(dyx) + ry * dimag(dyy) + rz * dimag(dyz))

            dshvec(1,k) = (dreal(dyx) - dreal(rdy) * rx) + im * (dimag(dyx) - dimag(rdy) * rx)
            dshvec(2,k) = (dreal(dyy) - dreal(rdy) * ry) + im * (dimag(dyy) - dimag(rdy) * ry)
            dshvec(3,k) = (dreal(dyz) - dreal(rdy) * rz) + im * (dimag(dyz) - dimag(rdy) * rz)
        enddo
    enddo
    
!     fill-in m<0  
    do l = 1, lmax
        do m = 1, l
            phase = dble((-1)**(m))

            if (trim(dof) == 'atomic') then
                k = shindex_atomic(l,-m)
                k1 = shindex_atomic(l,m)
            elseif (trim(dof) == 'magnetic') then
                k = shindex_magnetic(l,-m)
                k1 = shindex_magnetic(l,m)
            endif            

            shvec(k) = conjg(shvec(k1)) * phase
            dshvec(1:3,k) = conjg(dshvec(1:3,k1))*phase


        enddo
    enddo

end subroutine

!--------------------------------------------------------------------------
