!
!  tools for setting up radial functions
!
!--------------------------------------------------------------------------

subroutine calcCheb(n,x,cheb1,dcheb1)

    implicit none
    integer n
    double precision x, cheb1(0:n), dcheb1(0:n)
    double precision cheb2(0:n)
    double precision tx
    double precision, parameter :: eps = 1.d-10
    integer m

    tx = 2.d0*x

    cheb1(0) = 1.d0
    dcheb1(0) = 0.d0
    cheb2(0) = 1.d0  
    if ( n > 0 ) then
        cheb1(1) = x
        cheb2(1) = tx
        do m = 1, n-1
            cheb1(m+1) = tx*cheb1(m) - cheb1(m-1)
            cheb2(m+1) = tx*cheb2(m) - cheb2(m-1)
        enddo
        do m = 1, n
            dcheb1(m) = dble(m)*cheb2(m-1)
        enddo
    endif

end subroutine calcCheb

!--------------------------------------------------------------------------

subroutine setuplookuptables
  
    implicit none
    integer ntot

    ! number of bins for lookup tables
    ntot = 10000

    ! set up radial functions
    call setuplookupRadspline_atomic(ntot)
    call setuplookupRadspline_magnetic(ntot)
 
end subroutine setuplookuptables

!--------------------------------------------------------------------------

subroutine lookupRadspline_atomic(r,nradbase,nradial,lmax,elei,elej,fr,dfr,gr,dgr,cr,dcr)

    use lookuptables, only: nlut, rscalelookup_atomic, lutgrs_atomic, lutfrs_atomic, luthcrs_atomic
    implicit none
    double precision r
    integer nradbase, nradial, lmax, elei, elej
    double precision gr(1:nradbase)
    double precision dgr(1:nradbase)
    double precision fr(1:nradial,0:lmax)
    double precision dfr(1:nradial,0:lmax)
    double precision cr, dcr

    integer nl, nr, l
    double precision x, wl, wl2, wl3, w2l1, w3l2, c(0:3)

    x = r*rscalelookup_atomic
    nl = floor(x)
  
    if ( nl <= 0 ) then
        print *,'Encountered very small distance.'
        print *,'Stopping.'
        stop
    endif
  
    if ( nl >= nlut ) then
        gr(1:nradbase) = 0.d0
        dgr(1:nradbase) = 0.d0
        fr(1:nradial,0:lmax) = 0.d0
        dfr(1:nradial,0:lmax) = 0.d0
        cr = 0.d0
        dcr = 0.d0
    else
        wl =  x - dble(nl)
        wl2 = wl*wl
        wl3 = wl2*wl
        w2l1 = 2.d0*wl
        w3l2 = 3.d0*wl2
        do nr = 1, nradbase
            c(0:3) =  lutgrs_atomic(0:3,nr,nl,elei,elej)
            gr(nr) = c(0) + c(1)*wl + c(2)*wl2  + c(3)*wl3
            dgr(nr) =     ( c(1)    + c(2)*w2l1 + c(3)*w3l2 )*rscalelookup_atomic
        enddo
        do l = 0, lmax
            do nr = 1, nradial
                c(0:3) =  lutfrs_atomic(0:3,nr,l,nl,elei,elej)
                fr(nr,l) = c(0) + c(1)*wl + c(2)*wl2  + c(3)*wl3
                dfr(nr,l) =     ( c(1)    + c(2)*w2l1 + c(3)*w3l2 )*rscalelookup_atomic
            enddo
        enddo
        c(0:3) =  luthcrs_atomic(0:3,nl,elei,elej)
        cr = c(0) + c(1)*wl + c(2)*wl2  + c(3)*wl3
        dcr =     ( c(1)    + c(2)*w2l1 + c(3)*w3l2 )*rscalelookup_atomic
    endif

end subroutine lookupRadspline_atomic

!--------------------------------------------------------------------------

subroutine lookupRadspline_magnetic(r,nradbase,nradial,lmax,elei,elej,fr,dfr,gr,dgr)

    use lookuptables, only: nlut, rscalelookup_magnetic, lutgrs_magnetic, lutfrs_magnetic
    !luthcrs_magnetic
    implicit none
    double precision r
    integer nradbase, nradial, lmax, elei, elej
    double precision gr(1:nradbase)
    double precision dgr(1:nradbase)
    double precision fr(1:nradial,0:lmax)
    double precision dfr(1:nradial,0:lmax)

    integer nl, nr, l
    double precision x, wl, wl2, wl3, w2l1, w3l2, c(0:3)

    x = r*rscalelookup_magnetic
    nl = floor(x)

    if ( nl < 0 ) then
        print *,'Encountered very small moment.'
        print *,'Stopping.'
        stop
    endif
  
    if ( nl .eq. 0.d0 ) then
        gr(1:nradbase) = 1.d0
        dgr(1:nradbase) = 0.d0
        fr(1:nradial,0:lmax) = 1.d0
        dfr(1:nradial,0:lmax) = 0.d0
    else
        wl =  x - dble(nl)
        wl2 = wl*wl
        wl3 = wl2*wl
        w2l1 = 2.d0*wl
        w3l2 = 3.d0*wl2
        do nr = 1, nradbase
            c(0:3) =  lutgrs_magnetic(0:3,nr,nl,elei,elej)
            gr(nr) = c(0) + c(1)*wl + c(2)*wl2  + c(3)*wl3
            dgr(nr) =     ( c(1)    + c(2)*w2l1 + c(3)*w3l2 )*rscalelookup_magnetic
        enddo
        do l = 0, lmax
            do nr = 1, nradial
                c(0:3) =  lutfrs_magnetic(0:3,nr,l,nl,elei,elej)
                fr(nr,l) = c(0) + c(1)*wl + c(2)*wl2  + c(3)*wl3
                dfr(nr,l) =     ( c(1)    + c(2)*w2l1 + c(3)*w3l2 )*rscalelookup_magnetic
            enddo
        enddo
    endif

end subroutine lookupRadspline_magnetic

!--------------------------------------------------------------------------

subroutine setuplookupRadspline_atomic(ntot)

    use lookuptables
    use tables
    use functionparameters
    use global
    use modneigh
    implicit none
    integer ntot, n, nr, l
    double precision x, y
    double precision cu, lam, dc
    integer elei, elej

    double precision f0, f1, f1d1, f0d1
    double precision c(0:3)

    double precision gr(1:nradbase_atomic)
    double precision dgr(1:nradbase_atomic)
    double precision fr(1:nradial3_atomic,0:l3max_atomic)
    double precision dfr(1:nradial3_atomic,0:l3max_atomic)

    double precision f1g(1:nradbase_atomic)
    double precision f1gd1(1:nradbase_atomic)
    double precision f1f(1:nradbase_atomic,0:l3max_atomic)
    double precision f1fd1(1:nradbase_atomic,0:l3max_atomic)
    double precision dfrprev(1:nradbase_atomic,0:l3max_atomic)

    double precision f1cr, f1crd1
    double precision cr, dcr
    double precision lc, pc

    nlut = ntot

    ! cutoff is atomic cutoff
    rscalelookup_atomic = dble(nlut)/cutoff
    invrscalelookup_atomic = 1.d0/rscalelookup_atomic

    if (.not.allocated(lutfrs_atomic)) then
        allocate( lutfrs_atomic(0:3,1:nradial3_atomic,0:l3max_atomic,1:nlut,1:nelements,1:nelements))
        allocate( lutgrs_atomic(0:3,1:nradbase_atomic,1:nlut,1:nelements,1:nelements))
        allocate( luthcrs_atomic(0:3,1:nlut,1:nelements,1:nelements))
        allocate(scoreRnl_atomic(1:3))
    endif

    ! at r = rcut + eps the function and its derivatives is zero

    scoreRnl_atomic = 0.d0

    do elei = 1, nelements
        do elej = 1, nelements

            f1g = 0.d0
            f1gd1 = 0.d0
            f1f = 0.d0
            f1fd1 = 0.d0
            f1cr = 0.d0
            f1crd1 = 0.d0
            dfrprev = 0.d0 
            
            do n = nlut, 1, -1          
            
                x = invrscalelookup_atomic*dble(n)
                lam = lambda(1,elei,elej)
                cu = cut(1,elei,elej)
                dc = dcut(1,elei,elej)

                ! set up radial functions
                call radbase(radtype_atomic,lam,nradbase_atomic,cu,dc,x,gr,dgr)
                call radfunc(nradbase_atomic,nelements,elei,elej,nradial3_atomic,l3max_atomic,crad_atomic,gr,dgr,fr,dfr)

                ! hard core repulsion
                pc = corepara(1,elei,elej)
                lc = corepara(2,elei,elej)
                call radcore(x,pc,lc,cu,cr,dcr)
            
                do nr = 1, nradbase_atomic
                    f0 = gr(nr)
                    f1 = f1g(nr)
                    f0d1 = dgr(nr)*invrscalelookup_atomic
                    f1d1 = f1gd1(nr)
                    ! evaluate coefficients 
                    c(0) = f0
                    c(1) = f0d1
                    c(2) =  3.d0*(f1 - f0) - f1d1 - 2.d0*f0d1
                    c(3) = -2.d0*(f1 - f0) + f1d1 + f0d1
                    ! store coefficients
                    lutgrs_atomic(0:3,nr,n,elei,elej) = c(0:3)
                    ! evalute function values and derivatives at current position
                    f1g(nr) = c(0)
                    f1gd1(nr) = c(1)
                enddo
                
                do l = 0, l3max_atomic
                    do nr = 1, nradial3_atomic
                        f0 = fr(nr,l)
                        f1 = f1f(nr,l)
                        f0d1 = dfr(nr,l)*invrscalelookup_atomic
                        f1d1 = f1fd1(nr,l)
                        ! evaluate coefficients 
                        c(0) = f0
                        c(1) = f0d1
                        c(2) =  3.d0*(f1 - f0) - f1d1 - 2.d0*f0d1
                        c(3) = -2.d0*(f1 - f0) + f1d1 + f0d1
                        ! store coefficients
                        lutfrs_atomic(0:3,nr,l,n,elei,elej) = c(0:3)
                        ! evalute and store function values and derivatives at current position
                        f1f(nr,l) = c(0)
                        f1fd1(nr,l) = c(1)
                        y = x/cutoff
                        y = y*y
                        scoreRnl_atomic(1) = scoreRnl_atomic(1) + (fr(nr,l)**2)*y
                        scoreRnl_atomic(2) = scoreRnl_atomic(2) + (dfr(nr,l)**2)*y
                        scoreRnl_atomic(3) = scoreRnl_atomic(3) + ( ( ( dfr(nr,l) - dfrprev(nr,l) )*invrscalelookup_atomic)**2 )*y
                        dfrprev(nr,l) = dfr(nr,l)
                    enddo
                enddo
                
                f0 = cr
                f1 = f1cr
                f0d1 = dcr*invrscalelookup_atomic
                f1d1 = f1crd1
                ! evaluate coefficients 
                c(0) = f0
                c(1) = f0d1
                c(2) =  3.d0*(f1 - f0) - f1d1 - 2.d0*f0d1
                c(3) = -2.d0*(f1 - f0) + f1d1 + f0d1
                ! store coefficients
                luthcrs_atomic(0:3,n,elei,elej) = c(0:3)
                ! evalute function values and derivatives at current position
                f1cr = c(0)
                f1crd1 = c(1)

            enddo

        enddo
    enddo

    scoreRnl_atomic(1:3) = scoreRnl_atomic(1:3)/dble(nlut*nelements*nelements)

    totscoreRnl_atomic = sum(wscoreRnl_atomic(1:3)*scoreRnl_atomic(1:3))

end subroutine setuplookupRadspline_atomic

!--------------------------------------------------------------------------

subroutine setuplookupRadspline_magnetic(ntot)

    use lookuptables
    use tables
    use functionparameters
    use global
    use modneigh
    implicit none
    integer ntot, n, nr, l
    double precision x, y
    double precision cu, lam, dc
    integer elei, elej

    double precision f0, f1, f1d1, f0d1
    double precision c(0:3)

    double precision gr(1:nradbase_magnetic)
    double precision dgr(1:nradbase_magnetic)
    double precision fr(1:nradial2_magnetic,0:l2max_magnetic)
    double precision dfr(1:nradial2_magnetic,0:l2max_magnetic)

    double precision f1g(1:nradbase_magnetic)
    double precision f1gd1(1:nradbase_magnetic)
    double precision f1f(1:nradbase_magnetic,0:l2max_magnetic)
    double precision f1fd1(1:nradbase_magnetic,0:l2max_magnetic)
    double precision dfrprev(1:nradbase_magnetic,0:l2max_magnetic)

    double precision lc, pc

    nlut = ntot

    ! cutoff is magnetic cutoff
    rscalelookup_magnetic = dble(nlut)/cutoff_magnetic
    invrscalelookup_magnetic = 1.d0/rscalelookup_magnetic

    if (.not.allocated(lutfrs_magnetic)) then
        allocate( lutfrs_magnetic(0:3,1:nradial2_magnetic,0:l2max_magnetic,1:nlut,1:nelements,1:nelements))
        allocate( lutgrs_magnetic(0:3,1:nradbase_magnetic,1:nlut,1:nelements,1:nelements))
        allocate(scoreRnl_magnetic(1:3))
    endif

    ! at r = rcut + eps the function and its derivatives is zero


    scoreRnl_magnetic = 0.d0

    do elei = 1, nelements
        do elej = 1, nelements

            f1g = 0.d0
            f1gd1 = 0.d0
            f1f = 0.d0
            f1fd1 = 0.d0
            dfrprev = 0.d0 
            
            do n = nlut, 1, -1
            
                x = invrscalelookup_magnetic*dble(n)
                lam = lambda(2,elei,elej)
                cu = cut(2,elei,elej)
                dc = dcut(2,elei,elej)
                
                ! set up radial functions
                call radbase(radtype_magnetic,lam,nradbase_magnetic,cu,dc,x,gr,dgr)
                call radfunc(nradbase_magnetic,nelements,elei,elej,nradial2_magnetic,l2max_magnetic,crad_magnetic,gr,dgr,fr,dfr)

                do nr = 1, nradbase_magnetic
                    f0 = gr(nr)
                    f1 = f1g(nr)
                    f0d1 = dgr(nr)*invrscalelookup_magnetic
                    f1d1 = f1gd1(nr)
                    ! evaluate coefficients 
                    c(0) = f0
                    c(1) = f0d1
                    c(2) =  3.d0*(f1 - f0) - f1d1 - 2.d0*f0d1
                    c(3) = -2.d0*(f1 - f0) + f1d1 + f0d1
                    ! store coefficients
                    lutgrs_magnetic(0:3,nr,n,elei,elej) = c(0:3)
                    ! evalute function values and derivatives at current position
                    f1g(nr) = c(0)
                    f1gd1(nr) = c(1)
                enddo
                
                do l = 0, l2max_magnetic
                    do nr = 1, nradial2_magnetic
                        f0 = fr(nr,l)
                        f1 = f1f(nr,l)
                        f0d1 = dfr(nr,l)*invrscalelookup_magnetic
                        f1d1 = f1fd1(nr,l)
                        ! evaluate coefficients 
                        c(0) = f0
                        c(1) = f0d1
                        c(2) =  3.d0*(f1 - f0) - f1d1 - 2.d0*f0d1
                        c(3) = -2.d0*(f1 - f0) + f1d1 + f0d1
                        ! store coefficients
                        lutfrs_magnetic(0:3,nr,l,n,elei,elej) = c(0:3)
                        ! evalute and store function values and derivatives at current position
                        f1f(nr,l) = c(0)
                        f1fd1(nr,l) = c(1)
                        y = x/cutoff_magnetic
                        y = (1-y)*(1-y)
                        scoreRnl_magnetic(1) = scoreRnl_magnetic(1) + (fr(nr,l)**2)*y
                        scoreRnl_magnetic(2) = scoreRnl_magnetic(2) + (dfr(nr,l)**2)*y
                        scoreRnl_magnetic(3) = scoreRnl_magnetic(3) + & 
                                    &       ( ( ( dfr(nr,l) - dfrprev(nr,l) )*invrscalelookup_magnetic)**2 )*y
                        dfrprev(nr,l) = dfr(nr,l)
                    enddo
                enddo
            enddo

        enddo
    enddo

    scoreRnl_magnetic(1:3) = scoreRnl_magnetic(1:3)/dble(nlut*nelements*nelements)

    totscoreRnl_magnetic = sum(wscoreRnl_magnetic(1:3)*scoreRnl_magnetic(1:3))

end subroutine setuplookupRadspline_magnetic

!--------------------------------------------------------------------------

subroutine radbase(radtype,radpara,nradbase,cutoff,dcutoff,r,gr,dgr)

    implicit none
    integer radtype, nradbase
    double precision radpara, r, cutoff, dcutoff
    double precision cheb(0:nradbase), dcheb(0:nradbase)
    double precision gr(1:nradbase), dgr(1:nradbase)

    integer n, ncheb, nrad, ll
    double precision x, dx, fcut, dfcut, y1, dy1, y2, x0, env, denv, y0
    double precision, parameter :: pi = 3.14159265358979323846264338327950288419d0
    double precision dshift

    if (r >= cutoff ) then
        gr(1:nradbase) = 0.d0
        dgr(1:nradbase) = 0.d0
    else
        select case (radtype)
        case(0)
            print *,'Not working.'
            stop
            x = r/cutoff
            y1 = exp(-radpara*x)
            dy1 = -radpara*y1/cutoff
            y2 = exp(-radpara)
            ! x = 1 if r = 0 and x = -1 if r = cutoff
            x = 1.d0 - 2.d0*( y1 - y2 )/( 1.d0 - y2 )
            dx = - 2.d0*y1/( 1.d0 - y2 )
            call calcCheb(nradbase-1,x,cheb,dcheb)
            gr(1:nradbase) = (0.5d0 - 0.5d0*cheb(1:nradbase))
            dgr(1:nradbase) = -0.5d0*dcheb(1:nradbase)*dx    
        case(1)
            print *,'Not working.'
            stop
            do n = 1, nradbase
                x = ( ( float(2*(nradbase+1))/float(n+1) )*r/cutoff)**2
                gr(n) = exp(-x)
            enddo
        case(2)
            print *,'Not working.'
            stop
            x = r/cutoff
            y1 = exp(-radpara*x**2)
            y2 = exp(-radpara)
            ! x = 1 if r = 0 and x = -1 if r = cutoff
            x = 1.d0 - 2.d0*( y1 - y2 )/( 1.d0 - y2 )
            call calcCheb(nradbase,x,cheb,dcheb)
            gr(1:nradbase) = 0.5d0 - 0.5d0*cheb(1:nradbase)    
        case(3)
            ! Cheb exp scaling including cos envelope
            x0 = r/cutoff
            y1 = exp(-radpara*x0)
            dy1 = -radpara*y1/cutoff
            y2 = exp(-radpara)
            ! x = 1 if r = 0 and x = -1 if r = cutoff
            x = 1.d0 - 2.d0*( y1 - y2 )/( 1.d0 - y2 )
            dx = - 2.d0*dy1/( 1.d0 - y2 )
            ncheb = nradbase - 1
            call calcCheb(ncheb,x,cheb,dcheb)
            gr(1) = cheb(0)
            dgr(1) = dcheb(0)*dx
            do n = 2, nradbase
                gr(n) = 0.5d0 - 0.5d0*cheb(n-1)
                dgr(n) = -0.5d0*dcheb(n-1)*dx
            enddo
            env  =  0.5d0*( 1.d0 + cos( pi*x0 ) )
            denv = -0.5d0*sin( pi*x0 )*pi/cutoff
            do n = 1, nradbase
                dgr(n) = gr(n)*denv + dgr(n)*env
                gr(n) = gr(n)*env
            enddo
        case(4)
            ! Cheb exp scaling without cos envelope 
            x0 = r/cutoff
            y1 = exp(-radpara*x0)
            dy1 = -radpara*y1/cutoff
            y2 = exp(-radpara)
            ! x = 1 if r = 0 and x = -1 if r = cutoff
            x = 1.d0 - 2.d0*( y1 - y2 )/( 1.d0 - y2 )
            dx = - 2.d0*dy1/( 1.d0 - y2 )
            ncheb = nradbase - 1
            call calcCheb(ncheb,x,cheb,dcheb)
            gr(1) = cheb(0)
            dgr(1) = dcheb(0)*dx
            do n = 2, nradbase
                gr(n) = 0.5d0 - 0.5d0*cheb(n-1)
                dgr(n) = -0.5d0*dcheb(n-1)*dx    
            enddo
        case(5)
            ! Cheb without cos envelope, power law variable

            x0 = r/cutoff
            y0 = (1.d0 - x0)

            dy1 = y0**(radpara-1.d0)
            y1 = 1.d0 - dy1*y0
            dy1 = radpara/cutoff*dy1
            
            ! x = 1 if r = 0 and x = -1 if r = cutoff
            x = 2.d0*y1 - 1.d0
            dx = 2.d0*dy1
            ncheb = nradbase
            call calcCheb(ncheb,x,cheb,dcheb)
            do n = 1, nradbase
                gr(n) = 0.5d0 - 0.5d0*cheb(n)
                dgr(n) = -0.5d0*dcheb(n)*dx    
            enddo
        case(6)
            ! Cheb with cos envelope, power law variable

            x0 = r/cutoff
            y0 = (1.d0 - x0)

            dy1 = y0**(radpara-1.d0)
            y1 = 1.d0 - dy1*y0
            dy1 = radpara/cutoff*dy1
            
            ! x = 1 if r = 0 and x = -1 if r = cutoff
            x = 2.d0*y1 - 1.d0
            dx = 2.d0*dy1
            ncheb = nradbase -1
            call calcCheb(ncheb,x,cheb,dcheb)
            gr(1) = cheb(0)
            dgr(1) = dcheb(0)*dx
            do n = 1, nradbase -1
                gr(n+1) = 0.5d0 - 0.5d0*cheb(n)
                dgr(n+1) = -0.5d0*dcheb(n)*dx    
            enddo
            env  =  0.5d0*( 1.d0 + cos( pi*x0 ) )
            denv = -0.5d0*sin( pi*x0 )*pi/cutoff
            do n = 1, nradbase
                dgr(n) = gr(n)*denv + dgr(n)*env
                gr(n)  = gr(n)*env
            enddo
        case (7)
            ! Cheb without envelope, power law
            x0 = r/cutoff
            x = 1.d0 - 2.d0*x0*x0
            dx = -( 4.d0/cutoff )*x0

            ncheb = nradbase
            call calcCheb(ncheb,x,cheb,dcheb)
            do n = 1, nradbase
                gr(n) = cheb(n)
                dgr(n) = dcheb(n)*dx
            enddo

        case (8)
            ! Cheb without envelope, power law
            x0 = r/cutoff
            x = 1.d0 - 2.d0*x0
            dx = -2.d0/cutoff 

            ncheb = nradbase
            call calcCheb(ncheb,x,cheb,dcheb)

            do n = 1, nradbase
                gr(n) = cheb(n) - 1.d0
                dgr(n) = dcheb(n)*dx
            enddo
        case default
            print *,'radtype unknown', radtype
            stop
        end select
  
        ! for radtype = 3 a smooth cutoff is already included in the basis function  
        dx = cutoff - dcutoff
        if (r > dx ) then
            fcut = 0.5d0*( 1.d0 + cos( pi*(r-dx)/dcutoff ) )
            dfcut = -0.5d0*sin( pi*(r-dx)/dcutoff )*pi/dcutoff
            dgr(1:nradbase) = gr(1:nradbase)*dfcut + dgr(1:nradbase)*fcut
            gr(1:nradbase) = gr(1:nradbase)*fcut
        endif

    endif

end subroutine radbase

!--------------------------------------------------------------------------

subroutine radfunc(nradbase,nelements,ele1,ele2,nradial,lmax,crad,gr,dgr,fr,dfr)

    implicit none
    integer nradbase, nelements, ele1, ele2, nradial, lmax
    double precision crad(1:nradbase,1:nradial,0:lmax,1:nelements,1:nelements)
    double precision gr(1:nradbase), dgr(1:nradbase)
    double precision fr(1:nradial,0:lmax), dfr(1:nradial, 0:lmax)
    integer n,l,nr

    do l = 0, lmax
        do nr = 1, nradial
            fr(nr,l) = sum(crad(1:nradbase,nr,l,ele1,ele2)*gr(1:nradbase))
            dfr(nr,l) = sum(crad(1:nradbase,nr,l,ele1,ele2)*dgr(1:nradbase))
        enddo
    enddo
     
end subroutine radfunc

!--------------------------------------------------------------------------

subroutine radcore(r,pre,lambda,cutoff,cr,dcr)

    implicit none
    integer coretype
    double precision lambda, pre
    double precision cutoff
    double precision r
    double precision cr
    double precision dcr
    double precision, parameter :: pi = 3.14159265358979323846264338327950288419d0

    double precision r2, lr2, y, x0, env, denv

    if (r < cutoff ) then
        ! repulsion strictly positive and decaying
        pre = dabs(pre)
        lambda = dabs(lambda)
        
        r2 = r*r
        lr2 = lambda*r2
        if ( lr2 < 50.d0 ) then
            y = exp(-lr2)
            cr = pre*y/r
            dcr = -pre*y*(2.d0*lr2 + 1.d0)/r2
            x0 = r/cutoff
            env  =  0.5d0*( 1.d0 + cos( pi*x0 ) )
            denv = -0.5d0*sin( pi*x0 )*pi/cutoff
            dcr = cr*denv + dcr*env
            cr = cr*env
        else
            cr = 0.d0
            dcr = 0.d0
        endif
    else
        cr = 0.d0
        dcr = 0.d0
    endif

end subroutine radcore

!--------------------------------------------------------------------------

