!
!
! code for setting up basis functions B and gradients
!
!---------------------------------------------------------------------
subroutine prepB
      
    use core
    use global
    use functionparameters
    use tables
    use modprepA
    use modneigh
    implicit none

    double precision b1(1:nvar1)
    double precision b2(1:nvar2)
    double precision b3(1:nvar3)
    double precision b4(1:nvar4)
    double precision b5(1:nvar5)
    double precision b6(1:nvar6)
    double precision b7(1:nvar7)

    double complex xb1, xb2, xb3, xb4, xb5, xb6     
    
    
    integer n, l, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12
    integer m, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12
    integer k, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12
    integer k_atomic, k_magnetic, m_atomic, m_magnetic
    integer lint1, lint2
    integer na, na1, na2, na3, nr1, nr2, nr3, nr4, nr5, nr6, nr7, nr8, nr9, nr10, nr11
    integer mat, atomi, nat, nei, nclu, nmult, nr
    integer pack1count, pack2count
    integer  ij, ijstart, ijstop, ji
    integer ele, ele1, ele2, ele3, ele4, ele5, ele6
    
    double precision y, x
    double complex xiy
    double precision fr(1:nradbase), dfr(1:nradbase)
    double precision printcount, dprintcount

    double complex db2(-l2max_magnetic:l2max_magnetic,1:nvar2)
    double complex db3(1:2,-l3max_atomic:l3max_atomic,-l3max_magnetic:l3max_magnetic,1:nvar3)
    double complex db4(1:3,-l4max_atomic:l4max_atomic,-l4max_magnetic:l4max_magnetic,1:nvar4)
    double complex db5(1:4,-l5max_atomic:l5max_atomic,-l5max_magnetic:l5max_magnetic,1:nvar5)
    double complex db6(1:5,-l6max_atomic:l6max_atomic,-l6max_magnetic:l6max_magnetic,1:nvar6)
    
    double complex db2_0(-l2max_magnetic:l2max_magnetic,1:nvar2)
    double complex db3_0(-l3max_atomic:l3max_atomic,1:nvar3)
    double complex db4_0(-l4max_atomic:l4max_atomic,1:nvar4)
    double complex db5_0(-l5max_atomic:l5max_atomic,1:nvar5)
    double complex db6_0(-l6max_atomic:l6max_atomic,1:nvar6)
    
    double complex ap0, ap1, ap2, ap3, ap4, ap5, ap6
    double precision cf, cf_atomic, cf_magnetic
    
    printcount = 0.d0
    dprintcount = 0.1d0

      
    do atomi = 1, n2atom
        b1 = 0.d0
        
        do nclu = 1, nvar1

            ele = b1occ(nclu)
            nr = b1index(nclu)  ! n0'
            
            xiy = (0.d0,0.d0)
                
            k = shindex_magnetic(0,0)
                
            ap0 = grpack_magnetic_0(nr,atomi) - 1.d0

            xiy = ap0 
                
            b1(nclu) = dreal(xiy)

        enddo
        b1pack(1:nvar1,atomi) = b1(1:nvar1)

    enddo

    if (nradial2_magnetic == 0) then
       goto 1000
    endif      
      
    !--------------------------
    ! 2-body contributions
    do atomi = 1, n2atom
        b2 = 0.d0
        db2 = 0.d0
        db2_0 = (0.d0,0.d0)
        do nclu = 1, nvar2
            ele = b2occ(nclu)
            nr1 = b2index(1,nclu)  ! n0'
            nr2 = b2index(2,nclu)  ! n1'
            nr3 = b2index(3,nclu)  ! n1
            l = b2index(4,nclu)    ! l0'  (=l1')
            k = b2map(nclu)
            
            xiy = (0.d0,0.d0)
            
            do m = 1, madd2(k)
                m1 = b2lm(k,m)
                m2 = -m1
                
                k1 = shindex_magnetic(l,m1)
                k2 = shindex_magnetic(l,m2)
                
                ap0 = abase2pack_magnetic(nr1,k1,atomi)
                ap1 = abase2pack(nr3,nr2,k2,ele,atomi)
                
                cf = cg01_magnetic(k,m1)
            
                xiy = xiy + cf*ap0*ap1
                db2(m2,nclu) = db2(m2,nclu) + cf*ap0
                db2_0(m1,nclu) = db2_0(m1,nclu) + cf*ap1

                
            enddo
            b2(nclu) = dreal(xiy)
        enddo
        b2pack(1:nvar2,atomi) = b2(1:nvar2)
        db2pack(-l2max_magnetic:l2max_magnetic,1:nvar2,atomi) = db2(-l2max_magnetic:l2max_magnetic,1:nvar2) 
        db2pack_0(-l2max_magnetic:l2max_magnetic,1:nvar2,atomi) = db2_0(-l2max_magnetic:l2max_magnetic,1:nvar2)
        b2hcpack(atomi) = abase2hcpack(atomi)
    enddo

    !--------------------------
    ! 3-body contributions

    if ((nradial3_atomic == 0).and.(nradial3_magnetic == 0)) then
       goto 1000
    endif
    
    printcount = 0.d0
    dprintcount = 0.1d0
    nmult = 0

    do atomi = 1, n2atom
       b3 = 0.d0
       db3 = (0.d0,0.d0)
       db3_0 = (0.d0,0.d0)
       do nclu = 1, nvar3
          ele1 = b3occ(1,nclu)
          ele2 = b3occ(2,nclu)
          nr1 = b3index(1,nclu)  ! n1
          nr2 = b3index(2,nclu)  ! n2
          nr3 = b3index(3,nclu)  ! n0'
          nr4 = b3index(4,nclu)  ! n1'
          nr5 = b3index(5,nclu)  ! n2'
          l1  = b3index(6,nclu)  ! l1 (=l2)
          l2  = b3index(7,nclu)  ! l0'
          l3  = b3index(8,nclu)  ! l1'
          l4  = b3index(9,nclu)  ! l2'
          k_atomic = b3map_atomic(nclu)
          k_magnetic = b3map_magnetic(nclu)
          
          xiy = (0.d0,0.d0)
          
          do m_atomic = 1, madd3_atomic(k_atomic)             
             m1 =  b3lm_atomic(k_atomic,m_atomic)   ! m1
             m2 = -m1                               ! -m1
             
             k1 = shindex_atomic(l1,m1)
             k2 = shindex_atomic(l1,m2)
                          
             cf_atomic = cg02_atomic(k_atomic,m1)
             
             do m_magnetic = 1, madd3_magnetic(k_magnetic)
                m3 = b3lm_magnetic(k_magnetic,m_magnetic,1)  ! m0'
                m4 = b3lm_magnetic(k_magnetic,m_magnetic,2)  ! m1'
                m5 = b3lm_magnetic(k_magnetic,m_magnetic,3)  ! -M01'
             
                k3 = shindex_magnetic(l2,m3)
                k4 = shindex_magnetic(l3,m4)
                k5 = shindex_magnetic(l4,m5)
                
                ap0 = abase2pack_magnetic(nr3,k3,atomi)
                ap1 = abasepack(nr1,k1,nr4,k4,ele1,atomi)
                ap2 = abasepack(nr2,k2,nr5,k5,ele2,atomi)
                
                ! coupling coefficient for two spins
                cf_magnetic = cg02_magnetic(k_magnetic,m3,m4,m5)
                cf = cf_atomic*cf_magnetic
                                
                ! B
                ! multiply real and imag part seperately
                xiy = xiy + cf*ap0*ap1*ap2
                
                ! gradients dB/dA
                db3(1,m1,m4,nclu) = db3(1,m1,m4,nclu) + cf*ap0*ap2
                db3(2,m2,m5,nclu) = db3(2,m2,m5,nclu) + cf*ap0*ap1
                                
                ! gradients dB/dA_0
                db3_0(m3,nclu) = db3_0(m3,nclu) + cf*ap1*ap2
                
             enddo
          enddo
                    
          b3(nclu) = dreal( xiy )

       enddo
       
       b3pack(1:nvar3,atomi) = b3(1:nvar3)
       db3pack(1:2,-l3max_atomic:l3max_atomic,-l3max_magnetic:l3max_magnetic,1:nvar3,atomi) =  &
         &  db3(1:2,-l3max_atomic:l3max_atomic,-l3max_magnetic:l3max_magnetic,1:nvar3)
       db3pack_0(-l3max_magnetic:l3max_magnetic,1:nvar3,atomi) = db3_0(-l3max_magnetic:l3max_magnetic,1:nvar3)
       
    enddo

    !--------------------------
    ! 4-body contributions

    if ((nradial4_atomic == 0).and.(nradial4_magnetic == 0)) then
       goto 1000
    endif
    
    printcount = 0.d0
    dprintcount = 0.1d0
    nmult = 0

    do atomi = 1, n2atom
       b4 = 0.d0
       db4 = (0.d0,0.d0)
       db4_0 = (0.d0,0.d0)
       do nclu = 1, nvar4
          ele1 = b4occ(1,nclu)
          ele2 = b4occ(2,nclu)
          ele3 = b4occ(3,nclu)          
          nr1 = b4index(1,nclu)  ! n1
          nr2 = b4index(2,nclu)  ! n2
          nr3 = b4index(3,nclu)  ! n3
          nr4 = b4index(4,nclu)  ! n0'
          nr5 = b4index(5,nclu)  ! n1'
          nr6 = b4index(6,nclu)  ! n2'
          nr7 = b4index(7,nclu)  ! n3'
          l1  = b4index(8,nclu)  ! l1
          l2  = b4index(9,nclu)  ! l2 
          l3  = b4index(10,nclu) ! l3
          l4  = b4index(11,nclu) ! l0'
          l5  = b4index(12,nclu) ! l1'
          l6  = b4index(13,nclu) ! l2'
          l7  = b4index(14,nclu) ! l3'
          l8  = b4index(15,nclu) ! L01'          
          k_atomic = b4map_atomic(nclu)
          k_magnetic = b4map_magnetic(nclu)
          
          xiy = (0.d0,0.d0)
          
          do m_atomic = 1, madd4_atomic(k_atomic)
             m1 = b4lm_atomic(k_atomic,m_atomic,1)   ! m1
             m2 = b4lm_atomic(k_atomic,m_atomic,2)   ! m2
             m3 = b4lm_atomic(k_atomic,m_atomic,3)   ! m3
             
             k1 = shindex_atomic(l1,m1)
             k2 = shindex_atomic(l2,m2)
             k3 = shindex_atomic(l3,m3)
             
             cf_atomic = cg03_atomic(k_atomic,m1,m2,m3)
             
             do m_magnetic = 1, madd4_magnetic(k_magnetic)
                m4 = b4lm_magnetic(k_magnetic,m_magnetic,1)
                m5 = b4lm_magnetic(k_magnetic,m_magnetic,2)
                m6 = b4lm_magnetic(k_magnetic,m_magnetic,3)
                m7 = b4lm_magnetic(k_magnetic,m_magnetic,4)
                
                k4 = shindex_magnetic(l4,m4)                
                k5 = shindex_magnetic(l5,m5)
                k6 = shindex_magnetic(l6,m6)
                k7 = shindex_magnetic(l7,m7)
                
                ap0 = abase2pack_magnetic(nr4,k4,atomi)
                ap1 = abasepack(nr1,k1,nr5,k5,ele1,atomi)
                ap2 = abasepack(nr2,k2,nr6,k6,ele2,atomi)
                ap3 = abasepack(nr3,k3,nr7,k7,ele3,atomi)
                
                cf_magnetic = cg03_magnetic(k_magnetic,m4,m5,m6,m7)
                cf = cf_atomic*cf_magnetic
                
                ! B
                xb3 = ap0*ap2*ap3
                xb4 = cf*ap1
                xiy = xiy + xb4*xb3
                
                ! gradients dB/dA
                db4(1,m1,m5,nclu) = db4(1,m1,m5,nclu) + cf*xb3
                db4(2,m2,m6,nclu) = db4(2,m2,m6,nclu) + ap0*xb4*ap3
                db4(3,m3,m7,nclu) = db4(3,m3,m7,nclu) + ap0*xb4*ap2
                
                ! gradient dB/dA_0
                db4_0(m4,nclu) = db4_0(m4,nclu) + xb4*ap2*ap3 
                
             enddo   
          enddo
                       
          ! from here we are back to real numbers for b4                   
          b4(nclu) = dreal( xiy )
       enddo
       b4pack(1:nvar4,atomi) = b4(1:nvar4)
       db4pack(1:3,-l4max_atomic:l4max_atomic,-l4max_magnetic:l4max_magnetic,1:nvar4,atomi) =  &
       &  db4(1:3,-l4max_atomic:l4max_atomic,-l4max_magnetic:l4max_magnetic,1:nvar4)
       db4pack_0(-l4max_magnetic:l4max_magnetic,1:nvar4,atomi) = db4_0(-l4max_magnetic:l4max_magnetic,1:nvar4)
       
    enddo

    !--------------------------
    ! 5-body contributions

    if ((nradial5_atomic == 0).and.(nradial5_magnetic == 0)) then
       goto 1000
    endif
    
    printcount = 0.d0
    dprintcount = 0.1d0
    nmult = 0
       
    do atomi = 1, n2atom
       b5 = 0.d0
       db5 = (0.d0,0.d0)
       db5_0 = (0.d0,0.d0)
       do nclu = 1, nvar5
          ele1 = b5occ(1,nclu)
          ele2 = b5occ(2,nclu)
          ele3 = b5occ(3,nclu)
          ele4 = b5occ(4,nclu)          
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
          k_atomic = b5map_atomic(nclu)
          k_magnetic = b5map_magnetic(nclu)

          xiy = (0.d0,0.d0)
          
          do m_atomic = 1, madd5_atomic(k_atomic)
             m1 = b5lm_atomic(k_atomic,m_atomic,1)
             m2 = b5lm_atomic(k_atomic,m_atomic,2)
             m3 = b5lm_atomic(k_atomic,m_atomic,3)
             m4 = b5lm_atomic(k_atomic,m_atomic,4)
            
             k1 = shindex_atomic(l1,m1)
             k2 = shindex_atomic(l2,m2)
             k3 = shindex_atomic(l3,m3)
             k4 = shindex_atomic(l4,m4)
             
             cf_atomic = cg04_atomic(k_atomic,m1,m2,m3,m4)
             
             do m_magnetic = 1, madd5_magnetic(k_magnetic)
                m5 = b5lm_magnetic(k_magnetic,m_magnetic,1)
                m6 = b5lm_magnetic(k_magnetic,m_magnetic,2)
                m7 = b5lm_magnetic(k_magnetic,m_magnetic,3)
                m8 = b5lm_magnetic(k_magnetic,m_magnetic,4)
                m9 = b5lm_magnetic(k_magnetic,m_magnetic,5)
                
                k5 = shindex_magnetic(l5,m5)                
                k6 = shindex_magnetic(l6,m6)
                k7 = shindex_magnetic(l7,m7)
                k8 = shindex_magnetic(l8,m8)
                k9 = shindex_magnetic(l9,m9)
                
                ap0 = abase2pack_magnetic(nr5,k5,atomi)
                ap1 = abasepack(nr1,k1,nr6,k6,ele1,atomi)
                ap2 = abasepack(nr2,k2,nr7,k7,ele2,atomi)
                ap3 = abasepack(nr3,k3,nr8,k8,ele3,atomi)
                ap4 = abasepack(nr4,k4,nr9,k9,ele4,atomi)
                
                cf_magnetic = cg04_magnetic(k_magnetic,m5,m6,m7,m8,m9)
                cf = cf_atomic*cf_magnetic

                ! B
                xb1 = ap0*cf*ap1*ap2
                xb2 = ap3*ap4
                
                xiy = xiy + xb1*xb2
                
                ! gradients dB/dA
                db5(1,m1,m6,nclu) = db5(1,m1,m6,nclu) + ap0*cf*ap2*xb2
                db5(2,m2,m7,nclu) = db5(2,m2,m7,nclu) + ap0*cf*ap1*xb2
                db5(3,m3,m8,nclu) = db5(3,m3,m8,nclu) + xb1*ap4
                db5(4,m4,m9,nclu) = db5(4,m4,m9,nclu) + xb1*ap3
                
                ! gradient dB/dA_0
                db5_0(m5,nclu) = db5_0(m5,nclu) + cf*ap1*ap2*xb2
                
                !               nmult = nmult + 10 
             enddo
          enddo
          
             ! from here we are back to real numbers                   
             b5(nclu) = dreal( xiy )
       enddo
       b5pack(1:nvar5,atomi) = b5(1:nvar5)
       db5pack(1:4,-l5max_atomic:l5max_atomic,-l5max_magnetic:l5max_magnetic,1:nvar5,atomi) =  &
       &  db5(1:4,-l5max_atomic:l5max_atomic,-l5max_magnetic:l5max_magnetic,1:nvar5)
       db5pack_0(-l5max_magnetic:l5max_magnetic,1:nvar5,atomi) = db5_0(-l5max_magnetic:l5max_magnetic,1:nvar5)
    enddo

    !--------------------------
    ! 6-body contributions

    if ((nradial6_atomic == 0).and.(nradial6_magnetic == 0)) then
       goto 1000
    endif
    
    printcount = 0.d0
    dprintcount = 0.1d0
    nmult = 0
       
    do atomi = 1, n2atom
       b6 = 0.d0
       db6 = (0.d0,0.d0)
       db6_0 = (0.d0,0.d0)
       do nclu = 1, nvar6
          ele1 = b6occ(1,nclu)
          ele2 = b6occ(2,nclu)
          ele3 = b6occ(3,nclu)
          ele4 = b6occ(4,nclu)
          ele5 = b6occ(5,nclu)          
          nr1 = b6index(1,nclu)     ! n1
          nr2 = b6index(2,nclu)     ! n2
          nr3 = b6index(3,nclu)     ! n3
          nr4 = b6index(4,nclu)     ! n4
          nr5 = b6index(5,nclu)     ! n5
          nr6 = b6index(6,nclu)     ! n0'
          nr7 = b6index(7,nclu)     ! n1'
          nr8 = b6index(8,nclu)     ! n2'
          nr9 = b6index(9,nclu)     ! n3'
          nr10 = b6index(10,nclu)   ! n4'
          nr11 = b6index(11,nclu)   ! n5'          
          l1 = b6index(12,nclu)     ! l1
          l2 = b6index(13,nclu)     ! l2
          l3 = b6index(14,nclu)     ! l3
          l4 = b6index(15,nclu)     ! l4
          l5 = b6index(16,nclu)     ! l5
          l6 = b6index(17,nclu)     ! l0'
          l7 = b6index(18,nclu)     ! l1'
          l8 = b6index(19,nclu)     ! l2'
          l9 = b6index(20,nclu)     ! l3'
          l10 = b6index(21,nclu)    ! l4'
          l11 = b6index(22,nclu)    ! l5'          
          k_atomic = b6map_atomic(nclu)
          k_magnetic = b6map_magnetic(nclu)
          
          xiy = (0.d0,0.d0)
          
          do m_atomic = 1, madd6_atomic(k_atomic)
             m1 = b6lm_atomic(k_atomic,m_atomic,1)
             m2 = b6lm_atomic(k_atomic,m_atomic,2)
             m3 = b6lm_atomic(k_atomic,m_atomic,3)
             m4 = b6lm_atomic(k_atomic,m_atomic,4)
             m5 = b6lm_atomic(k_atomic,m_atomic,5)
             
             k1 = shindex_atomic(l1,m1)
             k2 = shindex_atomic(l2,m2)
             k3 = shindex_atomic(l3,m3)
             k4 = shindex_atomic(l4,m4)
             k5 = shindex_atomic(l5,m5)
             
             cf_atomic = cg05_atomic(k_atomic,m1,m2,m3,m4,m5)
          
             do m_magnetic = 1, madd6_magnetic(k_magnetic)
             
                m6 = b6lm_magnetic(k_magnetic,m_magnetic,1)
                m7 = b6lm_magnetic(k_magnetic,m_magnetic,2)
                m8 = b6lm_magnetic(k_magnetic,m_magnetic,3)
                m9 = b6lm_magnetic(k_magnetic,m_magnetic,4)
                m10 = b6lm_magnetic(k_magnetic,m_magnetic,5)
                m11 = b6lm_magnetic(k_magnetic,m_magnetic,6)
                
                k6 = shindex_magnetic(l6,m6)                
                k7 = shindex_magnetic(l7,m7)
                k8 = shindex_magnetic(l8,m8)
                k9 = shindex_magnetic(l9,m9)
                k10 = shindex_magnetic(l10,m10)
                k11 = shindex_magnetic(l11,m11)
                
                
                ap0 = abase2pack_magnetic(nr6,k6,atomi)
                ap1 = abasepack(nr1,k1,nr7,k7,ele1,atomi)
                ap2 = abasepack(nr2,k2,nr8,k8,ele2,atomi)
                ap3 = abasepack(nr3,k3,nr9,k9,ele3,atomi)
                ap4 = abasepack(nr4,k4,nr10,k10,ele4,atomi)
                ap5 = abasepack(nr5,k5,nr11,k11,ele5,atomi)
                
                cf_magnetic = cg05_magnetic(k_magnetic,m6,m7,m8,m9,m10,m11)
                cf = cf_atomic*cf_magnetic
                
                ! B
                xb1 = ap0*cf*ap1*ap2
                xb2 = ap3*ap4*ap5
                
                xiy = xiy + xb1*xb2
                
                ! gradients dB/dA
                db6(1,m1,m7,nclu) = db6(1,m1,m7,nclu) + ap0*cf*ap2*xb2
                db6(2,m2,m8,nclu) = db6(2,m2,m8,nclu) + ap0*cf*ap1*xb2
                db6(3,m3,m9,nclu) = db6(3,m3,m9,nclu) + xb1*ap4*ap5
                db6(4,m4,m10,nclu) = db6(4,m4,m10,nclu) + xb1*ap3*ap5
                db6(5,m5,m11,nclu) = db6(5,m5,m11,nclu) + xb1*ap3*ap4             
                !               nmult = nmult + 15    
                
                ! gradient dB/dA_0
                db6_0(m6,nclu) = db6_0(m6,nclu) + cf*ap1*ap2*xb2
                
             enddo   
          enddo
          
          b6(nclu) = dreal( xiy )
       enddo
       b6pack(1:nvar6,atomi) = b6(1:nvar6)
       db6pack(1:5,-l6max_atomic:l6max_atomic,-l6max_magnetic:l6max_magnetic,1:nvar6,atomi) =  &
       &   db6(1:5,-l6max_atomic:l6max_atomic,-l6max_magnetic:l6max_magnetic,1:nvar6)
       db6pack_0(-l6max_magnetic:l6max_magnetic,1:nvar6,atomi) = db6_0(-l6max_magnetic:l6max_magnetic,1:nvar6)
       
    enddo
       
    1000 continue
    
end subroutine prepB

!-------------------------------------------------------------------------------

  
