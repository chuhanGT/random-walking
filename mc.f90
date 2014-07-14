!-------------------------------------------------------
      function distance(box,x,y)! z-direction is not periodic
!-------------------------------------------------------
      real*8 x(3),y(3),d(3),r,distance,box(3)
      integer i
!
      r=0.d0
      do i=1,2
         d(i)=abs(x(i)-y(i))
         d(i)=mod(d(i),box(i))
         d(i)=d(i)-box(i)*int(d(i)*2.d0/box(i))
         r=r+d(i)**2
      enddo
      d(3)=abs(x(3)-y(3))
      r=r+d(3)**2
      distance=sqrt(r)
      return 
      end
!--------------------------------------------------------
! Every atom's neighbours has noxy capacity, store the atom
! number for those within the cutoff, all others store as 0
!--------------------------------------------------------
      subroutine neighbours(iat)
!--------------------------------------------------------
      use constants
      use system
      use surface
      use sf_mc
!
      integer i,n,iat
      real*8 origin(3),r,distance
!
      if (id(iat) .ne. ' O') then
         goto 100
      endif
!
      deallocate(nb)
      allocate(nb(noxy))
!
      n=0
      do i=1,noxy
          r=distance(box,x(1,iat),x(1,oxy(i)))
          if (r .lt. rcut/aul .and. r .gt. 0.d0) then 
               n=n+1
               nb(n)=oxy(i)
          endif
      enddo
      maxnbrs=n
!
      do i=n+1,noxy
          nb(i) = 0 
      enddo
!
100      end
!--------------------------------------------------------
      subroutine writeneighbours(iat)
!--------------------------------------------------------
      use constants
      use system
      use sf_mc
      integer iat,i,j
!
      maxnbrs=0
      call neighbours(iat)
      if (maxnbrs .gt. 0) then
      write(onb,*) maxnbrs+1
      write(onb,*) maxnbrs+1    
      write(onb,"(a2,3x,3(1x,f10.5))") 'N ',(x(i,iat)*aul,i=1,3)
      do i=1,maxnbrs
         write(onb,"(a2,3x,3(1x,f10.5))")'O ', (x(j,nb(i))*aul,j=1,3)
      enddo
      endif
!
      end      
!----------------------------------------------------------
       subroutine allneighbours
!----------------------------------------------------------
       use system
       integer i
!
       do i=1,natom
          call writeneighbours(i)
       enddo
       end
!---------------------------------------------------------
! isvac is an array of natom elements
!---------------------------------------------------------
      subroutine init_vacancy
!---------------------------------------------------------
      use system
      use surface
      use sf_mc
!
      integer i,j
!
      allocate(nb(noxy))
      allocate(isvac(natom))
      do i=1,natom
           isvac(i) = 0 
      enddo
!
      do i=1,nvac
          do j=1,natom 
             if (j .eq. vac(i)) then
                 isvac(j) = 1
             endif
          enddo
       enddo   
!     
       end
!--------------------------------------------------------
      function SlabEn(vaca,n,iszro) 
!--------------------------------------------------------
       use surface
       integer n,i,layer
       integer, dimension(n):: vaca
       real*8 SlabEn, e
       logical iszro
!
       SlabEn=0.d0
       do i=1,n
          if (vaca(i) .eq. 1) then
              layer=il(i)
            if (layer .gt. laylmt) then
                 e=0.d0
            else
              if (iszro) then
                    e=zr(layer)
              else
                    e=ba(layer)
              endif
             endif
                 SlabEn=SlabEn+e
           endif
        enddo
        return 
        end
!-------------------------------------------------------
        function LayerEn(iat,iszro)
!-------------------------------------------------------
        use system
        use surface
        integer iat, layer
        logical iszro
        real*8 LayerEn
!
        layer=il(iat)
        if (layer .gt. laylmt) then
            LayerEn=0.d0
        else
           if (iszro) then
              layerEn=zr(layer)
           else
              layerEn=ba(layer)
           endif
        endif
!
        return
        end
!-------------------------------------------------------
        subroutine move
!-------------------------------------------------------
        use constants
        use system
        use surface
        use sf_mc
!
        integer i,j,k,myoxy,n,choose,ivac
        real*8 Ei,Ef,deltaE,beta,ex,r,LayerEn,SlabEn
        real*8, dimension(natom):: temp
!
           beta=-aue/(kB*T)
           temp=isvac
!
           call random_number(r)
           ivac=int(nvac*r)+1
           ivac=max(1,ivac)
           ivac=min(ivac,nvac)
           call neighbours(vac(ivac)) 
           call random_number(r)
           k=int(maxnbrs*r)+1
           k=max(1,k)
           k=min(k,maxnbrs)
           choose=nb(k)
           if (isvac(choose) .eq. 0) then
                     isvac(choose) = 1
                     isvac(vac(ivac)) = 0
                     Ei=LayerEn(vac(ivac),zro)
                     Ef=LayerEn(choose,zro)
           else 
!              write(*,*) 'not move'
           endif
!
          deltaE=Ef-Ei
          if (deltaE .ge. zero) then
              ex=exp(beta*deltaE)
              call random_number(r)
              r=max(zero, r)
              r=min(r, one)
              if (ex .lt. r) then
                  isvac=temp
              endif
          endif
!
          k=0
          do i=1,natom
                if (isvac(i) .eq.1) then
                   k=k+1
                   vac(k)=i
                endif
           enddo
!
           if  (k .ne. nvac) then
                write(*,*) 'stuck!',k,nvac 
                stop 
           endif
!
       end
!---------------------------------------------------------
       subroutine mcstep
!---------------------------------------------------------
       use constants
       use system
       use surface
       use sf_mc
       integer i
       real*8 SlabEn
!
       do i=1,nvac
          call move
       enddo
!
       write(oen,"(f10.5)") SlabEn(isvac,natom,zro)*aue
       call check_equil
       end
!----------------------------------------------------------
      subroutine outxyz
!---------------------------------------------------------
      use constants
      use surface
      use system
      use sf_mc
!
      integer i,k
      character*2, dimension(natom):: names
      do i=1,natom
         if (isvac(i) .eq. 1) then
             names(i)=' N'
         else
             names(i)=id(i)
         endif
      enddo
!
      write(omc,*) nvac
      write(omc,*) nvac
      do i=1,natom
         if ( names(i) .eq. ' N') then
         write(omc,"(a2,3x,3(3x,f10.5))") names(i),(x(k,i)*aul,k=1,3)
         endif
      enddo
      end
!--------------------------------------------------------------
!count number of vacancies on each layer, compare with prob to 
!derermine whether reach equilibrium
!--------------------------------------------------------------
       subroutine reach_equil
!--------------------------------------------------------------
       use system
       use surface
       use sf_mc
       integer i,k,bin(nl+1),l
       real*8 ratio(nl+1)
!count vacancy numbers of each layer
       do i=1,nl+1
          bin(i)=0
       enddo
       do i=1,nvac
          k=vac(i)
          l=il(k)
          bin(l) = bin(l)+1
       enddo
!
       do i=1,nl+1
          ratio(i) = dble(bin(i))/dble(nvac)
       enddo
       write(14,"(f10.5)") (ratio(i),i=1,nl+1) 
       write(14,*) '----------------------------'
!now compare the real partition with the prob
       do i=1,nl+1
          if (ratio(i) .ge. prob(i)-offset .and.ratio(i) .le. prob(i)+offset) then
               isequil=.true.
          else
               isequil = .false.
               goto 200
          endif
        enddo
200     end
!-----------------------------------------------------------------
! Suppose chech_equil returning true last certain amount of steps
! then we accept it reach equilibrium
         subroutine check_equil
!-----------------------------------------------------------------
         use system
         use surface
         use sf_mc
         integer i,j,k    
!
         call reach_equil
         if (isequil) then
               eqcounter=eqcounter+1
         else
               eqcounter=0
         endif
!
        if (eqcounter .eq. eqstep) then
             eqcounter = 0
             write(*,*) 'reach equilibrium! Temperature is ', T,' and ',runstep,' th step'
        else
        endif
!
        end
!-----------------------------------------------------------------
! one runstep depends one how many mcsteps to print output     
        subroutine run_mc_step
!-----------------------------------------------------------------
        use sf_mc
        integer j
!
         j=0
         do runstep=1,totstep
              j=j+1
              call mcstep
              if (j.eq.nstep) then
                     j=0
                     call outxyz
              endif
         enddo
!
         end
!-----------------------------------------------------------------
