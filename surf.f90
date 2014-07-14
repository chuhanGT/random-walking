!--------------------------------------------------------------
      module constants
!--------------------------------------------------------------
        real*8, parameter :: zero=0.d0,half=0.5d0,one=1.d0
        real*8, parameter :: two=2.d0,three=3.d0,four=4.d0
        real*8, parameter :: kB=8.6173324e-5
!
        integer, parameter :: oxyz=3, ovac=8,omc=7,onb=9,oen=10
!
      end module constants
!---------------------------------------------------------------
      module system
!---------------------------------------------------------------
        integer nx,ny,nz
        integer natom
        real*8, parameter :: a0=4.28d0
        real*8, parameter :: aue = 27.2117d0
        real*8, parameter :: aul = 0.5291772083d0
        real*8  T      
        real*8 box(3)
        character*2, allocatable :: id(:)
        real*8,  allocatable :: x(:,:)
        integer switch

!
      end module system
!---------------------------------------------------------------
      module surface
!---------------------------------------------------------------
        use system,only:aue,aul
        integer nl,noxy,nvac
        integer, allocatable :: il(:), vac(:),oxy(:)
        real*8, allocatable :: lyen(:),prob(:)
        real*8 xOv
        logical zro
        integer, parameter :: laylmt=6
        real*8, parameter,dimension(6) :: ba=[-0.44/aue,0.02/aue, -0.33/aue,0.03/aue,-0.29/aue,0.00/aue]
        real*8, parameter,dimension(6) :: zr=[-0.39/aue,0.66/aue,0.16/aue,0.00/aue,0.d0,0.d0]
!
      end module surface
!-----------------------------------------------------------------------
      module sf_mc
!------------------------------------------------------------------
      use surface, only: noxy
      real*8 rcut,offset !offset is less than 1  
      integer, allocatable :: isvac(:),nb(:)
      integer maxnbrs,eqcounter,runstep
      integer  totstep,nstep,eqstep
      logical isequil
      end module
!-----------------------------------------------------------------------
      program mc___code
!-----------------------------------------------------------------------
      use system
      use surface
      use sf_mc
      implicit none
      integer i,j
      real*8 SlabEn
!
      read(*,*) nx,ny,nl,zro,xOv,rcut
      read(*,*) totstep,nstep,eqstep
      read(*,*) T
      read(*,*) switch
      read(*,*) offset
      call random_seed()
      call buildSlab
      call vacancy
      if (switch .eq. 1) then
         call rand_vacancy
      else
         call man_vacancy
      endif
      call sf_energy
      call init_vacancy
      call allneighbours
!      do i=T,300,-200
!      T=i
      call lyr_prob
      call run_mc_step
!      enddo
!
      stop
      end
!-----------------------------------------------------------------------
      subroutine sf_energy
!-----------------------------------------------------------------------
      use constants
      use system
      use surface
      implicit none
      integer i,j,k,num
      real*8 Etot,le
!
      write(*,*) 'Vacancy coordinates:'
      write(*,*) '------------------------------------------------'
      Etot=0.d0
      do i=1,nvac
         num=vac(i)
         write(*,"(i6,3x,3(3x,f10.5),i5)") vac(i),(x(k,vac(i))*aul,k=1,3),il(num)
         if (il(num) .gt. laylmt) then
            le=0.d0
         else if (zro) then
             le=zr(il(num))
         else
             le=ba(il(num))  
         endif
         Etot=Etot+le
      enddo
!
      write(*,*) '----------------------------------'
      write(*,"(a20,i5)") 'Total vacancies:  ',nvac
      write(*,*) '----------------------------------'
      write(*,"(a30,f5.2)") 'Initial total Energy (eV):    ',Etot*aue
      write(*,*) '----------------------------------'
!    
      return
      end
!-------------------------------------------------------------------------
      subroutine buildSlab
!-------------------------------------------------------------------------
      use constants
      use system
      use surface
      implicit none
      real*8 basis(3,5),z
      character*2 isp(5)
      integer i,j,k,num,iu,layerNum
     
!
      if (zro) then
         natom=3*nx*ny*(nl+1)+2*nx*ny*nl
         z=a0*nl/aul
      else
         natom=3*nx*ny*nl+2*nx*ny*(nl+1)
         z=-a0/(two*aul)
      endif
!
      allocate(x(3,natom))
      allocate(id(natom))
      allocate(il(natom))
!
      box(1)=a0*nx/aul
      box(2)=a0*ny/aul
      box(3)=a0*nl/aul
!      
      isp(1)='Ba'
      do i=1,3
          basis(i,1)=a0/(two*aul)
      enddo
      isp(2)='Zr'
      do i=1,3
         basis(i,2)=zero
      enddo
      isp(3)=' O'
      basis(1,3)=a0/(two*aul)
      basis(2,3)=zero
      basis(3,3)=zero
      isp(4)=' O'
      basis(1,4)=zero
      basis(2,4)=a0/(two*aul)
      basis(3,4)=zero
      isp(5)=' O'
      basis(1,5)=zero
      basis(2,5)=zero
      basis(3,5)=a0/(two*aul)
!
      num=0
      do iu=1,5
        do i=0,nx-1
           do j=0,ny-1
            do k=0,nl-1
               num=num+1
               x(1,num)=basis(1,iu)+a0*i/aul
               x(2,num)=basis(2,iu)+a0*j/aul
               x(3,num)=basis(3,iu)+a0*k/aul
               id(num)=isp(iu)
               il(num)=layerNum(x(3,num))
            enddo
         enddo
       enddo
      enddo
!
      if (zro) then
       do iu=2,4
          do i=0,nx-1
             do j=0,ny-1
                num=num+1
                x(1,num)=basis(1,iu)+a0*i/aul
                x(2,num)=basis(2,iu)+a0*j/aul
                x(3,num)=z
                id(num)=isp(iu)
                il(num)=1
              enddo
           enddo
        enddo
      else
        do i=0,nx-1
           do j=0,ny-1
              num=num+1
              x(1,num)=basis(1,1)+a0*i/aul
              x(2,num)=basis(2,1)+a0*j/aul
              x(3,num)=z
              id(num)=isp(1)
              il(num)=1
              num=num+1
              x(1,num)=basis(1,5)+a0*i/aul
              x(2,num)=basis(2,5)+a0*j/aul
              x(3,num)=z
              id(num)=isp(5)
              il(num)=1
           enddo
        enddo
      endif


      if (num .ne. natom) then
         write(*,*) 'Error! Number of atom is not correct'
         write(*,*)num, natom
         stop
      endif
!
      write(oxyz,*) natom
      write(oxyz,*) natom
      do i=1,natom
         write(oxyz,"(a2,3x,3(1x,f12.5),a3,i5)") id(i),(x(k,i)*aul,k=1,3),' # ',il(i)
      enddo
!
      call lyr_energy
!
      return 
      end
!--------------------------------------------------------------
      function layerNum(z)
!--------------------------------------------------------------
      use system
      use surface
      use constants
      real*8 z,halfa0
      integer layerNum,start,height
!
      if (zro) then 
         start = 1
     else
         start=2
      endif
!
      halfa0=half*a0/aul
      height=z/halfa0
      height=height+start
      if (height .le. nl+1) then
          layerNum=height
      else
          layerNum=2*nl+2-height
      endif
!
      return 
      end
!-------------------------------------------------------------------
      subroutine vacancy
!-------------------------------------------------------------------
     use system 
     use surface
!
      integer i,j,k,c,l
!
      if (zro) then
         noxy=3*nx*ny*nl+2*nx*ny
      else
         noxy=3*nx*ny*nl+nx*ny
      endif
!
      allocate(oxy(noxy))
!
      nvac=int(xOv*noxy)
      allocate(vac(nvac))
!
      c=0
      j=0
      do i=1,natom
       c=c+1
       if (id(i) .eq. ' O') then
          j=j+1 
          oxy(j)=c
       endif
      enddo
!
      if (j .ne. noxy) then
        write(*,*) 'Number of oxygens invalid'
        write(*,*) j, noxy
        stop
      endif
!
      end
!------------------------------------------------------------
        subroutine rand_vacancy
!------------------------------------------------------------
         use constants
         use system
          use surface
!
         integer i,j,k,l
         real*8 r
         logical isvac
!
     call random_number(r)
      k=int(noxy*r)+1
      k=max(1,k)
      k=min(k,noxy)
      vac(1)=oxy(k)
!
      i=2
      do while (i .le. nvac)
       call random_number(r)
       k=int(noxy*r)+1
       k=max(1,k)
       k=min(k,noxy)
       l=0
       do j=1,i-1
          if (oxy(k) .ne. vac(j)) then
                 l=l+1
          endif
       enddo
       if ( i-1 .eq.l ) then 
            vac(i)=oxy(k)
            i=i+1
       endif
      enddo
!   
     write(ovac,*) noxy
     write(ovac,*) 'layer #, atom index'
      do i =1,noxy
      isvac=.false.
      do j=1,nvac
        if (oxy(i).eq.vac(j)) then
           isvac=.true.
        endif
       enddo
       if (isvac) then
          write(ovac,"(a2,3(3X,f10.5),a3,2(3x,i5))") ' N',(x(k,oxy(i))*aul,k=1,3),' # ',il(oxy(i)),oxy(i)
         else
           write(ovac,"(a2,3(3X,f10.5),a3,2(3x,i5))") id(oxy(i)),(x(k,oxy(i))*aul,k=1,3),' # ',il(oxy(i)),oxy(i)
       endif

      enddo
      end               
!--------------------------------------------------------------------
      subroutine man_vacancy !let as more vac. on first layer as possible
!--------------------------------------------------------------------
      use constants
      use system
      use surface
      integer i,j,k
      logical vv
!
      k=0
      do i=1,natom
           if(id(i) .eq. ' O' .and. il(i) .eq. 1) then
                 k=k+1
                 vac(k)=i
                 if (k .eq. nvac) goto 100
           endif
      enddo
100   continue
      if (k .lt. nvac) then
                do i=1,natom
                    if (id(i) .eq. ' O' .and. il(i) .ne. 1) then
                       k=k+1
                       if (k .eq. nvac+1) goto 200 
                       vac(k)=i
                     endif
                enddo
       endif
200    continue
!
       write(ovac,*) noxy
       write(ovac,*) noxy
       do i =1,noxy
           vv=.false.
           do j=1,nvac
              if (oxy(i).eq.vac(j)) then
                 vv=.true.
              endif
           enddo
           if (vv) then
                write(ovac,"(a2,3(3X,f10.5),a3,2(3x,i5))") ' N',(x(k,oxy(i))*aul,k=1,3),' # ',il(oxy(i)),oxy(i)
           else
                write(ovac,"(a2,3(3X,f10.5),a3,2(3x,i5))") id(oxy(i)),(x(k,oxy(i))*aul,k=1,3),' # ',il(oxy(i)),oxy(i)
           endif
      enddo
!
      end
!----------------------------------------------------------------------------
!we need to know each layer's energy
!----------------------------------------------------------------------------
      subroutine lyr_energy
!----------------------------------------------------------------------------
      use surface
      integer i
!
      allocate(lyen(nl+1))
      allocate(prob(nl+1))
!
      do i=1,nl+1
         if (i .le. laylmt) then
             if (zro) then
                lyen(i) = zr(i)
             else
                lyen(i) = ba(i)
             endif
         else
            lyen(i) = 0.d0
         endif
      enddo
! 
      write(*,*) '----------------------------------'
      write(*,*) 'Layer #      Energy(eV)'
      write(*,*) '----------------------------------'
      do i=1,nl+1
      write(*,"(i5,3x,f10.5)") i,lyen(i)*aue
      enddo
      write(*,*) '----------------------------------'
      end
!------------------------------------------------------------------------
! we need the probability for vacancy staying on each layer
!------------------------------------------------------------------------
      subroutine lyr_prob
!------------------------------------------------------------------------
      use constants, only: kB
      use system , only: T
      use surface
      integer i,n
      real*8 z,beta,bolt(nl+1),p,e
!
      beta=-aue/(kB*T)
!
      do i=1,nl+1
         bolt(i)=exp(beta*lyen(i))
      enddo
!
      z=0.d0
      do i=1,nl
         z=z+2*bolt(i)
      enddo
      z=z+bolt(nl+1)
!
      do i=1,nl
        prob(i) = 2*bolt(i)/z
      enddo
      prob(nl+1)=bolt(nl+1)/z
!
      write(*,*) 'Layer #     Prob.     Particle #'
      write(*,*) '----------------------------------------------------'
      n=0
      p=0.d0
      e=0.d0
      do i=1,nl+1
      n=n+int(prob(i)*nvac+0.5d0)
      p=p+prob(i)
      e=e+lyen(i)*aue*int(prob(i)*nvac+0.5d0)
      write(*,"(i5,3x,f10.5,3x,i5)") i,prob(i),int(prob(i)*nvac+0.5d0)
      enddo
      write(*,*) '-----------------------------------------------------'
      write(*,"(a7,3x,f10.5,3x,i5)") 'Total: ',p,n
      write(*,*) '-----------------------------------------------------'
      write(*,"(a40,f10.5)") 'Predicted equilibrium total energy(eV): ',e
      write(*,*) '-----------------------------------------------------'
      

      end 
