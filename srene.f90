!
      module constants
!
        real*8, parameter :: zero=0.d0,half=0.5d0,one=1.d0
        real*8, parameter :: two=2.d0,three=3.d0,four=4.d0
!
      end module constants
!
      module system
!
        integer nx,ny,nz
        integer natom
        real*8, parameter :: a0=4.28d0
        real*8, parameter :: aue = 27.2117d0
        real*8, parameter :: aul = 0.5291772083d0
        real*8 box(3)
        character*2, allocatable :: id(:)
        real*8,  allocatable :: x(:,:)
!
      end module system
!
      module sr_parameters
!     
        use system,only:aue,aul
	real*8, parameter :: Aab=1000/aue, Aao=2000/aue,Abo=3000/aue
	real*8, parameter :: bab=aul/0.5d0, bao=aul/0.3d0, bbo=aul/0.2d0
        real*8, parameter :: cutoff=3.71d0/aul
!   a=Ba    b=Zr    o=O
      end module sr_parameters 
!-----------------------------------------------------------------------
      program mc___code
!-----------------------------------------------------------------------
      implicit none
!
!

      call sr_energy
!
!
      stop
      end
!-----------------------------------------------------------------------
      subroutine sr_energy
!-----------------------------------------------------------------------
      use constants
      use system
      use sr_parameters
      implicit none
      integer ia,ja,k
      real*8 d(3),sru,Etot,u
!
      read(*,*) nx,ny,nz
      call buildSystem
!
      Etot=0.d0
      do ia=1,natom
        do ja=ia+1,natom
!
               u=sru(box,x(1,ia),x(1,ja),d,id(ia),id(ja))
               write(*,*) ia,id(ia),ja,id(ja),u
               Etot=Etot+u
!
        enddo
      enddo
!
      write(*,*) Etot
      return
      end
!-----------------------------------------------------------------------
      function pbcr(box,x,y,d)
!-----------------------------------------------------------------------
      implicit none
      real*8 , parameter :: zero=0.d0,two=2.d0
      integer k
      real*8 box(3),x(3),y(3),d(3),pbcr,r
!
      r=zero
      do k=1,3
        d(k)=x(k)-y(k)
        d(k)=mod(d(k),box(k))
        d(k)=d(k)-box(k)*int(d(k)*two/box(k))
        r=r+d(k)*d(k)
      enddo
      pbcr=sqrt(r)
      return
      end
!-----------------------------------------------------------------------
      function sru(box,x,y,d,isp,jsp)
!-----------------------------------------------------------------------
      use sr_parameters
      implicit none
      real*8 x(3),y(3),box(3),d(3),r,sru,pbcr
      character*2 isp,jsp
!
      r=pbcr(box,x,y,d)
!
      if (r .lt. cutoff) then
!
      if (isp .eq. jsp) then
         sru=0.d0
      else if ((isp .eq. 'Ba' .and. jsp .eq. 'Zr') .or.(isp .eq. 'Zr' .and. jsp .eq. 'Ba')) then
         sru=Aab*exp(-bab*r)
      else if ((isp .eq. 'Ba' .and. jsp .eq. ' O') .or.(isp .eq. ' O' .and. jsp .eq. 'Ba')) then
         sru=Aao*exp(-bao*r)
      else if ((isp .eq. 'Zr' .and. jsp .eq. ' O') .or.(isp .eq. ' O' .and. jsp .eq. 'Zr')) then
         sru=Abo*exp(-bbo*r)
      endif
!
      else
        sru=0.d0
      endif
!
      return
      end

!-------------------------------------------------------------------------
      subroutine buildSystem
!-------------------------------------------------------------------------
      use constants
      use system
      implicit none
      real*8 basis(3,5)
      character*2 isp(5)
      integer i,j,k,num,iu
!
      natom = 5*nx*ny*nz
      allocate(x(3,natom))
      allocate(id(natom))
!
      box(1)=a0*nx/aul
      box(2)=a0*ny/aul
      box(3)=a0*nz/aul
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
            do k=0,nz-1
               num=num+1
               x(1,num)=basis(1,iu)+a0*i/aul
               x(2,num)=basis(2,iu)+a0*j/aul
               x(3,num)=basis(3,iu)+a0*k/aul
               id(num)=isp(iu)
            enddo
         enddo
       enddo
     enddo
!
     if (num .ne. natom) then
         write(*,*) 'Error! Number of atom is not correct'
         write(*,*)num, natom
         stop
     endif
!
     return 
     end
               
