	program test
        implicit none
        integer natom,i,k
        real*8, allocatable:: x(:,:)
        integer, allocatable::nl(:)
        character*3 aa
        character*2 bb
!     
        open(unit=100,file='slab.xyz') 
        read(100,*) natom
        read(100,*) 
        allocate(x(3,natom))
        allocate(nl(natom))
        i=1
        read(100,"(a2,3x,3(1x,f12.5),a3,i5)") bb,(x(k,i),k=1,3),aa,nl(i) 
        write(*,*) natom
        write(*,*) i,(x(k,i),k=1,3),nl(i)
       end
