	module execute
!----------------------------------------------------------
	integer conti
        character, parameter::readx='slab.xyz'
        end module
!-------------------------------------------------------
	subroutine readX
!-------------------------------------------------------
        use system
        use surface
        use constants
        use execute
        integer i,k
        character*3 
!
        open(oxyz,file = readx)  
        read(oxyz,*) natom
        read(oxyz,*)
        allocate(x(3,natom))
        allocate(nl(natom)) 
        allocate(id(natom))
        do i=1,natom
            read(oxyz,format=100) id(i),(x(k,i),k=1,3),      
