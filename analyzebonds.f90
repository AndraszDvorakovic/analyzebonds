program AnalyzeBonds
implicit none
real,allocatable,dimension(:,:) :: array
integer :: Nlines,i,j,Maxsteps,step,Natoms
character :: junk
real :: time,timestep,attofs
attofs=0.0241888425
Natoms=15
timestep=4
Maxsteps=6000

allocate(array(Natoms*Maxsteps,3))

open(10,file='movie.xyz')
open(20,file='movieanalyzed.xyz',status='replace')

do i=1,Maxsteps
read(10,*,end=100) 
read(10,*,end=100) junk,junk,step
write(20,*) step
time=timestep*step*attofs
do j=1,Natoms
read(10,*) junk,array(j*step,1),array(j*step,2),array(j*step,3) 
write(20,*) array(j*step,1),array(j*step,2),array(j*step,3),time,'fs'
end do
end do


100 close(10)
close(20)
end program AnalyzeBonds

