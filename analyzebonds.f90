program AnalyzeBonds
implicit none
real,allocatable,dimension(:,:) :: array,bonds,arrayhelp
integer,allocatable,dimension(:,:) :: existbonds
integer :: i,j,k,m,n,o,Maxsteps,step,Natoms,Nbonds,choice, Maxbonds,Ntraj
character*30 :: junk, mainfolder, moviesource, trajfolder, xn, form, outputfile1,outputfile2, finalfolder
real :: time,timestep,attofs
attofs=0.0241888425

Natoms=15                            
timestep=4			    !timestep of the simulation
Maxsteps=6000                       !max number of geometries in movie.xyz (you can use higher number)
Maxbonds=14                         !maximum of real bonds in molecule (matters only if you want to scan more distances)
Ntraj=100                           !number of trajectories to scan
mainfolder='ION1TRAJ_0ST'           !there are located folders with TRAJ.n folders
trajfolder='/TRAJ.'

!----------------------------------------------------------------------------------------------
!Which bonds exist? Which bonds we want?
open(30,file='existbonds.txt',status='replace')
allocate(existbonds(Natoms,2))

print *, 'How many bonds? (All bonds = 100 Personal option)'
read *, choice

if (choice == 100) then				!This option works only 
Nbonds=14					!with predefined bonds 						
existbonds(1,1)=1				!in this array
existbonds(1,2)=2				!every first row means the first atom
existbonds(2,1)=1				!every second row means the second atom
existbonds(2,2)=3					
existbonds(3,1)=1
existbonds(3,2)=4
existbonds(4,1)=1
existbonds(4,2)=5
existbonds(5,1)=5
existbonds(5,2)=6
existbonds(6,1)=2
existbonds(6,2)=10
existbonds(7,1)=2
existbonds(7,2)=11
existbonds(8,1)=2
existbonds(8,2)=12
existbonds(9,1)=4
existbonds(9,2)=13
existbonds(10,1)=4
existbonds(10,2)=14
existbonds(11,1)=4
existbonds(11,2)=15
existbonds(12,1)=3
existbonds(12,2)=7
existbonds(13,1)=3
existbonds(13,2)=8
existbonds(14,1)=3
existbonds(14,2)=9
	do k=1,Nbonds
		write(30,*) existbonds(k,1),existbonds(k,2)
	end do
end if

if (choice > Maxbonds .AND. choice /= 100) then									
  print *, 'Too many bonds for this molecule'
  stop
endif

if (choice /= 100) then								!This works for any 
print *, 'Type pairs of atom numbers to scan bond lenght between them.'		!distance between any atoms
print *, 'The order of atoms is the same as in movie.xyz' 			!
print *, 'Separate pairs by enter.'   
  Nbonds=choice
  	do o=1,Nbonds
  		read *, existbonds(o,1),existbonds(o,2)
        if (existbonds(o,1) > Natoms .OR. existbonds(o,2) > Natoms) then
          print *, 'There are only', Natoms, ' atoms'
          stop
        end if
  		write(30,*) existbonds(o,1),existbonds(o,2)
  		print *, Nbonds-o,' left.'
  	end do
end if

!-----------------------------------------------------------------------------------------------------
!Analyzing bond lenghts
!-----------------------------------------------------------------------------------------------------
allocate(array(Natoms,3))
allocate(arrayhelp(1,3))
allocate(bonds(Maxsteps,Nbonds+1))

finalfolder=TRIM(mainfolder)//'analyzed'
call system('mkdir '//TRIM(finalfolder))

do n=1,Ntraj
	if (n < 10) then
		form = '(I1)'
	end if
	if (n >= 10 .AND. n <100) then
		form = '(I2)'
	end if
	if (n >= 100 .AND. n <1000) then			!					
		form = '(I3)'					!
	end if							!	
	write (xn,form) n					!
	moviesource=TRIM(mainfolder)//TRIM(trajfolder)//xn      !Here I am preparing names 
	outputfile1=xn						!of input and output files 
	outputfile2=xn					        !to change with n - number of trajectory

	open(10,file=TRIM(moviesource)//'/movie.xyz',status='old')                         
	open(20,file=TRIM(outputfile1)//'_bonds1to7.txt',status='replace')                 

	if (Nbonds >= 8) then                                                   !In case of Nbonds > 14 it is necessary  
  	  open(21,file=TRIM(outputfile2)//'_bonds8to14.txt',status='replace')   !to create another outputfile.
          write(21,*) '#Time (fs)', '        #Bond lenght (A)'		        !There are max 7 scans per outputfile
	endif
	write(20,*) '#Time (fs)', '        #Bond lenght (A)'

	do i=1,Maxsteps						!Here I am getting data from movie
		read(10,*,end=100)                              !counting bond lenghts and writing them
    		read(10,*,end=100) junk,junk,step
    	    	time=timestep*step*attofs
    
    		do j=1,Natoms
	       		read(10,*) junk,array(j,1),array(j,2),array(j,3)
        	end do
        
        	do m=1,Nbonds
      			arrayhelp(1,1)=(array(existbonds(m,1),1)-array(existbonds(m,2),1))**2.0
   				arrayhelp(1,2)=(array(existbonds(m,1),2)-array(existbonds(m,2),2))**2.0
           		arrayhelp(1,3)=(array(existbonds(m,1),3)-array(existbonds(m,2),3))**2.0
       			bonds(i,1)=time
   				bonds(i,m+1)=sqrt(arrayhelp(1,1)+arrayhelp(1,2)+arrayhelp(1,3))
       		end do
        
write(20,*) bonds(i,1), bonds(i,2), bonds(i,3), bonds(i,4), bonds(i,5), bonds(i,6), bonds(i,7), bonds(i,8)
		if (Nbonds >= 8) then
write(21,*) bonds(i,1), bonds(i,9), bonds(i,10),bonds(i,11),bonds(i,12),bonds(i,13), bonds(i,14), bonds(i,15)
   		endif  
	end do



100    call system('cp '//TRIM(outputfile1)//'_bonds1to7.txt '//TRIM(finalfolder))
       call system('rm '//TRIM(outputfile1)//'_bonds1to7.txt')
       if (Nbonds >= 8) then
       call system('cp '//TRIM(outputfile2)//'_bonds8to14.txt '//TRIM(finalfolder))
       call system('rm '//TRIM(outputfile2)//'_bonds8to14.txt')
       endif
       call system('echo -n '//xn//'_')
       
  close(10)
	close(20)
	close(21)

end do
print *, 'Done'
print *, 'To plot more files in xmgrace type   xmgrace -nxy file1 -nxy file2'
print *, 'Thank god, it works'

close(30)

end program AnalyzeBonds


