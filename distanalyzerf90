program DistanceAnalyzer
implicit none
real,allocatable,dimension(:,:) :: array,bonds,arrayhelp
integer,allocatable,dimension(:,:) :: existbonds
integer :: i,j,k,m,n,o,q,Maxsteps,step,Natoms,Nbonds,choice,Ntraj
character*30 :: junk, mainfolder, moviesource, trajfolder, xn, form, outputfile1,outputfile2,outputfile3, finalfolder
real :: time,timestep,attofs
attofs=0.0241888425

open (40,file='inanalyze.in', status='old')
read (40,*) junk, Natoms
read (40,*) junk, timestep
read (40,*) junk, Maxsteps
read (40,*) junk, Ntraj
read (40,*) junk, mainfolder
read (40,*) junk, trajfolder
!print *, mainfolder
                     

!----------------------------------------------------------------------------------------------
!Which bonds exist? Which bonds we want?
open(30,file='existbonds.txt',status='replace')
allocate(existbonds(Natoms,2))

print *, 'How many bonds?'
read *, choice

if (choice == 1000) then			
Nbonds=14											
existbonds(1,1)=1			
existbonds(1,2)=2				
existbonds(2,1)=1				
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


if (choice /= 1000) then								
print *, 'Write pairs of atom numbers to scan distance between them.'		
print *, 'The order of atoms is the same as in movie.xyz' 		
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
	outputfile2=xn
    outputfile3=xn					        !
!print *,TRIM(moviesource)//'/movie.xyz'
	open(10,file=TRIM(moviesource)//'/movie.xyz',status='old')                         
		open(20,file=TRIM(outputfile1)//'_bonds1to6.txt',status='replace')                 
    		write(20,*) '#Time (fs)', '        #Bond lenght (A)'
	if (Nbonds >= 7) then                                                   !In case of Nbonds > 14 it is necessary  
  	  open(21,file=TRIM(outputfile2)//'_bonds7to12.txt',status='replace')   !to create another outputfile.
          write(21,*) '#Time (fs)', '        #Bond lenght (A)'		        !There are max 7 scans per outputfile
	endif
    if (Nbonds >= 13) then
		open(22,file=TRIM(outputfile3)//'_bonds13to18.txt',status='replace')
			write(22,*) '#Time (fs)', '        #Bond lenght (A)'
    end if
    
	do i=1,Maxsteps									!Here I am getting data from movie
		read(10,*,end=100)                              !computing bond lenghts and writing them
    		read(10,*,end=100) junk,junk,step
    
    		do j=1,Natoms
	       		read(10,*) junk,array(j,1),array(j,2),array(j,3)
        	end do

        	do m=1,Nbonds
      			arrayhelp(1,1)=(array(existbonds(m,1),1)-array(existbonds(m,2),1))**2.0
   				arrayhelp(1,2)=(array(existbonds(m,1),2)-array(existbonds(m,2),2))**2.0
           		arrayhelp(1,3)=(array(existbonds(m,1),3)-array(existbonds(m,2),3))**2.0
   				bonds(i,m+1)=sqrt(arrayhelp(1,1)+arrayhelp(1,2)+arrayhelp(1,3))
       		end do
            
         time=timestep*step*attofs
         bonds(i,1)=time
         
			if (Nbonds >= 7 .AND. Nbonds <13) then
  		write(20,*) bonds(i,1),(bonds(i,q),q=2,7)
   		write(21,*) bonds(i,1),(bonds(i,q),q=8,Nbonds+1)   
   			endif 
            if (Nbonds >= 13) then
   		write(20,*) bonds(i,1),(bonds(i,q),q=2,7)
   		write(21,*) bonds(i,1),(bonds(i,q),q=8,13)
   		write(22,*) bonds(i,1),(bonds(i,q),q=14,Nbonds+1)
   			else
   		write(20,*) bonds(i,1),(bonds(i,q),q=2,Nbonds+1) 
            endif  
	end do

	   call system('cp '//TRIM(outputfile1)//'_bonds1to6.txt '//TRIM(finalfolder))
       call system('rm '//TRIM(outputfile1)//'_bonds1to6.txt')
       if (Nbonds >= 7) then
       call system('cp '//TRIM(outputfile2)//'_bonds7to12.txt '//TRIM(finalfolder))
       call system('rm '//TRIM(outputfile2)//'_bonds7to12.txt')
       endif
       if (Nbonds >= 13) then
       call system('cp '//TRIM(outputfile3)//'_bonds13to18.txt '//TRIM(finalfolder))
       call system('rm '//TRIM(outputfile3)//'_bonds13to18.txt')
       end if
       call system('echo -n '//xn//'_')
       
100 close(10)
	close(20)
	close(21)
    close(22)
    
end do

print *, 'Done'
print *, 'To plot more files in xmgrace write   xmgrace -nxy file1 -nxy file2'


close(30)
close(40)

end program DistanceAnalyzer
