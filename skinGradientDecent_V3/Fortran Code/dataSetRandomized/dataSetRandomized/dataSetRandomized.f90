program dataSetRandomized
use ifport
implicit none

integer, parameter      :: L = 67280, C = 4
integer,dimension(L,C)  :: X
integer,dimension(L)    :: R, G, B, Z
integer,dimension(L)    :: F
integer,dimension(2)    :: s
integer                 :: i, j, k
real                    :: q
real                    :: t,t_(2),t1,t2,ta(2)

call CPU_TIME(t1); t = etime(t_);

open( unit = 10, file = 'R.csv', action = 'readwrite' )
open( unit = 20, file = 'G.csv', action = 'readwrite' )
open( unit = 30, file = 'B.csv', action = 'readwrite' )
open( unit = 40, file = 'Z.csv', action = 'readwrite' )
open( unit = 100, file = 'dataSetRandomized.csv', action = 'readwrite')

do i = 1,L,1
    
    read(10,*) R(i)
    read(20,*) G(i)
    read(30,*) B(i)
    read(40,*) Z(i)
    F(i) = 0;
enddo

do j = 1,C,1
    do i = 1,L,1
        
        if(j .eq. 1) then
            
            X(i,j) = R(i);            
        elseif(j .eq. 2) then
            
             X(i,j) = G(i); 
        elseif(j .eq. 3) then
            
             X(i,j) = B(i); 
        elseif(j .eq. 4) then
            
             X(i,j) = Z(i);            
        endif     
        
    enddo
enddo

F(1) = 1;
write(100,*) X(F(1),:)
do k = 2,L,1   
    
    call random_number(q);
    
    F(k) = int(floor(q*L+1));
    
    do while( findloc(F(1:k-1),F(k),1) )
        
        call random_number(q);
        
        F(k) = int(floor(q*L+1));      
    enddo
    
    write(100,*) X(F(k),:);
enddo

print*," "
write(*,*) 'Program has used', t, 'seconds of CPU time.'
write(*,*) 'User Time: ',t_(1),' System Time: ',t_(2)

call CPU_TIME(t2)
write(*,*) 'Program has used', t2,'seconds.'
pause

end program dataSetRandomized