program simpleNN
use ifport
implicit none

integer, parameter        :: N = 3, M = 10
integer                   :: j
integer                   :: aa = 1, bb = 1, cc = 1, dd = 1
integer                   :: ee = 1, ff = 1, gg = 1, hh = 1
integer                   :: ii = 1, jj = 1
integer(8)                :: rr = 0
integer(8)                :: W
integer, dimension(N,M)   :: V
integer, dimension(M)     :: A

real	                  :: t,t_(2),t1,t2,ta(2)

call CPU_TIME(t1); t = etime(t_);


do while( W < N*M ) 

     if( aa <= size(V,1) ) then 

        V(aa,1) = 1;

        aa = aa + 1;
        
        rr = rr + 1;

    else if( aa > size(V,1) .and. bb <= size(V,1) ) then
        
        V(bb,2) = 1;
        
        bb = bb + 1; aa = 1; V(:,1) = 0;

        rr = rr + 1; 

    else if( bb > size(V,1) .and. cc <= size(V,1) ) then

        V(cc,3) = 1;

        cc = cc + 1; bb = 1; V(:,2) = 0;
        
        rr = rr + 1;

	else if( cc > size(V,1) .and. dd <= size(V,1) ) then

        V(dd,4) = 1;

        dd = dd + 1; cc = 1; V(:,3) = 0;
        
        rr = rr + 1;

    else if( dd > size(V,1) .and. ee <= size(V,1) ) then

        V(ee,5) = 1;

        ee = ee + 1; dd = 1; V(:,4) = 0;
        
        rr = rr + 1;

    else if( ee > size(V,1) .and. ff <= size(V,1) ) then

        V(ff,6) = 1;

        ff = ff + 1; ee = 1; V(:,5) = 0;
        
        rr = rr + 1;

    else if( ff > size(V,1) .and. gg <= size(V,1) ) then

        V(gg,7) = 1;

        gg = gg + 1; ff = 1; V(:,6) = 0;
        
        rr = rr + 1;

    else if( gg > size(V,1) .and. hh <= size(V,1) ) then

        V(hh,8) = 1;

        hh = hh + 1; gg = 1; V(:,7) = 0;
        
        rr = rr + 1;
        
     else if( hh > size(V,1) .and. ii <= size(V,1) ) then

        V(ii,9) = 1;

        ii = ii + 1; hh = 1; V(:,8) = 0;
        
        rr = rr + 1;
     else if( ii > size(V,1) .and. jj <= size(V,1) ) then

        V(jj,10) = 1;

        jj = jj + 1; ii = 1; V(:,9) = 0;
        
        rr = rr + 1;
        
    end if
    
    do j = 1,size(V,2),1
        
        A(j) = sum(V(:,j));         
    enddo
    
    W = sum(A);
    
end do

print*, rr

print*," "
write(*,*) 'Program has used', t, 'seconds of CPU time.'
write(*,*) 'User Time: ',t_(1),' System Time: ',t_(2)

call CPU_TIME(t2)
write(*,*) 'Program has used', t2,'seconds.'
pause

end program simpleNN