clear all; close all; clc;

N = 32; M = 8;

V = zeros(32,8);

rr = 0;

aa = 1; bb = 1; 
cc = 1; dd = 1; 
ee = 1; ff = 1;
gg = 1; hh = 1;

while( sum(sum(V,1),2) < M*N )

    if( aa <= size(V,1) )   

        V(aa,1) = 1;

        aa = aa + 1;
        
        rr = rr + 1;

    elseif( aa > size(V,1) && bb <= size(V,1) )
        
        V(bb,2) = 1;
        
        bb = bb + 1; aa = 1; V(:,1) = 0;

        rr = rr + 1; 

    elseif( bb > size(V,1) && cc <= size(V,1) )

        V(cc,3) = 1;

        cc = cc + 1; bb = 1; V(:,2) = 0;
        
        rr = rr + 1;

     elseif( cc > size(V,1) && dd <= size(V,1) )

        V(dd,4) = 1;

        dd = dd + 1; cc = 1; V(:,3) = 0;
        
        rr = rr + 1;

    elseif( dd > size(V,1) && ee <= size(V,1) )

        V(ee,5) = 1;

        ee = ee + 1; dd = 1; V(:,4) = 0;
        
        rr = rr + 1;

    elseif( ee > size(V,1) && ff <= size(V,1) )

        V(ff,6) = 1;

        ff = ff + 1; ee = 1; V(:,5) = 0;
        
        rr = rr + 1;

    elseif( ff > size(V,1) && gg <= size(V,1) )

        V(gg,7) = 1;

        gg = gg + 1; ff = 1; V(:,6) = 0;
        
        rr = rr + 1;

    elseif( gg > size(V,1) && hh <= size(V,1) )

        V(hh,8) = 1;

        hh = hh + 1; gg = 1; V(:,7) = 0;
        
        rr = rr + 1;
    end

end


