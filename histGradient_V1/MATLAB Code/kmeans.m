function [ Y ] = kmeans( S, V )

global imageLength classGroups classType RA A

    p = 2;

    ii = 1; jj = 1;  

    W = [ 1e0 1e1 1e2 ]; 
    
    while( jj < N*imageLength*A ) % Modified convergence criterion.

        if ( ii == imageLength*A )

            ii = 1;
        end

        n = V( ii, 1 );
            
        for k = 1:1:size(RA,3)
            for j = 1:1:size(S,2)-1
                for i = 1:1:size( S, 1 )
            
                    D( i, j, k ) = ( ( S( i, j ) - RA( i, j, k ) )^p )^( 1 / p );        
           
                end
            end

            Ci(k,1) = W(k,1) * mean(D(:,1,k)); 
            
            Cj(k,1) = W(k,1) * mean(D(:,2,k)); 
            
            Ck(k,1) = W(k,1) * mean(D(:,3,k)); 

            Cn(k,:) = [ Ci(k,1) Cj(k,1) Ck(k,1) ];            
        end       
        
        for k = 1:1:size(RA,3)
            for j = 1:1:size(S,2)-1
                for i = 1:1:size(S,1)
    
                    Y( i, j, k ) = ( ( S( i, j ) - Cn( k, j ) )^p )^(1/p); 
                end
            end   
        end

        ii = ii + 1; jj = jj + 1;
    end

end