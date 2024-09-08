function [ Y ] = kmeans( X, V )

global imageLength C RA

    p = 2;

    ii = 1; jj = 1; N = 1;
    
    while( jj < N * imageLength ) % Modified convergence criterion.

        if ( ii  == N * imageLength )
            ii = 1;
        end

        n = V( ii );

        for j = 1:C-3
            for i = 1:size( X, 1 )
                D( i, j ) = ( ( X( i, j ) - RA( n, j ) )^p )^( 1 / p );
            end
        end
    
        % Ci = mean( D( :, 1 ) ); Cj = mean( D( :, 2 ) ) ; Ck = mean( D( :, 3 ) ); Cn = [ Ci Cj Ck ]; 
        
        Ci = mean( D( :, 1 ) ); Cn = [ Ci ]; % Dimensionality reduction...
        
    
        for j = 1:C-3
            for i = 1:imageLength
                Y( i, j ) = ( ( X( i, j ) - Cn( j ) )^p )^(1/p); 
            end
        end   

        ii = ii + 1; jj = jj + 1;
    end

    for i = 1:imageLength
        Y( i, 4 ) = V( i );
    end
end