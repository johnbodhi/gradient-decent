function [ Y ] = kmeans( X, V )

global vectorLength C RA

    p = 2;

    ii = 1; jj = 1; N = 1;
    
    while( jj < N * vectorLength ) % Modified convergence criterion.

        if ( ii  == N * vectorLength )
            ii = 1;
        end

        n = V( ii );

        for j = 1:C-1
            for i = 1:size( X, 1 )
                D( i, j ) = ( ( X( i, j ) - RA( n, j ) )^p )^( 1 / p );
            end
        end
        
        Ci = mean( D( :, 1 ) ); Cn = Ci; % Dimensionality reduction...
        
        for j = 1:C-1
            for i = 1:vectorLength
                Y( i, j ) = ( ( X( i, j ) - Cn( j ) )^p )^(1/p); 
            end
        end   

        ii = ii + 1; jj = jj + 1;
    end

    for i = 1:vectorLength
        Y( i, 4 ) = V( i );
    end
end