function [ Y ] = randomizedPhotos( dataSet, Nr, Mr, L, C, N )

    M = 0;
    for j = 3:1:size(N,2)

        M = M + N(1,j);
    end
    N = M;

    A = zeros( Nr * Mr, C, N ); Z = zeros( L, C );

    ii = 1;
    for k = 1:N
        for i = 1:( Nr * Mr )
            A( i,1:C, k ) = dataSet( ii,1:C ); ii = ii + 1;
        end
    end     

    if ( isfile( 'randomizedPhotos.csv' ) )
    
        Y = readmatrix('randomizedPhotos.csv');
    
    else
    
        i = zeros( N,1 );

        for k = 1:N

            i( k ) = randi( [ 1, N ] );

            while ( find( i( 1:k-1, 1 ) == i( k ) ) )
                i( k ) = randi( [ 1, N ] );      
            end

            B( : , :, k)  = A( :, :, i( k ) );
        end        

        for k = 2:2:N         
            Z((k-2)*Nr*Mr+1:k*Nr*Mr,1:C) = cat( 1, B( :, :, k-1 ), B( :, :, k ) );
        end

        randomizedPhotos = Z; Y = randomizedPhotos;   
                
        writematrix( randomizedPhotos, 'randomizedPhotos.csv' );               
    end    
end