function [ Y, i ] = randomizeAll( dataSet )

    global frameLength
    
    N = size(dataset,1);

    Y = zeros(size(dataSet,1),size(dataSet,2));

    M = 0;
    for i = 1:1:size(N,1)

        M = M + N(i,1);
    end
    N = M;

    A = zeros( frameLength, size(dataSet,2), N );

    ii = 1;
    for k = 1:N
        for i = 1:1:imageLength
            A( i,1:size(dataSet,2),k) = dataSet(ii,1:size(dataSet,2)); ii = ii + 1;
        end
    end     

    if ( isfile( 'randomizedFrames.csv' ) )
    
%         Y = readmatrix('randomizedFrames.csv');    
    else
    
        i = zeros( N,1 );

        for k = 1:N

            i( k ) = randi( [ 1, N ] );

            while ( find( i( 1:k-1, 1 ) == i( k ) ) )
                i( k ) = randi( [ 1, N ] );      
            end

            B( : , :, k)  = A( :, :, i( k ) );
        end        

        for k = 2:1:N         
            Y((k-2)*frameLength+1:k*frameLength,1:size(dataSet,2)) =...
                cat( 1, B( :, :, k-1 ), B( :, :, k ) );
        end 
                
%          writematrix( Y, 'randomizedFrames.csv' );               
    end    
end