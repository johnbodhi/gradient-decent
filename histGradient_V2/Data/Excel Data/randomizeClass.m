function [ Y ] = randomizeClass( dataSet, N )

    global imageLength

    A = zeros( imageLength, size(dataSet,2), N ); 

    ii = 1;
    for k = 1:N
        for i = 1:( imageLength )
            A( i,1:size(dataSet,2), k ) = dataSet( ii,1:size(dataSet,2) ); ii = ii + 1;
        end
    end     
    
    i = zeros( N,1 );

    for k = 1:N

        i( k ) = randi( [ 1, N ] );

        while ( find( i( 1:k-1, 1 ) == i( k ) ) )
            i( k ) = randi( [ 1, N ] );      
        end

        B( : , :, k)  = A( :, :, i( k ) );
    end        

    for k = 2:1:N         
        Y((k-2)*imageLength+1:k*imageLength,1:size(dataSet,2)) = cat( 1, B( :, :, k-1 ), B( :, :, k ) );
    end  
end