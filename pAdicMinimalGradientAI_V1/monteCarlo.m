function [ Y ] = monteCarlo( dataSet )
    
    N = size(dataSet,1);
    
    A = dataSet; 
    
    % Casino stack...
    
    % A = cat(1,A,A,A);
    
    Y = zeros( size(dataSet,1), size(dataSet,2), size(dataSet,3) );

    if ( isfile( 'randomizedFrames.csv' ) )
    
        % Y = readmatrix('randomizedFrames.csv');
    
    else
    
        i = zeros( size(dataSet, 1), 1 );

        for k = 1:N

            i( k ) = randi( [ 1, N ] );

            while ( find( i( 1:k-1, 1 ) == i( k ) ) )
                
                i( k ) = randi( [ 1, N ] );      
            end

            B( k, : , : )  = A( i( k ), :, : );
        end        

        for k = 2:1:N     
            
            Y((k-2)+1:k,:,:) = cat( 1, B( k-1, :, : ), B( k, :, : ) );
        end 
                
        % writematrix( Y, 'randomizedFrames.csv' );        
    end   
    
end