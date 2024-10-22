function [ H ] = kmeans( dataSet )

    global RA BPI BFLAG

    S = dataSet(:,1:end-1,:);

    W = zeros(size(RA,1),size(RA,2));

    for i = 1:1:size(W,2)
    
        if ( i == 1 )

            W(i,:) = 1e0;

        elseif ( i == 2 )

            if( BFLAG )

                W(i,:) = 1e0;

                BPI(1,1) = 0; 

            elseif( ~BFLAG ) % Backprop...

                W(i,:) = 1e2; 

                BPI(1,2) = 0; 
                
            end
        end
    end

    for i = 1:1:size(S,1)
        for j = 1:1:size(S,2)
    
            if ( BPI(1,1) )
    
                D( i, j ) = ( ( S( i, j ) - RA( BPI(1,1), j, 2 ) )^2 )^( 0.5 );  
                
            elseif ( BPI(1,2) )
    
                D( i, j ) = ( ( S( i, j) - RA( BPI(1,2), j, 2 ) )^2 )^( 0.5 );
                
            end
        end

        Cn(i,:) = W(i,1) * mean(D(i,:));
    end
         
    for i = 2:size(S,1)
        for j = 2:size(S,2)

            H( i, j, 2 ) = ( ( S( i, j, 2 ) - Cn( i-1, 1 ) )^2 )^( 0.5 ); 
        end
    end    
end