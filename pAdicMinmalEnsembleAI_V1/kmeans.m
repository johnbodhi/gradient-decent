function [ Y ] = kmeans( dataSet, l_ )

global RA BPI BFLAG

    S = dataSet(:,1:end-1,:);

    V = size(BPI,2)*size(dataSet,2); N = 1;

    p = 2;

    W = zeros(size(RA,1),size(RA,2));

    % We need to weight our centroids, and backpropagate the
    % following gradient for every other class in a group.

    for ii = 1:1:size(W,2)
    
        if ( ii == 1 )

            W(ii,:) = 1e0;

        elseif ( ii == 2 )

            if( BFLAG )

                W(ii,:) = 1e0;

                BPI(1,1) = 0; 

            elseif( ~BFLAG ) % Backprop...

                W(ii,:) = 1e2; 

                BPI(1,2) = 0; 
                
            end
        end
        
    end
  
    for ii = 1:1:(N*V)  % Modified convergence criterion.

        for k = 1:1:size(BPI,2)
            for j = 1:1:size(S,2)
                for i = 1:1:size(S,1)
            
                    if ( BPI(1,1) )
            
                        D( i, j, k ) = ( ( S( i, j, k ) -...
                            RA( BPI(1,1), j, 2 ) )^p )^( 1 / p );  
                        
                    elseif ( BPI(1,2) )
            
                        D( i, j, k ) = ( ( S( i, j, k ) -...
                            RA( BPI(1,2), j, 2 ) )^p )^( 1 / p );
                        
                    end                    
                end
            end
        
            Ci(k,1) = W(k,1,1) * mean(D(1,:,k)); 
            
            Cn(k,:) = [ Ci(k,1) ];        
        end       
         
        for i = 2:size(S,1)
            for j = 2:size(S,2)

                H( i, j, 2 ) = ( ( S( i, j, 2 ) - ...
                    Cn( i-1, 1 ) )^p )^(1/p); 
            end
        end
        Y = H(:,:,:); 

    end
    
    for i = 2:1:size(l_,1)

        Y(i,end,2) = l_(i, 1, 2); % Re-append labels for completeness...
    end

end