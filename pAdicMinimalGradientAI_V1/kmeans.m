function [ Y ] = kmeans( dataSet, l_ )

global classGroups classType RA BPI BFLAG

    S = dataSet(:,1:end-1,:);

    V = size(RA,1)*size(dataSet,2);

    p = 2;

    W = zeros(size(RA,1),size(RA,2));

    % We need to weight our centroids, and backpropagate the
    % following gradient for every other class in a group.

    for kk = 1:1:size(W,2)
        for ii = 1:1:size(W,2)
        
            if ( ii == 1 )
    
                W(ii,:,kk) = 1e0;

            elseif ( ii == 2 )

                if( BFLAG )
    
                    W(ii,:,kk) = 1e0;

                elseif( ~BFLAG )

                    W(ii,:,kk) = 1e2; % Backprop...
                    
                end
            end
            
        end
    end

    ii = 2; jj = 1; N = 1;
    
    while( jj < N*V ) % Modified convergence criterion.

        if ( ii == N*V )

            ii = 2;
        end

        n = l_( ii, 1, 2 );
       
        % We need to mask the batch indexes.

        B = zeros(1,1);
        for j = 0:1:size(RA,1)-1
            
            B(1,j+1) = j;
            
            % B(1,j+1) = 2*j+1;

            % B(2,j+1) = 2*(j+1);
        end

        for j = 1:1:size(B,2)
            % for i = 1:1:size(B,1)

                if( n ~= B(1,j) )
                    
                    B(1,j) = 0;
                end
                
            % end
        end

       TA = find( B(1,:) == n ); % TB = find( B(2,:) == n );

       % if( TA )

           % TB = 0;
       % elseif( TB )

           % TA = 0;
       % end

       % We can find all 2norms with the mask...
        
        for k = 1:1:size(BPI,2)
            for j = 1:1:size(S,2)
                for i = 1:1:size(S,1)
            
                    if ( TA )
            
                        D( i, j, k ) = ( ( S( i, j, k ) -...
                            RA( TA , j, 2 ) )^p )^( 1 / p );  
                        
                    % elseif ( TB )
            
                        % D( i, j, k ) = ( ( S( i, j, k ) -...
                           % RA( TB, j, 2 ) )^p )^( 1 / p );
                        
                    end            
                end
            end
        
            Ci(k,1) = W(k,1,1) * mean(D(1,:,k)); 
            
            % Cj(k,1) = W(k,2,1) * mean(D(2,:,k)); 
            % 
            % Ck(k,1) = W(k,3,1) * mean(D(3,:,k)); 
            % 
            % Cn(k,:) = [ Ci(k,1) Cj(k,1) Ck(k,1) ];
            
            Cn(k,:) = [ Ci(k,1) ];        
        end       
        
        for bp = 1:1:size(BPI,2)
            for i = 1:size(S,1)
                for j = 1:size(S,2)
    
                    H( i, j, bp ) = ( ( S( i, j ) - ...
                        Cn( bp, 1 ) )^p )^(1/p); 
                end
            end   
        end

        Y = H(:,:,1); Y = cat(1,Y,H(:,:,2));

        ii = ii + 1; jj = jj + 1;
    end
    
    for i = 2:1:size(l,1)

        Y(i,end,2) = l_(i, 1, 2); % Re-append labels for completeness...
    end

end