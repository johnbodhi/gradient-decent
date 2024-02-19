function [ Y ] = kmeans( dataSet, l )

global classGroups classType GFLAG RA BPI KFLAG

    S = dataSet;

    V = size(RA,1)*size(dataSet,2);

    p = 2;

    ii = 1; jj = 1; N = 1;

    W  = zeros( size(classType,2), 3, size(classGroups,2) );

    % We need to weight our centroids, and backpropagate the
    % following gradient for every other class in a group.

    for kk = 1:1:size(classGroups, 2)
        for ii = 1:1:size(classType, 2)
        
            if ( ii == 1 )
    
                W(ii,:,kk) = 1e0;

            elseif ( ii == 2 )

                if( KFLAG )
    
                    W(ii,:,kk) = 1e0;

                elseif( ~KFLAG )

                    W(ii,:,kk) = 1e2; % Backprop...

                end

            end
            
        end
    end
    
    while( jj < N*V ) % Modified convergence criterion.

        if ( ii == N*V )

            ii = 1;
        end

        n = l( ii, 1 );
       
        % We need to mask the batch indexes.

        B = zeros(2,1);
        for j = 0:1:size(classGroups, 2)-1
            
            B(1,j+1) = 2*j+1;

            B(2,j+1) = 2*(j+1);
        end

        for j = 1:1:size(B,2)
            for i = 1:1:size(B,1)

                if( n ~= B(i,j) )
                    B(i,j) = 0;
                end
                
            end
        end

       TA = find(B(1,:) == n ); TB = find(B(2,:) == n );

       if( TA )

           TB = 0;
       elseif( TB )

           TA = 0;
       end

       % We can find all 2norms with the mask...

        if( l ) 
            
            for k = 1:1:size(BPI,2)
                for j = 1:size(RA,2)
                    for i = 1:size(S,1)
                
                        if ( TA == BPI(1,1) )
                
                            D( i, j, k ) = ( ( S( i, j, k ) - RA( BPI(1,1), j, k ) )^p )^( 1 / p );  
                            
                        elseif ( TB == BPI(1,2) )
                
                            D( i, j, k ) = ( ( S( i, j, k ) - RA( BPI(1,2), j, k ) )^p )^( 1 / p );
                            
                        end            
                    end
                end

                Ci(k,1) = W(k,1,1) * mean(D(1,:,k)); 
                
                Cj(k,1) = W(k,2,1) * mean(D(2,:,k)); 
                 
                Ck(k,1) = W(k,3,1) * mean(D(3,:,k)); 
    
                Cn(k,:) = [ Ci(k,1) Cj(k,1) Ck(k,1) ];
                
                % Cn(k,:) = [ Ci(k,1) ];

            end       
            
        else 

        % We can use k-means as a classifier...

%             for k = 1:1:size(RA,3)
%                 for ii = 1:1:size(RA,1)
%                     for j = 1:C-1                    
%                         for i = 1:size( X, 1 )
%         
%                             D( i, j, ii, k ) = ( ( X( i, j ) - RA( ii, j, k ) )^p )^( 1 / p );  
%                         end
%                     end    
% 
%                     Z_(ii,1) = sum(sum(D(:,:,ii,k)));
%                 end
% 
%                 Z(k,1) = min(Z_);
%             end
% 
%             [~,W] = min(Z);
%        
%             Ci = mean(D(:,1,W)); Cj = mean(D(:,2,W)); Ck = mean(D(:,3,W)); 
%     
%             Cn = [ Ci Cj Ck ];  
        end
        
        for k = 1:1:size(BPI,2)
            for j = 1:size(S,3)
                for i = 1:size(S,1)
    
                    H( i, j, k ) = ( ( S( i, j, k ) - Cn( k, 1 ) )^p )^(1/p); 
                end
            end   
        end

        Y = H(:,:,:); Y = cat(1,Y,H(:,:,2));

        ii = ii + 1; jj = jj + 1;
    end

    for i = 1:1:V

        Y(i,end,1) = l(i, 1); % Re-append labels for completeness...
    end
end