function [ Z ] = gradientDecentUS( F )

    global W RA C imageLength classType classGroups

    Y = F(:,1:C-1); % We can remove all labels from the data.

    eps = 1e1; ii = 0;

    % We need to generate the learning rate, and find the gradient. This
    % gradient is modified for infinite convolutional locii...

    N = size( RA, 3 ) * size( RA, 3 );

    RA_ = RA; 

    kk = 1; ll = 1;

    while ( ll < N - 1 )
        
        for mm = 2:1:size(classGroups,2)
            for k = 2:1:size(classType,2)
                for j = 1:1:C-1
                    for i = 2:imageLength
    
%                         gamma( i, j, mm ) = abs(...
%                                             ( Y( i, j ) - Y( i-1, j ) ) *...
%                                             ( RA( k, j, mm ) - RA( k-1, j, mm-1 ) ) + eps ) /...
%                                             ( abs( RA( k, j, mm ) - RA( k-1, j, mm-1 ) ) + eps )^2;

                        gamma( i, j, mm ) = 1;
                    end
                end
    
                Y = F(:,1:C-1);
            end
        end
    
        for mm = 1:1:size(classGroups,2)
            for j = 1:1:C-1
                for i = 1:1:size(classType,2)
                    
                    eps( i, j, mm ) = RA( i, j, mm ) / W( i, j, mm );

                end
            end
        end
        
        for mm = 1:1:size(classGroups,2)
            for k = 1:1:size(classType,2)
                for j = 1:1:C-1
                    for i = 1:1:imageLength
            
                        if ( RA( k, j, mm ) - Y( i, j ) > 0 )
                            while( Y( i, j ) < RA( k, j, mm ) )
                                Y( i, j ) = Y( i, j ) + ( gamma( i, j, k ) * RA( k, j, mm ) +...
                                    eps( k, j, mm ) );
                                ii = ii + 1;
                            end 
    
                            iiVec( i, j, k, mm ) = ii; ii = 0;
                        end
        
                        Y = F(:,1:C-1); 
                        
                        if ( RA( k, j, mm ) - Y( i, j ) < 0 )
                            while( Y( i, j ) > RA( k, j, mm ) )
                                Y( i, j ) = Y( i, j ) - ( gamma( i, j, k ) * RA( k, j, mm ) +...
                                    eps( k, j, mm ) );
                                ii = ii + 1;
                            end
    
                            iiVec( i, j, k, mm ) = ii; ii = 0;
                        end      
        
                    end   
                end
    
                Y = F(:,1:C-1);        
            end
    
            Y = F(:,1:C-1);
        end
    
        cc = 1;
        for mm = 1:1:size(iiVec,4)
            for k = 1:1:size(iiVec,3)
        
                S( cc, 1 ) = sum( sum( iiVec( :, :, k, mm ) ) );
    
                cc = cc + 1;
            end
        end
    
        [ X( kk, 1 ), P( kk, 1 ) ] = min( S ); 
        
        kk = kk + 1;

        ll = ll + 1;

        RA( 2, :, : ) = circshift( RA( 2, :, : ), 1, 3 );

        Y = F(:,1:C-1); 
    end
    
    RA = RA_;
end