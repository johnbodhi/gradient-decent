function [ Z ] = gradientDecent( F )

    global RA classType classGroups

    Y = F(:,:,1:size(F,2)-1); % We can remove all labels from the data.

    eps = 1e1; ii = 0;

    % We need to generate the learning rate, and find the gradient. This
    % gradient is modified for infinite convolutional locii...

   
    for k = 2:1:size(RA,3)
        for j = 1:1:size(Y,2)
            for i = 2:size(Y,1) % We process the entire split frame at once.

                gamma( i, j ) = abs( ...
                                    ( Y( i, 1 ) - Y( i-1, 1 ) ) *...
                                    ( RA( k, j ) - RA( k-1, j ) ) + eps ) /...
                                    ( abs( RA( k, j ) - RA( k-1, j) ) + eps )^2;
            end
        end

        Y = F(:,:,1:size(F,2)-1);
    end


    for mm = 1:1:size(classGroups,2)
        for j = 1:1:size(Y,2)
            for i = 1:1:size(classType,2)
                
                eps( i, j, mm ) = RA( i, j, mm ); % Constant step size.
            end
        end
    end
    
    for mm = 1:1:size(classGroups,2)
        for k = 1:1:size(classType,2)
            for j = 1:1:size(Y,2)
                for i = 1:1:size(Y,1)
        
                    if ( RA( k, j, mm ) - Y( i, j ) > 0 )
                        while( Y( i, j ) < RA( k, j, mm ) )
                            Y( i, j ) = Y( i, j ) + ( gamma( i, j, k ) * RA( k, j, mm ) +...
                                eps( k, j, mm ) );
                            ii = ii + 1;
                        end 

                        iiVec( i, j, k, mm ) = ii; ii = 0;
                    end
    
                    Y = F(:,:,1:size(F,2)-1); 
                    
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

            Y = F(:,:,1:size(F,2)-1);        
        end

        Y = F(:,:,1:size(F,2)-1);
    end

    cc = 1;
    for mm = 1:1:size(iiVec,4)
        for k = 1:1:size(iiVec,3)
    
            S( cc, 1 ) = sum(sum(iiVec(:,:,k,mm)));

            cc = cc + 1;
        end
    end

    [ ~, Z ] = min( S( :, 1 ) ); % Decisions are not constrained to groups.
end