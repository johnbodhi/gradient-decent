function [ Z ] = gradientDecent( F )

    global RA

    Y = F(:,:,1:size(F,2)-1); % We can remove all labels from the data.

    eps = 1e1; ii = 0; M = 1.00;

    % We need to generate the learning rate, and find the gradient. This
    % gradient is modified for infinite convolutional locii...

   
    for k = 2:1:size(RA,3)
        for j = 1:1:ceil(M*size(Y,2))
            for i = 2:size(Y,1) % We process the entire split frame at once.

                gamma( i, j ) = abs( ...
                                    ( Y( i, j ) - Y( i-1, j ) ) *...
                                    ( RA( i, j, k ) - RA( i-1, j, k-1 ) ) + eps ) /...
                                    ( abs( RA( i, j, k ) - RA( i-1, j, k-1 ) ) + eps )^2;
            end
        end

        Y = F(:,:,1:size(F,2)-1);
    end

    for k = 1:1:size(RA,3)
        for j = 1:1:ceil(M*size(Y,2))
            for i = 1:1:size(Y,1)
                
                eps( i, j, k ) = min(RA( i, :, k )); % Constant step size.
            end
        end
    end
    
    for k = 1:1:size(RA,3)
        for j = 1:1:ceil(M*size(Y,2))
            for i = 1:1:size(Y,1)
    
                if ( RA( i, j, k ) - Y( i, j ) > 0 )
                    while( Y( i, j ) < RA( i, j, k ) )
                        Y( i, j ) = Y( i, j ) + ( gamma( i, j ) * RA( i, j, k ) +...
                            eps( i, j, k ) );
                        ii = ii + 1;
                    end 

                    iiVec( i, j, k ) = ii; ii = 0;
                end

                Y = F(:,:,1:size(F,2)-1); 
                
                if ( RA( i, j, k ) - Y( i, j ) < 0 )
                    while( Y( i, j ) > RA( i, j, k ) )
                        Y( i, j ) = Y( i, j ) - ( gamma( i, j ) * RA( i, j, k ) +...
                            eps( i, j ,k ) );
                        ii = ii + 1;
                    end

                    iiVec( i, j, k ) = ii; ii = 0;
                end      

            end   
        end

        Y = F(:,:,1:size(F,2)-1);        
    end

    cc = 1;
    for k = 1:1:size(iiVec,3)

        S( cc, 1 ) = sum(sum(iiVec(:,:,k )));

        cc = cc + 1;
    end

    [ ~, Z ] = min( S( :, 1 ) ); % Decisions are not constrained to groups.
end