function [ Z ] = gradientDecent( F )

    global RA W

    Y = F(:,1:size(F,2)-1,:); % We can remove all labels from the data.

    eps = 1e1; ii = 0;

    % We need to generate the learning rate, and find the gradient. This
    % gradient is modified for infinite convolutional locii...

    for k = 2:1:size(RA,3)
        for j = 1:1:size(Y,2)
            for i = 2:size(Y,1)

                gamma( i, j, k ) = abs( ...
                                    ( Y( i, j, k ) - Y( i-1, j, k-1 ) ) *...
                                    ( RA( i, j, k ) - RA( i-1, j, k-1 ) ) + eps ) /...
                                    ( abs( RA( i, j, k ) - RA( i-1, j, k-1 ) ) + eps )^2;
            end
        end

        Y = F(:,1:size(F,2)-1,:);
    end

    for k = 1:1:size(RA,3)
        for j = 1:1:size(RA,2)
            for i = 1:1:size(RA,1)
                
                eps( i, j, k ) = max(RA( i, :, k )) / W( 1, k ); % Constant step size.
            end
        end
    end
    
    for k = 1:1:size(RA,3)
        for j = 1:1:size(Y,2)
            for i = 1:1:size(Y,1)
    
                if ( RA( i, j, k ) - Y( i, j, k ) > 0 )
                    while( Y( i, j, k ) < RA( i, j, k ) )
                        Y( i, j, k ) = Y( i, j, k ) + ( gamma( i, j, k ) * RA( i, j, k ) +...
                            eps( i, j, k ) );
                        ii = ii + 1;
                    end 

                    WALK( i, j, k ) = ii; ii = 0;
                end

                Y = F(:,1:size(F,2)-1,:); 
                
                if ( RA( i, j, k ) - Y( i, j, k ) < 0 )
                    while( Y( i, j, k ) > RA( i, j, k ) )
                        Y( i, j, k ) = Y( i, j, k ) - ( gamma( i, j, k ) * RA( i, j, k ) +...
                            eps( i, j ,k ) );
                        ii = ii + 1;
                    end

                    WALK( i, j, k ) = ii; ii = 0;
                end      

            end   
        end

        Y = F(:,1:size(F,2)-1,:);        
    end

    for i = 1:1:size(WALK,1)

        S(i,1) = sum(sum(WALK(i,:,:),2),3);
    end

    [ ~, Z ] = min(S(:,1)); % Decisions are not constrained to groups.
end