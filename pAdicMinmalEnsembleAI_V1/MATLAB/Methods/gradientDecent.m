function [ Z ] = gradientDecent( F )

    global RA BPI

    Y = F(:,:,:);

    eps = 1e1; ii = 0;

    % We need to generate the learning rate, and find the gradient. 
    % This gradient is modified for infinite convolutional locii...

    for k = 2:1:size(Y,3)
        for j = 2:1:size(Y,2)
            for i = 2:size(Y,1)

                % We can try multiple partial finite differences...

                gamma( i, j, k ) = abs( ...
                                    ( Y( i, j, k ) - Y( i-1, j-1, k-1 ) ) *...
                                    ( RA( i, j, k ) - RA( i-1, j-1, k-1 ) ) + eps ) /...
                                    ( abs( RA( i, j, k ) - RA( i-1, j-1, k-1 ) ) + eps );
            end
        end
        Y = F(:,:,:);
    end
    eps = 0;

    eps_ = 0.0001;
    for k = 2:1:size(Y,3)
        for j = 2:1:size(Y,2)
            for i = 2:1:size(Y,1)
                
                % eps( i, j, k ) = Y( i, j, k ) + eps_; % Non-constant step size.
                eps( i, j, k ) = eps_;                  % Constant step size.
            end
        end
    end
    
    % We need to pass updated addresses in RA 
    % forward after optimization, and index only
    % those corresponding adresses in the sample...

    for k = 2:1:size(RA,3)
        for j = 2:1:size(Y,2)
            for i = 2:1:size(Y,1)
    
                if ( RA( i, j, k ) - Y( i, j, k ) > 0 )
                    while( Y( i, j, k ) < RA( i, j, k ) )
                        Y( i, j, k ) = Y( i, j, k ) + ( gamma( i, j, k ) * RA( i, j, k ) +...
                            eps( i, j, k ) );
                        ii = ii + 1;
                    end 

                    WALK( i, j, k ) = ii; ii = 0;
                end

                Y = F(:,:,:); 
                
                if ( RA( i, j, k ) - Y( i, j, k ) < 0 )
                    while( Y( i, j, k ) > RA( i, j, k ) )
                        Y( i, j, k ) = Y( i, j, k ) - ( gamma( i, j, k ) * RA( i, j, k ) +...
                            eps( i, j, k ) );
                        ii = ii + 1;
                    end

                    WALK( i, j, k ) = ii; ii = 0;
                end      

            end   
        end
        Y = F(:,:,:);        
    end
    
    ii = 1;
    for i = 2:1:size(WALK,1)

        S(ii,1) = sum(sum(WALK(i,2:end,2:end),2),3);

        ii = ii + 1;
    end

    [ ~, Z ] = min(S(BPI(1,1):BPI(1,2),1)); % Decisions are not constrained to groups. 

    if( Z == 1 )

        Z = BPI(1,1);

    elseif ( Z == 2 )

        Z = BPI(1,2);
        
    end

end