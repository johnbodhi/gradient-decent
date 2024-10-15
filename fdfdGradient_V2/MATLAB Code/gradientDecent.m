function [ Z ] = gradientDecent( F )

    global W RA C vectorLength

    Y = F(:,1:C-1); % We can remove all labels from the data.

    eps = 1e1; fieldClasses = 2; ii = 0;

    for k = 2:fieldClasses
        for j = 1:C-1
            for i = 2:vectorLength
                gamma( i, j, k ) = abs( ...
                    ( Y( i, j ) - Y( i - 1, j ) ) *...
                    ( RA( k, j ) - RA( k-1, j ) ) + eps ) /...
                    ( abs( RA( k, j ) - RA( k-1, j ) ) + eps )^2;
            end
        end
        Y = F(:,1:C-1);
    end

    for k = 1:fieldClasses
        for j = 1:C-1
            eps( k, j ) = RA( k, j ) / W( k, j ); % Efficieny problem...
        end
    end
    
    for k = 1:fieldClasses
        for j = 1:C-1
            for i = 1:vectorLength
    
                if ( RA( k, j ) - Y( i, j ) > 0 )
                    while( Y( i, j ) < RA( k, j ) )
                        Y( i, j ) = Y( i, j ) + ( gamma( i, j, k ) * RA( k, j ) + eps( k, j ) );
                        ii = ii + 1;
                    end 
                    iiVec( i, j, k ) = ii; ii = 0;
                    % runningAverage( RA( k, j ) );
                end

                Y = F(:,1:C-1);
                
                if ( RA( k, j ) - Y( i, j ) < 0 )
                    while( Y( i, j ) > RA( k, j ) )
                        Y( i, j ) = Y( i, j ) - ( gamma( i, j, k ) * RA( k, j ) + eps( k, j ) );
                        ii = ii + 1;
                    end
                    iiVec( i, j, k ) = ii; ii = 0;
                    % runningAverage( RA( k, j ) );
                end      

            end   
        end
        Y = F(:,1:C-1);        
    end        

    for k = 1:1:size(iiVec,3)

        S(k) = sum(sum(iiVec(:,:,k)));
    end
    
    [ ~, Z ] = min(S(:));
end