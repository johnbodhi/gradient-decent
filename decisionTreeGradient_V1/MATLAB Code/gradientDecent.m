function [ Z ] = gradientDecent( F )

    global W C RA imageLength classType classGroups A uu vv

    Y = F(:,1:C-3); % We can remove all labels from the data.

    eps = 1e1; ii = 0;

    % We need to generate the learning rate, and find the gradient. This
    % gradient is modified for infinite convolutional locii...

    for mm = 2:1:size(classGroups,2)
        for k = 2:1:size(classType,2)
            for j = 1:1:C-3
                for i = 2:imageLength*A % We process the entire split frame at once.

                    gamma( i, j, mm ) = abs( ...
                                        ( Y( i, 1 ) - Y( i-1, 1 ) ) *...
                                        ( RA( k, j, mm ) - RA( k-1, j, mm-1 ) ) + eps ) /...
                                        ( abs( RA( k, j, mm ) - RA( k-1, j, mm-1 ) ) + eps )^2;
                end
            end

            Y = F(:,1:C-3);
        end
    end

    for mm = 1:1:size(classGroups,2)
        for j = 1:1:C-3
            for i = 1:1:size(classType,2)
                
                eps( i, j, mm ) = RA( i, j, mm ) / W( i, j, mm ); % Constant step size.
            end
        end
    end
    
    for mm = 1:1:size(classGroups,2)
        for k = 1:1:size(classType,2)
            for j = 1:1:C-3
                for i = 1:1:imageLength*A
        
                    if ( RA( k, j, mm ) - Y( i, j ) > 0 )
                        while( Y( i, j ) < RA( k, j, mm ) )
                            Y( i, j ) = Y( i, j ) + ( gamma( i, j, k ) * RA( k, j, mm ) +...
                                eps( k, j, mm ) );
                            ii = ii + 1;
                        end 

                        iiVec( i, j, k, mm ) = ii; ii = 0;
                    end
    
                    Y = F(:,1:C-3); 
                    
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

            Y = F(:,1:C-3);        
        end

        Y = F(:,1:C-3);
    end

    cc = 1;
    for mm = 1:1:size(iiVec,4)
        for k = 1:1:size(iiVec,3)
    
            S( cc, 1 ) = sum(sum(iiVec(:,:,k,mm)));

            cc = cc + 1;
        end
    end

    % We can implement segmentation in the classification process by
    % classifying only single groups at a time.

%     [ ~, Z ] = min( S( uu:vv, 1 ) ); % Decisions are constrained to groups.
% 
%     if( Z == 1 )
% 
%         Z = uu;
%     elseif ( Z == 2 )
% 
%         Z = vv;        
%     end

    [ ~, Z ] = min( S( :, 1 ) ); % Decisions are not constrained to groups.
end