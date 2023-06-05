function [ Z ] = gradientDecent( F )

    global W RA C imageLength

    Y = F(:,1:C-1); % We can remove all labels from the data.

    eps = 1e1; skinClasses = 2; ii = 0; ALPHA = 10;

    for k = 2:skinClasses
        for j = 1:C-1
            for i = 2:imageLength
                gamma( i, j, k ) = abs( ...
                    ( Y( i, j ) - Y( i - 1, j ) ) *...
                    ( RA( k, j ) - RA( k - 1, j ) ) + eps ) /...
                    ( norm( RA( k, j ) - RA( k-1, j ) ) + eps )^2;
            end
        end
        Y = F(:,1:C-1);
    end

    for k = 1:skinClasses
        for j = 1:C-1
            eps( k, j ) = RA( k, j ) / W( k, j );
        end
    end
    
    for k = 1:skinClasses
        for j = 1:C-1
            for i = 1:imageLength
    
                if ( RA( k, j ) - Y( i, j ) > 0 )
                    while( Y( i, j ) < RA( k, j ) )
                        Y( i, j ) = Y( i, j ) + ( gamma( i, j, k ) * RA( k, j ) + eps( k, j ) );
                        ii = ii + 1;
                    end 
                    iiVec( i, j, k ) = ii; ii = 0;
                end

                Y = F(:,1:C-1);
                
                if ( RA( k, j ) - Y( i, j ) < 0 )
                    while( Y( i, j ) > RA( k, j ) )
                        Y( i, j ) = Y( i, j ) - ( gamma( i, j, k ) * RA( k, j ) + eps( k, j ) );
                        ii = ii + 1;
                    end
                    iiVec( i, j, k ) = ii; ii = 0;
                end      

            end   
        end
        Y = F(:,1:C-1);        
    end        

    for k = 1:skinClasses
        for j = 1:C-1
            for i = 1:skinClasses
                if( size( find( iiVec( :, j, k ) <= ALPHA ), 1 ) )
                    T( i, j, k ) = size( find( iiVec( :, j, k ) <= ALPHA ), 1 );
                else 
                    T( i, j, k ) = 0;
                end  
            end
        end
    end

    for k = 1:skinClasses
        for j = 1:C-1
            for i = 1:skinClasses
                if( T( i, j, k ) < imageLength )
                    T( i, j, k ) = 0;
                end  
            end
        end
    end

   for k = 1:skinClasses
       for i = 1:skinClasses                
            M( i, k ) = size( find( T( i, :, k ) ~= 0 ), 2 );
       end
   end

    for j = 1:skinClasses          
        if( find( M( :, j ) ~= 0 ) )
            Z = j;     
        end
    end
end
