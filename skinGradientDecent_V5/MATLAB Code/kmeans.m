function [ Y ] = kmeans( X, V )

global imageLength classGroups C RA

    p = 2;

    ii = 1; jj = 1; N = 1;
    
    while( jj < N * imageLength ) % Modified convergence criterion.

        if ( ii  == N * imageLength )
            ii = 1;
        end

        n = V( ii );
       
        % We need to mask our observations to index RA (the way it is shaped)...

        B = zeros(2,1);
        for j = 0:1:size(classGroups,2)
            
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

       TA = find(B(1,:) == n); TB = find(B(2,:) == n);

       if( TA )

           TB = 0;
       elseif( TB )

           TA = 0;
       end

       % We can find all 2norms with the mask...

        for j = 1:C-1
            for i = 1:size( X, 1 )
        
                if ( TA )
        
                    D( i, j ) = ( ( X( i, j ) - RA( 1, j, TA ) )^p )^( 1 / p );        
                elseif ( TB )
        
                    D( i, j ) = ( ( X( i, j ) - RA( 2, j, TB ) )^p )^( 1 / p );
                end                
        
            end
        end

        Ci = mean(D(:,1)); Cj = mean(D(:,2)); Ck = mean(D(:,3)); 

        Cn = [ Ci Cj Ck ]; 

        % Ci = mean(D(:,1));  Cn = [ Ci ]; 

        for j = 1:C-1
            for i = 1:imageLength

                Y( i, j ) = ( ( X( i, j ) - Cn( j ) )^p )^(1/p); 
            end
        end   

        ii = ii + 1; jj = jj + 1;
    end

    for i = 1:imageLength

        Y( i, 4 ) = V( i );
    end

end