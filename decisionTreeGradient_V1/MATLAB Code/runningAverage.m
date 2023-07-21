function [ ZA ] = runningAverage( RGB )

global i n ZA R G B W classType classGroups imgDecision CFLAG

    if( CFLAG )

        % During classification we can update RA with the decisions made
        % by the gradient.

        R( imgDecision, i ) = RGB( 1, 1 ); 
        G( imgDecision, i ) = RGB( 1, 2 ); 
        B( imgDecision, i ) = RGB( 1, 3 );
    else

        R( n, i ) = RGB( 1, 1 ); 
        G( n, i ) = RGB( 1, 2 ); 
        B( n, i ) = RGB( 1, 3 );
    end

    uu = 1;
    for kk = 1:1:size(classGroups,2)
        for jj = 1:1:size(ZA,2)/3
            for ii = 1:1:size(classType,2)  
                
                if ( uu == n )
    
                    ZA(ii,1*jj*M,kk) = W(ii,1*jj*M,kk) * mean( R( n, : ), 2 );
        
                    ZA(ii,2*jj*M,kk) = W(ii,2*jj*M,kk) * mean( G( n, : ), 2 );
        
                    ZA(ii,3*jj*M,kk) = W(ii,3*jj*M,kk) * mean( B( n, : ), 2 );
                end
            end
        end

        uu = uu + 1;
    end
end