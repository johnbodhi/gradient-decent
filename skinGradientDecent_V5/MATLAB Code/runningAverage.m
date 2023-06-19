function [ RA ] = runningAverage( RGB )

global i n RA R G B W classType classGroups imgDecision CFLAG

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
        for ii = 1:1:size(classType,2)  
            
            if ( uu == n )

                RA(ii,1,kk) = W(ii,1,kk) * mean( R( n, : ), 2 );
    
                RA(ii,2,kk) = W(ii,2,kk) * mean( G( n, : ), 2 );
    
                RA(ii,3,kk) = W(ii,3,kk) * mean( B( n, : ), 2 );
            end
        end
        
        uu = uu + 1;
    end
    
end