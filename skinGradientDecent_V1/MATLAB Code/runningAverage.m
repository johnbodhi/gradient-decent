function [ RA ] = runningAverage( RGB )

global i n R G B RA W imgDecision CFLAG
    
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

    W = [ 1e1 1e1 1e1;...
            1e2 1e2 1e2;...
            1e3 1e3 1e3;...
            1e4 1e4 1e4;...
            1e5 1e5 1e5;...
            1e6 1e6 1e6     ];
    
    RA = [  W(1,1) * mean( R( 1, : ) ) W(1,2) * mean( G( 1, : ) ) W(1,3) * mean( B( 1, : ) ) 
            W(2,1) * mean( R( 2, : ) ) W(2,2) * mean( G( 2, : ) ) W(2,3) * mean( B( 2, : ) ) 
            W(3,1) * mean( R( 3, : ) ) W(3,2) * mean( G( 3, : ) ) W(3,3) * mean( B( 3, : ) ) 
            W(4,1) * mean( R( 4, : ) ) W(4,2) * mean( G( 4, : ) ) W(4,3) * mean( B( 4, : ) ) 
            W(5,1) * mean( R( 5, : ) ) W(5,2) * mean( G( 5, : ) ) W(5,3) * mean( B( 5, : ) )
            W(6,1) * mean( R( 6, : ) ) W(6,2) * mean( G( 6, : ) ) W(6,3) * mean( B( 6, : ) ) ];
end
