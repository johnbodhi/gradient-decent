function [ RA ] = runningAverage( S_ )

global i n R RA W imgDecision CFLAG

    if( CFLAG )

        % During classification we can update RA with the decisions made
        % by the gradient.

        R( imgDecision, i ) = S_( 1, 1 ); 
    else

        R( n, i ) = S_( 1, 1 ); 
    end

    W = [ 1e0; 1e1; 1e0; 1e0 ];

    RA = [  W(1,1) * mean( R( 1, : ) ) 
            W(2,1) * mean( R( 2, : ) )
            W(3,1) * mean( R( 3, : ) ) 
            W(4,1) * mean( R( 4, : ) )  ];
end