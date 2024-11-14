function [ Z ] = frameDecision( Y ) 

    global GD BiCGSTAB

    % We can utilize two different gradients to make a decision about a
    % frame.

    if ( GD )

        Z = gradientDecent( Y );

    elseif ( BiCGSTAB )

        Z = BiCGSTAB_( Y );

    end
 
end

