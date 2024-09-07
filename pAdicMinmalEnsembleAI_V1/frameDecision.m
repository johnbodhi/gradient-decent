function [ Z ] = frameDecision( Y ) 

    Z = gradientDecent( Y );
    
    % Z = GMRES( Y );
end

