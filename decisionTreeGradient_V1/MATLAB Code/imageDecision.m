function [ Z ] = imageDecision( Y ) 

    global ZA

    % We can implement decision trees on RA for an unknown input signal.

    % Z = gradientDecent( Y );

    % We can implement dyadic decision trees...

    M = size(ZA,2); N = M / 3;

    for i = 1:1:N

        X = ZA(:,(i-1)*M/N+1:i*M/N,:);

        Z_(i,1) = gradientDecent( Y, X );
    end
    
    Z = min(Z_);
end