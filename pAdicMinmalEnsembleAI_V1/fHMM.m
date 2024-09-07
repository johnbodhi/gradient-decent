function [ S ] = fHMM( A_ )
    
    % We can use a factorial hiddem markov model to converge as a non-convex
    % minimizer. The cost minimization is the discriminator ie. the 
    % output length.

    N = size(A_,1); M = size(A_,2);

    R = zeros(N,M); DELTA = zeros(1,2); E = 1.0; ii = 1;
    
    T = 1e-4;
    
    while( E > T )
        
        DELTA(1,1) = DELTA(1,2);

        for i = 1:1:N
            for j = 1:1:M
                    
                R(i,j) = A_(i,j) - DELTA(1,1);
            end
        end

        Rd = abs(det(R)); Ri = inv(R);
        
        % Y(:,1) = mean(A_,1); X(:,1) = mean(A_,2);
         
        Y(:,1) = mean(R,1); X(:,1) = mean(R,2);

        DELTA(1,2) = (2*pi*Rd)^(-0.5)*exp(-0.5*(Y-X)'*Ri*(Y-X));
         
        E = abs(DELTA(1,2) - DELTA(1,1));
         
        S(ii,1) = DELTA(1,2); ii = ii + 1;
    end
    F = length(S); disp(F);
  
end 