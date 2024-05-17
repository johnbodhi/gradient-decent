function [ DELTA ] = fHMM( A_ )
    
    % We can use a factorial hiddem markov model to converge as a 
    % non-convex minimizer. The cost minimization is the discriminator 
    % ie.the output length of the covergence.

    N = size(A_,1); M = size(A_,2);

    R = zeros(N,M); DELTA = zeros(1,2); 
    
    E = 1.0; ii = 1; jj = 1; kk = 1;
    
    T = 1e-4;
    
    for k = 1:1:N
        
        S_(:,1) = A_(k,:)
    
        for i = 1:1:N
            for j = 1:1:M
        
            D_(i,j) = S_S(ii,1) - S_(j,1);
            end
            ii = ii + 1;
        end
        W(1,1)  = 1;
        
        while( E > T )
        
            DELTA(1,1) = DELTA(1,2);
        
            for i = 1:1:N
                for j = 1:1:M
                    
                    R(i,j) = D_(i,j) - W(1,1)*DELTA(1,1);
                end
            end
            Rd = abs(det(R)); Ri = inv(R);
        
            % Y(:,1) = mean(D_,1); X(:,1) = mean(D_,2);
         
            Y(:,1) = mean(R,1); X(:,1) = mean(R,2);

            DELTA(1,2) = (2*pi*Rd)^(-0.5)*exp(-0.5*(Y-X)'*Ri*(Y-X));
         
            E = abs(DELTA(1,2) - DELTA(1,1));
         
            S(jj,kk) = DELTA(1,2); jj = jj + 1;
        end
        F(kk,1) = length(S(:,kk)); disp(F); kk = kk + 1;
    end
    
end 
