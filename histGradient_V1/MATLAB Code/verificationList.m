function [ Y ] = verificationList( N )

    global pp

    T = sum(N);

    ii = 1; jj = 1; kk = 1; pp = 1;
    for i = 1:N(ii,1):2*T
    
        while ( jj <= N(ii,1) )
    
            Y(kk,1) = ii;
    
            jj = jj + 1; kk = kk + 1;
        end
        ii = ii + 1; jj = 1;
    
        if ( size(Y,1) >= T ) % Unsupervised observations.
            break;
        end
    end
    
end