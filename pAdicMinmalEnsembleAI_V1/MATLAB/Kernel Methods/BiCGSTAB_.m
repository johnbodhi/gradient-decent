function [ S_ ] = BiCGSTAB_( A, ~, B )
    
    global TOL LIMIT
    
    S_ = 0;
    
    % R_0(:,1) = B(:,1) - A(:,1).*X(:,1);
    
    R_0(:,1) = B(:,1) - A(:,1);
            
    R(:,1)   = R_0(:,1);
    
    RHO_0    = dot(R_0(:,1), R_0(:,1));
            
    RHO(1,1) = RHO_0;
    
    P(:,1)   = R_0(:,1);
    
    while( TOL >= LIMIT )
                
        C(:,1) = A(:,1).*P(:,1);
                
        ALPHA  = RHO(1,1) / dot( R(:,1), C(:,1) );
    
        % H(:,1) = X(:,1) + ALPHA.*P(:,1);
    
        S(:,1) = R(:,1) - ALPHA.*C(:,1);           
    
        T(:,1) = A(:,1).*S(:,1);
    
    
        OMEGA  = dot( T(:,1), S(:,1) ) / dot( T(:,1), T(:,1) );
    
        % X(:,1) = H(:,1) + OMEGA(:,1).*S(:,1);
    
        R(:,1) = S(:,1) - OMEGA(:,1).*T(:,1);
    
    
        RHO(1,2) = dot( R_0(:,1), R(:,1) );
    
    
        TOL      = abs( RHO(1,2) / RHO_0 );
    
    
        BETA     = ( RHO(1,2) / RHO(1,1) ) * ( ALPHA / OMEGA );
    
        P(:,1)   = R(:,1) + BETA.*( P(:,1) - OMEGA.*C(:,1) );
    
        RHO(1,1) = RHO(1,2);
                    
        S_ = S_ + 1;
    end 
    TOL = 1;
    
end