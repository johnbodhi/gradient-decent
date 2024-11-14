function [ S_ ] = BiCGSTAB_( A, ~, B )
    
    global TOL LIMIT RA BPI

    TOL = 1; LIMIT = 1e-2; S_ = 0; 

    for j = 1:1:size(RA,1)
        
        ii = 1;

        % R_0(:,1) = B(:,1) - A(:,1).*X(:,1);
        
        R_0(:,1) = B(:,j) - A(:,1);
                
        R(:,1)   = R_0(:,1);
        
        RHO_0    = dot(R_0(:,1), R_0(:,1));
                
        RHO(1,1) = RHO_0;
        
        P(:,1)   = R_0(:,1);

        for i = 1:1:size(RA,2)
        
            while( TOL >= LIMIT )
     
                C(i,1) = A(i,1) * P(i,1);
                        
                ALPHA  = RHO(1,1) / dot( R(i,1), C(i,1) );
            
                % H(i,1) = X(i,1) + ALPHA*P(i,1);
            
                S(i,1) = R(i,1) - ALPHA * C(i,1);           
            
                T(i,1) = A(i,1) * S(i,1);
            
            
                OMEGA  = dot( T(i,1), S(i,1) ) / dot( T(i,1), T(i,1) );
            
                % X(i,1) = H(i,1) + OMEGA(i,1) * S(i,1);
            
                R(i,1) = S(i,1) - OMEGA(i,1) * T(i,1);
            
            
                RHO(1,2) = dot( R_0(i,1), R(i,1) );
            
            
                TOL      = abs( RHO(1,2) / RHO_0 );
            
            
                BETA     = ( RHO(1,2) / RHO(1,1) ) * ( ALPHA / OMEGA );
            
                P(i,1)   = R(i,1) + BETA * ( P(i,1) - OMEGA * C(i,1) );
            
                RHO(1,1) = RHO(1,2);
                            
                S_ = S_ + 1;
            end 
            TOL = 1; W_(ii,1) = S_; 

            S_ = 0; ii = ii + 1;

        end

        WALK_(j,1) = sum(W_,1);
    end

    Z_ = min(WALK_);    

    [ ~, Z ] = min(S(BPI(1,1):BPI(1,2),1)); % Decisions are not constrained to groups. 

    if( Z == 1 )

        Z = BPI(1,1);

    elseif ( Z == 2 )

        Z = BPI(1,2);        
    end
end