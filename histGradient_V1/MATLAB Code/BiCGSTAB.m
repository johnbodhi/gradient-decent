function [ Z_ ] = BiCGSTAB( X_, Y_ )

    global classType classGroups frameLength

    N   = size(Y_,1); 
    
    M   = size(classType,2);

    SUM = size(classGroups,2) * N * M; 

    % SUP = simpleNN(N,M); 
    
    UB  = 34912;

    SUP = size(classGroups,2) * UB;
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    LIMIT = 1e-5;

    while( sum(sum(sum(V,1),2),3) < SUM )

        K = ceil( ii / SUP ) + 1;

        if( aa <= size(V,2) )

            V(aa,1,K) = 1;

            aa = aa + 1;    
            
            A(:,1) = Y_(aa,:,K);            
            
        elseif( aa > size(V,3) && bb <= size(V,3) )

            V(bb,2,K) = 1; V(:,1,K) = 0;
            
            bb = bb + 1; aa = 1;
            
            B(:,1) = frameLength.*X_(bb,:,K);           

        elseif( bb > size(V,1) && cc <= size(V,1) )

            V(cc,3,K) = 1; V(:,2,K) = 0;

            cc = cc + 1; bb = 1;
            
            X(:,1) = Y_(cc,:,K);

        end
        
        R_0(:,1) = B(:,1) - A(:,1).*X(:,1); 
        
        R(:,1)   = R_0(:,1);

        RHO_0    = dot(R_0(:,1), R_0(:,1)); 
        
        RHO(1,1) = RHO_0;

        P(:,1)   = R_0(:,1);

        while( TOL <= LIMIT )
            
            V(:,1) = A(:,1).*P(:,1);


            ALPHA  = RHO(1,1) / dot( R(:,1), V(:,1) );

            H(:,1) = X(:,1) + ALPHA.*P(:,1);

            S(:,1) = R(:,1) - ALPHA.*V(:,1);            

            T(:,1) = A(:,1).*S(:,1);


            OMEGA  = dot( T(:,1), S(:,1) ) / dot( T(:,1), T(:,1) );

            X(:,1) = H(:,1) + OMEGA(:,1).*S(:,1);

            R(:,1) = S(:,1) - OMEGA(:,1).*T(:,1);


            RHO(1,2) = dot( R_0(:,1), R(:,1) );


            TOL      = RHO(1,2) / RHO_0;


            BETA     = ( RHO(1,2) / RHO(1,1) ) * ( ALPHA / OMEGA );

            P(:,1)   = R(:,1) + BETA.*( P(:,1) - OMEGA.*V(:,1) );

            RHO(1,1) = RHO(1,2);

            kk = kk + 1;
        end 

        WALK_(aa,bb,cc) = kk; kk = 1;
        
        ii = ii + 1;
    end

    WALK = sort(WALK_, 1,'descend');
    WALK = sort(WALK,  3,'descend');

    Z_ = 
end