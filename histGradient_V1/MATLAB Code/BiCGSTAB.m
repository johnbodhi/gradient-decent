function [ RA ] = BiCGSTAB( X_, Y_ )

    global classType classGroups BINS frameLength

    N   = size(Y_,1); 
    
    M   = size(classType,2);

    O   = size(classGroups,2);

    SUM = N * M * O; 

    V   = zeros( N, M, O );

    % UB = simpleNN( N, M ); 
    
    UB  = 34912;

    SUP = size(classGroups,2) * UB;
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    LIMIT = 1e-5;

    if( N == 1 )

        X_ = frameLength.*X_; A = Y_;

        for ii = 1:1:size(X_,1)

            R_0(:,1) = X_(ii,:,1) - Y_(1,:,1).*Y_(1,:,1); 
            
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
    
            WALK_(1,jj) = kk; jj = jj + 1; kk = 1;

            [ ~, Z_ ] = min(WALK_);
        end
        
        RA(Z_,:,:) = ( RA(Z_,:,:) + Y_ ) ./ 2;
    else

        B(:,1) = frameLength.*X_(1,:,1);
        
        X(:,1) = Y_(1,:,1);
        
        while( sum(sum(sum(V,1),2),3) < SUM )
    
            K = ceil( ii / SUP ) + 1;
    
            if( aa <= size(V,1) )
    
                V(aa,1,K) = 1;
    
                aa = aa + 1;    
                
                A(:,1) = Y_(aa,:,K);            
                
            elseif( aa > size(V,1) && bb <= size(V,1) )
                
                V(bb,2,K) = 1; V(:,1,K) = 0;
    
                bb = bb + 1; aa = 1;
                
                X(:,1) = Y_(bb,:,K);
    
            elseif( bb > size(V,1) && cc <= size(V,2) )
    
                V(cc,3,K) = 1; V(:,2,K) = 0;
                
                cc = cc + 1; bb = 1;
                
                B(:,1) = frameLength.*X_(cc,:,K); 
                
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
    
        WALK  = sort(WALK_, 1,'descend');
        WALK  = sort(WALK,  3,'descend');
    
        EP_MU = 1; 
    
        for j = 1:1:size(WALK,2)    
            for k = 2:1:size(WALK,1)
    
                Z_(k,j) = WALK(1,j,k);
        
                if( ( WALK(1,j,k) - WALK(1,j,k-1) ) >= EP_MU )
                    
                    break;                
                end
            end    
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RA = zeros(M,BINS,O); ii = 1; ll = 1;

    for k = 1:1:size(classGroups,2)
        for j = 1:1:size(Z_,2)
            for i = 1:1:size(Z_,1)
        
                if( Z_(i,j) ~= 0 )

                    RA(ii,:,k) = RA(ii,:,k) + Y_(Z_(i,j),:,k); 

                    ll = ll + 1;
                else

                    break;
                end
                
            end
            RA(ii,:,k) = RA(ii,:,k) / ll; ll = 1;
            ii = ii + 1; 
        end
        ii = 1;
    end

endfunction [ RA ] = BiCGSTAB( X_, Y_ )

    global classType classGroups BINS frameLength

    N   = size(Y_,1); 
    
    M   = size(classType,2);

    O   = size(classGroups,2);

    SUM = N * M * O; 

    V   = zeros( N, M, O );

    % UB = simpleNN( N, M ); 
    
    UB  = 34912;

    SUP = size(classGroups,2) * UB;
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    LIMIT = 1e-5;

    if( N == 1 )

        X_ = frameLength.*X_; A = Y_;

        for ii = 1:1:size(X_,1)

            R_0(:,1) = X_(ii,:,1) - Y_(1,:,1).*Y_(1,:,1); 
            
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
    
            WALK_(1,jj) = kk; jj = jj + 1; kk = 1;

            [ ~, Z_ ] = min(WALK_);
        end
        
        RA(Z_,:,:) = ( RA(Z_,:,:) + Y_ ) ./ 2;
    else

        B(:,1) = frameLength.*X_(1,:,1);
        
        X(:,1) = Y_(1,:,1);
        
        while( sum(sum(sum(V,1),2),3) < SUM )
    
            K = ceil( ii / SUP ) + 1;
    
            if( aa <= size(V,1) )
    
                V(aa,1,K) = 1;
    
                aa = aa + 1;    
                
                A(:,1) = Y_(aa,:,K);            
                
            elseif( aa > size(V,1) && bb <= size(V,1) )
                
                V(bb,2,K) = 1; V(:,1,K) = 0;
    
                bb = bb + 1; aa = 1;
                
                X(:,1) = Y_(bb,:,K);
    
            elseif( bb > size(V,1) && cc <= size(V,2) )
    
                V(cc,3,K) = 1; V(:,2,K) = 0;
                
                cc = cc + 1; bb = 1;
                
                B(:,1) = frameLength.*X_(cc,:,K); 
                
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
    
        WALK  = sort(WALK_, 1,'descend');
        WALK  = sort(WALK,  3,'descend');
    
        EP_MU = 1; 
    
        for j = 1:1:size(WALK,2)    
            for k = 2:1:size(WALK,1)
    
                Z_(k,j) = WALK(1,j,k);
        
                if( ( WALK(1,j,k) - WALK(1,j,k-1) ) >= EP_MU )
                    
                    break;                
                end
            end    
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RA = zeros(M,BINS,O); ii = 1; ll = 1;

    for k = 1:1:size(classGroups,2)
        for j = 1:1:size(Z_,2)
            for i = 1:1:size(Z_,1)
        
                if( Z_(i,j) ~= 0 )

                    RA(ii,:,k) = RA(ii,:,k) + Y_(Z_(i,j),:,k); 

                    ll = ll + 1;
                else

                    break;
                end
                
            end
            RA(ii,:,k) = RA(ii,:,k) / ll; ll = 1;
            ii = ii + 1; 
        end
        ii = 1;
    end

end
