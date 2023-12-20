function [ RA ] = BiCGSTAB( X_, Y_ )

    global classType classGroups BINS frameLength

    N   = size(Y_,1); 
    
    M   = size(classType,2);

    O   = size(classGroups,2);
    
    % IT = simpleNN(N,M);
    
    IT = O*(N^(M-1)+M);
    
    % CONV = O*(N^O) - ( O*(N^O) - IT );
    
    % SUP = O*N^M

    V  = zeros( N, M, O );
    
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
                
                C(:,1) = A(:,1).*P(:,1);
                
    
                ALPHA  = RHO(1,1) / dot( R(:,1), C(:,1) );
    
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
        
        while( sum(sum(sum(V,1),2),3) < IT )
    
            K = floor( ii / IT ) + 1;
    
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
                
                C(:,1) = A(:,1).*P(:,1);
                
    
                ALPHA  = RHO(1,1) / dot( R(:,1), C(:,1) );
    
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

            WALK_(aa,bb,cc,K) = kk; kk = 1;
            
            ii = ii + 1;
        end
    
        [WALK_, L] = sort(WALK_(:,:,:,1), 2); W_ = WALK_(:,:,:,1);
        
        for k = 1:1:size(WALK_,3)
            for i = 1:1:size(WALK_,1)
                
                WALKA(i,k) = mean(WALK_(i,:,k),2);
            end
        end
        
        jj = 1;
        for k = 1:1:size(WALK_,3)
            for i = 1:1:size(WALK_,1)
                
                while( WALK_(i,:,k) )
                
                     [WALKB(i,jj,k), LA(:,1)] = mode(W_(i,:,k),2);
                     
                     jj = jj + 1;
                     
                     for j = 1:1:size(L,1)
                         
                          W_(i,LA(j,1),k) = NaN;
                     end
                     L = 0;
                end 
                jj = 1:
            end
        end
        
        W_ = WALK_; ii = 1;
        for k = 1:1:size(WALK_,3)
            for j = 1:1:size(WALK_,2)
                
                while( WALK_(:,j,k) )
                
                     [WALKB(i,ii,k), LB(:,1)] = mode(W_(:,j,k),1);
                     
                     jj = jj + 1;
                     
                     for i = 1:1:size(L,1)
                         
                          W_(j,LB(i,1),k) = NaN;
                     end
                     L = 0;
                end 
                ii = 1:
            end
        end
    
        EP_MUA = 1;
        EP_MUB = 1;
        
        ii = 1;
        
        for k = 1:1:size(WALK_,3)
            for j = 1:1:size(WALK_,2)   
                for i = 2:1:size(WALK_,1)
                
                    Z_(i-1,j,k,1) = WALKA(i,j,k) - WALKA(i-1,j,k);
                
                    Z_(i-1,j,k,2) = WALKB(j,i,k) - WALKB(j,i-1,k);
                
                    if( ( WALKA(i,j) - WALKA(i-1,j) ) >= EP_MUA ||...
                            WALKB(j,i) - WALKB(j,i-1) ) >= EP_MUB )
                    
                        break;                
                    end
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

                    RA(ii,:,k) = RA(ii,:,k) + Y_(L,:,k); 

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