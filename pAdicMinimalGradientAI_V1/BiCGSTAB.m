function [ Z_ ] = BiCGSTAB( X_, Y_ )

    global classType classGroups BINS frameLength

    N   = size(Y_,1); 
    
    M   = size(classType,2);

    O   = size(classGroups,2);
    
    CONTAINMENT  = 4; % Average containment limit... (0.25N)
    
    IT  = O*N*M;
    
    V   = zeros( N, M, O );
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    TOL = 1; LIMIT = 1e-3; 

    if( N == 1 )

        X_ = frameLength.*X_; A = Y_;

        for ii = 1:1:size(X_,1)

            R_0(:,1) = X_(ii,:,1) - Y_(1,:,1).*Y_(1,:,1);
            
            R(:,1)   = R_0(:,1);
    
            RHO_0    = dot(R_0(:,1), R_0(:,1)); 
            
            RHO(1,1) = RHO_0;
    
            P(:,1)   = R_0(:,1);          
    
            while( TOL >= LIMIT )
                
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
    
            WALK_(ii,1) = kk; kk = 1;        
        end
        
        [ ~, Z_ ] = min(WALK_);
        
        RA(Z_,:,:) = ( RA(Z_,:,:) + Y_ ) ./ 2;
        
    elseif( N > 1 )

        B(:,1) = frameLength.*X_(1,:,1);
        
        X(:,1) = Y_(1,:,1);
            
        ii = 1;
        
        while( sum(sum(sum(V,1),2),3) < IT )
    
            Ks = floor( ii / IT )+1;
    
            if( aa <= size(V,1) )
    
                V(aa,1,Ks) = 1;
                    
                A(:,1) = Y_(aa,:,Ks);
    
                aa = aa + 1;                          
                
            elseif( aa > size(V,1) && bb <= size(V,1) )
                
                V(bb,2,Ks) = 1; V(:,1,Ks) = 0;
    
                X(:,1) = Y_(bb,:,Ks);
                    
                bb = bb + 1; aa = 1;
    
            elseif( bb > size(V,1) && cc <= CONTAINMENT )
    
                V(cc,3,Ks) = 1; V(:,2,Ks) = 0;
                    
                Y_ = monteCarlo(Y_);
                    
                B(:,1) = frameLength.*mean(Y_(1:cc,:,Ks),1);
                    
                cc = cc + 1; bb = 1;
                
            end
            
            R_0(:,1) = B(:,1) - A(:,1).*X(:,1);
            
            R(:,1)   = R_0(:,1);
    
            RHO_0    = dot(R_0(:,1), R_0(:,1)); 
            
            RHO(1,1) = RHO_0;
    
            P(:,1)   = R_0(:,1);
    
            while( TOL >= LIMIT )
                
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

            WALK_(aa,bb,cc,Ks) = kk; kk = 1;
            
            ii = ii + 1;
        end
        
        for K = 1:1:size(classGroups,2)
        
            W_ = WALK_(:,:,:,K);
        
            jj = 1;
            for k = 1:1:size(WALK_,3)
                for i = 1:1:size(WALK_,1)
                
                    while( WALK_(i,:,k) )                                                   
                         
                         [WALKA(i,jj,k), LA(i,jj,k)] = ...
                             find(W_(i,:,k) == min(W_(i,:,k),2));
                     
                         jj = jj + 1;
                      
                         for j = 1:1:size(LA,2)
                         
                             W_(i,LA(i,j,k),k) = NaN;
                         end
                         L = 0;
                     
                    end 
                    jj = 1;
                        
                end
            end
         
            W_ = WALK_(:,:,:,K); 
            
            ii = 1;
            for k = 1:1:size(WALK_,3)
                for j = 1:1:size(WALK_,2)
                
                    while( WALK_(:,j,k) )
                
                         [WALKB(ii,j,k), LB(j,ii,k)] = ...
                             find(W_(:,j,k) == min(W_(:,j,k),1));
                     
                         ii = ii + 1;
                     
                         for i = 1:1:size(L,1)
                         
                             W_(j,LB(i,j,k),k) = NaN;
                         end
                         L = 0;
                     
                    end 
                    ii = 1;
                
                end
            end
            
            EP_MU = 0.9;
        
            RA = zeros(M,BINS,O);
        
            ii = 1; ll = 1;
            for k = 1:1:size(WALK_,3)
                for i = 2:1:size(WALK_,1)
                        
                    Z_(i-1,j,k,1) = WALKA(j,i,k) - WALKA(j,i-1,k);
                
                    Z_(i-1,j,k,2) = WALKB(i,j,k) - WALKB(i-1,j,k);
                        
                    E(1,1) = Z_(i-1,j,k,1)/Z_(i,j,k,1);
                    
                    E(1,2) = Z_(i-1,j,k,2)/Z_(i,j,k,2);
                    
                    
                    if( E(1,1) >= EP_MU && E(1,2) >= EP_MU )
                    
                        RA(k,:,K) = RA(k,:,K) + Y_(LA(1,i-1,k),:,K);
                        
                        RB(k,:,K) = RB(k,:,K) + Y_(LB(1,i-1,k),:,K);
                       
                        ll = ll + 1;
                    end  
                 end
            
                 RAA(k,:,K) = RAA(k,:,K)/ll;
                
                 RAB(k,:,K) = RAB(k,:,K)/ll;
                
                 ll = 1;
            end
         end
    end 
    Z_ = ( RAA + RAB ) ./ 2; 
    
end