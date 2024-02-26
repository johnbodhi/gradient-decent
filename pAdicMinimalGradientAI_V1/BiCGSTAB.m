function [ Z_ ] = BiCGSTAB( Y_ )

    global classType classGroups BINS frameLength

    N   = size(Y_,1); 
    
    M   = size(classType,2)+1;

    O   = size(classGroups,2);
    
    CONTAINMENT = 0.25*N; % Objective function containment limit.
    
    IT  = floor(N^((M-1)+CONTAINMENT/N));
    
    V   = zeros( N, M );
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    TOL = 1; LIMIT = 1e-2;

    if( N == 1 )
        
        %{
        
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
        
        %}
        
    elseif( N > 1 )

        B(:,1) = frameLength.*Y_(1,:,1);
        
        X(:,1) = Y_(1,:,1);
        
        while( sum(sum(V,1),2) <= IT )
    
            if( aa <= size(V,1) )
    
                V(aa,1) = 1
                    
                A(:,1) = Y_(aa,:);
                
                aa = aa + 1;            
                
            elseif( aa > size(V,1) && bb <= size(V,1) )
                
                V(bb,2) = 1; V(:,1) = 0;
    
                X(:,1) = Y_(bb,:);
                    
                bb = bb + 1; aa = 1;
    
            elseif( bb > size(V,1) && cc <= CONTAINMENT )
    
                V(cc,3) = 1; V(:,2) = 0;
                    
                Y_ = monteCarlo(Y_);
                    
                B(:,1) = frameLength.*mean(Y_(1:cc,:),1);
                    
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
    
                S(:,1) = R(:,1) - ALPHA.*C(:,1);            
    
                T(:,1) = A(:,1).*S(:,1);
    
    
                OMEGA  = dot( T(:,1), S(:,1) ) / dot( T(:,1), T(:,1) );
    
                X(:,1) = H(:,1) + OMEGA(:,1).*S(:,1);
    
                R(:,1) = S(:,1) - OMEGA(:,1).*T(:,1);
    
    
                RHO(1,2) = dot( R_0(:,1), R(:,1) );
    
    
                TOL      = abs( RHO(1,2) / RHO_0 );
    
    
                BETA     = ( RHO(1,2) / RHO(1,1) ) * ( ALPHA / OMEGA );
    
                P(:,1)   = R(:,1) + BETA.*( P(:,1) - OMEGA.*C(:,1) );
    
                RHO(1,1) = RHO(1,2);
                    
                kk = kk + 1;
            end 
            WALK_(aa,bb,cc) = kk; kk = 1; 
                
            TOL = 1;
        end
    end
            
    % We need to sort the minimum walk accumulations in the 
    % form of a chart, and a transpose chart.
        
    W_ = WALK_(:,:,:);
        
    jj = 1;
    for k = 1:1:size(WALK_,3)
        for i = 1:1:size(WALK_,1)
                
             while( W_(i,:,k) )                                                   
                         
                 [ WA, LA ] = find(W_(i,:,k) == min(W_(i,:,k),2));
                         
                 WALKA(i,jj,k) = size(find(WA),2);
                     
                 jj = jj + 1;
                      
                 for j = 1:1:size(LA,2)
                         
                     W_(i,LA(1,j),k) = NaN;
                 end
                 LA = 0;
             end 
             jj = 1;
        end
    end
         
    W_ = WALK_(:,:,:); 
            
    ii = 1;
    for k = 1:1:size(WALK_,3)
        for j = 1:1:size(WALK_,2)
                
            while( W_(:,j,k) )
                
               [ WB, LB ] = find(W_(:,j,k) == min(W_(:,j,k),1));
                         
               WALKB(i,ii,k) = size(find(WB),2);
                     
               ii = ii + 1;
                     
               for i = 1:1:size(LB,2)
                         
                   W_(j,LB(1,i),k) = NaN;
               end
               LB = 0;
            end 
            ii = 1;
        end
    end
            
    % We need to separate the walk accumulation deltas from eachother
    % according to an average edge amplitude. This is the discrimination
    % value for a class.
            
    EP_MU = 0.9; % We can  step through an epsilon vector as well once 
                 % we have experience...
        
    RA = zeros(M,BINS);
        
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
            
        RAA(k,:) = RAA(k,:) / ll;
           
        RAB(k,:) = RAB(k,:) / ll;
                
        ll = 1;
    end
    Z_ = ( RAA + RAB ) ./ 2;
end
     
    
