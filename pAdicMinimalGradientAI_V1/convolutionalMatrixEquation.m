function [ Z_ ] = convolutionalMatrixEquation( Y_ )

    global classType BINS frameLength TOL LIMIT
    
        
    N   = 0.5*size(Y_,1); 
    
    M   = size(classType,2)+1;
    
    
    CONTAINMENT = 0.5*N; % Objective function containment limit.
    
    IT  = N*M-CONTAINMENT;
    
    V   = zeros( N, M );
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    TOL = 1; LIMIT = 1e-2;
    
    B(:,1) = frameLength.*Y_(1,:);
        
    X(:,1) = Y_(1,:);
    
    
    F = sum(sum(V,1),2);
    
    while( F < IT )
    
        if( aa <= N )
    
            V(aa,1) = 1
                    
            A(:,1) = Y_(aa,:);
                
            WALK_(aa,bb,cc) = BiCGSTAB_(A,X,B);
                
            aa = aa + 1;
        
        elseif( aa >= N && bb <= N )
                
            V(bb,2) = 1; V(:,1) = 0;
    
            X(:,1) = Y_(bb,:);
                
            WALK_(aa,bb,cc) = BiCGSTAB_(A,X,B);
                
            bb = bb + 1; aa = 1;
        
        elseif( bb >= N && cc <= CONTAINMENT )
                
            V(cc,3) = 1; V(:,2) = 0;
                    
            Y_ = monteCarlo(Y_);
                    
            B(:,1) = frameLength.*mean(Y_(1:cc,:),1);
                
            WALK_(aa,bb,cc) = BiCGSTAB_(A,X,B);
                
            cc = cc + 1; bb = 1;
            
        end
        F = sum(sum(V,1),2)
        
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
        
    RA = zeros(M,BINS); RB = zeros(M,BINS);
        
    ll = 1;
    for k = 1:1:size(WALK_,3)
        for j = 1:1:size(WALK_,2)
            for i = 2:1:size(WALK_,1)
                        
                Z_(i-1,j,k,1) = WALKA(j,i,k) - WALKA(j,i-1,k);
            
                Z_(i-1,j,k,2) = WALKB(i,j,k) - WALKB(i-1,j,k);
            
                        
                E(1,1) = Z_(i-1,j,k,1) / Z_(i,j,k,1);
                    
                E(1,2) = Z_(i-1,j,k,2) / Z_(i,j,k,2);
            
                    
                if( E(1,1) >= EP_MU && E(1,2) >= EP_MU )
            
                    RA(k,:) = RA(k,:) + Y_(LA(1,i-1,k),:);
               
                    RB(k,:) = RB(k,:) + Y_(LB(1,i-1,k),:);
             
                    ll = ll + 1;
                end
            end
        end     
        RA(k,:) = RA(k,:) / ll;
           
        RB(k,:) = RB(k,:) / ll;
                
        ll = 1;
    end
    Z_ = ( RA + RB ) ./ 2;
    
end