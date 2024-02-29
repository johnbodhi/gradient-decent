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
    
            V(aa,1) = 1;
                    
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
        F = sum(sum(V,1),2);
        
    end
            
    % We need to sort the minimum walk accumulations in the 
    % form of a chart, and a transpose chart.
        
    W_ = WALK_(:,:,:); 
    
    [WA,LA] = sort(W_,2);
    
    [WB,LB] = sort(W_,1);
        
    % We need to separate the walk accumulation deltas from eachother
    % according to an average edge amplitude. This is the discrimination
    % value for a class.
            
    EP_MU = 0.1; % We can  step through an epsilon vector as well once 
                 % we have experience...
        
    RA = zeros(M,BINS); RB = zeros(M,BINS);
    
    Z_ = zeros(size(W_,1),size(W_,2),size(W_,3),2);
        
    ll = 0;
    for k = 1:1:size(WALK_,3)
        for j = 1:1:size(WALK_,2)
            for i = 2:1:size(WALK_,1)
                        
                Z_(i-1,j,k,1) = WA(j,i,k) - WA(j,i-1,k);
            
                Z_(i-1,j,k,2) = WB(i,j,k) - WB(i-1,j,k);
            
                        
                E(1,1) = Z_(i-1,j,k,1) / Z_(i,j,k,1);
                    
                E(1,2) = Z_(i-1,j,k,2) / Z_(i,j,k,2);
            
                    
                if( E(1,1) <= EP_MU && E(1,2) <= EP_MU )
            
                    RA(k,:) = RA(k,:) + Y_(LA(j,i-1,k),:);
               
                    RB(k,:) = RB(k,:) + Y_(LB(i-1,j,k),:);
             
                    ll = ll + 1;
                end
                
            end
        end
        
        if( ll )
            
            RA(k,:) = RA(k,:) / ll;
           
            RB(k,:) = RB(k,:) / ll;
            
            ll = 0;
        end
        
    end
    Z_ = ( RA + RB ) ./ 2;
    
end
