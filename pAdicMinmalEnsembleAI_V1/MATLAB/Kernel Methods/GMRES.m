function [ Z_ ] = GMRES( A_ )

    global classType frameLength TOL LIMIT M
    
        
    N   = size(A_,1); 
    
    M   = size(classType,2)+1;
    
    
    CONTAINMENT = 1.0*N; % Objective function containment limit.
    
    IT  = N^CONTAINMENT;
    
    V   = zeros( N, M );

    ii = 1; jj = 0; kk = 1; ll = 0;
    
    aa = 1; bb = 1;
    

    TOL = 1; LIMIT = 1e-2;
    
    B(:,1) = frameLength.*mean(A_(1,:),1);
        
    % X(:,1) = Y_(1,:);
    
    
    F = sum(sum(V,1),2);
    
    while( F < IT )

        if( aa <= N )
    
            V(aa,1) = 1;

            A_ = monteCarlo(A_);
                    
            A(:,1) = A_(aa,:);

            B(:,1) = frameLength.*mean(A_(1:bb,:),1);
                
            WALK_(:,:) = BiCGSTAB_(A,[],B);
                
            aa = aa + 1; 
        
        elseif( aa >= N && bb <= CONTAINMENT )

            V(bb,2) = 1; V(:,1) = 0;
                    
            A_ = monteCarlo(A_);
                    
            B(:,1) = frameLength.*mean(A_(1:bb,:),1);
                
            WALK_(:,:) = BiCGSTAB_(A,[],B); 

            bb = bb + 1;

        end
        F = sum(sum(V,1),2); aa = 1; bb = 1;
        
    end
            
    % We need to sort the minimum walk accumulations in the 
    % form of a chart, and a transpose chart.
        
    W_ = WALK_( :, : ); 
    
    [W,~] = sort(W_,1);
    
    [ ~, Z_ ] = min(W(1,:));
   
    % We can use GMRES to return a decision, or a filter.

end