function [ Z_ ] = GMRES( Y_ )

    global classType BINS frameLength TOL LIMIT M
    
        
    N   = 0.5*size(Y_,1); 
    
    M   = size(classType,2)+1;
    
    
    CONTAINMENT = 0.5*N; % Objective function containment limit.
    
    IT  = N^(M-CONTAINMENT);
    
    V   = zeros( N, M );
    
    
    ii = 1; kk = 1;

    aa = 1; bb = 1; cc = 1;

    TOL = 1; LIMIT = 1e-2;
    
    
    B(:,1) = frameLength.*Y_(1,:);
        
    % X(:,1) = Y_(1,:);
    
    
    F = sum(sum(V,1),2);
    
    while( F < IT )
    
        if( aa <= N )
    
            V(aa,1) = 1;
                    
            A(:,1) = Y_(aa,:);
                
            WALK_(aa,bb) = BiCGSTAB_(A,[],B);
                
            aa = aa + 1;
        
        %{ 
            
            elseif( aa >= N && bb <= N )
                
            V(bb,2) = 1; V(:,1) = 0;
    
            X(:,1) = Y_(bb,:);
                
            WALK_(aa,bb,cc) = BiCGSTAB_(A,X,B);
                
            bb = bb + 1; aa = 1;
            
        %}
        
        elseif( aa >= N && bb <= CONTAINMENT )
                
            V(bb,2) = 1; V(:,1) = 0;
                    
            Y_ = monteCarlo(Y_);
                    
            B(:,1) = frameLength.*mean(Y_(1:cc,:),1);
                
            WALK_(aa,bb) = BiCGSTAB_(A,[],B);
                
            bb = bb + 1; aa = 1;
            
        end
        F = sum(sum(V,1),2);
        
    end
            
    % We need to sort the minimum walk accumulations in the 
    % form of a chart, and a transpose chart.
        
    W_ = WALK_( :, :, : ); 
    
    [ Z_ ] = averageMapping( W_ );
    
end