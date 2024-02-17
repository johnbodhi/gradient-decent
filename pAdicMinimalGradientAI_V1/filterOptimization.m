function [ RA ] = filterOptimization( RA, A, l)
    
    global BINS

    % Allocate state-space to find all minimal distributions...

    % We need to generate all combinations of deleted delta trains 
    % and store them in state-space. Then, find the combinations 
    % of minimal target distributions among every class that are 
    % most efficient with maxmimal f-measures.
        
    S  = permn([0 1], BINS);

    V_ = zeros(size(S,1),size(RA,2),size(RA,1),size(RA,3));

    for k = 1:1:size(V_,4)
        for m = 1:1:size(V_,3)
            for i = 1:1:size(V_,1)
                for j = 1:1:size(V_,2)

                    if ( S(i,j) )

                        V_(i,j,m,k) = RA(m,j,k);             
                    else
                         
                        V_(i,j,m,k) = 0;
                    end

                end
            end
        end
    end

    % We need to reshape for a more efficient mapping...

    V = cat(1,V_(:,:,:,1),V_(:,:,:,2),V_(:,:,:,3));

    % We need to apply a sub-gradient within the SVM to optimize 
    % the filter efficiency.

    % We need to convolve! This is the shave, or keying.
    
    N = size(RA,1);
    
    M = size(V,1);
    
    O = size(RA,3);
    
    IT  = N*M*O; 
    
    B   = zeros(N,M,O);
    
    T = 1.0;
    
    for i = 2:1:size(RA,1)
      
        ii = 1;

        aa = 1; bb = 1; 

        while( sum(sum(sum(B,1),2),3) < IT )
        
            K  = ceil( ii / IT )+1;
            Ka = ceil( aa / size(S,1) )+1;
            Kb = ceil( bb / size(S,1) )+1;

            if( aa <= size(V,1) )
            
                B(aa,1,K) = 1;

                RA(i,2:end,Ka) = V(aa,:,1);

                aa = aa + 1;
            
                [ D, E ] = classifier( A, l );
            
                [ ~, ~, ACC(K), ~ ] = fMeasure( D, E ); 
            
                if( ACC(K) >= T )
            
                    [ ~, M(K) ] = max(ACC);
                
                    X(:,:,K) = RA;         
                end
            elseif( aa > size(V,1) && bb <= size(V,1) )
            
                B(bb,2,K) = 1; B(:,1,K) = 0;
             
                RA(i,2:end,Kb) = V(bb,:,2); 

                bb = bb + 1; aa = 1;
            
                [ D, E ] = classifier( A, l );
            
                [ ~, ~, ACC(K), ~ ] = fMeasure( D, E ); 
            
                if( ACC(K) >= T )
            
                    [ ~, M(K) ] = max(ACC);
                
                    X(:,:,K) = RA;         
                end
            end
            ii = ii + 1;
            
        end
    end
    
    for k = 2:1:size(RA,3)
        
        RA(2:end,2:end,k) = X(:,:,k-1); 
    end
    Q = RA;
    
end
    