function [ RA ] = buildManifold( A, l )

    global classType classGroups BINS RA BPI ii ij KFLAG
    
    KFLAG = 1;
    
    for j = 1:1:size(classType,2)
       for k = 1:1:size(classGroups,2)
           
            TL_(:,j,k) = (1:1:size(A,1)); % Sample tracking labels.
        end
    end
    
    F   = A;

    N   = size(A,1);
    
    M   = size(classType,2); 

    O   = size(classGroups,2);
    
    V   = zeros(N-M,1); % Stencil limit...
    
    IT  = O*N*M; 
    
    B   = zeros(N,M,O);
    
    X   = zeros(N,BINS+1,O,1);

    % Convoltution with a sub-gradient!!

    ii = 1;

    aa = 1; bb = 1; 
    
    T  = 1.0;
    
    MM = zeros(1,size(classGroups,2));

    while( sum(sum(sum(B,1),2),3) < IT )

        K = ceil( ii / IT ) + 1;

        if( aa <= size(V,1) )

            B(aa,1,K) = 1;

            for q = 1:1:size(A,1)
                
                % We can utilize monte carlo methods in the convergence.
                
                % A = cat(2,A,TL_);
                
                % F = monteCarlo(A); 
                
                % TL_(:,1,:) = F(:,end,:); F = F(:,1:end-1,:);
                
                % We need to track the samples used in the 
                % construction of the filter to enssure they
                % are unique or exclusive, in that every sample is
                % used in thefilters construction, and used no more 
                % than once. To accomplish this we
                % append temporary labels to the sample space, and 
                % and store them for the convergence criterion.
                
                RA(BPI(1,1),2:end,K) = mean(F(1:aa,:,K-1),1);
    
                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC(K), ~ ] = fMeasure( D, E ); 

                if( ACC(K) >= T )
                    
                    [ ~, MM(K)] = max(ACC);
                
                    X(:,1:end-1,K) = RA; 
                    
                    X(:,end,K) = TL_(:,1,K);  
                end
                F  = circshift(A,1); TL_ = circshift(TL_,1);
                
                l  = circshift(l,1);
            end
            aa = aa + 1;

        elseif( aa > size(V,1) && bb <= size(V,1) )

            B(bb,2,K) = 1; B(:,1,K) = 0;

            for q = 1:1:size(A,1)
                
                % A = cat(2,A,TL_);
                
                % F = monteCarlo(A); 
                
                % TL_(:,2:) = F(:,end,:); F = F(:,1:end-1,:);
                 
                RA(BPI(1,2),2:end,K) = mean(F(1:bb,:,K-1),1);
               
                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC(K), ~ ] = fMeasure( D, E ); 

                if( ACC(K) >= T )
                    
                    [ ~, MM(K)] = max(ACC);
                
                    X(:,1:end-1,K) = RA; 
                    
                    X(:,end,K) = TL_(:,2,K);      
                end
                F  = circshift(A,1); TL_ = circshift(TL_,1);
                
                l  = circshift(l,1);
                
            end     
            bb = bb + 1; aa = 1;
            
        end  
        ii = ii + 1;
        
    end
    
    for k = 1:1:size(X,3)
        RA = X(:,:,k); 
    end
    Q = RA;
    
    ii = ii + 2; ij = ij + 2;
end