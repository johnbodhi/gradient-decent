function [ RA ] = SVM( dataSet, N, Observation )

    global RA classType Supervision Randomized Optimized Train

    A = histogramization( dataSet, N, Observation ); 

    EIGEN_FRAMES = size(A,1);
    
    SEGMENTS     = EIGEN_FRAMES / size(classType,2);

    if ( ~Supervision && Train && ~Randomized )
   
        for k = 1:1:size(A,3)      
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B(i,j,k) = norm( ( A(j,:,k) - A(i,:,k) ), Inf );
                    I(i,1,1) = i;
                end
            end
        end        
        % [ ~, I ] = sort(B,2);

        [ I, ~ ] = combinations( I, SEGMENTS ); % Permutation windows about A.
        
        ii = 1; jj = 1; NN = SEGMENTS*size(I,1);   
        for k = 1:1:size(A,3)
            for i = 1:1:size(classType,2)
                
                while ( jj <= SEGMENTS )
                
                    RA(i+1,2:end,k+1) = RA(i+1,2:end,k+1) + A(I(ii,jj,i),1:end-1,k);
    
                    if( ii >= size(I,1) )

                        ii = 0; jj = jj + 1;
                    end
                    ii = ii + 1;
                end
                jj = 1;        
            end
        end

        RA = RA ./ NN; 

    elseif( ~Supervision && Train && Randomized )
    
    
    
    
        
        
    else

        if ( Train )

            for k = 1:1:size(A,3)
                for i = 1:1:size(A,1)
    
                    A(i,size(A,2),k) = Observation(i,1);
                end
            end 
        else
            
            A = dataSet;
        end

        EIGEN_FRAMES = size(A,1);
    
        SEGMENTS     = EIGEN_FRAMES / size(classType,2);
 
        ii = 1; jj = 1; kk = 1;
        for k = 1:1:size(A,3)
            while ( jj <= size(classType,2) )
                for i = 1:1:size(A,1)

                    if( A(i,size(A,2),k) == Observation(kk,1) && ii <= SEGMENTS )
    
                        RA(jj+1,2:end,k+1) = RA(jj+1,2:end,k+1) + A(i,1:end-1,k); 

                        ii = ii + 1; kk = kk + 1;
                    end                        
                end   
                jj = jj + 1; ii = 1;
            end   
            jj = 1; kk = 1;
        end
 
        RA = RA ./ SEGMENTS; 

    end

    if ( Optimized )        
    
        [ RA ] = filterOptimization( dataSet, N, Observation );
    end

    % plotDist();
    % dope();
    % plotDist();
end