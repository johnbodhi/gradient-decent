function [ RA ] = SVM( dataSet, N, Observation )

    global RA classType Supervision Randomized Train

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
                
                    RA(i,:,k) = RA(i,:,k) + A(I(ii,jj,i),1:end-1,k);
    
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

        for k = 1:1:size(A,3)      
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B(i,j,k) = norm( ( A(j,1:end-1,k) - A(i,1:end-1,k) ), Inf );
                    I(i,1,1) = i;
                end
            end
        end  

        [ B, I ] = sort(B,2);

        [ ~, W ] = combinations( I, SEGMENTS ); % Permutation windows about A.

        ii = 1; jj = 1; NN = SEGMENTS*size(W,1);
        for k = 1:1:size(A,3)            
            for i = 1:1:size(RA,1)
                
                while ( ii <= size(W,1) )
                    for kk = 1:1:size(W,3)
                    
                        % Convergence criterion to avoid a bricked RA...

                        while ( jj <= SEGMENTS )
                        
                            RA(i,:,k) = RA(i,:,k) + A(W(ii,jj,kk),1:end-1,k);
            
                            jj = jj + 1;                      
                        end
                        jj = 1;
                        
                    end
                    ii = ii + 1;
                end
                ii = 1;
            end                 
        end
    
        RA = RA ./ NN;
        
        % The per element histogram magnitude ratios appear to be identical 
        % bewteen the supervised and unsupervised cases for each group.
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
    
                        RA(jj,:,k) = RA(jj,:,k) + A(i,1:end-1,k); 

                        ii = ii + 1; kk = kk + 1;
                    end                        
                end   
                jj = jj + 1; ii = 1;
            end   
            jj = 1; kk = 1;
        end
 
        RA = RA ./ SEGMENTS; 
    end

    plotDist(RA);

end   