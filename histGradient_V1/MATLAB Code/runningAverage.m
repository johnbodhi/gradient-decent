function [ RA ] = runningAverage( dataSet, Observation )

    global classType imageLength RA BINS Supervision Train

    if ( ~Supervision && Train )
    
    % Unsupervised training...
    
        for k = 1:1:size(dataSet,2)-1
            for i = 1:1:size(dataSet,1)/(imageLength)
    
                A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k),BINS);
            end
        end
    
        EIGEN_FRAMES = size(A,1);
    
        SEGMENTS     = EIGEN_FRAMES / size(classType,2);
    
        for k = 1:1:size(A,3)      
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B(i,j,k) = norm( ( A(j,:,k) - A(i,:,k) ), Inf );
                    I(i,1,1) = i;
                end
            end
        end        
        % [ ~, I ] = sort(B,2);

        I = combinations( I, SEGMENTS );
        
        RA = zeros(size(classType,2),size(A,2),size(A,3)); 
        
        ii = 1; jj = 1;    
        for k = 1:1:size(A,3)
            for i = 1:1:size(classType,2)  
                
                while ( jj < SEGMENTS )
                
                    RA(i,:,k) = RA(i,:,k) + A(I(ii,jj,i),:,k);
    
                    if( ii >= size(I,1) )

                        ii = 0; jj = jj + 1;
                    end
                    ii = ii + 1;
                end
                jj = 1;        
            end
        end
    
        RA = RA ./ size(I,1); 
        
        % The per element histogram magnitude ratios appear to be identical 
        % bewteen the supervised and unsupervised cases for each group.

    else

        if ( Train )
        
            for k = 1:1:size(dataSet,2)-1
                for i = 1:1:size(dataSet,1)/(imageLength)
        
                    A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k),BINS);
                end
            end

            A = cat(2,A,zeros(size(A,1),1,size(A,3)));

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

        if ( Train )

            RA = RA ./ SEGMENTS; 
        else

            RA = RA ./ size(RA,1); 
        end
    end

end   
