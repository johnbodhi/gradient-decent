function [ RA ] = runningAverage( dataSet, Sup, verObservation )

    global classType classGroups imageLength

    N = size(classType,2)*size(classGroups,2); M = 0;

    if ( ~Sup )
    
        % Unsupervised training...
    
        for k = 1:1:size(dataSet,2)-1
            for i = 1:1:size(dataSet,1)/(imageLength)
    
                A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k),20);
            end
        end
    
        EIGEN_FRAMES = size(A,1);
    
        SEGMENTS     = EIGEN_FRAMES / size(classType,2) + M;
    
        for k = 1:1:size(A,3)      
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B(i,j,k) = norm( ( A(j,:,k) - A(i,:,k) ), Inf );
                end
            end
        end
        
    %     [ B, I ] = sort(B,1);
        [ ~, I ] = sort(B,2);
    
        RA = zeros(size(classType,2),size(A,2),size(A,3)); jj = 1;
    
        for k = 1:1:size(A,3)
            for i = 1:1:size(classType,2)
                
                while ( jj < SEGMENTS )
                
                    RA(i,:,k) = RA(i,:,k) + A(I(i,jj),:,k);
    
                    jj = jj + 1;
                end
                jj = 1;
    
            end
        end
    
        RA = RA ./ SEGMENTS;

    elseif( Sup )
        
        for k = 1:1:size(dataSet,2)-1
            for i = 1:1:size(dataSet,1)/(imageLength)
    
                A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k),20);
            end
        end

        A = cat(2,A,zeros(size(A,1),1,size(A,3)));

        for k = 1:1:size(A,3)
            for i = 1:1:size(A,1)

                A(i,size(A,2),k) = verObservation(i,1);
            end
        end

        for k = 1:1:size(A,3)
            for j = 1:1:size(A,2)
                
                if( A() == )

                    RA(,j,k) = mean(A)

                end
            end
        end

        
    end



end
   
