function [ RA ] = runningAverage( dataSet )

    global classType classGroups imageLength

    N = size(classType,2)*size(classGroups,2);
    
    % Unsupervised training.

    for k = 1:1:size(dataSet,2)-1
        for i = 1:1:size(dataSet,1)/(imageLength)

            A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k));
        end
    end

    EIGEN_FRAMES = size(A,1);

    SEGMENTS     = EIGEN_FRAMES / size(classType,2);

    for k = 1:1:size(A,3)      
        for j = 1:1:size(A,1)
            for i = 1:1:size(A,1)

                B(i,j,k) = norm( ( A(j,:,k) - A(i,:,k) ), Inf );
            end
        end
    end
    
    [ ~, I ] = sort(B,2);

    RA = zeros(8,10,3); jj = 1;

    for k = 1:1:size(A,3)
        for i = 1:1:size(classType,2)
            
            while ( jj < SEGMENTS )
            
                RA(i,:,k) = RA(i,:,k) + A(I(i,jj),:,k);

                jj = jj + 1;
            end
            jj = 1;
            
        end
    end

end
   
