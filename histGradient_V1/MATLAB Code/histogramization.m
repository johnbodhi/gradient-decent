function [ Y ] = histogramization( X )

    global imageLength

    for k = 1:1:size(dataSet,2)-1
        for i = 1:imageLength:size(dataSet,1)
    
            X_(i,:,k) = hist(X((i-1)*imageLength+1:i*imageLength,k));
        end
    end
    
    Y = X_(:,:,1);
    for k = 2:1:size(A,3)
            
        Y = cat( 1, Y, X_(:,:,k) );
    end
    
end