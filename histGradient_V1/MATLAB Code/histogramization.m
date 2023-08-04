function [ Y ] = histogramization( X )

    global imageLength

    for k = 1:1:size(X,2)-1
        for i = 1:1:size(X,1)/(imageLength)

            Y(i,:,k) = hist(X((i-1)*imageLength+1:i*imageLength,k),20);
        end
    end
    
end