function [ F ] = histogramization( X, Y )

    global imageLength

    for k = 1:1:size(X,2)-1
        for i = 1:1:size(X,1)/(imageLength)

            F(i,:,k) = hist(X((i-1)*imageLength+1:i*imageLength,k),20);
        end
    end

    F = cat(2,F,zeros(size(F,1),1,size(F,3)));

    for k = 1:1:size(F,3)
        for i = 1:1:size(F,1)

            F(i,size(F,2),k) = Y(i,1);
        end
    end
    
end