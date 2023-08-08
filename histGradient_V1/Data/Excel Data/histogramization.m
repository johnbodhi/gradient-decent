function [ F ] = histogramization( X, Y, Z )

    global imageLength BINS Randomized


    if ( Randomized )

        [ X, L ] = randomizeAll( X, Y ); % Randomize all photos.

        for k = 1:1:size(X,2)-1
            for i = 1:1:size(X,1)/(imageLength)
    
                F(i,:,k) = hist(X((i-1)*imageLength+1:i*imageLength,k),BINS);
            end
        end

        F = cat(2,F,zeros(size(F,1),1,size(F,3)));

    else

        for k = 1:1:size(X,2)-1
            for i = 1:1:size(X,1)/(imageLength)
    
                F(i,:,k) = hist(X((i-1)*imageLength+1:i*imageLength,k),BINS);
            end
        end
    
        F = cat(2,F,zeros(size(F,1),1,size(F,3)));
    end

    if( Randomized )

        for k = 1:1:size(F,3)
            for j = 1:1:size(F,1)
    
                F(j,size(F,2),k) = L(j,1);
            end
        end

    else        

        for k = 1:1:size(F,3)
            for i = 1:1:size(F,1)
    
                F(i,size(F,2),k) = Z(i,1);
            end
        end
    end
end