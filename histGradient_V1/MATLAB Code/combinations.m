function [ F ] = combinations( X, Y )

    global classType

    F(:,:,1) = uniqueperms(X(1:Y,1)); % Shout out...

    for k = 2:1:size(classType,2)
        for j = 1:1:size(F,2)
            for i = 1:1:size(F,1)
                
                F(i,j,k) = F(i,j,k-1)+Y;
            end
        end 
    end
    
end