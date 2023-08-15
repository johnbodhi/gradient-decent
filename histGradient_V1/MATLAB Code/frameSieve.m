function [ Y ] = frameSieve( S )

    global RA

    % S = sort(S,"ascend");

    for k = 1:1:size(RA,3)
        for i = 1:1:size(RA,1)
                
            Y(i,:,k) = S(1,1:end-1,k); 
            L(i,1,k) = i;
        end
    end

    Y = cat(2,Y,L);    
end