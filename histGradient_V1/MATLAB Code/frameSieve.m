function [ Y ] = frameSieve( S )

    global RA

    for i = 1:1:size(RA,1)
            
        Y(i,:) = S(1,1:end-1); L(i,1) = i;
    end

    Y = cat(2,Y,L);    
end