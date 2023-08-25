function [ Y ] = frameSieve( S )

    global RA

    SPLITS = 8;

    Y = zeros(size(RA,1),size(RA,2),size(RA,3));

    for k = 1:1:size(RA,3)-1
       for i = 1:1:SPLITS
                
            Y(i+1,2:end,k+1) = S(1,1:end-1,k); 
            L(i+1,1,k+1) = i;
        end
    end

    Y = cat(2,Y,L);    
end