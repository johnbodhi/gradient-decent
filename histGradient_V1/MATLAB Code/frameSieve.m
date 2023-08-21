function [ Y ] = frameSieve( S )

    global RA Optimized

    % S = sort(S,"ascend");

    SPLITS = 8;

    if( Optimized )

        for k = 1:1:size(RA,4)
           for i = 1:1:SPLITS
                    
                Y(i,:,k) = S(1,1:end-1,k); 
                L(i,1,k) = i;
            end
        end
    else

        for k = 1:1:size(RA,3)
           for i = 1:1:SPLITS
                    
                Y(i,:,k) = S(1,1:end-1,k); 
                L(i,1,k) = i;
            end
        end
    end

    Y = cat(2,Y,L);    
end