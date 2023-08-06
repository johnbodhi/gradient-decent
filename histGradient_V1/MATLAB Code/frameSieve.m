function [ Y ] = frameSieve( S )

    global RA

    Y = S(:,1:end-1); Y = sort( Y,'descend' );

    for i = 1:1:size(RA,1)
        for j = 1:1:size(RA,2)-1
            
            Y(i,j) = S(1,j);      
        end      
        L(i,1) = i;
    end

    Y = cat(2,Y,L);    
end