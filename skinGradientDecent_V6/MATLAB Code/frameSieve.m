function [ Y ] = frameSieve( S )

    global A RA X imageLength

    A = size(RA,1); M = imageLength;

    Y = S; Y = sort(Y);

    for i = 1:1:A-1
        
        Y = cat( 1, Y, Y );
    end

    cc = 1;
    for i = 1:M:size(Y,1)

        if( i <= imageLength )

            Y((cc-1)*M+1:cc*M,4) = X(1,1); 
        elseif( i > imageLength )

            Y((cc-1)*M+1:cc*M,4) = X(1,2); 
        end
        
        cc = cc + 1;
    end    
end