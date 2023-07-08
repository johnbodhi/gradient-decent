function [ Y ] = frameSieve( X )

    global A RA imageLength

    A = size(RA,1); M = imageLength;

    Y = X;

    for i = 1:1:A-1
        
        Y = cat( 1, Y, X );
    end

    uu = 1;
    for i = 1:M:size(Y,1)

        Y((uu-1)*M+1:uu*M,4) = uu; 
        
        uu = uu + 1;
    end    
end