function [ Y ] = addNoise( X )

    Y_ = reconstruct( X ); R = 0;
    
    for i= 1:1:size(Y,1)
        
        while( A_(i,1) - R < DATARANGE-1 || A_(i,1) + R > DATARANGE-1 )
            
            R = randi([-255 255],1);
        end
        
        Y(i,1) = Y_(i,1) + R; R = 0;
    end
        
    % AWGN = rand([size(A_,1),size(A_,2)]);
        
    % A_ = A_ + AWGN;
end