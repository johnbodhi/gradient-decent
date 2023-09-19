function [ Y ] = filterUpdate( X )

    Y = krylovSubspace( X );   

end

