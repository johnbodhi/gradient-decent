function [ Y ] = krylovSubspace( X )
    
    Y = BiCGSTAB( X );    
    
end