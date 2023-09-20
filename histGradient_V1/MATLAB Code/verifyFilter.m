function [ Y ] = verifyFilter( X )
    
    Y = krylovSubspace( X );   
    
end