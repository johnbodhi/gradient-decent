function [ Y ] = initializeFilter( X )
    
    Y = krylovSubspace( X );   
    
end