function [ Y ] = randomizeClassType( X )

    Z = 1;
    randArr(1,1) = randi( [ 1, size( X, 2 ) ] );
    for j = 2:1:size(X,2)
        while( Z )
            randArr( 1, j ) = randi( [ 1, size( X, 2 ) ] );
    
            Z = size(find(randArr(1,1:j-1) == randArr(1,j)),2);
        end
        Z = 1;
        A = size(find(randArr),2);
        if( ~A )
            break;
        end
    end
    
    Y = randArr;
end
