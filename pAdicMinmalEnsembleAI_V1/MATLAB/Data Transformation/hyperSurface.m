function [ Y_, Z_ ] = hyperSurface( X_ )

    C_ = zeros(size(X_,1),size(X_,2),size(X_,3));

    D_ = zeros(size(X_,1),size(X_,2),size(X_,3));

    for k = 1:1:size(X_,3)
       for i = 1:1:size(X_,1)

           [C_(i,:,k), D_(i,:,k)] = sort(X_(i,:,k),2,'descend');
       end
    end

    for k = 1:1:size(X_,3)

       Y_(:,:,k)  = C_(:,:,k);

       Z_(:,:,k) = D_(:,:,k);
    end
end