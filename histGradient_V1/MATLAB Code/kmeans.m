function [ Y ] = kmeans( S, V )

global imageLength RA W

    ii = 1; jj = 1;  

    W = [ 0 0 1e0 ]; 

    while( jj < imageLength ) % Modified convergence criterion.
            
        for k = 1:1:size(RA,3)
            for j = 1:1:size(S,2)-1
                for i = 1:1:size(S,1)
            
                    D(i,j,k) = norm((S(i,j) - RA(V(i,1),j,k)), Inf);        
           
                end
            end

            for ii = 1:1:size(D,1)

                Ci(ii,k) = W(1,k) * mean(D(ii,:,k)); 
                
                Cj(ii,k) = W(1,k) * mean(D(ii,:,k)); 
                
                Ck(ii,k) = W(1,k) * mean(D(ii,:,k)); 
            end

            Cn(:,:,k) = [ Ci(:,k) Cj(:,k) Ck(:,k) ];            
        end       
        
        for k = 1:1:size(RA,3)
            for j = 1:1:size(S,2)-1
                for i = 1:1:size(S,1)
    
                    Y(i,j,k) = norm((S(i,j) - Cn(i,k)), Inf); 
                end
            end   
        end

        ii = ii + 1; jj = jj + 1;
    end

    Y = cat(2,Y,S(:,size(S,2),:));
end