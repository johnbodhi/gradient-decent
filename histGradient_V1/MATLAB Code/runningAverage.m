function [ RA ] = runningAverage( dataSet )

    global imageLength
    
    % Unsupervised training.

    for k = 1:1:size(dataSet,2)
        for i = 1:imageLength:size(dataSet,1)

            A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k));
        end
    end

    while ( A ~= 0 )

        for k = 1:1:size(A,3)
            for j = 1:1:size(dataSet,2)
                for ii = 1:1:size(A,1)
                    for i = 1:1:size(A,1)
        
                        B(i,ii,j,k) = norm( ( A(ii,:,k) - A(i,:,k) ), 2 );
                    end
                end
            end
        end
    
        E = zeros(size(A,1),size(A,2),size(A,3));
    
        for k = 1:1:size(A,3)
            for j = 1:1:size(dataSet,2)
                for ii = 1:1:size(A,1)
                    for i = 1:1:size(A,1)
        
                        if( i ~= ii )
        
                            [ ~, E(i,j,k) ] = min(B(i,:,j,k),[],2); 
                        end
                        
                    end
                end
            end
            V(j,k) = min(E(:,j,k),[],1);
        end
    
        for k = 1:1:size(A,3)
            for j = 1:1:size(dataSet,2)
                for i = 1:1:size(A,1)
        
                    F(i,:,k) = ( A(i,:,k) + A(V(j,k),:,k) ) / 2; 
                    
                    A(V(j,k),:,k) = 0; V(j,k) = 0;   
                end
            end
        end
    
        % We need a recusion on F to reduce the number of histograms to the
        % number of classes being evaluated...
    
    
        ll = 0; % Resize accumulator.
    
        for k = 1:1:size(dataSet,2)
            for ii = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B_(i,ii,k) = norm( ( F(ii,:,k) - F(i,:,k) ), 2 );
                end
            end
        end
        
        for k = 1:1:size(dataSet,2)
            for ii = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    if( i ~= ii )
    
                        [ ~, E_(i,k) ] = min(B_(i,:,k),[],2); 
                    end
                    
                end
            end
            V_(k) = min(E(:,k),[],1);
        end
    
         for k = 1:1:size(F,3)
        
            F(i,:,k) = ( F(i,:,k) + F_(V_(k),:,k) ) / 2; 
            
            A(V_(k),:,k) = 0; V_(k) = 0; 
         end
    
         for k = 1:1:size(F,3)
    
            F(V_(k),:,k) = 0;
         end
    
         ll = ll + 1;
         for k = 1:1:size(F,3)
             for j = 1:1:size(F,2)
                for i = 1:1:size(F,1)-ll
    
                    if ( F(i,j,k) ~= 0 )
                        X(i,j,k) = F(i,j,k);
                    end
                end
             end
         end
         F = 0; F = X;




         


    end
    
    RA = F;
end
   
