function [ RA ] = runningAverage( dataSet )

    global classType classGroups imageLength

    N = size(classType,2)*size(classGroups,2);

    Np = 40; Mp = 40;
    
    % Unsupervised training.

    for k = 1:1:size(dataSet,2)-1
        for i = 1:1:size(dataSet,1)/(Np*Mp)

            A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k));
        end
    end

    while ( size(A,1) > N )

        for k = 1:1:size(A,3)        
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B(i,j,k) = norm( ( A(j,:,k) - A(i,:,k) ), 2 );
                end
            end
        end
    
        for k = 1:1:size(A,3) 
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    if( i ~= j )
    
                        [ ~, E(j,k) ] = min(B(:,j,k),[],1); 
                    end                    
                end
            end  
            W(k) = min(E(:,k),[],1);
        end
    
        ii = 1;
        for k = 1:1:size(A,3)
            for j = 1:1:size(A,1)
    
                V  = [ A(ii,:,k); A(W(k),:,k) ];
                V_ = sum(V,1)/2;
    
                F(ii,:,k) = V_; 
                
                A(W(k),:,k) = NaN; W(k) = NaN;   
    
                ii = ii + 1;
            end
            ii = 1;
        end
    
        % We need a recusion on F to reduce the number of histograms to the
        % number of classes being evaluated until we reach the convergence criterion...
    
        ll = 0; % Resize accumulator.
    
        for k = 1:1:size(A,3)
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B_(i,j,k) = norm( ( F(j,:,k) - F(i,:,k) ), 2 );
                end
            end
        end
        
        for k = 1:1:size(A,3) 
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    if( i ~= j )
    
                        [ ~, E_(j,k) ] = min(B_(:,j,k),[],1); 
                    end                    
                end
            end  
            W_(k) = min(E_(:,k),[],1);
        end
    
        ii = 1;
        for k = 1:1:size(F,3)
            for j = 1:1:size(F,1)
    
                V  = [ F(ii,:,k); F(W(j,k),:,k) ];
                V_ = sum(V,1)/2;
            
                F(ii,:,k) = V_;
                
                A(W(k),:,k) = NaN; F(W(k),:,k) = NaN;
    
                ii = ii + 1;
            end
            ii = 1;
        end
    
         ll = ll + 1;
         for k = 1:1:size(F,3)
             for j = 1:1:size(F,2)
                for i = 1:1:size(F,1)-ll
    
                    if ( ~isnan(F(i,j,k)) && ~isnan(A(i,j,k)) )
    
                        Q(i,j,k) = F(i,j,k);
    
                        P(i,j,k) = A(i,j,k);
                    end
                end
             end
         end
         F = 0; F = Q; A = 0; A = P;    

    end
    
    RA = F;
end
   
