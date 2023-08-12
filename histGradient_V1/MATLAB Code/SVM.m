function [ RA ] = SVM( dataSet, N, Observation )

    global classType classGroups imageLength RA BINS Supervision Train

    if ( ~Supervision && Train )
    
    % Unsupervised training with an SVM!

        A = zeros(sum(N),BINS,size(classGroups,2)); Z = 1;
    
        for k = 1:1:size(dataSet,2)-1
            for ii = 1:1:sum(N)
                X_(:,1) = dataSet((ii-1)*imageLength+1:ii*imageLength,k);
                for j = 0:1:BINS-1
                    for i = 1:1:size(X_,1)
                    
                        if( j ==  X_(i,1) )
    
                            A(ii,X_(i,1)+1,k) = Z;
                            Z = Z + 1;                  
                        end
                    end 
                    Z = 1;
                end
            end
        end    

        % Hist is deprecated!
    
%         for k = 1:1:size(dataSet,2)-1
%             for i = 1:1:size(dataSet,1)/(imageLength)
%     
%                 A_(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k),BINS);
%             end
%         end
     
        EIGEN_FRAMES = size(A,1);
    
        SEGMENTS     = EIGEN_FRAMES / size(classType,2);
    
        for k = 1:1:size(A,3)      
            for j = 1:1:size(A,1)
                for i = 1:1:size(A,1)
    
                    B(i,j,k) = norm( ( A(j,:,k) - A(i,:,k) ), Inf );
                    I(i,1,1) = i;
                end
            end
        end        
        % [ ~, I ] = sort(B,2);

        [ I, W ] = combinations( I, SEGMENTS ); % Permutation windows abuut A.
        
        RA = zeros(size(classType,2),size(A,2),size(A,3)); 
        
        ii = 1; jj = 1; NN = 0;   
        for k = 1:1:size(A,3)
            for i = 1:1:size(classType,2)
                
                while ( jj < SEGMENTS )
                
                    RA(i,:,k) = RA(i,:,k) + A(I(ii,jj,i),:,k); NN = NN + 1;
    
                    if( ii >= size(I,1) )

                        ii = 0; jj = jj + 1;
                    end
                    ii = ii + 1;
                end
                jj = 1;        
            end
        end

        RA = RA ./ NN; 
        
        %RA = sort(RA,2);  

%         RA = zeros(size(classType,2),size(A,2),size(A,3)); 
% 
%         ii = 1; jj = 1; NN = 0;
%         for k = 1:1:size(A,3)            
%             for i = 1:1:size(RA,1)
%                 while ( ii <= size(W,1) )
%                     for kk = 1:1:size(W,3)
%                     
%                         % Convergence criterion...
% 
%                         while ( jj <= SEGMENTS )
%                         
%                             RA(i,:,k) = RA(i,:,k) + A(W(ii,jj,kk),:,k);
%             
%                             jj = jj + 1; NN = NN + 1;                       
%                         end
%                         jj = 1;
%                         
%                     end
%                     ii = ii + 1;
%                 end
%                 ii = 1;
%             end                 
%         end
%     
%         RA = RA ./ NN; RA = sort(RA,2);
        
        % The per element histogram magnitude ratios appear to be identical 
        % bewteen the supervised and unsupervised cases for each group.
    else

        if ( Train )
        
            for k = 1:1:size(dataSet,2)-1
                for i = 1:1:size(dataSet,1)/(imageLength)
        
                    A(i,:,k) = hist(dataSet((i-1)*imageLength+1:i*imageLength,k),BINS);
                end
            end

            A = cat(2,A,zeros(size(A,1),1,size(A,3)));

            for k = 1:1:size(A,3)
                for i = 1:1:size(A,1)
    
                    A(i,size(A,2),k) = Observation(i,1);
                end
            end

        else
            
            A = dataSet;
        end

        EIGEN_FRAMES = size(A,1);
    
        SEGMENTS     = EIGEN_FRAMES / size(classType,2);
 
        ii = 1; jj = 1; kk = 1; NN = 0;
        for k = 1:1:size(A,3)

            while ( jj <= size(classType,2) )

                for i = 1:1:size(A,1)
                
                    if( A(i,size(A,2),k) == Observation(kk,1) && ii <= SEGMENTS )
    
                        RA(jj,:,k) = RA(jj,:,k) + A(i,1:end-1,k); 

                        ii = ii + 1; kk = kk + 1; NN = NN + 1;
                    end    
                    
                end   
                jj = jj + 1; ii = 1;
            end   
            jj = 1; kk = 1;
        end
 
        RA = RA ./ NN; 
        
        %RA = sort(RA,2);
    end

end   