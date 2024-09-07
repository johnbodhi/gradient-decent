function [ RA, TL_ ] = CNN(A_, l)

    global classType classGroups BINS RA Q BPI ii ij BFLAG
    
    ii = 2; ij = 3;

    BPI = [ ii ij ];
    
    for j = 1:1:size(classType,2)
           
        TL_(:,j) = (1:1:size(A,1)); % Sample tracking labels.
    end
    
    F   = A_;

    N   = size(A_,1);
    
    M   = size(classType,2); 

    % O   = size(classGroups,2);
    
    
    CONTAINMENT = 0.25*N;
    
    V   = zeros(N-CONTAINMENT,M); % Stencil limit...
    
    IT  = N^M;
    
    
    B   = zeros(N,M);
    
    X   = zeros(N,BINS+1,2);
    

    % Convoltution with a sub-gradient!!

    aa = 1; bb = 1; 
    
    T  = 1.0; ii = 1;
    
    % MM = zeros(1,size(classGroups,2)); % Random forest expansion...
   
    F = sum(sum(B,1),2);

    while( sum(sum(B,1),2) < IT )

        if( aa <= size(V,1) )

            B(aa,1) = 1;

            for q = 1:1:size(A_,1)
                
                % We can utilize Monte Carlo methods
                % in the convergence.
                
                A = cat(2,A_,TL_);
                
                F = monteCarlo(A); 
                
                TL_(:,1) = F(:,end); F = F(:,1:end-2);
                
                
                RA(BPI(1,1),2:end,2) = mean(F(1:aa,:),1); 
                
                Q = RA; BFLAG = 0;

    
                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 


                if( ACC >= T )
                    
                    % [ ~, MM ] = max(ACC);
                
                    X(BPI(1,1),1:end-1,1) = RA; 
                    
                    X(BPI(1,1),end,1) = TL_(:,1); aa_ = aa; 
                end
                F  = circshift(A,1); TL_ = circshift(TL_,1);
                
            end
            aa = aa + 1;

        elseif( aa >= size(V,1) && bb < size(V,1) )

            B(bb,2) = 1; B(:,1) = 0;

            for q = 1:1:size(A_,1)
                
                A = cat(2,A_,TL_);
                
                F = monteCarlo(A); 
                
                TL_(:,2) = F(:,end); F = F(:,1:end-1);
                
                 
                RA(BPI(1,2),2:end,2) = mean(F(1:bb,:),1);

                Q = RA; BFLAG = 1;
               

                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC >= T )
                    
                    WALK_(ii,1) =...
                        BiCGSTAB_(RA(BPI(1,1),:,:),[],RA(BPI(1,2),:,:));
                    
                        ii = ii + 1;
                        
                    if( ~size(find( TL_(1:aa_,1) == TL_(1:bb,2)),2) )
                        
                        X_(ii,1:end-1,1) = RA;
                    end
         
                    if( q == N )
                        
                        % [ ~, MM ] = max(ACC);
                        
                        [ ~, I ]  = min(WALK);
                        
                        X(BPI(1,2),1:end-1,1) = X_(I,:,:);
                    
                        X(BPI(1,2),end,1) = TL_(1:bb,2);
                    end
                    
                end
                F  = circshift(A,1); TL_ = circshift(TL_,1);
                
            end     
            bb = bb + 1; aa = 1;
            
        end

        F = sum(sum(B,1),2);
    end
    RA = X(:,:,:); Q = RA;
    
end