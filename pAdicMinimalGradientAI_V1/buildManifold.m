function [ RA, TL_ ] = buildManifold(A, l)

    global classType classGroups BINS RA Q BPI ii ij BFLAG
    
    for j = 1:1:size(classType,2)
           
        TL_(:,j) = (1:1:size(A,1)); % Sample tracking labels.
    end
    
    F   = A;

    N   = size(A,1);
    
    M   = size(classType,2); 

    O   = size(classGroups,2);
    
    
    CONTAINMENT = 0.25*N;
    
    V   = zeros(N-CONTAINMENT,M); % Stencil limit...
    
    IT  = N^M;
    
    
    B   = zeros(N,M);
    
    X   = zeros(N,BINS+1,2);
    

    % Convoltution with a sub-gradient!!
    
    ii = 1;

    aa = 1; bb = 1; 
    
    T  = 1.0;
    
    % MM = zeros(1,size(classGroups,2));
   
    F = sum(sum(B,1),2);

    while( sum(sum(B,1),2) < IT )

        if( aa <= size(V,1) )

            B(aa,1) = 1;

            for q = 1:1:size(A,1)
                
                % We can utilize Monte Carlo methods
                % in the convergence.
                
                A = cat(2,A,TL_);
                
                F = monteCarlo(A); 
                
                TL_(:,1) = F(:,end); F = F(:,1:end-2);
                
                
                RA(BPI(1,1),2:end,2) = mean(F(1:aa,:),1); 
                
                Q = RA; BFLAG = 1;
    
                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC >= T )
                    
                    [ ~, MM ] = max(ACC);
                
                    X(:,1:end-1,:) = RA; 
                    
                    X(:,end,:) = TL_(:,1);  
                end
                F  = circshift(A,1); TL_ = circshift(TL_,1);
                
                l  = circshift(l,1);
            end
            aa = aa + 1;

        elseif( aa >= size(V,1) && bb < size(V,1) )

            B(bb,2) = 1; B(:,1) = 0;

            for q = 1:1:size(A,1)
                
                A = cat(2,A,TL_);
                
                F = monteCarlo(A); 
                
                TL_(:,2) = F(:,end); F = F(:,1:end-1);
                
                 
                RA(BPI(1,2),2:end,2) = mean(F(1:bb,:),1);

                Q = RA; BFLAG = 0;
               
                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC >= T )
                    
                    [ ~, MM ] = max(ACC);
                
                    X(:,1:end-1,:) = RA; 
                    
                    X(:,end,:) = TL_(:,2);      
                end
                F  = circshift(A,1); TL_ = circshift(TL_,1);
                
                l  = circshift(l,1);
            end     
            bb = bb + 1; aa = 1;
            
        end  
        ii = ii + 1;

        F = sum(sum(B,1),2); 
    end
    RA = X(:,:,:); Q = RA;
    
    ii = ii + 2; ij = ij + 2;
end