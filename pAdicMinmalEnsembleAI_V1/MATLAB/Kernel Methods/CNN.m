function [ RA ] = CNN( A_ )

    global classType classGroups BINS RA W Q BPI ii ij BFLAG
    
    ii = 2; ij = 3;

    BPI = [ ii ij ];
    
    
    N   = size(A_,1);
    
    M   = size(classType,2); 

    O   = size(classGroups,2);


    for j = 1:1:M
           
        TL_(:,j) = (1:1:N); % Sample tracking labels (supervision).
    end


    CONTAINMENT = floor(0.50 * N);
    
    V   = zeros(N-CONTAINMENT,M); % Stencil limit...
    
    IT  = size(V,1)^M;
    
    
    B   = zeros(N,M);
    
    X   = zeros(N,BINS+1,2);
    

    % Convolution with a sub-gradient!!

    aa = 1; bb = 1; ss = 0;
    
    T  = 1.0; jj = 1;


    while( sum(sum(B,1),2) < IT )

        if( aa <= size(V,1) )

            B(aa,1) = 1;

            for q = 1:1:N 
                
                % We can utilize Monte Carlo methods
                % in a stencil to hasten convergence.
                
                A = cat(2,A_,TL_);
                
                F = monteCarlo(A); A = 0;
                
                F = F(:,1:end-1); TL_(:,1) = F(:,end); 
                
                
                RA(BPI(1,1),2:end,2) = W(1,:,1).*mean(F(1:aa,:),1); 
                
                BFLAG = 1; % Back propagation index for k-means.

    
                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 


                if( ACC >= T )
                
                    X(BPI(1,1),1:end-1,1) = RA(:,:,:); aa_ = aa; 

                    Q = X; % Filter reset.
                end
                A_  = circshift(A_,1); TL_ = circshift(TL_,1);
                
            end
            aa = aa + 1;

        elseif( aa >= size(V,1) && bb < size(V,1) )

            B(bb,2) = 1; B(:,1) = 0;

            for q = 1:1:N
                
                A = cat(2,A_,TL_);
                
                F = monteCarlo(A); A = 0;
                
                F = F(:,1:end-1); TL_(:,2) = F(:,end); 
                
                 
                RA(BPI(1,2),2:end,2) = W(2,:,1).*mean(F(1:bb,:),1);

                BFLAG = 0;
               

                [ D, E ] = classifier( F, l );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC >= T )

                    % We neeed to find the closest non-zero length(s) from the
                    % prior classification segment contained in the filter
                    % manifold...
                    
                    WALK_( jj, 1 ) = BiCGSTAB_( RA(BPI(1,1),:,:), [], RA(BPI(1,2),:,:) ); jj = jj + 1;

                    % We need to reinforce exlusivity of the samples 
                    % in the averages comprising the filter.

                    if ( aa_ > bb )

                        K = zeros(aa_,2);

                        K(:,1) = TL_(1:aa,1); K(:,2) = TL_(1:bb,2);

                        for i = 1:1:aa_

                            if( K(i,1) == K(i,2) && K(i,1) ~= 0 && K(i,2) ~= 0 )

                                ss = ss + 1;
                            end

                        end

                    elseif ( aa_ < bb )

                        K = zeros(bb,2);

                        K(:,1) = TL_(1:aa,1); K(:,2) = TL_(1:bb,2);

                        for i = 1:1:bb

                            if( K(i,1) == K(i,2) && ( K(i,1) ~= 0 && K(i,2) ~= 0 ) )

                                ss = ss + 1;
                            end

                        end
                           
                    end
                    K = 0;

                    % If the gradient exists for some subset of mutually exclusive
                    % labels, we assign a new filter to the array of new filters that
                    % have at least met the accuracy threshold criteria. We
                    % can then choose which filter is most minimum, ie. the
                    % prior classes dual subset.

                    if( ~ss )

                        X_(jj,1:end-1,1) = RA(:,:,:);
                    end

          
                    if( q == N )
                        
                        [ ~, I ]  = min( WALK_(:,1) );
                        
                        X(BPI(1,2),1:end-1,1) = X_(I,:,:);

                        % We need to remove the exclusive values of the
                        % sample set from the pool to accelerate 
                        % convergence of the convolution...

                        uu = 1;
                        for i = 1:1:size(A_,1)
                            
                            if( A_(i,end,1) ~= K(i,1) || A_(i,end,1) ~= K(i,2) )

                                A(uu,:,1) = A_(i,:,1); uu = uu + 1;
                            end                           
                        end                        
                   
                    end
                    Q = X; % Filter reset.
                    
                end
                A_  = circshift(A_,1); TL_ = circshift(TL_,1);
                
            end     
            bb = bb + 1; aa = 1; 
            
        end
        % We need to increment the backpropagation indexes by two for the
        % next set of comparisons, and reset the dual space walk index.
        ii = ii + 2; ij = ij + 2; jj = 1;

    end
    RA = Q;
    
end