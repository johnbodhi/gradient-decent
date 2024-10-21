function [ F ] = histogramization( X, Y, Z )

    global frameLength classGroups DATARANGE BINS Supervision Randomized...
        Noise

    if ( Noise )
    
        [ X ] = addNoise( X );
    end 

    if ( Randomized )

        [ X ] = monteCarlo( X ); % Randomize all samples...
    end

    % Generate histograms...

    SCALE = floor( DATARANGE / BINS );

    F_ = zeros(sum(Y),DATARANGE,size(classGroups,2)); H = 1;

    for k = 1:1:size(X,2)-1
        for ii = 1:1:sum(Y)
            X_(:,1) = X((ii-1)*frameLength+1:ii*frameLength,k);
            for j = 0:1:DATARANGE
                for i = 1:1:size(X_,1)
                
                    if( j ==  X_(i,1) )

                        F_(ii,X_(i,1)+1,k) = H;
                        H = H + 1;
                    end
                end 
                H = 1;
            end
        end
    end

    if ( BINS < DATARANGE )

        F = zeros(sum(Y),BINS,size(classGroups,2));
    
        for k = 1:1:size(F_,3)
            for i = 1:1:size(F_,1)
                for j = 1:1:BINS
                
                    F(i,j,k) = sum(F_(i,(j-1)*SCALE+1:j*SCALE,k),2);
                end
            end
        end
        
    else

        F = F_;
    end

    F = F ./ frameLength;
    
    if( Noise )
        
        AWGN = rand([size(F,1),size(F,2),size(F,3)]);
        
        F = F + AWGN;
    end

    % F = cat(2,F,zeros(size(F,1),1,size(F,3)));

    % Append labels...

    if( Randomized && Supervision )

        for k = 1:1:size(F,3)
            for j = 1:1:size(F,1)
    
                F(j,size(F,2),k) = L(j,1);
            end
        end
    elseif( ~Randomized && Supervision )        

        for k = 1:1:size(F,3)
            for i = 1:1:size(F,1)
    
                F(i,size(F,2),k) = Z(i,1);
            end
        end
    elseif( ( Randomized && ~Supervision ) ||...
            ( ~Randomized && ~Supervision ) )
        
    end
    
end