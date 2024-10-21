function [ F ] = complexHistogramization( X,Y,Z )

    global frameLength classGroups DATARANGE...
        Supervision Randomized

   if ( Noise )
    
        [ X ] = addNoise( X );
    end 

    if ( Randomized )

        [ X ] = monteCarlo( X ); % Randomize all samples...
    end

    % Generate histograms...
    
    % SCALE = floor( DATARANGE / BINS );

    C_    = zeros(frameLength*size(classGroups,2)^2,2);

    F_    = zeros(sum(Y),DATARANGE^2,1); H = 1; jj = 1;

    aa = 1; bb = 1;
    cc = 1; dd = 1;
    uu = 1; vv = 1;    
   
    for ii = 1:1:sum(Y)

        X_(:,:) = X((ii-1)*frameLength+1:ii*frameLength,:);
        
        while( dd <= 2*size(C_,1) )

            if( uu <= frameLength*size(classGroups,2)^2 )
            
                C_(uu,vv) = X_(aa,cc); 
                
                aa = aa + 1;
                uu = uu + 1;                     
                dd = dd + 1;
            end

            if( vv <= 1 )
            
                if( aa > size(X_,1) )
    
                    aa = 1; bb = bb + 1;
                elseif( bb > size(classGroups,2) )
                    
                    bb = 1; cc = cc + 1;
                elseif( cc > size(classGroups,2) )

                    cc = 1; uu = 1;
                    vv = vv + 1;                              
                end
                
            elseif( vv > 1 )
                
                if( aa > size(X_,1) )
    
                    aa = 1; cc = cc + 1;
                elseif( cc > size(classGroups,2) )
                    
                    cc = 1;    
                end
            end

        end
        
        for i = 0:1:DATARANGE-1
            for j = 0:1:DATARANGE-1
                for k = 1:1:size(C_,1)

                    if( i == C_(k,1) && j == C_(k,2) ) 
                        
                        H = H + 1;
                    end
                end
                F_(ii,jj,1) = H; H = 0;
                
                jj = jj + 1; 
            end                
        end
        jj = 1;

        aa = 1; bb = 1;
        cc = 1; dd = 1;
        uu = 1; vv = 1;
    end

    F = F_;

    F = F ./ size(C_,1);

    F = cat(2,F,zeros(size(F,1),1,size(F,3)));

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