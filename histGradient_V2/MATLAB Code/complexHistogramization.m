function [ F ] = complexHistogramization( X, Y, Z )

    global frameLength classGroups DATARANGE BINS Supervision Randomized

    if ( Randomized )

        [ X, L ] = randomizeAll( X, Y ); % Randomize all photos.
    end

    % Generate histograms...
    
    % SCALE = floor( DATARANGE / BINS );

    F_    = zeros(sum(Y),DATARANGE^2,1); H = 1; jj = 1;



   
    for ii = 1:1:sum(Y)

        X_(:,:) = X((ii-1)*frameLength+1:ii*frameLength,:);
        
        while( jj <= size(classGroups,2) )
            
            C_(uu,vv) = X_(ii,jj);

            
            
            


        end
        
        
        
        
         
        
        for i = 0:1:DATARANGE
            for j = 0:1:DATARANGE
                for k = 1:1:size(C_,1)

                    if( i == C_(k,1) && j == C_(k,2) ) 
                        
                        H = H + 1;
                    end
                end
                F_(ii,jj,1) = H; H = 1;
                
                jj = jj + 1; 
            end                
        end

    end

    % if ( BINS < DATARANGE )
    % 
    %     % F = zeros(sum(Y),BINS,size(classGroups,2));
    %     % 
    %     % for k = 1:1:size(F_,3)
    %     %     for i = 1:1:size(F_,1)
    %     %         for j = 1:1:BINS
    %     % 
    %     %             F(i,j,k) = sum(F_(i,(j-1)*SCALE+1:j*SCALE,k),2);
    %     %         end
    %     %     end
    %     % end
    % 
    % else
    % 
    % F = F_;
    % end

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