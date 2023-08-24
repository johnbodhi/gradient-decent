function [ V ] = filterOptimization( dataSet, N, Observation )

    global classType BINS RA

    % Allocate state-space to find all minimal distributions...

    % We need to generate all combinations of deleted delta trains 
    % and store them in state-space. Then, find the combinations 
    % of minimal target distributions 
    % among every class that are most efficient with maxmimal f-measures.
        
    S = permn( [ 0 1 ], BINS );

    V_ = zeros(size(S,1),size(RA,2),size(RA,1),size(RA,3));

    for k = 1:1:size(V_,4)
        for m = 1:1:size(V_,3)
            for i = 1:1:size(V_,1)
                for j = 1:1:size(V_,2)-1

                    if ( S(i,j) )

                        V_(i,j+1,m,k) = RA(m,j+1,k);                        
                    else
                         
                        V_(i,j+1,m,k) = 0;
                    end
                end
            end
        end
    end

    % We need to reshape for a more efficient mapping...

    V = cat(1,V_(:,:,:,2),V_(:,:,:,3),V_(:,:,:,4));

    % We need to apply a sub-gradient within the SVM to optimize the filter
    % efficiency.

    dataSet = histogramization( dataSet, N, Observation );

    SUP     = size(V,1)^size(classType,2);

    ii = 1;

    aa = 1; bb = 1; 
    cc = 1; dd = 1; 
    ee = 1; ff = 1; 
    hh = 1;

    while( ii <= SUP )

        K = ceil( aa / size(S,1) );

        if( aa <= size(V,1) )

             RA(1,:,K) = V(aa,:,1); 
             aa = aa + 1;
        elseif( aa > size(V,1) && bb <= size(V,1) )
             
             RA(2,:,K) = V(bb,:,2); 
             bb = bb + 1; aa = 1;
        elseif( bb > size(V,1) && cc <= size(V,1) )

             RA(3,:,K) = V(cc,:,3); 
             cc = cc + 1; bb = 1;
        elseif( cc > size(V,1) && dd <= size(V,1) )

             RA(4,:,K) = V(dd,:,4); 
             dd = dd + 1; cc = 1;
        elseif( dd > size(V,1) && ee <= size(V,1) )

             RA(5,:,K) = V(ee,:,5); 
             ee = ee + 1; dd = 1;
        elseif( ee > size(V,1) && ff <= size(V,1) )

             RA(6,:,K) = V(ff,:,6); 
             ff = ff + 1; ee = 1;
        elseif( ff > size(V,1) && gg <= size(V,1) )

             RA(7,:,K) = V(gg,:,7); 
             gg = gg + 1; ff = 1;
        elseif( hh > size(V,1) )

             RA(8,:,K) = V(hh,:,8); 
             hh = hh + 1; gg = 1;
        end
        ii = ii + 1;        

        % [ D, E ]  = classifier( dataSet, Observation );
        % 
        % [ PREC REC ACC F1 ] = fMeasure( D, E ); 
        % 
        % if( F1 >= 0.95 )
        % 
        %     A(rr,:) = F1;
        % 
        %     X(:,:,:,rr) = RA; rr = rr + 1;
        % end
    end

end