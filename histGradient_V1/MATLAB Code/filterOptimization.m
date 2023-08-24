function [ V ] = filterOptimization( dataSet, N, Observation )

    global classType BINS RA

    % Allocate state-space to find all minimal distributions...

    % We need to generate all combinations of deleted delta trains 
    % and store them in state-space. Then, find the combinations 
    % of minimal target distributions 
    % among every class that are most efficient with maxmimal f-measures.
        
    S = permn( [ 0 1 ], BINS );

    V_ = zeros(size(S,1),size(S,2),size(RA,1),size(RA,3));

    for k = 1:1:size(V_,4)
        for m = 1:1:size(V_,3)
            for i = 1:1:size(V_,1)
                for j = 1:1:size(V_,2)

                    if ( S(i,j) )

                        V_(i,j,m,k) = RA(m,j,k);                        
                    else
                         
                        V_(i,j,m,k) = 0;
                    end
                end
            end
        end
    end

    % Reshape...

    V = cat(1,V_(:,:,:,2),V_(:,:,:,3),V_(:,:,:,4));

    % We need to apply a sub-gradient within the SVM to optimize the filter
    % efficiency.

    dataSet = histogramization( dataSet, N, Observation );

    SUP     = size(V,1)^size(classType,1);

    ii = 1;

    aa = 1; bb = 1; 
    cc = 1; dd = 1; 
    ee = 1; ff = 1; 
    hh = 1;

    while( ii <= SUP )


        if(  )

             RA(1,:,k) = V(aa,:,1,k); 
             aa = aa + 1;

        elseif(  )

             RA(2,:,k) = V(bb,:,2,k); 
             bb = bb + 1; aa = 1;

        elseif(  )

             RA(3,:,k) = V(cc,:,3,k); 
             cc = cc + 1; bb = 1;

        elseif(  )

             RA(4,:,k) = V(dd,:,4,k); 
             dd = dd + 1; cc = 1;

        elseif(  )

             RA(5,:,k) = V(ee,:,5,k); 
             ee = ee + 1; dd = 1;

        elseif(  )

             RA(6,:,k) = V(ff,:,6,k); 
             ff = ff + 1; ee = 1;

        elseif(  )

             RA(7,:,k) = V(gg,:,7,k); 
             gg = gg + 1; ff = 1;

        elseif(  )

             RA(8,:,k) = V(hh,:,8,k); 
             hh = hh + 1; gg = 1;
        end
        ii = ii + 1;        

        [ D, E ]  = classifier( dataSet, Observation );

        [ PREC REC ACC F1 ] = fMeasure( D, E ); 

        if( F1 >= 0.95 )

            A(rr,:) = [ F1 i j k ]; disp(A)

            X(:,:,:,rr) = RA; rr = rr + 1;
        end

    end

end