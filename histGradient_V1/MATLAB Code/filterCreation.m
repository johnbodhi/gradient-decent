function [ RA ] = filterCreation( dataSet )

    global classType classGroups frameLength

    [ RA ] = initializeFilter( dataSet );

    SUP = size(classGroups,2)*size(dataSet,1)^size(classType,2);

    V   = (size(dataSet,1) - size(classType,2));

    Z(:,1) = (1:1:size(dataSet,1)); L(:,1) = Z;

    for j = 2:1:size(classType,2)

        L(:,j) = circshift(Z,j-1);
    end

    ii = 1;

    aa = 1; bb = 1; 
    cc = 1; dd = 1; 
    ee = 1; ff = 1; 
    gg = 1; hh = 1;

    rr = 1; 
    
    while( ii <= SUP )

        K = ceil( ii / SUP )+1;

        if( aa <= size(V,1) )

            for q = 1:1:size(dataSet,1)
            
                C = sum(dataSet(1:aa,:),1);
    
                RA(2,:,K) = RA(2,:,K) + C;
                RA        = RA / (aa*frameLength);
    
                aa = aa + 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1;
            end

        elseif( aa > size(V,1) && bb <= size(V,1) )

            for q = 1:1:size(dataSet,1)

                C =  sum(dataSet(1:bb,:),1);
                 
                RA(3,:,K)  = RA(3,:,K) + C; 
                RA         = RA / (bb*frameLength);
    
                bb = bb + 1; aa = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; 
            end

        elseif( bb > size(V,1) && cc <= size(V,1) )

            for q = 1:1:size(dataSet,1)
    
                C = sum(dataSet(1:cc,:),1);
                 
                RA(4,:,K) = RA(4,:,K) + C; 
                RA        = RA / (cc*frameLength);
    
                cc = cc + 1; bb = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1;
            end

        elseif( cc > size(V,1) && dd <= size(V,1) )

            for q = 1:1:size(dataSet,1)

                C = sum(dataSet(1:dd,:),1);
                 
                RA(5,:,K) = RA(5,:,K) + C; 
                RA        = RA / (dd*frameLength);
                dd = dd + 1; cc = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; 
            end

        elseif( dd > size(V,1) && ee <= size(V,1) )

            for q = 1:1:size(dataSet,1)

                C = sum(dataSet(1:ee,:),1);
                 
                RA(6,:,K)  = RA(6,:,K) + C; 
                RA         = RA / (ee*frameLength); 
    
                ee = ee + 1; dd = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; 
            end

        elseif( ee > size(V,1) && ff <= size(V,1) )

            for q = 1:1:size(dataSet,1)

                C = sum(dataSet(1:ff,:),1);
                 
                RA(7,:,K) = RA(7,:,K) + C; 
                RA         = RA / (ff*frameLength); 
    
                ff = ff + 1; ee = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; 
            end

        elseif( ff > size(V,1) && gg <= size(V,1) )

            for q = 1:1:size(dataSet,1)

                C = sum(dataSet(1:gg,:),1);
                 
                RA(8,:,K) = RA(8,:,K) + C; 
                RA        = RA / (gg*frameLength); 
    
                gg = gg + 1; ff = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; 
            end

        elseif( hh > size(V,1) )

            for q = 1:1:size(dataSet,1)

                C = sum(dataSet(1:hh,:),1);
                 
                RA(9,:,K) = RA(9,:,K) + C; 
                RA        = RA / (hh*frameLength); 
                
                hh = hh + 1; gg = 1;
                
                [ D, E ] = classifier( dataSet, Observation );
                
                [ PREC REC ACC F1 ] = fMeasure( D, E ); 
    
                dataSet = circshift(dataSet,1);
                
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; 
            end
        end
        ii = ii + 1;
    end
    RA = X(:,:,:,M);   

end