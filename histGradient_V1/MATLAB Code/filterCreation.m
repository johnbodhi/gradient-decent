function [ Y ] = filterCreation( dataSet )

    global classType classGroups frameLength

    SUP = size(classGroups,2)*size(dataSet,1)^size(classType,2);

    V   = (size(dataSet,1) - size(classType,2));

    rr = 1; ii = 1;

    aa = 1; bb = 1; 
    cc = 1; dd = 1; 
    ee = 1; ff = 1; 
    gg = 1; hh = 1;

    % We need to convolve! This is the shave, or keying.

    while( ii <= SUP )

        K = ceil( ii / SUP )+1;

        if( aa <= size(V,1) )
            
            C = sum(dataSet(1:aa,:),2);

            Y(2,:,K) = Y(2,:,K) + C;
            Y        = Y / (aa*frameLength);
            aa = aa + 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;            
            end

        elseif( aa > size(V,1) && bb <= size(V,1) )

            C =  sum(dataSet(1:bb,:),2);
             
            Y(3,:,K)  = Y(3,:,K) + C; 
            Y         = Y / (bb*frameLength);

            bb = bb + 1; aa = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break          
            end
        elseif( bb > size(V,1) && cc <= size(V,1) )

            C = sum(dataSet(1:cc,:),2);
             
            Y(4,:,K) = Y(4,:,K) + C; 
            Y        = Y / (cc*frameLength);

            cc = cc + 1; bb = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;             
            end
        elseif( cc > size(V,1) && dd <= size(V,1) )

            C = sum(dataSet(1:dd,:),2);
             
            Y(5,:,K) = Y(5,:,K) + C; 
            Y        = Y / (dd*frameLength);
            dd = dd + 1; cc = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;                
            end
        elseif( dd > size(V,1) && ee <= size(V,1) )

            C = sum(dataSet(1:ee,:),2);
             
            Y(6,:,K)  = Y(6,:,K) + C; 
            Y         = Y / (ee*frameLength); 

            ee = ee + 1; dd = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;              
            end
        elseif( ee > size(V,1) && ff <= size(V,1) )

            C = sum(dataSet(1:ff,:),2);
             
            Y(7,:,K) = Y(7,:,K) + C; 
            Y         = Y / (ff*frameLength); 

            ff = ff + 1; ee = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;                
            end
        elseif( ff > size(V,1) && gg <= size(V,1) )

            C = sum(dataSet(1:gg,:),2);
             
            Y(8,:,K) = Y(8,:,K) + C; 
            Y        = Y / (gg*frameLength); 

            gg = gg + 1; ff = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;                
            end
        elseif( hh > size(V,1) )

            C = sum(dataSet(1:hh,:),2);
             
            Y(9,:,K) = Y(9,:,K) + C; 
            Y        = Y / (hh*frameLength); 
            
            hh = hh + 1; gg = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 

            dataSet = circshift(dataSet,1);
            
            if( F1 >= 0.99 )
            
                A(rr,1) = ACC; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = Y; rr = rr + 1; break;               
            end
        end
        ii = ii + 1;

    end
    Y = X(:,:,:,M);   
end