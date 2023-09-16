function [ V ] = filterOptimization( dataSet, Observation )

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

    SUP     = size(V,1)^size(classType,2);

    rr = 1; ii = 1;

    aa = 1; bb = 1; 
    cc = 1; dd = 1; 
    ee = 1; ff = 1; 
    gg = 1; hh = 1;

    % We need to convolve! This is the shave, or keying.

    while( ii <= SUP )

        Ka = ceil( aa / size(S,1) )+1;
        Kb = ceil( bb / size(S,1) )+1;
        Kc = ceil( cc / size(S,1) )+1;
        Kd = ceil( dd / size(S,1) )+1;
        Ke = ceil( ee / size(S,1) )+1;
        Kf = ceil( ff / size(S,1) )+1;
        Kg = ceil( gg / size(S,1) )+1;
        Kh = ceil( hh / size(S,1) )+1;

        if( aa <= size(V,1) )

            RA(2,:,Ka) = V(aa,:,2);
            aa = aa + 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;            
            end

        elseif( aa > size(V,1) && bb <= size(V,1) )
             
            RA(3,:,Kb) = V(bb,:,3); 
            bb = bb + 1; aa = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break          
            end
        elseif( bb > size(V,1) && cc <= size(V,1) )

            RA(4,:,Kc) = V(cc,:,4); 
            cc = cc + 1; bb = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;             
            end
        elseif( cc > size(V,1) && dd <= size(V,1) )

            RA(5,:,Kd) = V(dd,:,5); 
            dd = dd + 1; cc = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;                
            end
        elseif( dd > size(V,1) && ee <= size(V,1) )

            RA(6,:,Ke) = V(ee,:,6); 
            ee = ee + 1; dd = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;              
            end
        elseif( ee > size(V,1) && ff <= size(V,1) )

            RA(7,:,Kf) = V(ff,:,7); 
            ff = ff + 1; ee = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;                
            end
        elseif( ff > size(V,1) && gg <= size(V,1) )

            RA(8,:,Kg) = V(gg,:,8); 
            gg = gg + 1; ff = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;                
            end
        elseif( hh > size(V,1) )

            RA(9,:,Kh) = V(hh,:,9); 
            hh = hh + 1; gg = 1;
            
            [ D, E ] = classifier( dataSet, Observation );
            
            [ PREC REC ACC F1 ] = fMeasure( D, E ); 
            
            if( F1 >= 0.99 )
            
                A(rr,1) = F1; [ ~, M ] = max(A);
                
                X(:,:,:,rr) = RA; rr = rr + 1; break;               
            end
        end
        ii = ii + 1;

    end
    RA = X(:,:,:,M);

end