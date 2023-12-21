function [ RA ] = filterCreation( A )

    global classType classGroups

    N   = size(A,1); 
    
    M   = size(classType,2); 

    O   = size(classGroups,2);

    SUM = O*N^(M+V); 
    
    B   = zeros(N,M,O);

    % SUP = simpleNN(N,M); 
    
    UB  = 1171591994624;

    SUP = size(classGroups,2)*UB;

    V   = (size(A,1) - size(classType,2)); % Pattern limit...

    % Shift initial entries for faster convergence...

    Z(:,1) = (1:1:size(A,1)); 
    
    L(:,1) = Z;

    for j = 2:1:size(classType,2)

        L(:,j) = circshift(Z,j-1);
    end

    % Convoltuion with a sub-gradient!!

    ii = 1;

    aa = 1; bb = 1; 
    cc = 1; dd = 1; 
    ee = 1; ff = 1; 
    gg = 1; hh = 1;

    rr = 1; xx = 1;
    
    T = 0.95;

    while( sum(sum(sum(B,1),2),3) < SUM )

        K = ceil( ii / SUP ) + 1;

        if( aa <= size(V,1) )

            B(aa,1,K) = 1;

            aa = aa + 1;

            for q = 1:1:size(A,1)
            
                C = sum(A(1:aa,:,K),1);
    
                RA(2,:,K) = RA(2,:,K) + C; RA = RA / aa;
    
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift( A, 1 );
                
                rr = rr + 1;
            end

        elseif( aa > size(V,1) && bb <= size(V,1) )

             B(bb,2,K) = 1; B(:,1,K) = 0;

             bb = bb + 1; aa = 1;

            for q = 1:1:size(A,1)

                C =  sum(A(1:bb,:,K),1);
                 
                RA(3,:,K)  = RA(3,:,K) + C; RA = RA / bb;
               
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1; 
            end

        elseif( bb > size(V,1) && cc <= size(V,1) )

            B(cc,3,K) = 1; B(:,2,K) = 0;

            cc = cc + 1; bb = 1;

            for q = 1:1:size(A,1)
    
                C = sum(A(1:cc,:,K),1);
                 
                RA(4,:,K) = RA(4,:,K) + C; RA = RA / cc;
                              
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1;
            end

        elseif( cc > size(V,1) && dd <= size(V,1) )

            B(dd,4,K) = 1; B(:,3,K) = 0;

            dd = dd + 1; cc = 1;

            for q = 1:1:size(A,1)

                C = sum(A(1:dd,:,K),1);
                 
                RA(5,:,K) = RA(5,:,K) + C; RA = RA / dd;
                
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1; 
            end

        elseif( dd > size(V,1) && ee <= size(V,1) )

            B(ee,5,K) = 1; B(:,4,K) = 0;

            ee = ee + 1; dd = 1;

            for q = 1:1:size(A,1)

                C = sum(A(1:ee,:,K),1);
                 
                RA(6,:,K)  = RA(6,:,K) + C; RA = RA / ee; 
    
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E );

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1; 
            end

        elseif( ee > size(V,1) && ff <= size(V,1) )

            B(ff,6,K) = 1; B(:,5,K) = 0;

            ff = ff + 1; ee = 1;

            for q = 1:1:size(A,1)

                C = sum(A(1:ff,:,K),1);
                 
                RA(7,:,K) = RA(7,:,K) + C; RA = RA / ff;
                
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E ); 

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1; 
            end

        elseif( ff > size(V,1) && gg <= size(V,1) )

            B(gg,7,K) = 1; B(:,5,K) = 0;

            gg = gg + 1; ff = 1;

            for q = 1:1:size(A,1)

                C = sum(A(1:gg,:,K),1);
                 
                RA(8,:,K) = RA(8,:,K) + C; RA = RA / gg; 
    
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E );

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1; 
            end

        elseif( gg > size(V,1) && hh <= size(V,1) )

            B(hh,8,K) = 1; B(:,7,K) = 0;

            hh = hh + 1; gg = 1;

            for q = 1:1:size(A,1)

                C = sum(A(1:hh,:,K),1);
                 
                RA(9,:,K) = RA(9,:,K) + C; RA = RA / hh; 
               
                [ D, E ] = classifier( A, Observation );
                
                [ ~, ~, ACC, ~ ] = fMeasure( D, E );

                if( ACC > T )

                    A(xx,1) = ACC; xx = xx + 1;
                end
    
                A = circshift(A,1);
                
                rr = rr + 1; 
            end
        end
        ii = ii + 1;
    end

    [ RA ] = verifyFilter( RA, A );

end