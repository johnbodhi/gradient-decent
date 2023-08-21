function [ V ] = filterOptimization( dataSet, N, Observation )

    global classType BINS RA Q

    % Allocate state-space to find all minimal distributions...

    % We need to generate all combinations of deleted delta trains 
    % and store them in state-space. Then, find the combinations 
    % of minimal target distributions 
    % among every class that are most efficient with maxmimal f-measures.
        
    S = permn( [ 0 1 ], BINS );

    V = zeros(size(S,1),size(S,2),size(classType,2),size(RA,3));

    for k = 1:1:size(V,4)
        for m = 1:1:size(V,3)
            for i = 1:1:size(V,1)
                for j = 1:1:size(V,2)

                    if ( S(i,j) )

                        V(i,j,m,k) = RA(m,j,k);                        
                    else

                        V(i,j,m,k) = 0;
                    end
                end
            end
        end
    end

    % We need to apply a sub-gradient within the SVM to optimize the filter
    % efficiency.

    dataSet = histogramization( dataSet, N, Observation );

    rr = 1;
    for k = 1:1:size(V,4)
        for m = 1:1:size(V,3)            
            for i = 1:1:size(V,1)

                RA(m,:,k) = V(i,:,m,k);            
             
                [ D, E ] = classifier( dataSet, Observation );
    
                [ PREC REC ACC F1 ] = fMeasure( D, E );
    
                AVE = [ PREC REC ACC F1 ]; 
    
                if( F1 == 1 )
    
                    A = [ i j ]; disp(A)

                    X(:,:,rr) = RA; rr = rr + 1;
                end

                RA = Q;
            end
        end
    end

end