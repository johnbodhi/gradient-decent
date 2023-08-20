function filterOptimization()

    global classType BINS RA

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
    
end