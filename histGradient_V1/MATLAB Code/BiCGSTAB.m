function [ Z_ ] = BiCGSTAB( X_, Y_ )

    global classType classGroups frameLength

    RA  = X;

    N   = size(Y,1); 
    
    M   = size(classType,2);

    SUM = size(classGroups,2) * N * M; 
    
    V   = zeros(N,size(classGroups,2),size(classGroups,2));

    % SUP = simpleNN(N,M); 
    
    UB  = 1004;

    SUP = size(classGroups,2) * UB;
    
    ii = 1;

    aa = 1; bb = 1; 

    while( sum(sum(sum(V,1),2),3) < SUM )

        K = ceil( ii / SUP ) + 1;

        if( aa <= size(V,2) )

            V(aa,1,K) = 1;

            aa = aa + 1;    
            
            A(:,1) = Y_(aa,:,K);            
            
        elseif( aa > size(V,3) && bb <= size(V,3) )

            V(bb,2,K) = 1; V(:,1,K) = 0;
            
            bb = bb + 1; aa = 1;
            
            B(:,1) = frameLength * X_(bb,:,K);           

        elseif( bb > size(V,1) && cc <= size(V,1) )

            V(cc,3,K) = 1; V(:,2,K) = 0;

            cc = cc + 1; bb = 1;
            
            X_(:,1) = Y_(cc,:,K);              
            
        end







    end

       
end