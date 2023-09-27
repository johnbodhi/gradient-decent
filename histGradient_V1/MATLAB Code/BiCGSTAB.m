function [ Z ] = BiCGSTAB( X, Y )

    global classType classGroups

    RA  = X;

    N   = size(Y,1); 
    
    M   = size(classType,2);

    SUM = size(classGroups,2) * N * M; 
    
    B   = zeros(N,size(classGroups,2),size(classGroups,2));

    R   = zeros(N,M,size(classGroups,2));

    % SUP = simpleNN(N,M); 
    
    UB  = 1004;

    SUP = size(classGroups,2)*UB;
    
    ii = 1;

    aa = 1; bb = 1; 

    while( sum(sum(sum(B,1),2),3) < SUM )

        K = ceil( ii / SUP ) + 1;

        if( aa <= size(B,2) )

            B(aa,1,K) = 1;

            aa = aa + 1;    
            
            

            
        elseif( aa > size(B,3) && bb <= size(B,3) )

            B(bb,2,K) = 1; B(:,1,K) = 0;
            
            bb = bb + 1; aa = 1;
            
             
            
            
        elseif( bb > size(B,1) && cc <= size(B,1) )

            B(cc,3,K) = 1; B(:,2,K) = 0;

            cc = cc + 1; bb = 1;
            
            
            
            
        end


        




    end

       
end