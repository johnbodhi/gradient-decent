function [ Z_ ] = averageMapping( W_ )
    
    global M BINS
    
    [WA,LA] = sort(W_,2);
    
    [WB,LB] = sort(W_,1);
        
    % We need to separate the walk accumulation deltas from each other
    % according to an average edge amplitude. This is the discrimination
    % value for a class.
            
    EP_MU = 0.1; 
    
    % We can  step through an epsilon vector as well once 
    % we have experience...
        
    RA_ = zeros(M,BINS); RB_ = zeros(M,BINS);
    
    Z_ = zeros(size(W_,1),size(W_,2),size(W_,3),2);
        
    ll = 0;
    for k = 1:1:size(W_,3)
        for j = 1:1:size(W_,2)
            for i = 2:1:size(W_,1)
                        
                Z_(i-1,j,k,1) = WA(j,i,k) - WA(j,i-1,k);
            
                Z_(i-1,j,k,2) = WB(i,j,k) - WB(i-1,j,k);

                        
                E(1,1) = Z_(i-1,j,k,1) / Z_(i,j,k,1);
                    
                E(1,2) = Z_(i-1,j,k,2) / Z_(i,j,k,2);
            
                    
                if( E(1,1) <= EP_MU && E(1,2) <= EP_MU &&...
                        LA(j,i-1,k) ~= LB(j,i-1,k))
            
                    RA_(k,:,j) = RA_(k,:,j) + Y_(LA(j,i-1,k),:);
               
                    RB_(k,:,j) = RB_(k,:,j) + Y_(LB(j,i-1,k),:);
             
                    ll = ll + 1;
               end
                
            end
        end
        
        if( ll )
            
            RA_(k,:,j) = RA_(k,:,j) / ll;
           
            RB_(k,:,j) = RB_(k,:,j) / ll;
            
            ll = 0;
        end
        
    end
    RA = ( RA_ + RB_ ) ./ 2;
    
end