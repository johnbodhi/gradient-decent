function [ RA ] = runningAverage( rgbData, Observation )

global RA R G B W X classType classGroups CFLAG

    if( CFLAG )

        jj = 1;
        for ii = 1:1:size(rgbData,1)

            R(Observation(ii,1),jj) = rgbData(ii,1); 
%             G(Observation(ii,1),jj) = rgbData(ii,2); 
%             B(Observation(ii,1),jj) = rgbData(ii,3); 

            jj = jj + 1;
        end

        uu = 1;        
        for kk = 1:1:size(classGroups,2)
            for jj = 1:1:size(classType,2) 
                for ii = 1:1:size(rgbData,1)
                
                    if ( uu == Observation(ii,1) && uu == X(1,jj) )
        
                        RA(jj,1,kk) = W(jj,1,kk) * mean( R( Observation(ii,1), : ), 2 );            
%                         RA(jj,2,kk) = W(jj,2,kk) * mean( G( Observation(ii,1), : ), 2 );            
%                         RA(jj,3,kk) = W(jj,3,kk) * mean( B( Observation(ii,1), : ), 2 ); 
                        
                    end
                end
                uu = uu + 1;
            end           
        end

    else

        jj = 1;
        for ii = 1:1:size(rgbData,1)

            R(Observation(ii,1),jj) = rgbData(ii,1); 
%             G(Observation(ii,1),jj) = rgbData(ii,2); 
%             B(Observation(ii,1),jj) = rgbData(ii,3); 

            jj = jj + 1;
        end

        uu = 1;
        for kk = 1:1:size(classGroups,2)
            for jj = 1:1:size(classType,2)  
                for ii = 1:1:size(rgbData,1)
                
                    if ( uu == Observation(ii,1) )
        
                        RA(jj,1,kk) = W(jj,1,kk) * mean( R( Observation(ii,1), : ), 2 );            
%                         RA(jj,2,kk) = W(jj,2,kk) * mean( G( Observation(ii,1), : ), 2 );            
%                         RA(jj,3,kk) = W(jj,3,kk) * mean( B( Observation(ii,1), : ), 2 );
                        
                    end
                end
            end
            
            uu = uu + 1;
        end    
    end
    
end