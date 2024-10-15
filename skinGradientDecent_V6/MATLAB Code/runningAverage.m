function [ RA ] = runningAverage( RGB, rgbData, skinObservation )

global RA R G B W X classType classGroups CFLAG

    if( CFLAG )

        % During classification we can update RA with the decisions made
        % by the gradient.

        jj = 1;
        for ii = 1:1:size(rgbData,1)

            R(skinObservation(ii,1),jj) = rgbData(ii,1); 
%             G(skinObservation(ii,1),jj) = rgbData(ii,2); 
%             B(skinObservation(ii,1),jj) = rgbData(ii,3); 

            jj = jj + 1;
        end

        uu = 1;        
        for kk = 1:1:size(classGroups,2)
            for jj = 1:1:size(classType,2) 
                for ii = 1:1:size(rgbData,1)
                
                    if ( uu == skinObservation(ii,1) && uu == X(1,jj) )
        
                        RA(jj,1,kk) = W(jj,1,kk) * mean( R( skinObservation(ii,1), : ), 2 );            
                        %RA(jj,2,kk) = W(jj,2,kk) * mean( G(skinObservation(ii,1), 2 );            
                        %RA(jj,3,kk) = W(jj,3,kk) * mean( B(skinObservation(ii,1), 2 ); 
                        
                    end
                end
                uu = uu + 1;
            end           
        end

    else

        jj = 1;
        for ii = 1:1:size(RGB,1)

            R(skinObservation(ii,1),jj) = RGB(ii,1); 
%             G(skinObservation(ii,1),jj) = RGB(ii,2); 
%             B(skinObservation(ii,1),jj) = RGB(ii,3); 

            jj = jj + 1;
        end

        uu = 1;
        for kk = 1:1:size(classGroups,2)
            for jj = 1:1:size(classType,2)  
                for ii = 1:1:size(RGB,1)
                
                    if ( uu == skinObservation(ii,1) )
        
                        RA(jj,1,kk) = W(jj,1,kk) * mean( R( skinObservation(ii,1), : ), 2 );            
                        %RA(jj,2,kk) = W(jj,2,kk) * mean( G( skinObservation(ii,1), : ), 2 );            
                        %RA(jj,3,kk) = W(jj,3,kk) * mean( B( skinObservation(ii,1), : ), 2 );
                        
                    end
                end
            end
            
            uu = uu + 1;
        end    
    end
    
end