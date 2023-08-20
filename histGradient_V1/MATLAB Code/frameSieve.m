function [ Y ] = frameSieve( S )

    global RA JX

    % S = sort(S,"ascend");

    SPLITS = 8;

    for k = 1:1:size(RA,3)
        for i = 1:1:SPLITS
                
            Y(i,:,k) = S(1,1:end-1,k); 
            L(i,1,k) = i;
        end
    end

    % T = 1;
    % for k = 1:1:size(Y,3)
    %     for i = 1:1:size(Y,1)
    %         [ ~, I ] = max(Y(i,:,k));
    %         for j = 1:1:size(Y,2)
    % 
    %             if( j > JX(i,k) + T || j < JX(i,k) - T )
    %                 Y(i,j,k) = 0;
    %             end
    %         end
    %     end
    % end

    Y = cat(2,Y,L);    
end