

dataSet = pixels;

for i = 1:1:size(dataSet,1)
    for j = 1:1:size(dataSet,2)
        if( dataSet(i,j) == 0 )
    
            dataSet(i,j) = 1;
        end
    end
end

% D = zeros(1,4);ii = 1;
% for i = 1:1:size(dataSet,1)    
%     if( dataSet(i,4) == 3 )
% 
%         D(ii,:) = dataSet(i,:);ii = ii + 1;
%     end
% end