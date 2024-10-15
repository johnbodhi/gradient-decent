clear all; close all; clc;  

% Training sequences...

% A1l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\Derm7pt\train\A1l"; 
% 
% A2l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\Derm7pt\train\A2l"; 
% 
% Cancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\breast_xrays\train\Cancerous"; 
% 
% nonCancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\breast_xrays\train\nonCancerous"; 
% 
% Pneumonia = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\chest_xray\train\PNEUMONIA";  
% 
% Normal = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\chest_xray\train\NORMAL"; 
% 
% Melanoma = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\MedNode\train\melanoma"; 
% 
% Naevus = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\MedNode\train\naevus"; 

% Test sequences...

A1l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\Derm7pt\test\A1l";

A2l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\Derm7pt\test\A2l";  

Cancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\breast_xrays\test\Cancerous"; 

nonCancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\breast_xrays\test\nonCancerous"; 

Pneumonia = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\chest_xray\test\PNEUMONIA";  

Normal = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\chest_xray\test\NORMAL"; 

Melanoma = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\MedNode\test\melanoma"; 

Naevus = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Image Data\MedNode\test\naevus"; 

directoryArr = [ A1l A2l...                 
                 Cancerous nonCancerous...
                 Pneumonia Normal...
                 Melanoma Naevus ];

str1_ = [ "A1l (" "A2l ("...          
          "Cancerous (" "nonCancerous ("...
          "Pneumonia (" "Normal ("...  
          "Melanoma (" "Naevus (" ];

for j = 1:size(str1_,2)
    cd( directoryArr( j ) );
    Images = dir('*.jpg'); numImages( j ) = 5;
    for i = 1:1:5
        str1 = str1_( j );  
        str2 = num2str( i ); 
        str3 = ").jpg";
    
        IMAGES( i, j ) = append( str1, str2, str3 );
    end    
end

% Nr = 40; Mr = 40;
% 
% NN = 0; sampleSize = Nr * Mr; l = 1;
% 
% for p = 1:size( IMAGES, 2 )
% 
%      cd( directoryArr( p ) )
%     
%     for k = 1:5
%     
% %         RGB = imread( IMAGES( k, p ) ); 
% %     
% %         RGB = imresize( RGB, [ Nr, Mr ] ); 
% % 
% %         I = rgb2gray(RGB);
% %     
% %         BW  = edge( I, 'Canny' );    
%      
%         h = imshow( IMAGES( k, p ) );
%         
%         im = imagemodel( h );    
%     
%         % Resized RGB vector...
%        
%         for j = 1:Nr
%             for i = 1:Mr          
%                 pixels( l, 1:3 ) = getPixelValue( im, i, j );
%                 pixels( l, 4:4 ) = p;
%                 l = l + 1;
%             end
%         end     
%     end          
% end
% dataSet = pixels;