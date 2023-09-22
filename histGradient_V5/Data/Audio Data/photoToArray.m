clear all; close all; clc;  

% Training sequences...

A1l = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\Derm7pt\train\A1l"; 

A2l = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\Derm7pt\train\A2l"; 

Cancerous = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\breast_xrays\train\Cancerous"; 

nonCancerous = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\breast_xrays\train\nonCancerous"; 

Pneumonia = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\chest_xray\train\PNEUMONIA";  

Normal = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\chest_xray\train\NORMAL"; 

Melanoma = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\MedNode\train\melanoma"; 

Naevus = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\MedNode\train\naevus"; 

% Test sequences...

% A1l = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\Derm7pt\test\A1l";
% 
% A2l = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\Derm7pt\test\A2l";  
% 
% Cancerous = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\breast_xrays\test\Cancerous"; 
% 
% nonCancerous = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\breast_xrays\test\nonCancerous"; 
% 
% Pneumonia = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\chest_xray\test\PNEUMONIA";  
% 
% Normal = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\chest_xray\test\NORMAL"; 
% 
% Melanoma = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\MedNode\test\melanoma"; 
% 
% Naevus = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data\MedNode\test\naevus"; 

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
    Images = dir('*.jpg'); numImages( j ) = numel( Images );
    for i = 1:1:numImages( j )
        str1 = str1_( j );  
        str2 = num2str( i ); 
        str3 = ").jpg";
    
        IMAGES( i, j ) = append( str1, str2, str3 );
    end    
end

% Np = 100; Mp = 100;
% 
% NN = 0; sampleSize = Np * Mp; l = 1;
% 
% for p = 1:size( IMAGES, 2 )
% 
%      cd( directoryArr( p ) )
%     
%     for k = 1:numImages( p )
%     
% %         RGB = imread( IMAGES( k, p ) ); 
% %     
% %         RGB = imresize( RGB, [ Np, Mp ] ); 
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
%         for j = 1:Np
%             for i = 1:Mp          
%                 pixels( l, 1:3 ) = getPixelValue( im, i, j );
%                 pixels( l, 4:4 ) = p;
%                 l = l + 1;
%             end
%         end     
%     end          
% end
% dataSet = pixels;