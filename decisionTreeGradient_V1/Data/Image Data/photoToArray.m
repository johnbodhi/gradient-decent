clear all; close all; clc;  

% Training sequences...

A1l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Data\Derm7pt\A1l"; 

A2l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Data\Derm7pt\A2l"; 

Cancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\breast_xrays\Cancerous"; 

nonCancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\breast_xrays\nonCancerous"; 

Pneumonia = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\chest_xray\train\PNEUMONIA";  

Normal = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\chest_xray\train\NORMAL"; 

Fractured = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\ "; 

nonFractured = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\ "; 

melanoma = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\MedNode\train\melanoma"; 

naevus = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\MedNode\train\naevus"; 

% Test sequences...

% A1l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Data\Derm7pt\test\A1l";
% 
% A2l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Data\Derm7pt\test\A2l";  
% 
% Cancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\breast_xrays\Cancerous"; 
% 
% nonCancerous = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\\breast_xrays\nonCancerous"; 
% 
% Pneumonia = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\chest_xray\test\PNEUMONIA";  
% 
% Normal = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\chest_xray\test\NORMAL"; 
% 
% Fractured = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Data\Derm7pt\A1l"; 
% 
% nonFractured = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\Data\Derm7pt\A2l"; 
% 
% melanoma = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\MedNode\test\melanoma"; 
% 
% naevus = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V1\Data\MedNode\test\naevus"; 

directoryArr = [ A1l A2l...
                 melanoma naevus...
                 Cancerous nonCancerous...
                 Pneumonia Normal...
                 Fractured nonFractured ];

str1_ = [ "A1l (" "A2l ("...
          "Cancerous (" "nonCancerous ("...
          "Pneumonia (" "Normal ("...
          "Fractured (" "nonFractured (" ];

for j = 1:size(str1_,2)
    cd( directoryArr( j ) );
    Images = dir('*.jpg'); numImages( j ) = numel( Images );
    for i = 1:1:numImages(j)
        str1 = str1_( j );  
        str2 = num2str( i ); 
        str3 = ").jpg";
    
        IMAGES( i, j ) = append( str1, str2, str3 );
    end    
end

Nr = 100; Mr = 100;

NN = 0; sampleSize = Nr * Mr; l = 1;

for p = 1:size( IMAGES, 2 )

     cd( directoryArr( p ) )
    
    for k = 1:numImages(p)
    
%         RGB = imread( IMAGES( k, p ) ); 
%     
%         RGB = imresize( RGB, [ Nr, Mr ] ); 
% 
%         I = rgb2gray(RGB);
%     
%         BW  = edge( I, 'Canny' );    
     
        h = imshow( IMAGES( k, p ) );
        
        im = imagemodel( h );    
    
        % Resized RGB vector...
       
        for j = 1:Nr
            for i = 1:Mr          
                pixels( l, 1:3 ) = getPixelValue( im, i, j );
                pixels( l, 4:4 ) = p;
                l = l + 1;
            end
        end     
    end          
end
dataSet = pixels;