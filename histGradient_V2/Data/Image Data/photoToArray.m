clear all; close all; clc;  

% Training sequences...

% Tumor    = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Brain MRI Images\train\yes"; 
% 
% NonTumor = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Brain MRI Images\train\no"; 
% 
% COVID    = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Covid-19 Images\train\Covid"; 
% 
% Normal = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Covid-19 Images\train\Normal"; 

% Test sequences...

Tumor    = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Brain MRI Images\test\yes"; 

NonTumor = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Brain MRI Images\test\no"; 

COVID    = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Covid-19 Images\test\Covid"; 

Normal   = "C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V2\Data\Image Data\Covid-19 Images\test\Normal"; 

directoryArr = [ Tumor NonTumor...                 
                 COVID Normal ];

str1_ = [ "Tumor (" "NonTumor ("...          
          "COVID (" "Normal (" ];

for j = 1:size(str1_,2)

    cd( directoryArr( j ) );    

    Images = dir('*.jpg'); 

    numImages( j ) = numel( Images );

    for i = 1:1:numImages( j )

        str1 = str1_( j );  
        str2 = num2str( i ); 
        str3 = ").jpg";
    
        IMAGES( i, j ) = append( str1, str2, str3 );
    end    
end

Np = 30; Mp = 30;

NN = 0; sampleSize = Np * Mp; l = 1;

for p = 1:size( IMAGES, 2 )

     cd( directoryArr( p ) )

    for k = 1:numImages( p )

        RGB = imread( IMAGES( k, p ) ); 

        RGB = imresize( RGB, [ Np, Mp ] ); 

%         I = rgb2gray(RGB);
%     
%         BW  = edge( I, 'Canny' );    

        h = imshow( RGB );

        im = imagemodel( h );    

        % Resized RGB vector...

        for j = 1:Np
            for i = 1:Mp          
                pixels( l, 1:3 ) = getPixelValue( im, i, j );
                pixels( l, 4:4 ) = p;
                l = l + 1;
            end
        end     
    end          
end
dataSet = pixels;