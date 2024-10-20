clear all; close all; clc;  

% Training sequences...

% train = "C:\Users\jmgar\OneDrive\Documents\GitHub\gradient-decent\pAdicMinmalEnsembleAI_V1\Data\Image Data\Derm7pt\train"; 

% Test sequences...

test = "C:\Users\jmgar\OneDrive\Documents\GitHub\gradient-decent\pAdicMinmalEnsembleAI_V1\Data\Image Data\Derm7pt\test";


directoryArr = [ test ];

str1_ = [ "test" ];

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


N = 50; M = 50; l = 1;

for p = 1:size( IMAGES, 2 )

     cd( directoryArr( p ) )

     for k = 1:numImages( p )

        RGB = imread( IMAGES( k, p ) ); 

        RGB = imresize( RGB, [ N, M ] ); 
% 
%         I = rgb2gray( RGB );
%     
%         BW  = edge( I, 'Canny' );    

         h = imshow( IMAGES( k, p ) );

         im = imagemodel( h );    

         % Resized RGB vector...

         for j = 1:N
             for i = 1:M          
                 pixels( l, 1:3 ) = getPixelValue( im, i, j );
                 pixels( l, 4:4 ) = p;
                 l = l + 1;
             end
         end     
    end          
end
dataSet = pixels;