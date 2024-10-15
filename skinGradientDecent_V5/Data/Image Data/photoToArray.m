clear all; close all; clc;  

A1l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\A1l";

A2l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\A2l";  

A3l = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\A3l";

Adl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Adl"; 

Ael = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Ael";

FAL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\FAL"; 
  
FBL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\FBL";

FCL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\FCL"; 

FDL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\FDL";

FEL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\FEL"; 

FFL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\FFL";

Fgl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Fgl";  
 
Fhl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Fhl";

Fil = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Fil"; 
 
Fll = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Fll";

% Fml = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Fml"; 

Gal = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Gal";

% Gbl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Gbl"; 

Gcl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Gcl";

Gdl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Gdl"; 

Ggl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Ggl";

Gzl = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\Gzl";  

NAL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NAL";

NBL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NBL"; 

NCL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NCL";

NDL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NDL"; 
  
NEL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NEL";

New = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\New"; 
  
NFL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NFL";

NGL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NGL"; 

NHL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NHL";

NIL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NIL"; 

NLL = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NLL";

NML = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test\NML"; 


% Test = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data\test";

% Unknown = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V3\Data\Image Data\Unknown";

directoryArr = [ A1l A2l A3l Adl Ael FAL ...
                 FBL FCL FDL FEL FFL Fgl ...
                 Fhl Fil Fll Gal Gcl Gdl ...
                 Ggl Gzl NAL NBL NCL NDL ...
                 NEL New NFL NGL NHL NIL ...
                 NLL NML ];

str1_ = [ "A1l (" "A2l (" "A3l (" "Adl (" "Ael (" "FAL (" ...
          "FBL (" "FCL (" "FDL (" "FEL (" "FFL (" "Fgl (" ...
          "Fhl (" "Fil (" "Fll (" "Gal (" "Gcl (" "Gdl (" ...
          "Ggl (" "Gzl (" "NAL (" "NBL (" "NCL (" "NDL (" ...
          "NEL (" "New (" "NFL (" "NGL (" "NHL (" "NIL (" ...
          "NLL (" "NML (" ];

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

% Nr = 25; Mr = 25;
% 
% NN = 0; sampleSize = Nr * Mr; l = 1;
% 
% for p = 1:size( IMAGES, 2 )
% 
%      cd( directoryArr( p ) )
%     
%     for k = 1:numImages(p)
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