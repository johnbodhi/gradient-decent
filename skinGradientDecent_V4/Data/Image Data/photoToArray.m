  
Melanoma = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V4\Data\Image Data\melanoma";

Naevus = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V4\Data\Image Data\naevus";  

Test = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V4\Data\Image Data\test";

% Unknown = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V3\Data\Image Data\Unknown";

directoryArr = [ Melanoma Naevus Test ];

str1_ = [ "Melanoma (" "Naevus (" "Test (" ];

for j = 1:size(str1_,2)
    cd( directoryArr( j ) );
    Images = dir('*.jpg'); numImages( j )= numel( Images );
    for i = 1:numImages( j )    
        str1 = str1_( j );  
        str2 = num2str( i ); 
        str3 = ").jpg";
    
        IMAGES( i, j ) = append( str1, str2, str3 );
    end    
end

Nr = 28; Mr = 28;

NN = 0; sampleSize = Nr * Mr; l = 1;

for p = 1:size( IMAGES, 2 )

     cd( directoryArr( p ) )
    
    for k = 1:numImages( p )
    
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