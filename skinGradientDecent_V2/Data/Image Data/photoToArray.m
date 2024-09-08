  
Benign = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V2\Data\Image Data\Benign";  

Malignant = "C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V2\Data\Image Data\Malignant";

directoryArr = [ Benign Malignant ];

str1_ = [ "Benign (" "Malignant (" ];

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

Nr = 58; Mr = 58;

NN = 0; sampleSize = Nr * Mr; l = 1;

for p = 1:size( IMAGES, 2 )

     cd( directoryArr( p ) )
    
    for k = 1:numImages( p )
    
        RGB = imread( IMAGES( k, p ) ); 
    
        RGB = imresize( RGB, [ Nr, Mr ] );
        
        I = rgb2gray(RGB);
    
        BW  = edge( I, 'Canny' );    
    
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