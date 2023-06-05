function [ X ] = dataSetRandomized(dataSet,R,C)

    if ( isfile( 'dataSetRandomized.csv' ) )

        X = readmatrix('dataSetRandomized.csv');

    else

        dataSetRandomized = zeros( R, C ); 

        i = zeros(R,1);

        for k = 1:R

            i( k ) = randi( [1, R ] );

            while ( find( i( 1:k-1, 1 ) == i( k ) ) )
                i( k ) = randi( [ 1, R ] );      
            end

            dataSetRandomized( k, 1:C )  = dataSet( i( k ), 1:C );
        end        

            writematrix(dataSetRandomized,'dataSetRandomized.csv');

            X = dataSetRandomized;
    end    
end