function [ RA ] = SVM( dataSet, N, Observation )

    global RA classType Supervision Randomized Optimized Train

    A = histogramization( dataSet, N, Observation ); 

    if ( ( ~Supervision && Train && ~Randomized ) ||...
         ( Supervision && Train && Randomized ) )

        [ RA ] = filterCreation( A );       
        
    elseif( ( ~Supervision && ~Train && ~Randomized ) ||...
            ( Supervision && ~Train && Randomized ) ) 
        
        [ RA ] = filterUpdate( dataSet );

    end

    if ( Optimized )        
    
        [ RA ] = filterOptimization( dataSet, N, Observation );
    end

end