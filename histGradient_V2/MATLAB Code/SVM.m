function [ RA ] = SVM( dataSet, N, Observation )

    global Supervision Randomized Optimized Train

    A = histogramization( dataSet, N, Observation );

    A = complexHistogramization( dataSet, N, Observation );

    if ( ( ~Supervision && Train && ~Randomized ) ||...
         ( Supervision && Train && Randomized ) )

        [ RA ] = filterCreation( A );       
        
    elseif( ( ~Supervision && ~Train && ~Randomized ) ||...
            ( Supervision && ~Train && Randomized ) ) 
        
        [ RA ] = filterUpdate( dataSet );

    end

    if ( Optimized )        
    
        [ RA ] = filterOptimization( A, Observation );
    end

end