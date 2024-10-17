function [ PREC, REC, ACC, F1 ] = fMeasure( D, E )    
    
    TP = ( D - E );
    
    FP = E;
    
    TN = 0;
    
    FN = 0;
    
    PREC = TP / ( TP + FP );
    
    REC  = TP / ( TP + FN );

    if( TP == 0 )

        REC = 0;
    end
    
    ACC  = ( TP + TN ) / ( TP + TN + FP + FN );
    
    F1 = TP / ( TP + 0.5 * ( FP + FN ) );
end