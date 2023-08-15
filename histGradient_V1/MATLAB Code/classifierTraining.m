function classifierTraining( dataSet, N, Observation )

    global RA Q classType BINS classGroups Train

    Train = 1;

    RA = zeros(size(classType,2),BINS,size(classGroups,2));
    
    [ RA ] =  SVM( dataSet, N, Observation );

    Q = RA; % Resetting RA during classification.
end