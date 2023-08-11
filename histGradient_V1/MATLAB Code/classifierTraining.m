function classifierTraining( dataSet, Observation )

    global RA Q classType BINS classGroups Train

    Train = 1;

    RA = zeros(size(classType,2),BINS,size(classGroups,2));
    
    [ RA ] =  SVM( dataSet, Observation );

    Q = RA; % Resetting RA during classification.
end