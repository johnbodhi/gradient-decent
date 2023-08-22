function classifierTraining( dataSet, N, Observation )

    global RA Q classType BINS classGroups Train

    Train = 1;

    RA = zeros(size(classType,2)+1,BINS+1,size(classGroups,2)+1);
    
    [ RA ] = SVM( dataSet, N, Observation );

    Q = RA; % Resetting RA during classification.
end