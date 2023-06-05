import numpy as np
import pandas as pd
import csv
import math
import statistics
import os

# Image Labels: Non-Cancerous = 1, Cancerous = 2, Unknown = 3.


def pixelClassifierTraining( data, obs, trainingN, totalN ):

    global i, R, RA, CFLAG

    CFLAG = 0;

    classes = 3;

    R = np.zeros((classes,totalN))

    for i in range(1,trainingN-1):

        RA = runningAverage( data[i], obs[i] );

    #print(RA)

########################################################################
            
def runningAverage( data, obs ):

    global i, R, RA, W, imgDecision, CFLAG
    

    if CFLAG and imgDecision == 1:

        R[0,i] = data;
    elif CFLAG and imgDecision == 2:

        R[1,i] = data;      
    elif CFLAG and imgDecision == 3:

        R[2,i] = data;

    ####################################################################
        
    if ~CFLAG and obs == 1:
        
        R[0,i] = data;
    elif ~CFLAG and obs == 2:

        R[1,i] = data;
    elif ~CFLAG and obs == 3:

        R[2,i] = data;
        
    ####################################################################

    W  = [ 1e1, 1e2, 1e3 ];

    RA = [ W[0] * statistics.mean( R[0,:] ),
           W[1] * statistics.mean( R[1,:] ),
           W[2] * statistics.mean( R[2,:] ) ];
    
    #print(RA)

    return RA

    ####################################################################

def imageClassification( data, obs, trainingN, testN ):

    global i, j, N, R, imageLength, imgDecision, CFLAG, predictions

    #print(N)

    CFLAG = 1;

    D = 0; E = 0;

    rgbData     = np.zeros(N[j]*Nr*Mr);

    labels      = np.zeros(sum(N));
    predictions = np.zeros(sum(N));

    A = np.zeros(sum(N))

    ii = -1; jj = 0; 

    for i in range(trainingN,trainingN+testN):

        ii += 1;

        rgbData[ii] = data[i];

        if ii == imageLength-1:

            #print(rgbData)

            rgbData = kmeans(rgbData, obs);

            imgDecision = gradientDecent(rgbData);

            #print(imgDecision)
                    
            labels[jj] = obs[ii];

            #print(labels)
            
            predictions[jj] = imgDecision;

            #print(predictions)

            ################################################################
            
            beta = 1;

            X = pfbeta( labels, predictions, beta )

            A = [ labels[jj], predictions[jj], X ];

            ii = -1; jj += 1;
                 
            # print(A)
       
        #runningAverage( rgbData[i], obs[i] );

        #print(RA)
        
    return 0

    #########################################################################

def kmeans( rgbData, obs ):

    global RA, imageLength
    
    p = 2;

    ii = 0; jj = 0;

    N = 1;

    D = np.zeros(imageLength);
    Y = np.zeros(imageLength);

    while jj < N * imageLength-1:

        if ii == N * imageLength-1:
            
            ii = 1;
        
        for i in range(len(D)):  
            D[i] = ( ( rgbData[i] - RA[int(obs[ii]-1)] )**p )**(1/p);
            #print(D[i])

        Ci = statistics.mean(D[:]);

        Cn = Ci;
   
        for i in range(0,imageLength-1):            
            Y[i] = ( ( rgbData[i] - Cn )**p )**(1/p)
            #print(Y[i])

        ii += 1;
        jj += 1;

    return Y

############################################################################

def gradientDecent( F ):

    global RA, W, imageLength

    classes = 3;
    
    ii = 0;

    eps = 1e1;

    Y     = np.zeros(imageLength);
    gamma = np.zeros((imageLength,classes));
    iiVec = np.zeros((imageLength,classes));
    S     = np.zeros(classes);

    for i in range(0,imageLength-1,1):
        Y[i] = F[i];

    for k in range(1,classes,1):        
        for i in range(1,imageLength-1,1):            
            gamma[i,k] = abs( ( F[i] - F[i-1] ) * ( RA[k] - RA[k-1] + eps ) ) \
            / (abs( RA[k] - RA[k-1] ) + eps )**2            
            #print(gamma[i,k])

    eps = np.zeros(classes);
    
    for k in range(0,classes,1):
        eps[k] = RA[k] / W[k]
        #print(eps[k])
            
    for k in range(0,classes,1):   
        for i in range(0,imageLength-1,1):
            
            if RA[k] - Y[i] > 0:
                while Y[i] < RA[k]:
                    Y[i] = Y[i] + ( gamma[i,k] * RA[k] + eps[k] );
                    ii += 1;
                    #print(ii)
                    
                iiVec[i,k] = ii;                    
                ii = 0;

            for i in range(0,imageLength-1,1):
                Y[i] = F[i];
                
            if RA[k] - Y[i] < 0:
                while Y[i] > RA[k]:
                    Y[i] = Y[i] - ( gamma[i,k] * RA[k] + eps[k] );
                    ii += 1;
                    #print(ii)
                    
                iiVec[i,k] = ii;                    
                ii = 0;
                
        for i in range(0,imageLength-1,1):
            Y[i] = F[i];

    for k in range(0,classes,1):
        S[k] = sum(iiVec[:,k]);
        #print(S[k])

    if S[0] > S[1]:
        
        Z = 2;        
    else:
        
        Z = 1;

    #print(Z)
    
    return Z

#########################################################################

def pfbeta(labels, predictions, beta):
    
    y_true_count = 0
    ctp          = 0
    cfp          = 0

    for idx in range(len(labels)):
        
        prediction = min( max(predictions[idx],0), 1)
        
        if (labels[idx]):
            
            y_true_count += 1
            ctp += prediction            
        else:
            
            cfp += prediction

    beta_squared = beta * beta
    c_precision  = ctp / (ctp + cfp)
    c_recall     = ctp / y_true_count
    
    if (c_precision > 0 and c_recall > 0):
        
        result = (1 + beta_squared) * (c_precision * c_recall) \
                 / (beta_squared * c_precision + c_recall)
        
        return result
    else:
        
        return 0

########################################################################

def outputFile( predictions ):
            
    # Shift the prediciton labels for test images.
    for jj in range(len(predictions)):

        if predictions[jj] == 1:

            predictions[jj] = 0      
        if predictions[jj] == 2:

            predictions[jj] = 1

    pred = [ '10008_R', '10008_L','10008_R', '10008_L' ]

    can  = [ predictions[0], predictions[1], predictions[2], predictions[3] ]

    arr  = { 'prediction_id': pred, 'cancer': can }

    df   = pd.DataFrame( arr )

    df.to_csv( 'submission.csv', index = False )

    print(df); df.info();

    print("Done!")
    
    return 0

#########################################################################

global N, predictions

classType = [ 1, 2, 3 ];

Nr = 58;

Mr = 58;

imageLength = Nr * Mr;

N = [ 0, 0, 4 ];

M = 1;

for j in range(len(classType)):

    data = np.zeros(imageLength);

    csv_  = open('trainR.csv','r').read()
    rows = csv_.split('\n');

    for i in range(1,imageLength-1):
        
        data[i] = int(rows[i]);

    obs  = np.zeros(imageLength);

    csv_  = open('trainObs.csv','r').read()
    rows = csv_.split('\n');

    for i in range(1,imageLength):        
        obs[i] = int(rows[i]);

    totalN    = np.size( data, 0 );
    trainingN = math.floor( 0.051 * totalN );
    testN     = math.floor( 0.00 * totalN );

    pixelClassifierTraining( data, obs, trainingN, totalN )

##########################################################################

    data = np.zeros(40368);
    M    = np.size(data,0);

    csv_ = open('testR.csv','r').read()
    rows = csv_.split('\n');

    for i in range(1,M):        
        data[i] = int(rows[i]);
        #print(data[i])

    obs  = np.zeros(40368);
    M    = np.size(obs,0);    

    csv_ = open('testObs.csv','r').read()
    rows = csv_.split('\n');

    for i in range(1,M):        
        obs[i] = int(rows[i]);
        #print(obs[i])

    totalN    = np.size( data, 0 );
    trainingN = math.floor( 0.0 * totalN );
    testN     = math.floor( 1.0 * totalN );

    data_ = np.zeros(N[j]*Mr*Nr);
    obs_  = np.zeros(N[j]*Mr*Nr);
    M     = np.size(data,0);

    ii = 0;
    for i in range(0,M):        
        if obs[i] == classType[j] and ii <= N[j]*Nr*Mr-1:          
            data_[ii] = data[i];
            obs_[ii]  = obs[i];
            ii += 1;
            #print(ii)

    del data
    del obs

    data = np.zeros(N[j]*Mr*Nr);
    obs  = np.zeros(N[j]*Mr*Nr);
    M    = np.size(data,0);
    
    for i in range(0,M):
        data[i] = data_[i];
        obs[i]  = obs_[i];
        #print(obs[i])

    del data_
    del obs_

    totalN    = N[j]*Mr*Nr;
    trainingN = math.floor( 0.0 * totalN );
    testN     = math.floor( 1.0 * totalN );

    #print(data)
    #print(obs)

    imageClassification( data, obs, trainingN, testN )

    if j == 2:

        outputFile(predictions)
        
# display(pd.read_csv('submission.csv').head())
# os.remove('submission.csv')
# !head ./submission.csv

