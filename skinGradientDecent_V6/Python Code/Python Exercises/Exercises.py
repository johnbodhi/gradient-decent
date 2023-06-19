import numpy as np
import pandas as pd
import csv
import math
import statistics
import os



def pfbeta(labels, predictions, beta):
    
    y_true_count = 0
    ctp = 0
    cfp = 0

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

#############################################################################

A = np.zeros((10,10,2));

N = np.size(A,0);

print(N);

print(" ")

#############################################################################

for j in range(0,N):

    print(" ")
    for i in range(0,N):

        A[i,j] = i;

        print(i)

print(" ")

#############################################################################

ii = 0;
while ii < 5:
    
    ii += 1;
    print("Hello!")

print(" ")

#############################################################################

for i in range(0,N):
    
    CFLAG = 1;

    if CFLAG:
        
        CFLAG = 0;
        
        print("Second sight is cool!")

print(" ")

#############################################################################
        
S = sum(sum(A[0:N,0:N,0],0),0);

print(S)

print(" ")

#############################################################################

M = statistics.mean(A[:,0,0])

print(M)

print(" ")

#############################################################################

##imageLength = 31;
##
##X    = np.zeros(imageLength);
##
##csv  = open('shortTrainR.csv', 'r').read()
##rows = csv.split('\n');
##
##for i in range(1,imageLength-1):
##    X[i] = int(rows[i]);
##    print(X[i])
##
##print(" ")

#############################################################################

##print(sum(X[:]))
##
##print(" ")

#############################################################################

C = np.arange(10);

print(C)

print(" ")

#############################################################################

print(np.linalg.norm(C));

print(" ")

#############################################################################

N = np.size(C,0);
D = np.zeros(N);

ii = -1;
for i in range(0,N-1):

    ii += 1;
    
    D[ii] = C[i];   
    
    print(D[ii])

print(" ")

#############################################################################

Y = np.zeros(C);
N = np.size(Y,0);

Cn = 5;
for i in range(0,N-1):

    Y[i] = math.floor(sum(abs( ( C[i] - Cn )**p )**(1/p)))

    print(Y[i])

print(" ")

##############################################################################

print(math.floor(12.34))

print(" ")

##############################################################################

labels      = [ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 ];
predictions = [ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 ];

beta        = 1;

print(pfbeta( labels, predictions, beta ))

labels      = [ 1, 1, 1, 1, 2, 2, 2, 2 ];
predictions = [ 1, 1, 1, 1, 2, 2, 2, 1 ];

beta        = 1;

print(pfbeta( labels, predictions, beta ))

labels      = [ 1, 1, 1, 1 ];
predictions = [ 1, 1, 1, 1 ];

beta        = 1;

print(pfbeta( labels, predictions, beta ))

labels      = [ 2, 2, 2, 2 ];
predictions = [ 2, 2, 2, 2 ];

beta        = 1;

print(pfbeta( labels, predictions, beta ))

labels      = [ 3, 3, 3, 3 ];
predictions = [ 3, 3, 3, 3 ];

beta        = 1;

print(pfbeta( labels, predictions, beta ))


##############################################################################

print(" ")

for i in range(1,3,1):

    print(i)

##############################################################################

#create a data in dictionary form
##data = {'Login email': ['laura@example.com', '2070', 'Laura', 'Grey'],
##       'Identifier': ['craig@example.com', '4081', 'Craig', 'Johnson'],
##        'First name': ['mary@example.coms', '9346', 'Mary', 'Jenkins'],
##       'Last name' : ['jamie@example.com', '5079', 'Jamie', 'Smith'] }
###convert the dictionary to a pandas dataframe
##df = pd.DataFrame(data)
###save the dataframe as a CSV file
##df.to_csv('second_written_data.csv')

##############################################################################
 
#open the CSV file
##with open('writen_data.csv', mode='w', newline = '') as file:
##    data = csv.writer(file)
## 
##    #write data to the file, row by row
##    data.writerow(['Login email', 'Identifier', 'First name', 'Last name'])
##    data.writerow(['laura@example.com', '2070', 'Laura', 'Grey'])
##    data.writerow(['craig@example.com', '4081', 'Craig', 'Johnson'])
##    data.writerow(['mary@example.com', '9346', 'Mary', 'Jenkins'])
##    data.writerow(['jamie@example.com', '5079', 'Jamie', 'Smith'])

##with open('submission.csv', mode = 'w', newline = '') as file:
##    data = csv.writer(file)
## 
##    data.writerow(['prediction_id', 'cancer'])
##    data.writerow(['10008_R', 1.0 ])
##    data.writerow(['10008_L', 1.0 ])
##    #data.writerow(['10008_L', min(max(predictions[2],0),1)])
##    #data.writerow(['10008_R', min(max(predictions[3],0),1)])
##
##print("Done!")

##############################################################################

nme = ["aparna", "pankaj", "sudhir", "Geeku"]
deg = ["MBA", "BCA", "M.Tech", "MBA"]
scr = [90, 40, 80, 98]

# dictionary of lists
dict = {'name': nme, 'degree': deg, 'score': scr}
	
df = pd.DataFrame(dict)

df.to_csv('submission.csv',index=False)

print(df)

