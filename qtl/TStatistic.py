'''
Created on Jan 13, 2015

@author: jiao
'''
import math
from __future__ import division

def twoSampleTTest(n1,n2,avg_x1,avg_x2,se1,se2):
    '''
    Statistical significance in the two sample problem (groups) for equal or unequal sample size, equal
    
    Algorithm: isites.harvard.edu/fs/docs/icb.topic1041076.files/lecture7_DiffExpr.ppt
    where n1 and n2 are the number of replicates, avg_x1 and avg_x2 are the mean of each sample group, se1 and se2 are the standard deviation,respectively.
    '''
    sp2  =  (n1-1)*(se1**se1) + (n2 -1)*(se2**se2) # Sp is an estimator of the common standard deviation of the two samples.
    t = (avg_x1 - avg_x2)/math.sqrt(sp2/n1+sp2/n2)
    return t