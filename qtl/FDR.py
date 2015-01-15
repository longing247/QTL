'''
Created on Dec 18, 2014

@author: jiao
'''
from numpy import array, empty   
import sys
import os
import time
import math
from ngslib import bh_rejected,bh_qvalues,get_pi0,storey_qvalues

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')


from qtl.models import Gene,Marker,LOD,ExperimentMarker,ExperimentGene
from django.db.models import Q


def significantLOD(exp,_gxe):
    '''
    LOD:

    The log of the odds (LOD) ratio provides a measure of the association between variation in a phenotype and genetic differences (alleles) at a particular chromosomal locus. 
    It also provides a measure of the strength of linkage between two markers and can be used to evaluate whether two or more markers to each other on the same chromosome.

    A LOD score is defined as the logarithm of the ratio of two likelihoods:
    (1) the likelihood for the alternative hypothesis (that there is a QTL) and (2) the likelihood of the null hypothesis (that there is no QTL). 
    Likelihoods are probabilities, but they are not Pr(hypothesis | data) but rather Pr(data | hypothesis). 
    That's why they are called likelihoods rather than probabilities. (The "|" symbol translates to "given the").

    In the two likelihoods, one has maximized over the various nuisance parameters (the mean phenotypes for each genotype group, or overall for the null hypothesis, and the residual variance). 
    Or one can say, one has plugged in the maximum likelihood estimates for these nuisance parameters.

    With complete data at a marker, the log likelihood for the normal model reduces to the (-n/2) times the log of the residual sum of squares.

    LOD values can be converted to LRS scores (likelihood ratio statistics) by multiplying by 4.61.
    
    [Williams RW, June 15, 2005, updated with text from Karl Broman, Oct 28, 2010]
    '''
    tic = time.time()
    #lod1 = -math.log10(lod_th)
    #lod2 = -lod1
    #print lod1,lod2
    #lod_list = LOD.objects.filter(Q(LOD_score__gte=lod_th) |Q(LOD_score__lte=-lod_th),experiment_name=exp).values('id','LOD_score').order_by('-LOD_score')
    lod_list = LOD.objects.filter(experiment_name=exp,gxe=_gxe).values_list('LOD_score',flat=True).order_by('-LOD_score')
    lod_list = list(lod_list)
    
    print 'Size of the QueryValueSet %d' % sys.getsizeof(lod_list)
    toc = time.time()
    print 'QuerySet was successfully finished in %f seconds'%(toc-tic)
    return lod_list   

def LODToPval(lod_score):
    '''
    The LOD is also roughly equivalent to the -log(P), where P is the probability of linkage (P = 0.001 => 3). 
    The LOD itself is not a precise measurement of the probability of linkage, but in general for F2 crosses and RI strains, values above 3.3 will usually be worth attention for simple interval maps.
    [Williams RW, June 15, 2005, updated with text from Karl Broman, Oct 28, 2010]
    '''
    pval = math.pow(10,-math.fabs(lod_score))
    return pval


def significantPval(lod_list):
    '''
    Add an extra attribute 'pval' into each entry resulted from significantLOD()
    '''
    tic = time.time()
    sorted_lod_list = sorted(map(abs,lod_list),reverse=True) #sort the list which take the absolute value; The bigger LOD socre, the smaller p valuue
    print 'Size of the sorted input lod list %d' % sys.getsizeof(sorted_lod_list)
    pval_list = []
    for lod in sorted_lod_list:
        pval_list.append(LODToPval(lod))
    del sorted_lod_list
    print 'Size of the sorted input lod list %d after deleting' % sys.getsizeof(pval_list)
    print 'Size of the sorted p value list %d' % sys.getsizeof(pval_list)
    toc = time.time()
    print 'Pval list length: %f' % len(pval_list)
    print 'pval list was converted successfully in %f seconds'%(toc-tic)
    return pval_list

def correct_pval_for_multiple_testing(pval, correction_type):                
    """                                                                                                   
    consistent with R, but not used yet.
    """
    tic = time.time()                                                                     
    pval = array(pval) 
    n = int(pval.shape[0])                                                                           
    adjust_pval = empty(n)
    if correction_type == "Bonferroni":                                                                   
        adjust_pval = n * pval
    elif correction_type == "Bonferroni-Holm":                                                            
        values = [ (pvalue, i) for i, pvalue in enumerate(pval) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            adjust_pval[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pval) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            pvalue, index = vals
            adjust_pval[index] = new_values[i]     
    toc = time.time()
    print 'adjust pval list was converted successfully in %f seconds'%(toc-tic)
    print 'Size of the adjust p value list %d' % sys.getsizeof(adjust_pval)                                                                                                             
    return adjust_pval

def adjustedPval(pval_list):
    '''
    return the adjusted P value list. Not used yet.
    '''
    tic = time.time()
    adjustPval = correct_pval_for_multiple_testing(pval_list,"Benjamini-Hochberg") 
    toc = time.time()
    print 'Adjusted p values in %f seconds'%(toc-tic)
    return adjustPval

def pvalThld(pval_list,fdr):
    '''
    max{i:Pi<=i*a/m} where i is the maximal index, a is the FDR contrl rate, and m is the size of the list
    which allows the corresponding p value is smaller than the product of i*a/m
    Benjamini and Hochberg FDR-controlling procedure
    '''
    thld = None
    for i in range(len(pval_list)-1,-1,-1):# Division of zero error
        if i!=0:
            if pval_list[i] <= ((i+1)*fdr/len(pval_list)):
                thld = pval_list[i]
                break
        else:
            if pval_list[i] > ((i+1)*fdr/len(pval_list)):
                raise ValueError('None of the p value will be rejected.')
            else:
                thld = pval_list[i]
    return thld   

def pvalToLOD(pval):
    '''
    Convert p value to LOD score.
    The LOD is also roughly equivalent to the -log(P), where P is the probability of linkage (P = 0.001 => 3). 
    '''
    lod = -math.log(pval,10)
    return lod 

if __name__=="__main__":
    #BH pval: 0.000139798657224 LOD: 3.854497    Storey:(i)5928 qval:0.0499958818396 pval:0.000160715957486
    #lod_list = significantLOD('Ligterink_2014',False) 
    #BH pval: 1.68594343678e-07 LOD: 6.773157    Storey: 6 0.0486991023789 1.68594343678e-07
    #lod_list = significantLOD('Ligterink_2014',True) 
    #BH pval: 0.000495268834226 LOD: 3.305159    Storey: 7554 0.0499964825558 0.000639533060296
    #lod_list = significantLOD('Keurentjes_2007',False) 
    #BH pval: 0.000816476426122 LOD: 3.00805635    Storey: 52956 0.0499968149156 0.00104801112293
    lod_list = significantLOD('Snoek_2012',True) 
    pval_list = significantPval(lod_list)
    #print pval_list[:10]

    #bh_rejected_list = bh_rejected(pval_list,0.05)
    #print bh_rejected_list[-1]

    #storey_qvalues_list =storey_qvalues(pval_list)
    #rejected_index = -1
    #for i in xrange(len(storey_qvalues_list)-1,-1,-1):
    #    if storey_qvalues_list[i] <= 0.05:
    #        rejected_index = i
    #        break
    
    #print rejected_index,storey_qvalues_list[rejected_index],pval_list[rejected_index]


'''
Size of the QueryValueSet 15593528
QuerySet was successfully finished in 40.243835 seconds
Size of the sorted input lod list 15593528
Size of the sorted input lod list 14122144 after deleting
Size of the sorted p value list 14122144
Pval list length: 3465216.000000
pval list was converted successfully in 33.959055 seconds
0.0 1.0
0.005 0.97883432427
0.01 0.971932065838
0.015 0.966585169854
0.02 0.962224965756
0.025 0.958343041538
0.03 0.954936380515
0.035 0.95175997283
0.04 0.948884752062
0.045 0.946202735598
0.05 0.943641942636
0.055 0.94124888256
0.06 0.939013644859
0.065 0.936903267362
0.07 0.934834430467
0.075 0.932887362928
0.08 0.931078165552
0.085 0.929277580671
0.09 0.927560908497
0.095 0.925875968769
0.1 0.924288061055
0.105 0.922661452865
0.11 0.92115761683
0.115 0.919663200914
0.12 0.918218701623
0.125 0.916863562247
0.13 0.915565157982
0.135 0.914244069977
0.14 0.913016342187
0.145 0.911837034315
0.15 0.910582401179
0.155 0.909422205917
0.16 0.90830900722
0.165 0.907180057549
0.17 0.906079576506
0.175 0.904986394307
0.18 0.903990386521
0.185 0.902976138273
0.19 0.901959700401
0.195 0.900907692782
0.2 0.899956164349
0.205 0.899047479491
0.21 0.898166013462
0.215 0.897246482252
0.22 0.896235247081
0.225 0.895197186237
0.23 0.894336033444
0.235 0.89340100328
0.24 0.892624920898
0.245 0.891844674931
0.25 0.89104171284
0.255 0.890188074778
0.26 0.889447303437
0.265 0.888662687482
0.27 0.887843208952
0.275 0.887050241502
0.28 0.886313596484
0.285 0.885552118722
0.29 0.884855516354
0.295 0.884112602144
0.3 0.88353840816
0.305 0.882987094367
0.31 0.882339542888
0.315 0.881721296602
0.32 0.881061705142
0.325 0.880395334663
0.33 0.879901212852
0.335 0.879321982094
0.34 0.878726979178
0.345 0.878119367578
0.35 0.877593866498
0.355 0.876997580096
0.36 0.87644879136
0.365 0.875918173439
0.37 0.875394707294
0.375 0.874929355053
0.38 0.874308016883
0.385 0.873798577926
0.39 0.873288356883
0.395 0.872759208532
0.4 0.872294348558
0.405 0.87177656967
0.41 0.87132338325
0.415 0.870733204648
0.42 0.870259229698
0.425 0.869794075702
0.43 0.869366326735
0.435 0.868791568109
0.44 0.868353413714
0.445 0.867950521998
0.45 0.867428544924
0.455 0.867027249455
0.46 0.86660943759
0.465 0.866063528583
0.47 0.865587904651
0.475 0.865181825538
0.48 0.86474740347
0.485 0.864325839473
0.49 0.864099714668
0.495 0.863737107171
0.5 0.863399568743
0.505 0.862988749978
0.51 0.862590160178
0.515 0.862175616835
0.52 0.861929193832
0.525 0.861678190517
0.53 0.861548331714l
0.535 0.861128339164
0.54 0.860820295436
0.545 0.860426200665
0.55 0.860122748802
0.555 0.859896134237
0.56 0.859461705863
0.565 0.859185796157
0.57 0.858771930107
0.575 0.858414190698
0.58 0.858194973237
0.585 0.857834179137
0.59 0.857494851129
0.595 0.857097266229
0.6 0.856736636331
0.605 0.856335461273
0.61 0.85604979196
0.615 0.855797928708
0.62 0.855566017423
0.625 0.855299448385
0.63 0.855126288606
0.635 0.854807651438
0.64 0.854577158955
0.645 0.854446664658
0.65 0.854140941616
0.655 0.85379457116
0.66 0.853432071937
0.665 0.853050998884
0.67 0.852768564132
0.675 0.852599087618
0.68 0.852299863847
0.685 0.852158793482
0.69 0.851845608544
0.695 0.851714228057
0.7 0.851663119028
0.705 0.851532996128
0.71 0.851313801765
0.715 0.851047426162
0.72 0.850907583085
0.725 0.850801482243
0.73 0.850645492336
0.735 0.850475993089
0.74 0.850281105798
0.745 0.850158926233
0.75 0.849933741504
0.755 0.849688764579
0.76 0.849712543172
0.765 0.849683301189
0.77 0.849523553134
0.775 0.84950035502
0.78 0.849691227434
0.785 0.849425219226
0.79 0.849373287195
0.795 0.849147080234
0.8 0.849307806497
0.805 0.849056480913
0.81 0.848960520483
0.815 0.848683103815
0.82 0.848810322679
0.825 0.848677666434
0.83 0.8486169914
0.835 0.848685561846
0.84 0.848186664266
0.845 0.848204817023
0.85 0.848387709934
0.855 0.84811352358
0.86 0.848297974581
0.865 0.848119862808
0.87 0.847899191899
0.875 0.847912511082
0.88 0.847982251804
0.885 0.847726814646
0.89 0.848083036986
0.895 0.84799771157
0.9 0.84725454344
0.905 0.847484194196
0.91 0.847380237063
0.915 0.847335346483
0.92 0.847681645242
0.925 0.847208370272
0.93 0.846329431041
0.935 0.846358606394
0.94 0.845507658589
0.945 0.844580697265
0.95 0.844738105792

'''        
            
            
            
    