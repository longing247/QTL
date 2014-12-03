'''
Created on Dec 4, 2014

@author: jiao
'''
import sys,os,time
import csv
import json
import math
import re
import itertools
from datetime import datetime

sys.path.append('/home/jiao/QTL')
#sys.path.append('/mnt/geninf15/prog/www/django/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Experiment,Gene,Marker,LOD,Parent,RIL,Metabolite,MParent,MRIL,MLOD,Genotype,ExperimentMarker,ExperimentGene,ExperimentRIL,ExperimentParent,ExperimentMetabolite
from qtl.models import TF_Family,TF_Family_Ref,Protein,Promoter_Binding_Site,TF


def getMarkersList(experiment):
    '''
    return all the markers.
    '''
    marker_list = ExperimentMarker.objects.filter(experiment_name = experiment).values_list('marker_name',flat=True) 
    return marker_list
    
def getMarkerNamesList(marker_list):
    marker_queryset_list = Marker.objects.filter(marker_name__in=marker_list).order_by('marker_chromosome','marker_cm').values_list('marker_name',flat=True)
    return marker_queryset_list

def getChrMarkers(chr_name,marker_list):
    '''
    return all the markers along on the chr_name.
    '''
    marker_chr_list = Marker.objects.filter(marker_chromosome = chr_name, marker_name__in=marker_list).order_by('marker_cm').values_list('marker_name',flat=True)
    return marker_chr_list

def geneUpload(f):    
    '''
    populating gene data into MySQL under table qtl_gene table.
    locus_identifier = models.CharField(max_length=15,primary_key=True) #AT1G01480
    gene_model_name = models.CharField(max_length=20,blank = True)#AT1G01480.1
    gene_model_description = models.TextField(blank = True)
    gene_model_type = models.CharField(max_length = 20,blank = True)
    primary_gene_symbol = models.TextField(blank = True)
    all_gene_symbols = models.TextField(blank = True)    
    '''
    
    tic = time.clock()
    csv_reader = csv.reader(f,delimiter='\t')
    i = 0
    for line in csv_reader:
        print line[0]
        if not 'Locus Identifier' in line[0]:
            i+=1    
            add_gene = Gene()
            add_gene.locus_identifier = line[0].upper()
            add_gene.gene_model_name = line[1]
            add_gene.gene_model_description = line[2]
            add_gene.gene_model_type = line[3]
            add_gene.primary_gene_symbol = line[4]
            add_gene.all_gene_symbols = line[5]
            add_gene.save()
        else:
            #print 'header %s' % line
            continue
    toc = time.clock()
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i
    

def markerUpload(f):
    '''
    populating makrer data into MySQL under table qtl_marker table. AT5G56840 
    marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
    marker_chromosome = models.IntegerField()#1
    marker_cm = models.DecimalField(max_digits = 3, decimal_places =1)#64.6
    marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639  
    '''
    tic = time.clock()
    exp_name = getExperimentName(f)
    csv_reader = csv.reader(f,delimiter='\t') 
    i=0    
    for line in csv_reader:      
        if not 'Marker' in line[0]:
            i+=1
            add_marker = Marker()
            add_marker.marker_name = line[0]
            add_marker.marker_chromosome = line[1]
            add_marker.marker_cm = line[2]
            add_marker.marker_phys_pos = line[3]
            add_marker.experiment_name = exp_name
            
            add_marker.save()               
        else:
            print 'header %s' % line
            continue
    toc = time.clock()
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i
def lodUpload(f):   
    '''
    Experiment: 
    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    LOD:
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gene_name = models.ForeignKey(Gene)
    marker_name = models.ForeignKey(Marker)
    
    Occurred: Segmentation fault (core dumped) because of insufficient memory
    After saving 10839 records
    Solution:,index
    delete from qtl_lod where locus_identifier_id = 'AT2G36270';
    split the data set from that row into another one.
    Anyhow this part of function needs to be optimized by either using mmap or generator.
    '''
    tic = time.clock()
    exp_name = getExperimentName(f)
    add_exp = Experiment()
    add_exp.experiment_name = exp_name
    add_exp.save()
    marker_list = []  
    i=0
    for line in getData(f):
        if 'QTL' in line[0]:
            for marker in range(1,len(line)):
                print line[marker]
                marker_list.append(line[marker])
        else:
            i+=1
            print i
            print line[0].upper()
            
            for gene in range(1,len(line)):
                
                if Gene.objects.filter(locus_identifier=line[0].upper()).exists():         
                    add_lod = LOD(experiment_name = Experiment.objects.get(pk=exp_name),
                                  locus_identifier = Gene.objects.get(pk= line[0].upper()),
                                  marker_name =  Marker.objects.get(pk=marker_list[gene-1]),
                                  LOD_score =  line[gene])
                    add_lod.save()
                else:
                    print 'obslete'
                    add_gene = Gene(locus_identifier = line[0].upper(), gene_model_name = 'obsolete') # gene is not registed in the database.
                    add_gene.save()
                    add_lod = LOD(experiment_name = Experiment.objects.get(pk=exp_name),
                                  locus_identifier = Gene.objects.get(pk= line[0].upper()),
                                  marker_name =  Marker.objects.get(pk=marker_list[gene-1]),
                                  LOD_score =  line[gene])
                    add_lod.save()
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i

def parentUpload(f):   
    '''
    Parent: 
    parent_type = models.CharField(max_length=20)rror
    expression = models.DecimalField(max_digits = 12, decimal_places = 10)
    locus_identifier = models.ForeignKey(Gene)
    '''
    tic = time.clock()   
    exp_name = getExperimentName(f)
    parent_list = []  
    i=0
    for line in getData(f):
        if 'Parents' in line[0]:
            for parent in range(1,len(line)):
                parent_list.append(line[parent])
        else:
            i+=1
            print i
            for exp in range(1,len(line)):
                if line[exp] == '':        #js_exp_list = simplejson.dumps(exp_list)
                    continue                 
                else:
                    add_parent = Parent(parent_type =  parent_list[exp-1],
                                        expression =  line[exp],
                                        locus_identifier = Gene.objects.get(pk= line[0].upper()),
                                        experiment_name = exp_name
                                        )
                    add_parent.save()
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i
    
def rilUpload(f):
    '''
    RIL: 
    locus_identifier = models.ForeignKey(Gene)
    ril_name = models.CharField(max_length=20)
    ril_type = models.CharField(max_len','))gth=20)
    ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    
    Suppose the dataset is like this:
    ;RIL_type_1; RIL_type_2;...
    gene; RIL3_6h;RIL4_RP;...
    AT1Gxxxxx;2.012;1.23;....
    AT1Gxxxxx;..............
    r_list:
    '''
    tic = time.time()   
    getAll(f,'RIL_type','gene')
    toc = time.time()            
    print 'in %f seconds' % (toc-tic)
        

def getAll(f,arg1,arg2):
    li1 = []
    li2 = []
    exp_name = getExperimentName(f)
    i = 0
    for line in getData(f):
        if arg1 in line[0]:
            print 'RIL type'
            for obj in range(1,len(line)):
                li1.append(line[obj])
        elif arg2 in line[0]:
            print 'RIL name'
            for obj in range(1,len(line)):
                li2.append(line[obj])        
        else:  
            i+=1    
            for exp in range(1,len(line)):
                if line[exp]: 
                    add_ril = RIL(ril_name = li2[exp-1],
                                  ril_type = li1[exp-1],
                                  ril_exp = line[exp],
                                  locus_identifier = Gene.objects.get(pk= line[0].upper()),
                                  experiment_name = exp_name)
                    add_ril.save()
                else:
                    continue
            print line[0]
            print '%d records added' % i
            
def getAllMetabolite(f,arg1,arg2):
    li1 = []
    li2 = []
    exp_name = getExperimentName(f)
    i = 0
    for line in getData(f):
        if arg1 in line[0]:
            print 'RIL type'
            for obj in range(1,len(line)):
                li1.append(line[obj])
        elif arg2 in line[0]:
            print 'RIL name'
            for obj in range(1,len(line)):
                li2.append(line[obj])        
        else:  
            i+=1    
            for exp in range(1,len(line)):
                if line[exp]: 
                    add_ril = MRIL(ril_name = li2[exp-1],
                                  ril_type = li1[exp-1],
                                  ril_exp = line[exp],
                                  metabolite_name = Metabolite.objects.get(pk= line[0].upper()),
                                  experiment_name = exp_name)
                    add_ril.save()
                else:
                    continue
            print line[0]
            print '%d records added' % i
    
def getData(f):
    csv_reader = csv.reader(f,delimiter='\t')
    for line in csv_reader:
        yield line
        

def metaboliteUpload(f):
    '''
    class metabolite model
    metabolite_name = models.CharField(max_length=50,primary_key=True)
    '''
    tic = time.clock()
    i=0
    for line in getData(f):
        if 'metabolitesID' in line[0]:
            continue
        else:
            i+=1
            add_metabolite = Metabolite()
            add_metabolite.metabolite_name = line[0].upper()
            

    toc = time.clock()
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i
    
def metaboliteParentUpload(f):   
    '''
    class MParent model
    parent_type = models.CharField(max_length=20)
    expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    metabolite_name = models.ForeignKey(Metabolite)
    '''
    tic = time.clock()  
    exp_name = getExperimentName(f)
    parent_list = []  
    i=0
    for line in getData(f):
        if 'Parents' in line[0]:
            for parent in range(1,len(line)):
                parent_list.append(line[parent])
                print line[parent]
        else:
            i+=1
            print i
            for exp in range(1,len(line)):
                if line[exp] == '':        #js_exp_list = simplejson.dumps(exp_list)
                    continue                 
                else:
                    add_parent = MParent(parent_type =  parent_list[exp-1],
                                        expression =  line[exp],
                                        metabolite_name = Metabolite.objects.get(pk= line[0].upper()),
                                        experiment_name = exp_name)
                    add_parent.save()
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i
    
def metaboliteRILUpload(f): 
    '''    
    class MRIL model
    metabolite_name = models.ForeignKey(Metabolite)
    ril_name = models.CharField(max_length=20)
    ril_type = models.CharField(max_length=20)
    ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    '''
    tic = time.time()   
    getAllMetabolite(f,'RIL_type','metabolites')
    toc = time.time()            
    print 'in %f seconds' % (toc-tic)
    

def metaboliteLODUpload(f):
    '''
    class Experiment: 
    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    class MLOD model
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gxe = models.BooleanField(default=False)
    metabolite_name = models.ForeignKey(Metabolite)
    marker_name = models.ForeignKey(Marker)
    
    NEED to be approved: Occurred: Segmentation fault (core dumped) because of insufficient memory
    After saving 10839 records
    Solution:,index
    delete from qtl_lod where locus_identifier_id = 'AT2G36270';
    split the data set from that row into another one.
    Anyhow this part of function needs to be optimized by either using mmap or generator.
    '''
    tic = time.clock()   
    add_exp = Experiment()
    exp_name = getExperimentName(f)
    add_exp.experiment_name = exp_name
    add_exp.save()
    marker_list = []  
    i=0
    for line in getData(f):
        if 'QTL' in line[0]:
            for marker in range(1,len(line)):
                print line[marker]
                marker_list.append(line[marker])
        else:
            i+=1
            print i
            for gene in range(1,len(line)):
                add_mlod = MLOD(experiment_name = Experiment.objects.get(pk=exp_name),
                              metabolite_name = Metabolite.objects.get(pk= line[0].upper()),
                              marker_name =  Marker.objects.get(pk=marker_list[gene-1]),
                              LOD_score =  line[gene])
                add_mlod.save()
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i
    
def envLODUpload(f):   
    '''
    Experiment: 
    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    LOD:
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gene_name = models.ForeignKey(Gene)
    marker_name = models.ForeignKey(Marker)
    
    Occurred: Segmentation fault (core dumped) because of insufficient memory
    After saving 10839 records
    Solution:,index
    delete from qtl_lod where locus_identifier_id = 'AT2G36270';
    split the data set from that row into another one.
    Anyhow this part of function needs to be optimized by either using mmap or generator.
    
    Segmentation fault (core dumped)
    1st stop at 10532 AT2G33390
    2nd stop at 10591 AT4G26170

    '''
    tic = time.clock()   
    add_exp = Experiment()
    exp_name = getExperimentName(f)
    add_exp.experiment_name = exp_name
    add_exp.save()
    marker_list = []  
    i=0
    for line in getData(f):
        if 'QTLEnv' in line[0]:
            for marker in range(1,len(line)):
                print line[marker]
                marker_list.append(line[marker])
        else:
            i+=1
            print i
            print line[0]
            for gene in range(1,len(line)):
                add_lod = LOD(experiment_name = Experiment.objects.get(pk=exp_name),
                              locus_identifier = Gene.objects.get(pk= line[0].upper()),
                              marker_name =  Marker.objects.get(pk=marker_list[gene-1]),
                              gxe = True,
                              LOD_score =  line[gene])
                add_lod.save()    
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i    

def envMLODUpload(f):                
    '''
    Experiment: 
    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    class MLOD model
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gxe = models.BooleanField(default=False)
    metabolite_name = models.ForeignKey(Metabolite)
    marker_name = models.ForeignKey(Marker)
    
    Occurred: Segmentation fault (core dumped) because of insufficient memory
    After saving 10839 records
    Solution:,index
    delete from qtl_lod where locus_identifier_id = 'AT2G36270';
    split the data set from that row into another one.
    Anyhow this part of function needs to be optimized by either using mmap or generator.
    '''
    tic = time.clock()   
    exp_name = getExperimentName(f)
    add_exp = Experiment()
    add_exp.experiment_name = exp_name
    add_exp.save()
    marker_list = []  
    i=0
    for line in getData(f):
        if 'QTLEnv' in line[0]:
            for marker in range(1,len(line)):
                print line[marker]
                marker_list.append(line[marker])
        else:
            i+=1
            print i
            print line[0]
            for meta in range(1,len(line)):
                add_lod = MLOD(experiment_name = Experiment.objects.get(pk=exp_name),
                              metabolite_name = Metabolite.objects.get(pk= line[0].upper()),
                              marker_name =  Marker.objects.get(pk=marker_list[meta-1]),
                              gxe = True,
                              LOD_score =  line[meta])
                add_lod.save()    
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i    


def findPeak(gene_name,gxe_,exp):
    '''find the peak marker and expression value
    argument2 gexpression refers LOD.gxe(BooleanField: Ture or False in Django, 1 or 0 in MySQL) enviorment interaction. 
    '''
    exp_list = LOD.objects.filter(locus_identifier = gene_name,gxe= gxe_,experiment_name= exp).order_by('-LOD_score')    
    marker = exp_list[0].marker_name.marker_name
    peak = exp_list[0].LOD_score 
    return marker,peak

def findMPeak(metabolite,gxe_,exp):
    '''find the peak marker and expression value of metabolite
    argument2 gexpression refers MLOD.gxe(BooleanField: Ture or False in Django, 1 or 0 in MySQL) enviorment interaction. 
    '''
    exp_list = MLOD.objects.filter(metabolite_name = metabolite,gxe= gxe_,experiment_name= exp).order_by('-LOD_score')    
    marker = exp_list[0].marker_name.marker_name
    peak = exp_list[0].LOD_score 
    return marker,peak
def yield_gene(gene_list):
    for gene in gene_list:
        yield gene
    

def marker_plot(gene_name,gxe_,experiment):
    '''
    plot expression (LOD scores) along the chromosome, gxe True: response to enviormental changes influencing by genotypes. 
    '''
    lod_list = []
    marker_list = []
    
    marker_query_list = ExperimentMarker.objects.filter(experiment_name = experiment).values_list('marker_name',flat=True)
    marker_query_list_set = Marker.objects.filter(marker_name__in=marker_query_list).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
    for marker in marker_query_list_set:
        marker_list.append(marker)
        exp = LOD.objects.get(locus_identifier = gene_name,experiment_name = experiment,marker_name = marker,gxe=gxe_).LOD_score
        lod_list.append(float('{0:.2f}'.format(exp)))
    '''
    #calculate number of markers in each chromosome.
    #it is needed for marker plot 
    chr_marker_list = []
    for group in Marker.objects.values('marker_chromosome').annotate(chr_group = Count('marker_chromosome')):
        nr_marker_per_chr={}
        nr_marker_per_chr[int(group['marker_chromosome'])]= group['chr_group']
        chr_marker_list.append(nr_marker_per_chr)
    '''
    return marker_list,lod_list

def m_marker_plot(metabolite,gxe_,experiment):
    '''
    plot expression (LOD scores) along the chromosome, gxe True: response to enviormental changes influencing by genotypes. 
    '''
    lod_list = []
    marker_list = []
    marker_query_list = ExperimentMarker.objects.filter(experiment_name = experiment).values_list('marker_name',flat=True)
    marker_query_list_set = Marker.objects.filter(marker_name__in=marker_query_list).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
    for marker in marker_query_list_set :
        marker_list.append(marker)
        exp = MLOD.objects.get(metabolite_name = metabolite,marker_name = marker,gxe=gxe_,experiment_name = experiment).LOD_score
        lod_list.append(float('{0:.2f}'.format(exp)))
    '''
    #calculate number of markers in each chromosome.
    #it is needed for marker plot 
    chr_marker_list = []
    for group in Marker.objects.values('marker_chromosome').annotate(chr_group = Count('marker_chromosome')):
        nr_marker_per_chr={}
        nr_marker_per_chr[int(group['marker_chromosome'])]= group['chr_group']
        chr_marker_list.append(nr_marker_per_chr)
    '''
    return marker_list,lod_list


def geneUpdate(f):
    '''
    Updating gene chromosome location to the database and writing the mismatched genes into a file named missing.txt
    '''   
    with open('missing.txt', 'a') as save_file:
        
        for line in getData(f):
            if 'AT' in line[0]:
                print line[0],line[1],line[2],line[3],type(line[3])
                if Gene.objects.filter(locus_identifier = line[0]).exists():
                    gene = Gene.objects.get(locus_identifier = line[0])
                    gene.start = line[1]
                    gene.end = line[2]
                    if '1' in line[3]:
                        gene.strand = True
                    elif '0' in line[3]:
                        gene.strand = False
                    else:
                        print 'WRONG'
                    gene.save()
                else:
                    save_file.write(line[0]+'\n')
        save_file.close()
        
def decimalFormat(number):
    return float('{0:.2f}'.format(number))

def genotypeUpload(f):
    
    '''
    class Genotype: 
    marker_name = models.ForeignKey(Marker)
    ril_name= models.ForeignKey(RIL)
    genotype = models.CharField(max_length=5)
    Function: add the genotype of certain genomic marker location of a RIL population
    Comment: Unfortunately because of the previously defined data strucute and consideration of data redundency, Genotype has only one relationship to Marker (many to one)
    Data integrity: Because the expression value of gene/trait in offspring RIL90 is missing, therefore the genotype of each marker in RIL90 will not be considered.   
    '''
    tic = time.clock()   
    exp_name = getExperimentName(f)
    ril_list = []
    i=0
    for line in getData(f):
        if 'Genotype' in line[0]:
            for ril in range(1,len(line)):
                print line[ril]
                ril_list.append(line[ril])
        else:
            i+=1
            print i
            for marker in range(1,len(line)):
                add_geno = Genotype(marker_name = Marker.objects.get(pk=line[0]),
                              ril_name =  ril_list[marker-1],
                              genotype =  line[marker],
                              experiment_name = exp_name)
                add_geno.save()
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i

def getExperimentName(f):
    
    fname = f.name
    match = re.findall(r'\D(\d{4}\D)',fname)
    if len(match)==1:
        ind = fname.index(match[0])
        return fname[:ind+4]
    else:
        raise NameError("Invalid file name!")
    




if __name__ == '__main__':
    pass