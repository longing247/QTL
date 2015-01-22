'''
Created on Oct 24, 2014

@author: jiao
'''

import sys
import os
import math
import time
import json

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,ExperimentMarker,ExperimentGene
from django.db.models import Q

class EQTLMapping:
    '''
    classdocs
    '''
    
    def __init__(self, **kwargs):
        '''
        EQTL class Constructor
        self.gene[string]           query gene
        self.thld[float]            default LOD threshold is 4.26(4.2644960, corresponding to a p value of 5.43881141028e-05 result from Statistic package.
        self.gxe [boolean]          set default value of gxe (experiment with environment permutation) to False
        Example of creating an instance: eqtls = EQTL(gene = 'AT1G01010',experiment = 'Ligterink_2014',thld = 10,gxe = True)
        
        '''
        if kwargs.get('gene'):
            if Gene.objects.filter(locus_identifier = kwargs.get('gene')).exists():#change it to validating from database
                self.gene = kwargs.get('gene')
                self.experiment = kwargs.get('experiment')
                self.thld = kwargs.get('thld',4.26)
                self.gxe = kwargs.get('gxe',False)
                self.lsi = kwargs.get('lsi',1.5) #lod support interval
                self.eqtls = self.getEQTL()
                self.eqtls_ci = self.getEQTLInterval()
                if self.eqtls:
                    #self.candidateTF = self.geneMapping()
                    self.cisTransProfile = self.cisTrans()
                
                
            else:
                raise ValueError("ERROR: GENE does not exist!")
        else:
            raise ValueError("ERROR: args ['gene'] is not assigned to the instance.")
        
    def getGeneName(self):
        return self.gene
    
    def getExperiment(self):
        return self.experiment
    
    def getThld(self):
        return self.thld
    
    def getGxe(self):
        return self.gxe
    
    def getGene(self):
        return Gene.objects.get(locus_identifier = self.gene)
    
    def getEQTL(self):
        '''
        detects eQTL from the marker-trait combination matrix of a given gene.
        caution: unicode needs to be converted to string.
        '''
        
        
        parent_1_eqtl = LOD.objects.filter(locus_identifier = self.gene,experiment_name = self.experiment, LOD_score__gte = self.thld, gxe = self.gxe).values('marker_name')
        parent_2_eqtl = LOD.objects.filter(locus_identifier = self.gene,experiment_name = self.experiment, LOD_score__lte = (-self.thld), gxe = self.gxe).values('marker_name')
        #peak = LOD.objects.filter(locus_identifier = self.gene,experiment_name = self.experiment, LOD_score__lte = (-self.thld), gxe = self.gxe).values('marker_name','LOD_score')
        eqtls = Marker.objects.filter(Q(marker_name__in=parent_1_eqtl) | Q(marker_name__in= parent_2_eqtl)).order_by('marker_chromosome','marker_cm').values_list('marker_name',flat=True)
        
        return list(eqtls)
    
    def getEQTLInterval(self):
        '''
        retrieves a list of  [confident interval] over the chromosome if any
        methods used: 1.5 LOD support interval for backcross or 1.8 LOD support interval for intercross
        marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
        marker_chromosome = models.IntegerField()#1
        marker_cm = models.DecimalField(max_digits = 3, decimal_places =1)#64.6
        marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639
        '''
        
        gene = self.getGene()
        chr = gene.chromosome
        #start = gene.start
        #end = gene.end
        #print self.gene,self.experiment,self.thld,self.gxe
        chr_start = {1:3631,2:1871,3:4342,4:1180,5:1251}
        chr_end = {1:30425192,2:19696821,3:23459800,4:18584524,5:26970641}
        
        marker_queryset_list = ExperimentMarker.objects.filter(experiment_name = self.experiment).values_list('marker_name',flat=True)
        #marker_queryset_list = LOD.objects.filter(experiment_name = self.experiment).values('marker_name').distinct()
        marker_list = Marker.objects.filter(marker_name__in=marker_queryset_list).order_by('marker_chromosome','marker_cm')
        
        if len(self.eqtls)!=0:
            chr_interval = {}
            for i in range(1,6):
                if Marker.objects.filter(marker_chromosome = i, marker_name__in=self.eqtls).exists():
                    chr_eqtls_unicode = list(Marker.objects.filter(marker_chromosome = i, marker_name__in=self.eqtls).values_list('marker_name',flat=True).order_by('marker_cm')) #unicode
                    chr_eqtls = [marker.encode('ascii','ignore') for marker in chr_eqtls_unicode] # convert to string
                    chr_all_markers_unicode = list(Marker.objects.filter(marker_chromosome=i,marker_name__in=marker_list).values_list('marker_name',flat=True).order_by('marker_cm'))
                    chr_all_markers = [marker.encode('ascii','ignore') for marker in chr_all_markers_unicode]
                    ref_list = []    
                    for j in range(len(chr_eqtls)):
                        ref_index = chr_all_markers.index(chr_eqtls[j])
                        ref_list.append(ref_index)
                    #group the continuous number in the list
                    interval_list = list(self.ranges(ref_list))
                    inv_li = []
                    
                    for eqtl in interval_list:
                        inv = {}
                        start,end = eqtl[0],eqtl[1]                    
                        if start==0:
                            inv['start'] = chr_start[i]
                        else:
                            prev = Marker.objects.get(marker_name = chr_all_markers[start-1]) 
                            prev_pos = int(prev.marker_phys_pos*1000000)
                                            
                            #define cis-trans eQTL by max{-log10P}
                            curr = Marker.objects.get(marker_name = chr_all_markers[start])                         
                            curr_pos = int(curr.marker_phys_pos*1000000)
                            prev_lod = float(math.fabs(LOD.objects.filter(locus_identifier_id = self.gene,marker_name = prev.marker_name,experiment_name = self.experiment,gxe = self.gxe).values_list('LOD_score',flat = True)[0]))
                            curr_lod = float(math.fabs(LOD.objects.filter(locus_identifier_id = self.gene,marker_name = curr.marker_name,experiment_name = self.experiment,gxe = self.gxe).values_list('LOD_score',flat = True)[0]))                                        
                            left_interval = int(curr_pos-self.lsi*float(curr_pos-prev_pos)/float(curr_lod-prev_lod))
                            if left_interval < prev_pos:
                                inv['start'] = prev_pos
                            else:
                                inv['start'] =  left_interval       
                        if end ==len(chr_all_markers)-1:
                            inv['end'] = chr_end[i]
                        else:    
                            next = Marker.objects.get(marker_name = chr_all_markers[end+1]) 
                            next_pos = int(next.marker_phys_pos*1000000)
                            curr = Marker.objects.get(marker_name = chr_all_markers[start])
                            curr_pos = int(curr.marker_phys_pos*1000000)
                            next_lod = float(math.fabs(LOD.objects.filter(locus_identifier_id = self.gene,marker_name = next.marker_name,experiment_name = self.experiment,gxe = self.gxe).values_list('LOD_score',flat = True)[0]))
                            curr_lod = float(math.fabs(LOD.objects.filter(locus_identifier_id = self.gene,marker_name = curr.marker_name,experiment_name = self.experiment,gxe = self.gxe).values_list('LOD_score',flat = True)[0]))
                            right_interval = int(curr_pos+self.lsi*float(next_pos-curr_pos)/float(curr_lod-next_lod))
                            if right_interval > next_pos: 
                                inv['end'] = next_pos
                            else:
                                inv['end'] = right_interval  
                        inv_li.append(inv)       
                    chr_interval[i] = inv_li     

            return chr_interval
    
    @staticmethod 
    def ranges(series_list):
        '''
        identify groups of continuous numbers in a list and group the continuous numbers as a sub-list
        for example,
        [6, 7, 8, 9, 10, 11]
        [(6, 11)]
        '''
        start,end = series_list[0],series_list[0]
        for n in series_list[1:]:
            if n-1 ==end: # Part of the group, bump the end
                end = n
            else: # Not part of the group, yild c8urrent goup and start a new
                yield start,end
                start = end = n
        yield start,end #Yield the last group
    
    def geneMapping(self):
        '''
        identify the genes(candidate TF) that are located in the confident interval
        '''
        tic1 = time.time()
        gene_list = []
        gene_queryset_list = LOD.objects.filter(experiment_name = self.experiment).values('locus_identifier').distinct()
        gene_all_list = Gene.objects.filter(locus_identifier__in=gene_queryset_list).order_by('chromosome','start')
        for chr in self.eqtls_ci:
            for interval in self.eqtls_ci[chr]:
                extend_gene = []
                int_start = interval['start']
                int_end = interval['end']
                extend_gene = Gene.objects.filter(locus_identifier__in=gene_all_list,chromosome = chr,start__gte=int_start,end__lte=int_end)
                gene_list.extend(extend_gene)
        
        toc1 = time.time()
        print 'geneMapping query from LOD in %f seconds %s' % (toc1-tic1,self.gene)
        return gene_list
    
    def cisTrans(self):
        '''
        Varifies whether a gene/trait has eQTL(s) under a certain threshold, and define whether the eQTL is cis or trans.
        Due to the limitation of the genetic marker resolution, at this stage we are not able to distinguish local-cis eQTL and local-trans eQTL. 
        Therefore local-trans eQTL are considered to be cis-eQTL.
        {1: [{'start': 9887359, 'end': 15926702}, {'start': 25827433, 'end': 28336537}], 2: [{'start': 6596380, 'end': 14491108}]}
        '''
        trait = self.getGene()
        chr = trait.chromosome
        start = trait.start
        end = trait.end
        cis = False
        trans = False
        nrOfTrans = 0 # individual
        if self.eqtls:
            for qtl_chr in self.eqtls_ci:
                if chr == qtl_chr:
                    for qtlPos in self.eqtls_ci[qtl_chr]:
                        if start in range(qtlPos['start'],qtlPos['end']) or end in range(qtlPos['start'],qtlPos['end']):
                            cis = True
                        else:
                            trans = True
                            nrOfTrans+=1
                else:
                    trans = True
                    nrOfTrans+=1
                    
            return cis,trans,nrOfTrans

def genomeWideEQTLMapping(_exp,_th,_gxe,_lsi):
    '''
    Genome wide EQTL mapping returns a list of dictionary of genes that have significant EQTL profile.
    '''
    tic = time.time()
    trait_list = LOD.objects.filter(experiment_name = _exp).values_list('locus_identifier',flat=True).distinct() # Which gene expression was measured.
    eQTL_profile_list = []
    for trait in trait_list:
        eQTL_profile = {}
        trait_encode = trait.encode('ascii','ignore')
        eQTL = EQTLMapping(gene = trait_encode,experiment = _exp,thld =_th,lsi=_lsi,gxe = _gxe)
        print eQTL.gene,eQTL.experiment,eQTL.thld,eQTL.gxe,eQTL.lsi
        if eQTL.eqtls:
            print eQTL.eqtls,eQTL.eqtls_ci
            #print eQTL.candidateTF
            print eQTL.cisTransProfile
            eQTL_profile['trait'] = trait_encode
            eQTL_profile['eQTL'] = eQTL.eqtls
            eQTL_profile['ci'] = eQTL.eqtls_ci
            #eQTL_profile['candidatesOfTF'] = eQTL.candidateTF
            eQTL_profile['cisTransProfile'] = eQTL.cisTransProfile
            eQTL_profile_list.append(eQTL_profile)
    toc = time.time()
    print 'Genome-wide eQTL mapping in %f seconds' % (toc-tic)
    return eQTL_profile_list     

if __name__=="__main__":
    #'AT1G03530'
    #eqtl_analysis = EQTLMapping(gene = 'AT1G01200',experiment = 'Ligterink_2014',thld = 3.85,gxe = False,ci=1.5)
    #eqtl_analysis = EQTLMapping(gene = 'AT3G01010',experiment = 'Ligterink_2014',thld = 2.3,gxe = False,ci=1.5)
    #eqtl_list = eqtl_analysis.eqtls
    #print eqtl_list 
    #eqtl_ci = eqtl_analysis.eqtls_ci
    #print eqtl_ci
    #gene_list = eqtl_analysis.geneMapping()
    #print gene_list
    #print len(gene_list)
    
    #Genome-wide eQTL mapping
    _exp = 'Ligterink_2014'
    _th = 3.85
    _gxe = False
    _lsi = 2
    genome_wide_eQTL_list = genomeWideEQTLMapping(_exp,_th,_gxe,_lsi)
    #with open('genome_wide_eQTL_mapping_Ligterink_2014_gxe0_3.85.txt','w') as outfile:
    with open('genome_wide_eQTL_mapping_Ligterink_2014_gxe0_3.85_2.txt','w') as outfile:
    #with open('genome_wide_eQTL_mapping_Ligterink_2014_gxe1_2.7.txt','w') as outfile:
    #with open('genome_wide_eQTL_mapping_Ligterink_2014_gxe1_2.7_2.txt','w') as outfile:
    #with open('genome_wide_eQTL_mapping_Keurentjes_2007_gxe0_3.3_1.txt','w') as outfile:
    #with open('genome_wide_eQTL_mapping_Keurentjes_2007_gxe0_3.3_2.txt','w') as outfile:
    #with open('genome_wide_eQTL_mapping_Snoek_2012_gxe1_3.01.txt','w') as outfile:
    #with open('genome_wide_eQTL_mapping_Snoek_2012_gxe1_3.01_2.txt','w') as outfile:
        json.dump(genome_wide_eQTL_list,outfile,sort_keys = True, indent = 4)
    
    
    
        