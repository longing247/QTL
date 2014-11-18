'''
Created on Oct 24, 2014

@author: jiao
'''

import sys
import os
import math

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,ExperimentMarker
from django.db.models import Q

class EQTLMapping:
    '''
    classdocs
    '''
    
    def __init__(self, **kwargs):
        '''
        EQTL class Constructor
        self.gene[string]           query gene
        self.thld[float]            default LOD threshold is 3.84, corresponding to 0.0001446
        self.gxe [boolean]          set default value of gxe (experiment with environment permutation) to False
        Example of creating an instance: eqtls = EQTL(gene = 'AT1G01010',experiment = 'Ligterink_2014',thld = 10,gxe = True)
        
        '''
        if kwargs.get('gene'):
            if Gene.objects.filter(locus_identifier = kwargs.get('gene')).exists():#change it to validating from database
                self.gene = kwargs.get('gene')
                self.experiment = kwargs.get('experiment')
                self.thld = kwargs.get('thld',3.84)
                self.gxe = kwargs.get('gxe',False)
                self.lsi = kwargs.get('lsi',1.5) #lod support interval
                self.eqtls = self.getEQTL()
                self.eqtls_ci = self.getEQTLInterval()
                
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
        eqtls = Marker.objects.filter(Q(marker_name__in=parent_1_eqtl) | Q(marker_name__in= parent_2_eqtl)).order_by('marker_chromosome','marker_cm')
        
        return eqtls
    
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
        print self.gene,self.experiment,self.thld,self.gxe
  
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
        gene_list = []
        gene_queryset_list = LOD.objects.filter(experiment_name = self.experiment).values('locus_identifier').distinct()
        gene_all_list = Gene.objects.filter(locus_identifier__in=gene_queryset_list).order_by('chromosome','start')
        for chr in self.eqtls_ci:
            for interval in self.eqtls_ci[chr]:
                extend_gene = []
                int_start = interval['start']
                int_end = interval['end']
                extend_gene = Gene.objects.filter(locus_identifier__in=gene_all_list,chromosome = chr,start__gte=int_start,end__lte=int_end)
                print chr,len(extend_gene)
                gene_list.extend(extend_gene)
        return gene_list
                
if __name__=="__main__":
    #'AT1G03530'
    eqtl_analysis = EQTLMapping(gene = 'AT1G03530',experiment = 'Ligterink_2014',thld = 2.3,gxe = False,ci=1.5)
    eqtl_list = eqtl_analysis.eqtls
    print eqtl_list
    eqtl_ci = eqtl_analysis.eqtls_ci
    print eqtl_ci
    gene_list = eqtl_analysis.geneMapping()
    #print gene_list
    print len(gene_list)


    
    
    
        
    
        