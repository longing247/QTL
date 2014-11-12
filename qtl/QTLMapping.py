'''
Created on Oct 24, 2014

@author: jiao
'''

import sys
import os

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD


class EQTL:
    '''
    classdocs
    '''
    
    def __init__(self, **kwargs):
        '''
        EQTL class Constructor
        self.gene[string]           query gene
        self.thld[float]            default LOD threshold is 3.84, corresponding to 0.0001446
        self.gxe [boolean]          set default value of gxe (experiment with environment permutation) to False
        Example of creating an instance: eqtls = EQTL(gene = 'AT1G01010') or eqtls = EQTL(gene = 'AT1G01010',thld = 10,gxe = True)
        
        '''
        if kwargs.get('gene'):
            if Gene.objects.filter(locus_identifier = kwargs.get('gene')).exists():#change it to validating from database
                self.gene = kwargs.get('gene')
                self.thld = kwargs.get('thld',3.84)
                self.gxe = kwargs.get('gxe',False)
                self.eqtls = self.getEQTL()
                if self.eqtls:
                    self.ci = self.getEQTLsConfidentInterval()
                else:
                    self.ci = None
            else:
                raise ValueError("ERROR: GENE does not exist!")
        else:
            raise ValueError("ERROR: args ['gene'] is not assigned to the instance.")
        
    def getGene(self):
        return self.gene
    
    def getThld(self):
        return self.thld
    
    def getGxe(self):
        return self.gxe
        
    def getEQTL(self):
        '''
        detects eQTL from the marker-trait combination matrix of a given gene.
        '''
        eqtls = []
        bay_eqtl = LOD.objects.filter(locus_identifier = self.gene, LOD_score__gte = self.thld, gxe = self.gxe).values_list('marker_name', flat = True)
        eqtls.extend(bay_eqtl)
        sha_eqtl = LOD.objects.filter(locus_identifier = self.gene, LOD_score__lte = (-self.thld), gxe = self.gxe).values_list('marker_name', flat = True)
        eqtls.extend(sha_eqtl)
        return eqtls
    
    def getEQTLsConfidentInterval(self):
        '''
        retrieves a list of genome regions [confident interval] over the chromosome if any
        marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
        marker_chromosome = models.IntegerField()#1
        marker_cm = models.DecimalField(max_digits = 3, decimal_places =1)#64.6
        marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639
        '''
        
        ref_marker_list = list(Marker.objects.values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm'))
        ref_marker_list= [i.encode('ascii','ignore') for i in ref_marker_list]
        ref_marker_dic = {}
    
        chr1_markers = Marker.objects.filter(marker_chromosome = 1).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
        chr2_markers = Marker.objects.filter(marker_chromosome = 2).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
        chr3_markers = Marker.objects.filter(marker_chromosome = 3).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
        chr4_markers = Marker.objects.filter(marker_chromosome = 4).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
        chr5_markers = Marker.objects.filter(marker_chromosome = 5).values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
        ref_marker_dic['1'] = list([i.encode('ascii','ignore') for i in chr1_markers])
        ref_marker_dic['2'] = list([i.encode('ascii','ignore') for i in chr2_markers])
        ref_marker_dic['3'] = list([i.encode('ascii','ignore') for i in chr3_markers])
        ref_marker_dic['4'] = list([i.encode('ascii','ignore') for i in chr4_markers])
        ref_marker_dic['5'] = list([i.encode('ascii','ignore') for i in chr5_markers])
        
        chr_start = {1:3631,2:1871,3:4342,4:1180,5:1251}
        chr_end = {1:30425192,2:19696821,3:23459800,4:18584524,5:26970641}
        
        
        eqtl_marker = [] 
        for eqtl in self.eqtls:
            marker = Marker.objects.get(marker_name = eqtl)
            marker_info = {}
            marker_info['marker'] = eqtl.encode('ascii','ignore')
            marker_info['chr'] = int(marker.marker_chromosome)
            marker_info['cm'] = float(marker.marker_cm)*1000000
            marker_info['pos'] = float(marker.marker_phys_pos)
            marker_info['index'] = ref_marker_list.index(eqtl)
            eqtl_marker.append(marker_info)
            
        eqtl_marker = sorted(eqtl_marker, key = lambda k:k['index'])
        ci = []
        start = 0
        end= 0
        chr = 0
        for m in eqtl_marker:
            frag = {}
            if len(ci)==0: # initialize the first QTL region
                if len(eqtl_marker) == 1:
                    if m['marker'] == ref_marker_dic[str(m['chr'])][0]:
                        frag['chr'] = m['chr']
                        frag['start'] = chr_start[m['chr']]
                        end = m['index']+1  
                        frag['end'] = float(Marker.objects.get(marker_name = ref_marker_list[end]).values_list('marker_phys_pos',flat=True)[0])*1000000      
                        ci.append(frag)
                    elif m['marker'] == ref_marker_dic[str(m['chr'])][-1]: #last marker in one of the chromosome
                        frag['chr'] = m['chr']
                        frag['start'] = float(Marker.objects.get(marker_name = ref_marker_list[m['index']-1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        frag['end'] = chr_end(m['chr'])
                        ci.append(frag)
                    else:
                        frag['chr'] = m['chr']
                        frag['start'] = float(Marker.objects.get(marker_name = ref_marker_list[m['index']-1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        frag['end'] = float(Marker.objects.get(marker_name = ref_marker_list[m['index']+1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        ci.append(frag)
                else:
                    if m['marker'] == ref_marker_dic[str(m['chr'])][0].encode('ascii','ignore'):
                        frag['chr'] = m['chr']
                        frag['start'] = chr_start[m['chr']] 
                        chr = m['chr']
                        end = m['index']+1     
                        
                    elif m['marker'] == ref_marker_dic[str(m['chr'])][-1].encode('ascii','ignore'): #last marker in one of the chromosome
                        frag['chr'] = m['chr']
                        frag['start'] = float(Marker.objects.filter(marker_name = ref_marker_list[m['index']-1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        frag['end'] = chr_end[m['chr']]               
                    else:
                        frag['start'] = float(Marker.objects.filter(marker_name = ref_marker_list[m['index']-1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        chr = m['chr']
                        end = m['index']+1   
                    
                
            
            if m['marker'] == eqtl_marker[-1]:  # last eqtl in the list
                if not frag['end']: 
                    frag['end'] = float(Marker.objects.get(marker_name = ref_marker_list[end]).values_list('marker_phys_pos',flat=True)[0])*1000000 
                ci.append(frag) 
                          
                    
            else: # if there is any other eQTLs detected. Either jump to the next eQTL or extend region with the next marker         
                chr_li = [i['chr'] for i in ci]
                chr_li = set(chr_li)
                if m['chr'] not in chr_li:# jump to another chromosome
                    if not frag['end']:
                        frag['end'] = float(Marker.objects.get(marker_name = ref_marker_list[end]).values_list('marker_phys_pos',flat=True)[0])*1000000 
                    ci.append(frag) # push the last eQTLs into the list
                    frag['chr'] = m['chr']
                    if m['marker'] == ref_marker_dic[str(m['chr'])][0].encode('ascii','ignore'):
                        frag['start'] = chr_start[m['chr']]
                        end = m['index']+1
                    elif m['marker'] == ref_marker_dic[str(m['chr'])][-1].encode('ascii','ignore'):
                        frag['start'] = float(Marker.objects.get(marker_name = eqtl_marker[m['index']-1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        frag['end'] = chr_end[m['chr']]
                    else:
                        frag['start'] = float(Marker.objects.get(marker_name = eqtl_marker[m['index']-1]).values_list('marker_phys_pos',flat=True)[0])*1000000
                        end = m['index']+1                
                    
                else: #extension 
                    if m['marker'] is not ref_marker_dic[str(m['chr'])][-1].encode('ascii','ignore'):
                        end = m['index']+1
                    else:
                        frag['end'] = chr_end[m['chr']]
            
        return ci
        
        

if __name__=="__main__":
    eqtl_analysis = EQTL(gene = 'AT3G50500')
    for marker in eqtl_analysis.ci:
        print marker
        
    
        