'''
Created on Oct 24, 2014

@author: jiao
'''

import sys
import os

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import LOD
from AGRIS import AGRISDB


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
        atcisdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtcisDB/atcisview.html?id='
        if kwargs.get('gene'):
            if AGRISDB.isATGene(atcisdb_pre,kwargs.get('gene')):
                self.gene = kwargs.get('gene')
                self.thld = kwargs.get('thld',3.84)
                self.gxe = kwargs.get('gxe',False)
            else:
                raise ValueError("ERROR: GENE does not exist!")
        else:
            raise ValueError("ERROR: args ['gene'] is not assigned to the instance.")
        
    def getGene(self):
        return self. gene
    
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
    
    

if __name__=="__main__":
    eqtl_analysis = EQTL(gene = 'AT3G50500')
    eqtls = eqtl_analysis.getEQTL()
    print len(eqtls)
    for eqtl in eqtls:
        print eqtl
        