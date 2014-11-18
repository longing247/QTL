'''
Created on Oct 24, 2014

@author: jiao
'''

import urllib2
from bs4 import BeautifulSoup
import sys
import os
import time
from lxml import etree
from io import StringIO

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD

class AGRISDB(object):
    '''
    a crawler fetches information from AtCisDB and AtTFDB databases of Gene Regulatory Information Server (AGRIS)
    '''
    def __init__(self, params):
        '''
        Constructor
        '''
        
    @staticmethod
    def isATGene(atcisdb_pre,gene):
        try: 
            response = urllib2.urlopen(atcisdb_pre+gene)
            soup = BeautifulSoup(response)
            if soup.find_all('table')[1].find_all('tr')[1].find_all('td')[1].string.strip(): 
                return True
            else:
                return False
        except urllib2.HTTPError, e:
            print 'ERROR code - %s' % e.code
         
    @staticmethod
    def isTF(atcisdb_pre,attfdb_pre,gene):
        '''
        returns True if the given gene is Transcription factor registered in the AtTFDB of AGRIS
        '''
        if AGRISDB.isATGene(atcisdb_pre,gene):
            try:
                response = urllib2.urlopen(attfdb_pre+gene)
                soup = BeautifulSoup(response)
                if soup.find_all('table')[1].find_all('tr')[0].find_all('td')[1].string:
                    response.close() 
                    return True
                else:
                    response.close() 
                    return False
            except urllib2.HTTPError, e:
                print 'ERROR code - %s' % e.code
        else:
            raise NameError("Gene does not exist")
        
    @staticmethod
    def hasTFBS(atcisdb_pre,gene):
        '''
        returns True if promoter region of the give gene has one or more transcription binding sites registered in the AtCisDB databases of AGRIS
        '''
        if AGRISDB.isATGene(atcisdb_pre,gene):
            try:
                response = urllib2.urlopen(atcisdb_pre+gene)
                soup = BeautifulSoup(response)
        
                if len(soup.find_all('table')[4].find_all('tr')) == 1:
                    response.close()  
                    return False
                else:
                    response.close()  
                    return True
            except urllib2.HTTPError, e:
                print 'ERROR code - %s' % e.code
        else:
            raise NameError("Gene does not exist")
        
    @staticmethod
    def getTFBS(atcisdb_pre,gene):
        '''
        returns list of Binding sites (type of dict) of the given gene
        for example AT1G10100: [{'bs': 'MyB1 binding site motif' , 'start': 3302233 , 'end': 3301951 ,'sequence': 'atccaacc','tff':'MYB',... }]
        '''
        if AGRISDB.isATGene(atcisdb_pre,gene):
            try:
                response = urllib2.urlopen(atcisdb_pre+gene)
                soup = BeautifulSoup(response)
                bs_family = []
                for tr in soup.find_all('table')[4].find_all('tr')[1:]:
                    bsf = {}
                    bsf['bs'] = tr.find_all('td')[0].string.encode('ascii','ignore').strip()
                    bsf['start'] = tr.find_all('td')[1].string.encode('ascii','ignore').strip()
                    bsf['end'] = tr.find_all('td')[2].string.encode('ascii','ignore').strip()
                    bsf['sequence'] = tr.find_all('td')[3].string.encode('ascii','ignore').strip()
                    if tr.find_all('td')[4].a:  # check if it exists
                        bsf['tff'] = tr.find_all('td')[4].a.string.encode('ascii','ignore').strip() 
                    else:
                        bsf['tff'] = None
                    bs_family.append(bsf)
            except urllib2.HTTPError, e:
                print 'ERROR code - %s' % e.code
        else:
            raise NameError("Gene does not exist")
    
    @staticmethod
    def getTF(attffamdb_pre,tff):
        '''
        returns list of genes encoding proteins belongs to the transcription factors family (tff)
       
        '''
        try:
            response = urllib2.urlopen(attffamdb_pre+tff)
            soup = BeautifulSoup(response)
            tf = []
            for tr in soup.find_all('table')[4].find_all('tr')[1:]:
                tf.append(tr.find_all('td')[0].a.string.encode('ascii','ignore').strip())
            return tf
        except urllib2.HTTPError, e:
                print 'ERROR code - %s' % e.code
    @staticmethod
    def getGenePosition(gene):
        '''
        fetch gene physical position from TAIR
        return datatype: tuple, (start,end,orientation)
        
        '''
        try:
            tic = time.clock()  
            tair_pre = 'http://arabidopsis.org/servlets/TairObject?type=locus&name='
            gene = gene.encode('ascii','ignore').strip()
            url = tair_pre+gene           
            response = urllib2.urlopen(tair_pre+gene)
            soup = BeautifulSoup(response)
            all_tables = soup.find('table')
            trace_tds = all_tables.find_all('td')
            start = None
            end = None
            orientation = None
            1
            for td in trace_tds:
                if 'nuc_sequence' in td:          
                    next_sib = td.findNextSibling('td')
                    next_next_sib = next_sib.findNextSibling('td')
                    pos = next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','')
                    ind = pos.index('-')
                    start = pos[:ind]
                    end = pos[ind+1:]
                    orientation = next_next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','')
            
            toc = time.clock()            
            print 'in %f seconds' % (toc-tic)
            return start,end,orientation
            
        except urllib2.HTTPError, e:
            print 'ERROR code - %s' % e.code
    @staticmethod
    def getPosition(tair_pre,input,output):
        #
        try:
            tic = time.clock()
            i=0
            gene_list = []
            with open(input,'r') as fi:
                for line in fi:
                    gene_list.append(line)
            with open(output,'w') as fo:   
                for gene in gene_list:
                    i+=1
                    gene = gene.encode('ascii','ignore').strip()
                    response = urllib2.urlopen(tair_pre+gene)
                    soup = BeautifulSoup(response)
                    #bs_family = []
                    all_tables = soup.find('table')
                    trace_tds = all_tables.find_all('td')
                    for td in trace_tds:
                        if 'nuc_sequence' in td:
                            
                            next_sib = td.findNextSibling('td')
                            next_next_sib = next_sib.findNextSibling('td')
                            pos = next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','').replace('\t','')
                            ind = pos.index('-')
                            start = pos[:ind].strip()
                            end = pos[ind+1:].strip()
                            orientation = next_next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','').replace('\t','')
                            fo.write(gene+'\t'+start+'\t'+end+'\t'+orientation+'\n')
                            print i,gene,start,end,orientation
                            break
            
            
            toc = time.clock()            
            print 'in %f seconds' % (toc-tic)
            print '%d records added' %i
                            
        except urllib2.HTTPError, e:
            print 'ERROR code - %s' % e.code

    @staticmethod
    def getPromotorSeq(gene):
        '''
        fetch gene promotor nucleotide sequence from TAIR
        return datatype: string
        
        '''
        try:
            tic = time.clock()  
            agris_prom_pre = 'http://arabidopsis.med.ohio-state.edu/AtcisDB/getpromseq.html?id='    
            gene = gene.encode('ascii','ignore').strip()    
            url = agris_prom_pre+gene
            response = urllib2.urlopen(url+gene)
   
            soup = BeautifulSoup(response)
            print soup
            
            toc = time.clock()            
            print 'in %f seconds' % (toc-tic)
            
        except urllib2.HTTPError, e:
            print 'ERROR code - %s' % e.code
    @staticmethod
    def updatePhysicalPos(f):
        i = 0
        with open(f,'r') as fi:
            for line in fi:
                i+=1
                linele = line.strip().split('\t')
                #print linele[0],linele[1],linele[2],linele[3]
                
                g = Gene.objects.get(locus_identifier = linele[0].strip())
                g.start = linele[1].strip()
                g.end = linele[2].strip()
                if linele[3]=='forward':
                    g.orientation = 1
                elif linele[3]=='reverse':
                    g.orientation = 0
                g.save()
        print i
                   
                    
                    
    @staticmethod
    def dismatch(f,f1):
        gene_list = [] 
        with open(f,'r') as fi:
            for line in fi:
                linele = line.strip().split('\t')
                for i in range(0,len(linele)):
                    gene_list.append(linele[0])
                    print linele[0]
        with open(f1,'r') as fi:
            for line in fi:
                linele = line.strip().split('\t')
                for i in range(0,len(linele)):
                    gene_list.append(linele[0])
                    print linele[0]
        gene_db = Gene.objects.values_list('locus_identifier',flat=True)
        gene_db = [gene.encode('ascii','ignore') for gene in gene_db]
        i = 0
        with open('dismatch.txt','w') as fo:
            for gene in gene_db:
                if gene not in gene_list:
                    i+=1
                    print i,gene
                    fo.write(gene+'\n')
                else:
                    continue
    @staticmethod        
    def findGene(fo):
        gene_db = Gene.objects.values_list('locus_identifier',flat=True)
        gene_db = [gene.encode('ascii','ignore') for gene in gene_db] 
        with open(fo,'w') as ff:
            for gene in gene_db:
                if len(gene) !=9:
                    print gene,len(gene)
                    ff.write(gene+'\n')
                else:
                    continue
    @staticmethod    
    def missingGene():
        '''
        fills up physical position of genes which are not given in TAIR9-genes with position.xlsx with fake data.
        The missing gene profile was saved in the parent directory named missing.txt.
        '''
        miss_genes = Gene.objects.filter(start__isnull=True)
        for gene in miss_genes:
            prefix = gene.locus_identifier.encode('ascii','ignore').strip()[:-3]
            suffix = int(gene.locus_identifier.encode('ascii','ignore').strip()[-3:])
            print gene.locus_identifier 
            for i in range(500):
                k=1
                suf = suffix+i+1
                next_gene = prefix+str(suf)
                if Gene.objects.filter(locus_identifier = next_gene).count() ==1:
                    next_=Gene.objects.get(locus_identifier = next_gene)
                    if not next_.start:
                        k+=1
                        continue
                    else:
                        g = Gene.objects.get(locus_identifier = gene.locus_identifier)
                        g.start = next_.start-(600*k)
                        g.end = next_.start-(300*k)
                        g.strand = next_.strand
                        g.save()
                        break

    @staticmethod    
    def getChromosome():
        '''
        populates the chromosome number for each Gene objects.
        '''
        gene_list = Gene.objects.all().values_list('locus_identifier',flat='True');
        for gene in gene_list:
            g = Gene.objects.get(locus_identifier = gene)
            g.chromosome =  int(gene[2])
    @staticmethod 
    def test(chr_eqtls,local_eqtls):
        ref_list = [] 
        for i in range(len(local_eqtls)):
            ref_index = chr_eqtls.index(local_eqtls[i])
            ref_list.append(ref_index)
        
        start,end = ref_list[0],ref_list[0]
        count = start
        for item in ref_list:
            if not count == item:
                yield start, end
                start,end = item,item
                count = item
            end = item
            count +=1
        yield start,end            


if __name__=="__main__":
    
    atcisdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtcisDB/atcisview.html?id='
    attfdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtTFDB/tfsummary.html?locusid='
    attffamdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtTFDB/tfbrowse.html?fam='  
    tair_pre = 'http://arabidopsis.org/servlets/TairObject?type=locus&name='
    #gene_list = AGRISDB.getGenes()
    #AGRISDB.getPosition(tair_pre,'dismatch.txt','output_dismatch.txt')  
    #AGRISDB.dismatch('physicalpos.txt','output_dismatch.txt') 
    #AGRISDB.updatePhysicalPos('physicalpos.txt')
    #AGRISDB.updatePhysicalPos('output_dismatch.txt')
    #AGRISDB.missingGene()
    #AGRISDB.getChromosome()
    #AGRISDB.getPromotorSeq('AT2G30420')  
    #AGRISDB.updatePhysicalPos('physicalpos.txt')
    #AGRISDB.dismatch('physicalpos.txt')
    #AGRISDB.findGene('wronglocus.txt')
    #AGRISDB.getPosition(tair_pre,'AT4G37710')
    #print AGRISDB.isTF(atcisdb_pre,attfdb_pre,'AT3G24650')
    #print AGRISDB.isTF(atcisdb_pre,attfdb_pre,'AT3G24650T')
    #print AGRISDB.isTF(1attfdb_pre,'AT3G50500')
    #print AGRISDB.hasTFBS(atcisdb_pre,'AT1G01010')
    #print AGRISDB.hasTFBS(atcisdb_pre,'AT3G24650')
    #AGRISDB.getTFBS(atcisdb_pre,'AT3G50500')
    #print AGRISDB.isATGene(atcisdb_pre,'AT1G10100T')
    #print AGRISDB.isATGene(atcisdb_pre,'AT1G10100')
    #AGRISDB.getTF(attffamdb_pre,'Homeobox')
    a = [1,2,3,4,5,6,7,8,9,0]
    b = [2,3,5,6,0]

    print (list(AGRISDB.test(a,b)))