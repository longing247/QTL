'''
Created on Oct 24, 2014

@author: jiao
'''

import urllib2
from bs4 import BeautifulSoup


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
        
        
if __name__=="__main__":
    
    atcisdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtcisDB/atcisview.html?id='
    attfdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtTFDB/tfsummary.html?locusid='
    attffamdb_pre = 'http://arabidopsis.med.ohio-state.edu/AtTFDB/tfbrowse.html?fam='       
    #print AGRISDB.isTF(atcisdb_pre,attfdb_pre,'AT3G24650')
    #print AGRISDB.isTF(atcisdb_pre,attfdb_pre,'AT3G24650T')
    #print AGRISDB.isTF(attfdb_pre,'AT3G50500')
    #print AGRISDB.hasTFBS(atcisdb_pre,'AT1G00010')
    #print AGRISDB.hasTFBS(atcisdb_pre,'AT3G24650')
    #AGRISDB.getTFBS(atcisdb_pre,'AT1G10100')
    #print AGRISDB.isATGene(atcisdb_pre,'AT1G10100T')
    #print AGRISDB.isATGene(atcisdb_pre,'AT1G10100')
    #AGRISDB.getTF(attffamdb_pre,'Homeobox')