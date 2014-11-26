'''
Created on Nov 26, 2014

@author: jiao
'''

import sys
import os
import urllib2
from bs4 import BeautifulSoup


sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Experiment,Gene,Marker,LOD,Parent,RIL,Metabolite,MParent,MRIL,MLOD,Genotype,ExperimentMarker,ExperimentGene,ExperimentRIL,ExperimentParent,ExperimentMetabolite
from qtl.models import TF_Family,TF_Family_Ref,Protein,Promoter_Binding_Site,TF


def promSeqUpload(f):
    i=0
    j = 0
    with open('PromoterSeqmissing.txt','w') as fo: # The gene exists in TAIR 9 (ATcisDB), but not yet been measured in any register experiment.
        
        with open(f,'r') as fi:
            for line in fi:
                
                content = line.split('\t')
                gene = content[0].upper()[:-2]
                seq = content[1].strip()
                
                if len(gene)==9 and Gene.objects.filter(locus_identifier=gene).exists():
                    i+=1
                    print gene
                    g = Gene.objects.get(locus_identifier = gene)
                    g.promoter_sequence = seq
                    g.save()
                else:
                    j+=1
                    g = Gene()
                    g.locus_identifier = gene
                    g.promoter_sequence = seq
                    g.save()
                    fo.write(gene+'\n')
    print i,j

def cdsSeqUpload(f): # NOT yet Used
    i = 0
    try:
        with open('CDSSeqmissing.txt','w') as fo: # The gene exists in TAIR 9 (ATcisDB), but not yet been measured in any register experiment.
            print 'openned output file'
            with open(f,'r') as fi:
                for line in fi:
                    i+=1
                    content = line.split('\t')
                    gene = content[0].upper()[:-2]
                    seq = content[1].strip()
                    if len(gene)==9 and Gene.objects.filter(locus_identifier=gene).exists():
                        print gene
                        g = Gene.objects.get(locus_identifier = gene)
                        g.coding_sequence = seq
                        g.save()
                    else:
                        fo.write(gene+'\n')
    except:
        print 'Didnot open.'

def promInfoUpload(fi):
    try:
        i = 0
        with open(fi,'r') as f:
            for line in f:
                i+=1
                content = line.strip().split('\t')
                gene = content[0].upper()[:-2]
                chr = int(content[1])
                ori = None
                if (content[2]) == '+':
                    ori = 1
                else:
                    ori = 0
                prom_start = int(content[3])
                prom_end = int(content[4])
                prom_type = content[5]
                prom_exam_org = content[6]
                if Gene.objects.filter(locus_identifier = gene).exists():
                    g = Gene.objects.get(locus_identifier = gene)
                    g.chromosome = chr
                    g.orientation = ori
                    g.promoter_start = prom_start
                    g.promoter_end = prom_end
                    g.promoter_type = prom_type
                    g.promoter_exam_org = prom_exam_org
                    g.save()
                else:
                    print 'NOT FOUND'
                print gene
                
        print i
    except IOError as e:
        print e.args
                

def getGenePosition(gene):
    tair_pre = 'http://arabidopsis.org/servlets/TairObject?type=locus&name='
    response = urllib2.urlopen(tair_pre+gene)
    soup = BeautifulSoup(response)
    all_tables = soup.find('table')
    trace_tds = all_tables.find_all('td')
    for td in trace_tds:
        if 'nuc_sequence' in td:
                            
            next_sib = td.findNextSibling('td')
            next_next_sib = next_sib.findNextSibling('td')
            pos = next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','').replace('\t','')
            ind = pos.index('-')
            start_ = pos[:ind].strip()
            end_ = pos[ind+1:].strip()
            orientation_ = next_next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','').replace('\t','')
            g = Gene.objects.get(locus_identifier = gene)
            g.start = start_
            g.end = end_
            g.orientation = orientation_
            g.save()
            break
        
def syncGene(gene):
    tair_pre = 'http://arabidopsis.org/servlets/TairObject?type=locus&name='
    response = urllib2.urlopen(tair_pre+gene)
    soup = BeautifulSoup(response)
    all_tables = soup.find('table')
    trace_tds = all_tables.find_all('td')
    for td in trace_tds:
        if 'nuc_sequence' in td:
                            
            next_sib = td.findNextSibling('td')
            next_next_sib = next_sib.findNextSibling('td')
            pos = next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','').replace('\t','')
            ind = pos.index('-')
            start_ = pos[:ind].strip()
            end_ = pos[ind+1:].strip()
            orientation_ = next_next_sib.text.encode('ascii','ignore').replace('bp','').replace(' ','').replace('\t','')
            g = Gene()
            g.locus_identifier = gene
            g.start = start_
            g.end = end_
            g.orientation = orientation_
            g.save()
            break
        
def updateGenePosition():
    gene_list = Gene.objects.filter(start__isnull = True)
    i = 0
    for gene in gene_list:
        i+=1
        gene_str = gene.locus_identifier.encode('ascii','ignore')
        getGenePosition(gene_str)
        print i,gene_str
    
def uploadTF_family(fi):  
    try:
        i = 0
        with open(fi,'r') as f:
            for line in f:
                i+=1
                content = line.strip().split('\t')
                tff = TF_Family()
                tff.tf_family_name = content[0].strip().encode('ascii','ignore')
                tff.tf_description = content[1].strip().encode('ascii','ignore')
                tff.save()
                print i,content[0]
    except IOError as e:
        print e.args

def uploadTF_family_ref(fi):
    '''
    TF_Family_Ref table
    tf_family_name = models.ForeignKey(TF_Family)
    reference = models.CharField(max_length = 200,blank = True)
    author = models.CharField(max_length = 200,blank = True)
    link = models.URLField(max_length = 200,blank = True)  
    '''
    try:
        i = 0
        with open(fi,'r') as f:
            for line in f:
                i+=1
                content = line.strip().split('\t')
                tff_name = content[0].strip().encode('ascii','ignore')
                tff_ref = content[1].strip().encode('ascii','ignore')
                tff_author = content[2].strip().replace('"','').decode('ascii','ignore')
                tff_link = content[3].strip()
                tff = TF_Family_Ref()
                tff.tf_family_name = TF_Family.objects.get(tf_family_name = tff_name)
                tff.reference = tff_ref
                tff.author = tff_author
                tff.link = tff_link
                tff.save()
                print i,content[0]
    except IOError as e:
        print e.args

def read_fasta(fp):
    name,seq = None,[]
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name,''.join(seq))
            name,seq = line,[]
        else:
            seq.append(line)
    if name: yield(name,''.join(seq))
    

def uploadProtein(TAIR_rel):
    '''
    Protein(models.Model):
    protein_name = models.CharField(max_length = 30,primary_key = True)
    locus_identifier = models.ForeignKey(Gene)
    protein_sequence = models.TextField(blank = True)
    protein_desc = models.TextField(blank = True)
    protein_symbol = models.CharField(max_length = 30,blank=True)
    '''
    try:
        with open(TAIR_rel,'r') as fp:
            i = 0
            for name,seq in read_fasta(fp):
                print name
                
                prot_name_list = name.split('|')
                prot_name = prot_name_list[0].strip()[1:]
                dot_index = prot_name.index('.')
                gene = prot_name[:dot_index]
                sym_ = prot_name_list[1].split(':')
                prot_symbol = sym_[1].strip()
                prot_desc = prot_name_list[2].strip()
                ###TAIR 6 has one extra alias attribute###
                #alias = prot_name_list[4].split(':','')
                #prot_alias = alias[1].strip()
                prot_seq = seq[:-1]
                if Gene.objects.filter(locus_identifier = gene).exists():
                    if not Protein.objects.filter(protein_name = prot_name).exists():
                        p = Protein()
                        p.protein_name = prot_name
                        p.locus_identifier = Gene.objects.get(locus_identifier = gene)
                        p.protein_sequence = prot_seq
                        p.protein_desc = prot_desc
                        if len(prot_symbol):
                            p.protein_symbol = prot_symbol
                        p.save()
                        print i,prot_name,prot_symbol
                else:
                    i+=1
                    syncGene(gene)
                    p = Protein()
                    p.protein_name = prot_name
                    p.locus_identifier = Gene.objects.get(locus_identifier = gene)
                    p.protein_sequence = prot_seq
                    p.protein_desc = prot_desc
                    if len(prot_symbol):
                        p.protein_symbol = prot_symbol
                    p.save()
                    print 'registered'
                    print i,prot_name,prot_symbol
                    
    except IOError as e:
        print e.args


    
def syncProtein(fi):
    try:
        with open(fi,'r') as f:
            i = 0
            for line in f:
                content = line.split('\t')
                prot = content[1].strip()
                if Protein.objects.filter(locus_identifier = prot).exists():
                    pass
                else:
                    print prot
            
    except IOError as e:
        print e.args       
        
def missProtein(fi): 
    
    try:
        with open(fi,'r') as f:
            i = 0
            for line in f:
                prot = line.strip()
                if Protein.objects.filter(locus_identifier = prot).exists():
                    pass
                else:
                    print prot
            
    except IOError as e:
        print e.args     

def uploadBS(fi):
    '''
    locus_identifier = models.ForeignKey(Gene)#
    binding_site_name = models.CharField(max_length = 30)
    binding_site_sequence = models.TextField()
    chr = models.IntegerField()
    start = models.IntegerField()
    end = models.IntegerField()
    tf_family_name = models.ForeignKey(TF_Family,blank = True)#
    motif = models.CharField(max_length = 20,blank = True)
    bs_atcisdb_color = models.CharField(max_length = 10,blank = True)
    fi: 'BindingSite.tbl' 

    '''
    try:
        with open(fi,'r') as f:
            i = 0
            for line in f:
                i+=1
                content = line.split('\t')
                tff = content[9].strip()  
                if not TF_Family.objects.filter(tf_family_name__iexact= tff).exists():                          
                    tff_ins = TF_Family()
                    tff_ins.tf_family_name = tff
                    tff_ins.save()            
                bs = Promoter_Binding_Site()
                gene_model = content[6].strip()
                dot_ind = gene_model.index('.')
                bs.locus_identifier = Gene.objects.get(locus_identifier = gene_model[:dot_ind].upper())
                bs.binding_site_name = content[1].strip()
                bs.binding_site_sequence = content[7].strip()
                bs.chr = int(content[2].strip())
                bs.start = int(content[4].strip())
                bs.end = int(content[5].strip())
                bs.tf_family_name = TF_Family.objects.get(tf_family_name__iexact = tff)#7 case insensitive
                if content[10].strip() is not 'NA':
                    bs.motif = content[10].strip()
                bs.bs_atcisdb_color = content[8].strip()
                bs.save()                                     
                print i,content[0].strip()
    except IOError as e:
        print e.args

def uploadTF(fi):
    '''
    upload ATtfDB predicted TF
    TF
    tf_family_name = models.ForeignKey(TF_Family)   
    locus_name = models.CharField(max_length=30)
    gene_name = models.CharField(max_length=20)
    description = models.CharField(max_length = 200,blank = True)
    motif = models.CharField(max_length = 100,blank = True)
    reference = models.CharField(max_length = 200,blank = True)
    author = models.CharField(max_length = 200,blank = True)
    link = models.URLField(max_length = 200,blank = True)  
    families_data.tbl
    '''
    
    try:
        with open(fi,'r') as f:
            i = 0
            for line in f:
                i+=1
                content = line.split('\t')
                tff = content[0].strip()  
                locus = content[1].strip().upper()
                gene = content[2].strip()
                desc = content[3].strip().replace('"','')
                if TF_Family.objects.filter(tf_family_name__iexact = tff).exists() and Gene.objects.filter(locus_identifier__iexact=locus).exists():
                    tf = TF()
                    tf.tf_family_name = TF_Family.objects.get(tf_family_name__iexact = tff)
                    tf.locus_name = locus
                    if gene is not 'NA':
                        tf.gene_name = gene
                    if desc is not 'NA':
                        tf.description = desc
                    tf.save()
                    print i,tff,locus
    except IOError as e:
        print e.args

def uploadTFRef(fi):  
    '''
    TF
    bs_name = models.CharField(max_length=100,blank = True)
    motif = models.CharField(max_length = 100,blank = True)
    reference = models.CharField(max_length = 200,blank = True)
    author = models.CharField(max_length = 200,blank = True)
    link = models.URLField(max_length = 200,blank = True)  
    'bindingsite_data.tbl'
    '''
    try:
        with open(fi,'r') as f:
            i = 0
            for line in f:
                i+=1
                content = line.split('\t')
                bs = content[0].strip()  
                locus = content[1].strip().upper()
                tff = content[2].strip()
                mot = content[3].strip()
                ref = content[4].strip()
                aut = content[5].strip().replace('"','')
                lin = content[6].strip()
                if TF_Family.objects.filter(tf_family_name__iexact = tff).exists() and Gene.objects.filter(locus_identifier__iexact=locus).exists():
                    tf_update =TF.objects.get(tf_family_name = tff,locus_name = locus)
                    tf_update.bs_name = bs
                    tf_update.motif = mot
                    tf_update.reference = ref
                    tf_update.author = aut
                    tf_update.link = lin
                    tf_update.save()
                print i, bs, locus,tff
                    
    except IOError as e:
        print e.args

if __name__ == '__main__':
    #promSeqUpload('PromoterSeq.tbl')
    #cdsSeqUpload('families_seq.tbl')
    #promInfoUpload('PromoterInfo.tbl')
    #updateGenePosition()
    #uploadTF_family('families_id.tbl')
    #uploadTF_family_ref('families_ref.tbl')
    #uploadProtein('TAIR9_pep_20090619.txt')
    #syncProtein('families_data.tbl')
    #uploadBS('BindingSite.tbl')
    #uploadBS('BindingSiteCut.tbl')
    #uploadTF('families_data.tbl')
    uploadTFRef('file/bindingsite_data.tbl')