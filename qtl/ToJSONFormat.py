'''
Created on Oct 15, 2014

@author: jiao
'''
import time,os,json

from django.core.serializers.json import DjangoJSONEncoder

from models import Gene,Marker,LOD
from views import getMarkersList,getChrMarkers



def toJSONFormat(lod_thld):
        '''
        Formulate JSON format object which will be later used to generate file for eQTL mapping plot
        '''

        # output_dic will be used to generate JSON file.
        tic = time.clock()
        output_dic = {}
        chr_nr = 5
        #define KEY chrnames
        # alternative way if it is used to another organism
        # for efficiency: define it derectly, or comment the following line
        # chr_size = Marker.objects.values('marker_chromosome').dictinct().count() # and use chr_size to generate ['1','2','3','4','5']
        output_dic["chrnames"]= ['1','2','3','4','5']
        
        # define KEY chr
        # it can be also achieved by calling
        # chr_start = Gene.objects.filter('chromosome'=1).orderby('start')[0].values('start') etc...
        output_dic["chr"] = {"1": {"start_bp": 3631,"end_bp": 30425192},
                             "2": {"start_bp": 1871,"end_bp": 19696821},
                             "3": {"start_bp": 4342,"end_bp": 23459800},
                             "4": {"start_bp": 1180,"end_bp": 18584524},
                             "5": {"start_bp": 1251,"end_bp": 26970641}}
        
        # define KEY pmarknames
        chr_marker_list_dic = {}
        for i in chr_nr:
            chr_marker_list_dic[i] = getChrMarkers(i)  #django.db.models.query.ValuesListQuerySet    
        output_dic["pmarknames"] = chr_marker_list_dic
        
        # define KEY markers
        output_dic["markers"] = getMarkersList()
        
        # define KEY pmark
        marker_list_dic = {}
        for marker in output_dic["markers"]:
            m = Marker.objects.filter(marker_name = marker)[0]
            marker_list_dic[marker] = {"chr":m.marker_chromosome,"pos_cM":m.marker_cm,"pos_Mbp":m.marker_phys_pos}
        output_dic["pmark"] = marker_list_dic
        
        # define KEY gene
        gene_list_dic = {}
        gene_queryset_list = Gene.objects.all().order_by('chromosome','start')
        for gene in gene_queryset_list:
            gene_list_dic[gene.locus_identifier]={"chr":gene.chromosome,"pos_Mbp":(gene.start/1000000)} 
        output_dic["gene"] = gene_list_dic
        
        # define KEY peaks
        peaks_list = []
        lod_list = LOD.objects.filter(LOD_score__gte=lod_thld,gxe = False)
        for lod in lod_list:
            peaks_list.append({"gene":lod.locus_identifier.locus_identifier,"marker":lod.marker_name.marker_name,"lod":lod.LOD_score})
        output_dic["peaks"] = peaks_list
        toc = time.clock()
        print 'in %f seconds' % (toc-tic)
        return output_dic
    
    
if __name__ == '__main__':
    pass