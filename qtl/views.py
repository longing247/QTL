import time,os
import numpy as np
import csv
import json
import pandas as pd
import logging


from django.core.serializers.json import DjangoJSONEncoder
from django.shortcuts import render,get_object_or_404,render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.http import HttpRequest as request
from django.template import RequestContext
from django.views.generic import FormView
from django.conf import settings 
from django.core.context_processors import csrf
from django.contrib import messages
from django.db.models import Count, Avg

from .models import Experiment,Gene,Marker,LOD,Parent,RIL,Metabolite,MParent,MRIL,MLOD
from .forms import GeneUploadFileForm,MarkerUploadFileForm, LODUploadFileForm,ParentUploadFileForm,RILUploadFileForm,MetaboliteUploadFileForm,MParentUploadFileForm,MRILUploadFileForm,MLODUploadFileForm,ENVLODUploadFileForm,ENVMLODUploadFileForm,GeneUpdateFileForm

from Arabidopsis import Arabidopsis
from MySQLCorrelation import mysqlCorrelationAll,mysqlCorrelationSingle

'''
WUR Seed Lab eQTL expression dataset was used for test purpose. 
The derived files: gene_test.txt, marker_test.txt, lod_test.txt, and etc are all from the original datadataset.
All the files were normalized use tab delimited text file format with extension (.txt) and were converted from Unicode to UTF8 by Notepad in order to be compatible with MySQL.
Python built-in module csv was invoked to parse the Django upoloaded files.
'''

# Get an instance of a logger
logger = logging.getLogger(__name__)

# Log an error message (for example)
#       logger.error('Error message!') or logger.debug('Debug message!')

def uploadView(request):
    '''handle parsing objects from txt file, and saving entries to MySQL DB
    '''
    args = {}
    args.update(csrf(request))
    if request.method=='POST':
        if request.FILES.get('genefile'): # instead of request.FILES['genefile'], using get('genefile') to avoid multi dictionary key value exception
            form = GeneUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                geneUpload(request.FILES.get('genefile'))# request.FILES['genefile'] is from <input name='*' /> and need to be identical to form attribute(genefile).
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')
                return render_to_response('qtl/upload.html',args)
                
        elif request.FILES.get('markerfile'):
            form = MarkerUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                markerUpload(request.FILES.get('markerfile'))# request.FILES['genefile'] is from <input name='*' /> and need to be identical to form attribute(genefile).
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')                  
                return render_to_response('qtl/upload.html',args)         
        elif request.FILES.get('expfile'):
            form = LODUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                lodUpload(request.FILES.get('expfile'))# 
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')  
                return render_to_response('qtl/upload.html',args)
        elif request.FILES.get('parentfile'):
            form = ParentUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                parentUpload(request.FILES.get('parentfile'))
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error') 
                return render_to_response('qtl/upload.html',args)
        elif request.FILES.get('rilfile'):
            form = RILUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                rilUpload(request.FILES['rilfile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args)
        elif request.FILES.get('metabolitefile'):
            form = MetaboliteUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                metaboliteUpload(request.FILES['metabolitefile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args)
        elif request.FILES.get('MParentFile'):
            form = MParentUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                metaboliteParentUpload(request.FILES['MParentFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args)   
        elif request.FILES.get('MRILFile'):
            form = MRILUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                metaboliteRILUpload(request.FILES['MRILFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args)   
        elif request.FILES.get('MLODFile'):
            form = MLODUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                metaboliteLODUpload(request.FILES['MLODFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args)
        elif request.FILES.get('envLODFile'):
            form = ENVLODUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                envLODUpload(request.FILES['envLODFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args)           
        
        elif request.FILES.get('envMLODFile'):
            form = ENVMLODUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                envMLODUpload(request.FILES['envMLODFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args) 
        elif request.FILES.get('envMLODFile'):
            form = ENVMLODUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                envMLODUpload(request.FILES['envMLODFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args) 
        elif request.FILES.get('geneUpDateFile'):
            form = GeneUpdateFileForm(request.POST, request.FILES)
            if form.is_valid():
                geneUpdate(request.FILES['geneUpDateFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args) 
        
        
        
    else: # request.method =='GET'
        return render_to_response('qtl/upload.html',args)

def indexView(request):
    return render_to_response('qtl/index.html',{})

def chromosomeView(request):
    
    chr_dic = {1:'Chr I',2:'Chr II',3:'Chr III', 4:'Chr IV',5:'Chr V'}
    color_list = {1:'Black',2:'Red',3:'Green',4:'blue',5:'cyan',6:'purple'}
    features_list = {}
    for i in range(1,6):
        li = []
        marker_list = Marker.objects.filter(marker_chromosome = i)
        for markers in marker_list:
            t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[i+1]
            li.append(t)
        features_list[chr_dic[i]] = li    
    at = Arabidopsis()
    at.draw(features_list)
    return HttpResponseRedirect('success/')


def searchGeneView(request):
    '''
    fetch query gene/trait description
    fetch expression value in both parent and RIL population and present those value in bar chart respectively
    plot LOD scores of gene expression value along with chromosome (against marker)
    '''
    if request.GET.get('gene'):
        search_gene = request.GET.get('gene')
        gene = Gene.objects.get(locus_identifier = search_gene)
        exp_list = Parent.objects.filter(locus_identifier = search_gene)
        
        parent_type_list = []
        expression_list = []
        for exp in exp_list:
            parent_type_list.append(exp.parent_type)
            expression_list.append(exp.expression)
        
        ril_list = RIL.objects.filter(locus_identifier = search_gene).values('ril_type').annotate(average = Avg('ril_exp'))
        ril_type_list = []
        ril_avg_exp_list = []
        
        for ril in ril_list:
            ril_type_list.append(ril['ril_type'])
            ril_avg_exp_list.append(ril['average'])       
        ril_avg_list = [ril_avg_exp_list]
        
        #Entire correlation calculation 
        #query_result = mysqlCorrelationAll(search_gene)
        
        #Back-end Correlation calculation 
        #gene_list = []
        #corr_list = []
        #for gene in query_result:
        #    gene_list.append(gene.locus_identifier)
        #    corr_list.append(gene.LOD_score)
      
        #cor_list = ril_correlation(search_gene)
        #cor_list_js = json.dumps(cor_list)
        gxe_ = False
        peak_marker,peak_lod = findPeak(search_gene,gxe_)
        js_search_gene = json.dumps(search_gene)
        js_parent_type_list = json.dumps(parent_type_list)
        js_express_list = json.dumps([expression_list],cls=DjangoJSONEncoder) 
        js_ril_type_list = json.dumps(ril_type_list)
        js_ril_avg_exp_list = json.dumps(ril_avg_list,cls=DjangoJSONEncoder) 
        js_peak_marker= json.dumps(peak_marker)
        js_peak_lod = json.dumps(peak_lod,cls=DjangoJSONEncoder)
        
        #js_gene_list = json.dump(gene_list)
        #js_corr_list = json.dumps(corr_list,cls=DjangoJSONEncoder)        
        marker_list,lod_list = marker_plot(search_gene,gxe_)
        js_marker_list= json.dumps(marker_list)
        js_lod_list = json.dumps([lod_list],cls=DjangoJSONEncoder)
        
        gxe_env = True
        peak_marker_env,peak_lod_env = findPeak(search_gene,gxe_env)
        js_peak_marker_env= json.dumps(peak_marker_env)
        js_peak_lod_env = json.dumps(peak_lod_env,cls=DjangoJSONEncoder)
        marker_env_list,lod_env_list = marker_plot(search_gene,gxe_env)
        js_marker_env_list= json.dumps(marker_env_list)
        js_lod_env_list = json.dumps([lod_env_list],cls=DjangoJSONEncoder)
        js_lod_all_list = json.dumps([lod_list,lod_env_list],cls=DjangoJSONEncoder)
        
        
        #it might be helpful to define request.session.set_expiry(value). 
        request.session['search_gene'] = js_search_gene.encode('ascii','ignore')[1:-1]   
        request.session['peak_marker'] = js_peak_marker.encode('ascii','ignore')[1:-1]
        request.session['peak_lod'] = float(js_peak_lod[1:-1]) 
        return render_to_response('qtl/gene.html',{'search_gene':js_search_gene,#searched gene
                                                        'gene':gene,#Gene instance returned from database 
                                                        'peak_marker_js': js_peak_marker, # highest peak marker
                                                        'peak_lod_js':js_peak_lod,# the LOD score of the highest peak marker
                                                        'js_marker_list':js_marker_list,# markers along the chromosome 
                                                        'js_lod_list':js_lod_list, # the corresponding LOD expression value of the searched gene/trait against the marker list. 
                                                        #'cor_list_js':cor_list_js, #calculate correlation between a certain trait and all the other genes.
                                                        #'js_gene_list':js_gene_list,
                                                        #'js_corr_list':js_corr_list,
                                                        #'js_query_result':query_result,
                                                        'js_peak_marker_env':js_peak_marker_env, # environment interaction peak marker
                                                        'js_peak_lod_env':js_peak_lod_env,# environment interaction peak marker expression
                                                        'js_marker_env_list':js_marker_env_list,#markers along the chromosome 
                                                        'js_lod_env_list': js_lod_env_list,# he corresponding LOD expression value of the searched gene/trait with environmental interaction against the marker list.
                                                        'js_lod_all_list': js_lod_all_list, # a list contain two lod expression (LOD G and LOD GxE) list elements
                                                        'parent_type_list':js_parent_type_list, # parent type 
                                                        'express_list': js_express_list, # the corresponding expression value in parents of the searched gene
                                                        'js_ril_type_list':js_ril_type_list, # RIL type
                                                        'js_ril_avg_exp_list':js_ril_avg_exp_list})# the corresponding expression value in RILs of the searched gene            
    else:
        return render_to_response('qtl/gene.html',{})
    
    
def searchOverlapQTLView(request):
    '''
    search candidate overlap gene/expression traits(s) for a query gene/trait.
    '''
    gxe_=False
    if request.session['search_gene'] and request.session['peak_marker'] and request.session['peak_lod']:
        
        search_gene = request.session['search_gene']
        peak_marker = request.session['peak_marker']
        peak_lod = request.session['peak_lod'] 
        lod_thld = 2.3
        target_traits = ''
        marker_list = []
        lod_list = []
        overlap_traits_exp_list = []
        traits_list = []
        
        
        if request.GET.get('trait') and request.GET.get('lod_thld'):    
            corr = {}   
            if request.GET.get('lod_thld').strip().isdigit():  
                lod_thld = request.GET.get('lod_thld').strip()         
                target_traits = request.GET.get('trait').strip().split(',')
                overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld,locus_identifier__in=target_traits,marker_name_id = peak_marker,gxe=gxe_).exclude(locus_identifier = search_gene)
                for trait in overlap_traits:
                    marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_)
                    traits_list.append(trait.locus_identifier.locus_identifier)
                    
                    overlap_traits_exp_list.append(lod_list) 
                    for t in mysqlCorrelationSingle(trait.locus_identifier.locus_identifier,search_gene):
                        corr[trait.locus_identifier.locus_identifier.encode('ascii','ignore')] = t.r 
                # unicode needs to be converted again, avoiding of coercing to Unicode: need string or buffer, Gene found error
                js_search_gene = json.dumps(search_gene)
                js_peak_marker = json.dumps(peak_marker)
                js_peak_lod = json.dumps(peak_lod,cls=DjangoJSONEncoder) 
                js_target_traits = json.dumps(target_traits)
                js_lod_thld = json.dumps(lod_thld,cls=DjangoJSONEncoder)
                js_marker_list = json.dumps(marker_list)
                js_traits_list = json.dumps(traits_list)
                js_overlap_traits_exp_list = json.dumps(overlap_traits_exp_list,cls=DjangoJSONEncoder)  
                
                return render_to_response('qtl/overlap.html',{'search_gene':js_search_gene,
                                                              'peak_marker':js_peak_marker,# request.session
                                                              'peak_lod':js_peak_lod,# request.session
                                                              'compare_traits':js_target_traits,#NULL
                                                              'lod_thld':lod_thld, # default 2.3
                                                              'overlap_traits': overlap_traits,
                                                              'marker_list':js_marker_list,# marker
                                                              'traits_list':js_traits_list, # traits
                                                              'overlap_traits_exp_list':js_overlap_traits_exp_list, # expression 
                                                              'corr':corr 
                                                          })          
        elif request.GET.get('trait'):   
            corr = {}   
            target_traits = request.GET.get('trait').strip().split(',')
            overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld,locus_identifier__in=target_traits,marker_name = peak_marker,gxe=gxe_).exclude(locus_identifier = search_gene).order_by('-LOD_score')
            for trait in overlap_traits:
                marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_) 
                traits_list.append(trait.locus_identifier.locus_identifier)
                overlap_traits_exp_list.append(lod_list)
                for t in mysqlCorrelationSingle(trait.locus_identifier.locus_identifier,search_gene):
                    corr[trait.locus_identifier.locus_identifier.encode('ascii','ignore')] = t.r 
            js_search_gene = json.dumps(search_gene)
            js_peak_marker = json.dumps(peak_marker)
            js_peak_lod = json.dumps(peak_lod,cls=DjangoJSONEncoder) 
            js_target_traits = json.dumps(target_traits)
            js_lod_thld = json.dumps(lod_thld,cls=DjangoJSONEncoder)
            js_marker_list = json.dumps(marker_list)
            js_traits_list = json.dumps(traits_list)
            js_overlap_traits_exp_list = json.dumps(overlap_traits_exp_list,cls=DjangoJSONEncoder)     
            return render_to_response('qtl/overlap.html',{'search_gene':js_search_gene,
                                                              'peak_marker':js_peak_marker,# request.session
                                                              'peak_lod':js_peak_lod,# request.session
                                                              'compare_traits':js_target_traits,#NULL
                                                              'lod_thld':lod_thld, # default 2.3
                                                              'overlap_traits': overlap_traits,
                                                              'marker_list':js_marker_list,# marker
                                                              'traits_list':js_traits_list, # traits
                                                              'overlap_traits_exp_list':js_overlap_traits_exp_list, # expression
                                                              'corr':corr  
                                                          })  
        
        elif request.GET.get('lod_thld'):
            corr = {}
            if request.GET.get('lod_thld').strip().isdigit():  
                lod_thld = request.GET.get('lod_thld').strip()
            overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld, marker_name = peak_marker,gxe=gxe_).exclude(locus_identifier = search_gene).order_by('-LOD_score')[:5]
            #counter = 1 # default: present the expression value of the top 5 co-regulated traits along the chromosome in the same figure
            for trait in overlap_traits:
                marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_)
                traits_list.append(trait.locus_identifier.locus_identifier)
                overlap_traits_exp_list.append(lod_list)
                for t in mysqlCorrelationSingle(trait.locus_identifier.locus_identifier,search_gene):
                    corr[trait.locus_identifier.locus_identifier.encode('ascii','ignore')] = t.r 
            js_search_gene = json.dumps(search_gene)
            js_peak_marker = json.dumps(peak_marker)
            js_peak_lod = json.dumps(peak_lod,cls=DjangoJSONEncoder) 
            js_target_traits = json.dumps(target_traits)
            js_lod_thld = json.dumps(lod_thld,cls=DjangoJSONEncoder)
            js_marker_list = json.dumps(marker_list)
            js_traits_list = json.dumps(traits_list)
            js_overlap_traits_exp_list = json.dumps(overlap_traits_exp_list,cls=DjangoJSONEncoder)   
            return render_to_response('qtl/overlap.html',{'search_gene':js_search_gene,
                                                          'peak_marker':js_peak_marker,# request.session
                                                          'peak_lod':js_peak_lod,# request.session
                                                          'compare_traits':js_target_traits,#NULL
                                                          'lod_thld':lod_thld, # default 2.3
                                                          'overlap_traits': overlap_traits,
                                                          'marker_list':js_marker_list,# marker
                                                          'traits_list':js_traits_list, # traits
                                                          'overlap_traits_exp_list':js_overlap_traits_exp_list, # expression 
                                                          'corr':corr 
                                                          }) 
        else:
            # __gte option needs to be placed in the first argument. Attention
            # tested in manage.py shell
            # traits = LOD.objects.filter(LOD_score__gte = 2.3, marker_name = 'MSAT318406')      
            # print overlap+traits gave error: coercing to Unicode: need string or buffer, Gene found 
            # exclude the query gene
            #.exclude(locus_identifier = search_gene)
            overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld, marker_name = peak_marker,gxe=gxe_).exclude(locus_identifier = search_gene).order_by('-LOD_score')[:5]#.annotate(correlation = mysqlCorrelationSingle(LOD.locus_identifier,search_gene))
            #counter = 1 # default: present the expression value of the top 5 co-regulated traits along the chromosome in the same figure
            corr = {}
            for trait in overlap_traits:

                marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_)
                traits_list.append(trait.locus_identifier.locus_identifier)
                overlap_traits_exp_list.append(lod_list)
                for t in mysqlCorrelationSingle(trait.locus_identifier.locus_identifier,search_gene):
                    corr[trait.locus_identifier.locus_identifier.encode('ascii','ignore')] = t.r 
            js_search_gene = json.dumps(search_gene)
            js_peak_marker = json.dumps(peak_marker)
            js_peak_lod = json.dumps(peak_lod,cls=DjangoJSONEncoder) 
            js_target_traits = json.dumps(target_traits)
            js_lod_thld = json.dumps(lod_thld,cls=DjangoJSONEncoder)
            js_marker_list = json.dumps(marker_list)
            js_traits_list = json.dumps(traits_list)
            js_overlap_traits_exp_list = json.dumps(overlap_traits_exp_list,cls=DjangoJSONEncoder)    
            return render_to_response('qtl/overlap.html',{'search_gene':js_search_gene,
                                                          'peak_marker':js_peak_marker,# request.session
                                                          'peak_lod':js_peak_lod,# request.session
                                                          'compare_traits':js_target_traits,#NULL
                                                          'lod_thld':lod_thld, # default 2.3
                                                          'overlap_traits': overlap_traits,
                                                          'marker_list':js_marker_list,# marker
                                                          'traits_list':js_traits_list, # traits
                                                          'overlap_traits_exp_list':js_overlap_traits_exp_list, # expression
                                                          'corr':corr 
                                                          })  
    
    else:
        return render_to_response('qtl/gene.html',{})

def searchMetaboliteView(request):
    '''
    fetch query metabolite information
    fetch expression value in both parent and RIL population and present those value in bar chart respectively
    plot LOD scores of gene expression value along with chromosome (against marker)
    '''
    if request.GET.get('metabolite'):
        search_metabolite = request.GET.get('metabolite')
        metabolite = Metabolite.objects.get(metabolite_name = search_metabolite)
        exp_list = MParent.objects.filter(metabolite_name = search_metabolite)
        
        parent_type_list = []
        expression_list = []
        for exp in exp_list:
            parent_type_list.append(exp.parent_type)
            expression_list.append(exp.expression)
        
        ril_list = MRIL.objects.filter(metabolite_name = search_metabolite).values('ril_type').annotate(average = Avg('ril_exp'))
        ril_type_list = []
        ril_avg_exp_list = []
        
        for ril in ril_list:
            ril_type_list.append(ril['ril_type'])
            ril_avg_exp_list.append(ril['average'])       
        ril_avg_list = [ril_avg_exp_list]
        
        #Entire correlation calculation 
        #query_result = mysqlCorrelationAll(search_gene)
        
        #Back-end Correlation calculation 
        #gene_list = []
        #corr_list = []
        #for gene in query_result:
        #    gene_list.append(gene.locus_identifier)
        #    corr_list.append(gene.LOD_score)
      
        #cor_list = ril_correlation(search_gene)
        #cor_list_js = json.dumps(cor_list)
        gxe_ = False
        peak_marker,peak_lod = findMPeak(search_metabolite,gxe_)
        js_search_metabolite = json.dumps(search_metabolite)
        js_parent_type_list = json.dumps(parent_type_list)
        js_express_list = json.dumps([expression_list],cls=DjangoJSONEncoder) 
        js_ril_type_list = json.dumps(ril_type_list)
        js_ril_avg_exp_list = json.dumps(ril_avg_list,cls=DjangoJSONEncoder) 
        js_peak_marker= json.dumps(peak_marker)
        js_peak_lod = json.dumps(peak_lod,cls=DjangoJSONEncoder)
        
        #js_gene_list = json.dump(gene_list)
        #js_corr_list = json.dumps(corr_list,cls=DjangoJSONEncoder)        
        marker_list,lod_list = m_marker_plot(search_metabolite,gxe_)
        js_marker_list= json.dumps(marker_list)
        js_lod_list = json.dumps([lod_list],cls=DjangoJSONEncoder)
        
        gxe_env = True
        peak_marker_env,peak_lod_env = findMPeak(search_metabolite,gxe_env)
        js_peak_marker_env= json.dumps(peak_marker_env)
        js_peak_lod_env = json.dumps(peak_lod_env,cls=DjangoJSONEncoder)
        marker_env_list,lod_env_list = m_marker_plot(search_metabolite,gxe_env)
        js_marker_env_list= json.dumps(marker_env_list)
        js_lod_env_list = json.dumps([lod_env_list],cls=DjangoJSONEncoder)
        js_lod_all_list = json.dumps([lod_list,lod_env_list],cls=DjangoJSONEncoder)
        
        
        #it might be helpful to define request.session.set_expiry(value). 
        request.session['search_metabolite'] = js_search_metabolite.encode('ascii','ignore')[1:-1]   
        request.session['peak_marker'] = js_peak_marker.encode('ascii','ignore')[1:-1]
        request.session['peak_lod'] = float(js_peak_lod[1:-1]) 
        return render_to_response('qtl/metabolite.html',{'search_metabolite':js_search_metabolite,#searched metabolite
                                                        'metabolite':metabolite,#Metabolite instance returned from database 
                                                        'peak_marker_js': js_peak_marker, # highest peak marker
                                                        'peak_lod_js':js_peak_lod,# the LOD score of the highest peak marker
                                                        'js_marker_list':js_marker_list,# markers along the chromosome 
                                                        'js_lod_list':js_lod_list, # the corresponding LOD expression value of the searched gene/trait against the marker list. 
                                                        #'cor_list_js':cor_list_js, #calculate correlation between a certain trait and all the other genes.
                                                        #'js_gene_list':js_gene_list,
                                                        #'js_corr_list':js_corr_list,
                                                        #'js_query_result':query_result,
                                                        'js_peak_marker_env':js_peak_marker_env, # environment interaction peak marker
                                                        'js_peak_lod_env':js_peak_lod_env,# environment interaction peak marker expression
                                                        'js_marker_env_list':js_marker_env_list,#markers along the chromosome 
                                                        'js_lod_env_list': js_lod_env_list,# he corresponding LOD expression value of the searched gene/trait with environmental interaction against the marker list.
                                                        'js_lod_all_list': js_lod_all_list, # a list contain two lod expression (LOD G and LOD GxE) list elements
                                                        'parent_type_list':js_parent_type_list, # parent type 
                                                        'express_list': js_express_list, # the corresponding expression value in parents of the searched gene
                                                        'js_ril_type_list':js_ril_type_list, # RIL type
                                                        'js_ril_avg_exp_list':js_ril_avg_exp_list})# the corresponding expression value in RILs of the searched gene            
    else:
        return render_to_response('qtl/metabolite.html',{})
def searchMarkerView(request):
    '''
    Query for a marker 
    '''
    if request.GET.get('marker'):
        query_marker_name = request.GET.get('marker')
        marker = Marker.objects.get(marker_name = query_marker_name)
        context_dict = {'marker':marker}
        
        return render_to_response('qtl/marker.html',context_dict)
    else:
        return render_to_response('qtl/marker.html',{})
    

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
            add_gene.locus_identifier = line[0]
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
    add_exp = Experiment()
    add_exp.experiment_name = f.name
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
                add_lod = LOD(experiment_name = Experiment.objects.get(pk=f.name),
                              locus_identifier = Gene.objects.get(pk= line[0]),
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
                                        locus_identifier = Gene.objects.get(pk= line[0]))
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
                                  locus_identifier = Gene.objects.get(pk= line[0]))
                    add_ril.save()
                else:
                    continue
            print line[0]
            print '%d records added' % i
            
def getAllMetabolite(f,arg1,arg2):
    li1 = []
    li2 = []
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
                                  metabolite_name = Metabolite.objects.get(pk= line[0]))
                    add_ril.save()
                else:
                    continue
            print line[0]
            print '%d records added' % i
    
def getData(f):
    csv_reader = csv.reader(f,delimiter='\t')
    for line in csv_reader:
        yield line
        

def upload_success(request):
    return render_to_response('qtl/success.html')

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
            add_metabolite.metabolite_name = line[0]
            

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
                                        metabolite_name = Metabolite.objects.get(pk= line[0]))
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
    add_exp.experiment_name = f.name
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
                add_mlod = MLOD(experiment_name = Experiment.objects.get(pk=f.name),
                              metabolite_name = Metabolite.objects.get(pk= line[0]),
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
    add_exp.experiment_name = f.name
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
                add_lod = LOD(experiment_name = Experiment.objects.get(pk=f.name),
                              locus_identifier = Gene.objects.get(pk= line[0]),
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
    add_exp = Experiment()
    add_exp.experiment_name = f.name
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
                add_lod = MLOD(experiment_name = Experiment.objects.get(pk=f.name),
                              metabolite_name = Metabolite.objects.get(pk= line[0]),
                              marker_name =  Marker.objects.get(pk=marker_list[meta-1]),
                              gxe = True,
                              LOD_score =  line[meta])
                add_lod.save()    
    toc = time.clock()            
    print 'in %f seconds' % (toc-tic)
    print '%d records added' %i    


def findPeak(gene_name,gexpression):
    '''find the peak marker and expression value
    argument2 gexpression refers LOD.gxe(BooleanField: Ture or False in Django, 1 or 0 in MySQL) enviorment interaction. 
    '''
    exp_list = LOD.objects.filter(locus_identifier = gene_name,gxe= gexpression)    
    if exp_list:
        peak = exp_list[0].LOD_score
        for exp in exp_list:
            if exp.LOD_score > peak:
                peak = exp.LOD_score

        marker = LOD.objects.filter(locus_identifier = gene_name,LOD_score = peak).values('marker_name')[0]['marker_name']   
    return marker,peak

def findMPeak(metabolite,gexpression):
    '''find the peak marker and expression value of metabolite
    argument2 gexpression refers MLOD.gxe(BooleanField: Ture or False in Django, 1 or 0 in MySQL) enviorment interaction. 
    '''
    exp_list = MLOD.objects.filter(metabolite_name = metabolite,gxe= gexpression)    
    if exp_list:
        peak = exp_list[0].LOD_score
        for exp in exp_list:
            if exp.LOD_score > peak:
                peak = exp.LOD_score

        marker = MLOD.objects.filter(metabolite_name = metabolite,LOD_score = peak).values('marker_name')[0]['marker_name']   
    return marker,peak

def yield_gene(gene_list):
    for gene in gene_list:
        yield gene
        
def exp_series(query_gene):
    '''accumulate expression value per gene or qtl 
    and cast to pandas Series object which will be later used in correlation calculation
    '''
    series_list = [] # expression list per gene or qtl
    
    #select ril_exp from RIL where locus_identifier = query_gene;
    exp_decimal = RIL.objects.filter(locus_identifier = query_gene).values_list('ril_exp', flat = True)# return Decimal object
    for exp in exp_decimal:
        exp_ = float('{0:.2f}'.format(exp)) # type 'float'      
        series_list.append(exp_) 
    
    return pd.Series(series_list) #pandas Series object which will be later used in correlation calculation

def marker_plot(gene_name,gxe_):
    '''
    plot expression (LOD scores) along the chromosome, gxe True: response to enviormental changes influencing by genotypes. 
    '''
    lod_list = []
    marker_list = []
    for marker in Marker.objects.values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm'):
        marker_list.append(marker)
        exp = LOD.objects.get(locus_identifier = gene_name,marker_name = marker,gxe=gxe_).LOD_score
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

def m_marker_plot(metabolite,gxe_):
    '''
    plot expression (LOD scores) along the chromosome, gxe True: response to enviormental changes influencing by genotypes. 
    '''
    lod_list = []
    marker_list = []
    for marker in Marker.objects.values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm'):
        marker_list.append(marker)
        exp = MLOD.objects.get(metabolite_name = metabolite,marker_name = marker,gxe=gxe_).LOD_score
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


def ril_correlation(query_gene):
    '''calculate Pearson correlation value of a query gene/qtl against the rest 29000 genes/qtls.
    It took about 8 mins to finish the entire proecess. 
    '''
    tic = time.time()
    gene_list = Gene.objects.values_list('locus_identifier',flat=True) #to retrieve all the genes/qtls
    #ril_list = RIL.objects.values_list('ril_name',flat = True) #to retrieve the names of RIL offsprings ***too redundant
    cor_dic ={} 
    pd_query_gene_exp_series = exp_series(query_gene)  
    i=0
    for gene in yield_gene(gene_list):
        i+=1
        #pd_target_gene_exp_series = exp_series(ril_list,gene.encode('ascii','ignore'))
        pd_target_gene_exp_series = exp_series(gene)
        cor_dic[gene.encode('ascii','ignore')] = pd_query_gene_exp_series.corr(pd_target_gene_exp_series) # calculate correlation of two series
        print i,cor_dic[gene.encode('ascii','ignore')]
    toc = time.time()
    print 'in %f seconds' % (toc-tic)
    return cor_dic

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
            
######################################
# main #
######################################
#if __name__ == "__main__":
    

