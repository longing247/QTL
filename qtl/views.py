import time,os
import numpy as np
import csv
import json
import pandas as pd
import logging
import math
import re
import itertools
from datetime import datetime

from django.core.serializers.json import DjangoJSONEncoder
from django.shortcuts import render,get_object_or_404,render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.http import HttpRequest as request
from django.template import RequestContext
from django.views.generic import FormView
from django.conf import settings 
from django.core.context_processors import csrf
from django.contrib import messages
from django.db.models import Count, Avg, Q
from django.core.serializers import serialize

from .models import Experiment,Gene,Marker,LOD,Parent,RIL,Metabolite,MParent,MRIL,MLOD,Genotype,ExperimentMarker,ExperimentGene,ExperimentRIL,ExperimentParent,ExperimentMetabolite
from .forms import GeneUploadFileForm,MarkerUploadFileForm, LODUploadFileForm,ParentUploadFileForm,RILUploadFileForm,MetaboliteUploadFileForm,MParentUploadFileForm,MRILUploadFileForm,MLODUploadFileForm,ENVLODUploadFileForm,ENVMLODUploadFileForm,GeneUpdateFileForm,GenotypeUploadFileForm

from Arabidopsis import Arabidopsis
from MySQLCorrelation import mysqlCorrelationAll,mysqlCorrelationSingle
#from OutputJson import outputJson
'''
WUR Seed Lab eQTL expression dataset was used for test purpose. 
The derived files: gene_test.txt, marker_test.txt, lod_test.txt, and etc are all from the original datadataset.
All the files were normalized use tab delimited text file format with extension (.txt) and were converted from Unicode to UTF8 by Notepad in order to be compatible with MySQL.
Python built-in module csv was invoked to parse the Django upoloaded files.
'''

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
        elif request.FILES.get('genotypeFile'):
            form = GenotypeUploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                genotypeUpload(request.FILES['genotypeFile'])
                return HttpResponseRedirect('success/')
            else:
                print 'form is not valid'
                messages.error(request,'Error')     
                return render_to_response('qtl/upload.html',args) 
        
        else: # No file was selected.
            return render_to_response('qtl/upload.html',args)
    else: # request.method =='GET'
        return render_to_response('qtl/upload.html',args)

def indexView(request):
    return render_to_response('qtl/index.html',{})

def chromosomeView(request):
    search_gene = request.session['search_gene']
    lod_list = request.session['js_lod_list_session'][1:-1].split(" ")
    experiment = request.session['experiment_name'].encode('ascii','ignore')    
    chr_dic = {1:'Chr I',2:'Chr II',3:'Chr III', 4:'Chr IV',5:'Chr V'}
    color_list = {1:'black',2:'blue',3:'Green',4:'purple',5:'cyan',6:'red'}
    features_list = {}
    j = 0

    marker_list_experiment = ExperimentMarker.objects.filter(experiment_name = experiment).values_list('marker_name',flat=True)
    for i in range(1,6):
        li = []        
        marker_list = Marker.objects.filter(marker_chromosome = i,marker_name__in=marker_list_experiment).order_by('marker_phys_pos')
        for markers in marker_list:
            print j,len(marker_list),len(lod_list)
            if float(lod_list[j]) <-3:
                t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[1]
            elif float(lod_list[j]) >=-3 and float(lod_list[j]) <0:
                t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[2]
            elif float(lod_list[j]) >= 0 and float(lod_list[j]) <3:
                t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[3]
            elif float(lod_list[j]) >=3 and float(lod_list[j]) < 5:
                t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[4]
            elif float(lod_list[j]) >=5 and float(lod_list[j]) <10:
                t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[5]
            else:
                t = int(markers.marker_phys_pos*1000000), (int(markers.marker_phys_pos*1000000)+100),None,markers.marker_name.encode('ascii','ignore'),color_list[6]
            
            li.append(t)
            j+=1
            print markers.marker_name
        features_list[chr_dic[i]] = li
        
    now = datetime.now()
    f_name = str(now).replace(' ','').replace(':','').replace('-','').replace('.','')+'.pdf' 
    at = Arabidopsis(f_name)
    at.draw(features_list)
    pdf_name = 'qtl/temp/'+at.f_name
    
    #args: 
    return render_to_response('qtl/chromosome.html',{'pdf_name':pdf_name
                                                     })


def searchGeneView(request):
    '''
    fetch query gene/trait description
    fetch expression value in both parent and RIL population and present those value in bar chart respectively
    plot LOD scores of gene expression value along with chromosome (against marker)
    '''
    if request.GET.get('gene') and request.GET.get('experiment'):     
        search_gene = request.GET.get('gene').encode('ascii','ignore').strip().upper()
       
        experiment = request.GET.get('experiment').encode('ascii','ignore').strip() #returns experiment name like Ligterink_2014 which are used as filter in the query
        
        
        if Gene.objects.filter(locus_identifier = search_gene).exists():
            render_dic = {}
            gene = Gene.objects.get(locus_identifier = search_gene)
            is_parent_exp = False # expression profile of both parent
            is_ril_exp = False # expression profile of RIL lines
            is_lod_exp = False # correlation profile of gene expression and genomic variance 
            is_gxp_exp = False # Gene expression phenotype
            is_gxe_exp = False # Environmental effects on gene expression phenotype
            
            
            js_search_gene = json.dumps(search_gene)
            render_dic['search_gene'] = js_search_gene 
            request.session['search_gene'] = js_search_gene.encode('ascii','ignore')[1:-1]   
            request.session['experiment_name'] = experiment 

            render_dic['gene'] = gene#Gene instance returned from database 
            # check whether parent expression data is asvailable.
            if Parent.objects.filter(locus_identifier = search_gene, experiment_name = experiment).exists():
                is_parent_exp = True
                render_dic['is_parent_exp'] = is_parent_exp
                
            # check whether RIL lines expression data is available.
            if RIL.objects.filter(locus_identifier = search_gene,experiment_name = experiment).exists():
                is_ril_exp = True
                render_dic['is_ril_exp'] = is_ril_exp        
            
            # check whether LOD lines expression data is available.
            if LOD.objects.filter(locus_identifier = search_gene,experiment_name = experiment).exists():
                is_lod_exp = True
                render_dic['is_lod_exp'] = is_lod_exp
            
            if LOD.objects.filter(locus_identifier = search_gene,experiment_name = experiment,gxe=False).exists():
                is_gxp_exp = True 
                render_dic['is_gxp_exp'] = is_gxp_exp
            
            if LOD.objects.filter(locus_identifier = search_gene,experiment_name = experiment,gxe=True).exists():
                is_gxe_exp = True
                render_dic['is_gxe_exp'] = is_gxe_exp
            
            if is_parent_exp:         
                exp_list = Parent.objects.filter(locus_identifier = search_gene, experiment_name = experiment)
                parent_type_list = []
                parent_expression_list = []
                for exp in exp_list:
                    parent_type_list.append(exp.parent_type)
                    parent_expression_list.append(decimalFormat(exp.expression))
                js_parent_type_list = json.dumps(parent_type_list)
                parent_exp_list = [parent_expression_list]
                js_parent_express_list = json.dumps(parent_exp_list,cls=DjangoJSONEncoder) 
                render_dic['parent_type_list']=js_parent_type_list # parent type 
                render_dic['parent_exp_list'] =js_parent_express_list # the corresponding expression value in parents of the searched gene
            if is_ril_exp:
                
                ril_list = RIL.objects.filter(locus_identifier = search_gene,experiment_name = experiment).values('ril_type').annotate(average = Avg('ril_exp'))
                ril_type_list = []
                ril_avg_exp_list = []
                
                for ril in ril_list:
                    ril_type_list.append(ril['ril_type'])
                    ril_avg_exp_list.append(decimalFormat(ril['average']))       
                ril_avg_list = [ril_avg_exp_list]
                
                js_ril_type_list = json.dumps(ril_type_list)
                js_ril_avg_exp_list = json.dumps(ril_avg_list,cls=DjangoJSONEncoder) 
                render_dic['js_ril_type_list']=js_ril_type_list # RIL type
                render_dic['js_ril_avg_exp_list']=js_ril_avg_exp_list# the corresponding expression value in RILs of the searched gene            

            if is_lod_exp:
                lod_list = []
                lod_env_list = []
                marker_list =  []
                if is_gxp_exp:
                    gxe_=False
                    peak_marker,peak_lod = findPeak(search_gene,gxe_,experiment)
                    js_peak_marker= json.dumps(peak_marker)
                    js_peak_lod = json.dumps(decimalFormat(peak_lod),cls=DjangoJSONEncoder)
                    marker_list,lod_list = marker_plot(search_gene,gxe_,experiment)
                    js_marker_list= json.dumps(marker_list)
                    js_lod_list = json.dumps([lod_list],cls=DjangoJSONEncoder)
                    render_dic['peak_marker_js'] =js_peak_marker # highest peak marker
                    render_dic['peak_lod_js']=js_peak_lod# the LOD score of the highest peak marker
                    render_dic['js_marker_list']=js_marker_list# markers along the chromosome 
                    render_dic['js_lod_list']=js_lod_list # the corresponding LOD expression value of the searched gene/trait against the marker list. 
                    request.session['peak_marker'] = js_peak_marker.encode('ascii','ignore')[1:-1]
                    request.session['peak_lod'] = js_peak_lod
                    
                if is_gxe_exp:
                    gxe_env = True
                    peak_marker_env,peak_lod_env = findPeak(search_gene,gxe_env,experiment)
                    js_peak_marker_env= json.dumps(peak_marker_env)
                    js_peak_lod_env = json.dumps(peak_lod_env,cls=DjangoJSONEncoder)
                    marker_env_list,lod_env_list = marker_plot(search_gene,gxe_env,experiment)
                    js_marker_env_list= json.dumps(marker_env_list)
                    js_lod_env_list = json.dumps([lod_env_list],cls=DjangoJSONEncoder)
                    js_lod_all_list = json.dumps([lod_list,lod_env_list],cls=DjangoJSONEncoder)
                         
                    render_dic['js_peak_marker_env']=js_peak_marker_env # environment interaction peak marker
                    render_dic['js_peak_lod_env']=js_peak_lod_env# environment interaction peak marker expression
                    render_dic['js_marker_env_list']=js_marker_env_list#markers along the chromosome 
                    render_dic['js_lod_env_list']=js_lod_env_list# he corresponding LOD expression value of the searched gene/trait with environmental interaction against the marker list.
                    
                    request.session['peak_marker_env'] = js_peak_marker_env.encode('ascii','ignore')[1:-1]
                    request.session['peak_lod_env'] = js_peak_lod_env
                
                if is_gxp_exp and is_gxe_exp:
                    js_marker_all_list= json.dumps(marker_list)
                    js_lod_all_list = json.dumps([lod_list,lod_env_list],cls=DjangoJSONEncoder)
                    render_dic['js_marker_all_list']=js_marker_all_list#markers along the chromosome 
                    render_dic['js_lod_all_list']=js_lod_all_list# a list contain two lod expression (LOD G and LOD GxE) list elements
                    
                       
            
            lod_str =  ''
            for lod in lod_list:
                lod_str +=' '+str(lod)
            request.session['js_lod_list_session'] = json.dumps(lod_str[1:],cls=DjangoJSONEncoder)
            
            return render_to_response('qtl/gene.html',render_dic)
                                                            
                                                                                        
        else:
            raise NameError('Query gene does not exist')
        
    else:
        exps = Experiment.objects.all().values_list('experiment_name',flat=True)
        return render_to_response('qtl/gene.html',{'exps':exps})
    
    
def searchOverlapQTLView(request):
    '''
    search candidate overlap gene/expression traits(s) for a query gene/trait.
    '''
    gxe_=False
    if request.session['search_gene'] and request.session['peak_marker'] and request.session['peak_lod']:
        
        search_gene = request.session['search_gene']
        peak_marker = request.session['peak_marker']
        peak_lod = request.session['peak_lod'] 
        experiment = request.session['experiment_name']
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
                overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld,locus_identifier__in=target_traits,marker_name_id = peak_marker,gxe=gxe_,experiment_name=experiment).exclude(locus_identifier = search_gene)
                for trait in overlap_traits:
                    marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_,experiment)
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
            overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld,locus_identifier__in=target_traits,marker_name = peak_marker,gxe=gxe_,experiment_name=experiment).exclude(locus_identifier = search_gene).order_by('-LOD_score')
            for trait in overlap_traits:
                marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_,experiment) 
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
            overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld, marker_name = peak_marker,gxe=gxe_,experiment_name=experiment).exclude(locus_identifier = search_gene).order_by('-LOD_score')[:5]
            #counter = 1 # default: present the expression value of the top 5 co-regulated traits along the chromosome in the same figure
            for trait in overlap_traits:
                marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_,experiment)
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
            overlap_traits = LOD.objects.filter(LOD_score__gte = lod_thld, marker_name = peak_marker,gxe=gxe_,experiment_name=experiment).exclude(locus_identifier = search_gene).order_by('-LOD_score')[:5]#.annotate(correlation = mysqlCorrelationSingle(LOD.locus_identifier,search_gene))
            #counter = 1 # default: present the expression value of the top 5 co-regulated traits along the chromosome in the same figure
            corr = {}
            for trait in overlap_traits:

                marker_list,lod_list = marker_plot(trait.locus_identifier,gxe_,experiment)
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
    if request.GET.get('metabolite') and request.GET.get('experiment'):     
        search_metabolite = request.GET.get('metabolite').encode('ascii','ignore').strip().upper()
       
        experiment = request.GET.get('experiment').encode('ascii','ignore').strip() #returns experiment name like Ligterink_2014 which are used as filter in the query
        
        
        if Metabolite.objects.filter(metabolite_name = search_metabolite).exists():
            render_dic = {}
            metabolite = Metabolite.objects.get(metabolite_name = search_metabolite)
            is_parent_exp = False # expression profile of both parent
            is_ril_exp = False # expression profile of RIL lines
            is_lod_exp = False # correlation profile of gene expression and genomic variance 
            is_gxp_exp = False # Gene expression phenotype
            is_gxe_exp = False # Environmental effects on gene expression phenotype
            
            
            js_search_metabolite = json.dumps(search_metabolite)
            render_dic['search_metabolite'] = js_search_metabolite 
            request.session['search_metabolite'] = js_search_metabolite.encode('ascii','ignore')[1:-1]   
            request.session['experiment_name'] = experiment 

            render_dic['metabolite'] = metabolite#Gene instance returned from database 
            # check whether parent expression data is asvailable.
            if MParent.objects.filter(metabolite_name = search_metabolite, experiment_name = experiment).exists():
                is_parent_exp = True
                render_dic['is_parent_exp'] = is_parent_exp
                
            # check whether RIL lines expression data is available.
            if MRIL.objects.filter(metabolite_name = search_metabolite,experiment_name = experiment).exists():
                is_ril_exp = True
                render_dic['is_ril_exp'] = is_ril_exp        
            
            # check whether LOD lines expression data is available.
            if MLOD.objects.filter(metabolite_name = search_metabolite,experiment_name = experiment).exists():
                is_lod_exp = True
                render_dic['is_lod_exp'] = is_lod_exp
            
            if MLOD.objects.filter(metabolite_name = search_metabolite,experiment_name = experiment,gxe=False).exists():
                is_gxp_exp = True 
                render_dic['is_gxp_exp'] = is_gxp_exp
            
            if MLOD.objects.filter(metabolite_name = search_metabolite,experiment_name = experiment,gxe=True).exists():
                is_gxe_exp = True
                render_dic['is_gxe_exp'] = is_gxe_exp
            
            if is_parent_exp:         
                exp_list = MParent.objects.filter(metabolite_name = search_metabolite, experiment_name = experiment)
                parent_type_list = []
                parent_expression_list = []
                for exp in exp_list:
                    parent_type_list.append(exp.parent_type)
                    parent_expression_list.append(decimalFormat(exp.expression))
                js_parent_type_list = json.dumps(parent_type_list)
                parent_exp_list = [parent_expression_list]
                js_parent_express_list = json.dumps(parent_exp_list,cls=DjangoJSONEncoder) 
                render_dic['parent_type_list']=js_parent_type_list # parent type 
                render_dic['parent_exp_list'] =js_parent_express_list # the corresponding expression value in parents of the searched gene
            if is_ril_exp:
                
                ril_list = MRIL.objects.filter(metabolite_name = search_metabolite,experiment_name = experiment).values('ril_type').annotate(average = Avg('ril_exp'))
                ril_type_list = []
                ril_avg_exp_list = []
                
                for ril in ril_list:
                    ril_type_list.append(ril['ril_type'])
                    ril_avg_exp_list.append(decimalFormat(ril['average']))       
                ril_avg_list = [ril_avg_exp_list]
                
                js_ril_type_list = json.dumps(ril_type_list)
                js_ril_avg_exp_list = json.dumps(ril_avg_list,cls=DjangoJSONEncoder) 
                render_dic['js_ril_type_list']=js_ril_type_list # RIL type
                render_dic['js_ril_avg_exp_list']=js_ril_avg_exp_list# the corresponding expression value in RILs of the searched gene            

            if is_lod_exp:
                lod_list = []
                lod_env_list = []
                marker_list =  []
                if is_gxp_exp:
                    gxe_=False
                    peak_marker,peak_lod = findMPeak(search_metabolite,gxe_,experiment)
                    js_peak_marker= json.dumps(peak_marker)
                    js_peak_lod = json.dumps(decimalFormat(peak_lod),cls=DjangoJSONEncoder)
                    marker_list,lod_list = m_marker_plot(search_metabolite,gxe_,experiment)
                    js_marker_list= json.dumps(marker_list)
                    js_lod_list = json.dumps([lod_list],cls=DjangoJSONEncoder)
                    render_dic['peak_marker_js'] =js_peak_marker # highest peak marker
                    render_dic['peak_lod_js']=js_peak_lod# the LOD score of the highest peak marker
                    render_dic['js_marker_list']=js_marker_list# markers along the chromosome 
                    render_dic['js_lod_list']=js_lod_list # the corresponding LOD expression value of the searched gene/trait against the marker list. 
                    request.session['peak_marker'] = js_peak_marker.encode('ascii','ignore')[1:-1]
                    request.session['peak_lod'] = js_peak_lod
                    
                if is_gxe_exp:
                    gxe_env = True
                    peak_marker_env,peak_lod_env = findMPeak(search_metabolite,gxe_env,experiment)
                    js_peak_marker_env= json.dumps(peak_marker_env)
                    js_peak_lod_env = json.dumps(peak_lod_env,cls=DjangoJSONEncoder)
                    marker_env_list,lod_env_list = m_marker_plot(search_metabolite,gxe_env,experiment)
                    js_marker_env_list= json.dumps(marker_env_list)
                    js_lod_env_list = json.dumps([lod_env_list],cls=DjangoJSONEncoder)
                    js_lod_all_list = json.dumps([lod_list,lod_env_list],cls=DjangoJSONEncoder)
                         
                    render_dic['js_peak_marker_env']=js_peak_marker_env # environment interaction peak marker
                    render_dic['js_peak_lod_env']=js_peak_lod_env# environment interaction peak marker expression
                    render_dic['js_marker_env_list']=js_marker_env_list#markers along the chromosome 
                    render_dic['js_lod_env_list']=js_lod_env_list# he corresponding LOD expression value of the searched gene/trait with environmental interaction against the marker list.
                    
                    request.session['peak_marker_env'] = js_peak_marker_env.encode('ascii','ignore')[1:-1]
                    request.session['peak_lod_env'] = js_peak_lod_env
                
                if is_gxp_exp and is_gxe_exp:
                    js_marker_all_list= json.dumps(marker_list)
                    js_lod_all_list = json.dumps([lod_list,lod_env_list],cls=DjangoJSONEncoder)
                    render_dic['js_marker_all_list']=js_marker_all_list#markers along the chromosome 
                    render_dic['js_lod_all_list']=js_lod_all_list# a list contain two lod expression (LOD G and LOD GxE) list elements
                    
            return render_to_response('qtl/metabolite.html',render_dic)
        else:
            raise NameError('Query gene does not exist')
        
    else:
        exps = Experiment.objects.all().values_list('experiment_name',flat=True)
        return render_to_response('qtl/metabolite.html',{'exps':exps})
def searchMarkerView(request):
    '''
    Query for a marker 
    '''
    if request.GET.get('marker'):
        query_marker_name = request.GET.get('marker').strip()
        marker = Marker.objects.get(marker_name = query_marker_name)
        context_dict = {'marker':marker}
        
        return render_to_response('qtl/marker.html',context_dict)
    else:
        return render_to_response('qtl/marker.html',{})

def eQTLPlotView(request):
    '''
    plot eQTL maping
    
    '''
    if request.session['search_gene'] and request.session['experiment_name']:
        tic = time.time()
        lod_thld = 2.3
        if request.GET.get('eQTL_lod_thld'):
            lod_thld = float(request.GET.get('eQTL_lod_thld').strip())
        search_gene = request.session['search_gene']
        experiment = request.session['experiment_name']
        # output_dic will be used to generate JSON file.
        output_dic = {}
        chr_nr = 5
        #define KEY chrnames
        # alternative way if it is used to another organism
        # for efficiency: define it derectly, or comment the following line
        # chr_size = Marker.objects.values('marker_chromosome').dictinct().count() # and use chr_size to generate ['1','2','3','4','5']
        output_dic['chrnames']= ['1','2','3','4','5']
        
        # define KEY chr
        # it can be also achieved by calling
        # chr_start = Gene.objects.filter('chromosome'=1).orderby('start')[0].values('start') etc...
        output_dic['chr'] = {'1': {'start_Mbp': 0.003631,'end_Mbp': 30.425192},
                             '2': {'start_Mbp': 0.001871,'end_Mbp': 19.696821},
                             '3': {'start_Mbp': 0.004342,'end_Mbp': 23.459800},
                             '4': {'start_Mbp': 0.001180,'end_Mbp': 18.584524},
                             '5': {'start_Mbp': 0.001251,'end_Mbp': 26.970641}}
        
        # define KEY pmarknames
        chr_marker_list_dic = {}
        marker_list = getMarkersList(experiment)
        for i in range(chr_nr):
            chr_marker_list_dic[i+1] = list(getChrMarkers(i+1,marker_list))  #django.db.models.query.ValuesListQuerySet
        output_dic['pmarknames'] = chr_marker_list_dic
        toc1 = time.time()
        print 'getChrMarker in %f seconds' % (toc1-tic)
        
        # define KEY markers
        output_dic['markers'] = list(getMarkerNamesList(marker_list))
        
        toc2 = time.time()
        print 'getMarkerNamesList in %f seconds' % (toc2-tic)
        
        # define KEY pmark
        marker_list_dic = {}
        marker_list = ExperimentMarker.objects.filter(experiment_name = experiment).values_list('marker_name',flat=True) 
        marker_queryset_list = Marker.objects.filter(marker_name__in=marker_list).order_by('marker_chromosome','marker_cm')
        for m in marker_queryset_list:
            m_info ={'chr':m.marker_chromosome,'pos_cM':float(m.marker_cm),'pos_Mbp':float(m.marker_phys_pos)}
            marker_list_dic[m.marker_name] = m_info
        output_dic['pmark'] = marker_list_dic
        
        toc3 = time.time()
        print 'pmark info in %f seconds' % (toc3-tic)
        # define KEY gene
        gene_list_dic = {}
        lod_gene_list = ExperimentGene.objects.filter(experiment_name = experiment).values_list('locus_identifier',flat=True)
        #gene_queryset_list = Gene.objects.filter(locus_identifier__in=lod_gene_list).order_by('chromosome','start')
        print len(lod_gene_list)
        counter = 0
        for gene in lod_gene_list:   
            gene_instance = Gene.objects.get(locus_identifier = gene)
            gene_list_dic[gene]={'chr':gene_instance.chromosome,'pos_Mbp':float(gene_instance.start)/1000000} 
        output_dic['gene'] = gene_list_dic
        
        
        toc4 = time.time()
        print 'gene list in %f seconds' % (toc4-tic)
        
        # define KEY peaks
        #sha_lod_thld = -lod_thld
        peaks_list = []
        neg_lod_thld = -lod_thld
        lod_bay_list = LOD.objects.filter(LOD_score__gte= lod_thld,gxe = False,experiment_name = experiment)
        lod_sha_list = LOD.objects.filter(LOD_score__lte = neg_lod_thld, gxe = False,experiment_name = experiment)
        for lod in lod_bay_list:
            lod_={}
            lod_['gene'] = lod.locus_identifier.locus_identifier
            lod_['marker'] = lod.marker_name.marker_name#.encode('ascii','ignore')
            lod_['lod'] = float(lod.LOD_score) 
            lod_['p'] = 'Bay'      
            peaks_list.append(lod_)
        for lod in lod_sha_list:
            lod_={}
            lod_['gene'] = lod.locus_identifier.locus_identifier
            lod_['marker'] = lod.marker_name.marker_name#.encode('ascii','ignore')
            lod_['lod'] = math.fabs(float(lod.LOD_score))
            lod_['p'] = 'Sha'       
            peaks_list.append(lod_)

        nr_eQTL = len(peaks_list)#number of dectected eQTL
    
        output_dic['peaks'] = peaks_list
        
        toc5 = time.time()
        print 'QTL detection in %f seconds' % (toc5-tic)
        peaks_gene_list = [g.next() for k,g in itertools.groupby(peaks_list,lambda x:x['gene'])]#sort peaks list by 'gene' 
        nr_gene = len(peaks_gene_list) # number of genes that the expression are likely to be regulated by eQTL
        pv = math.pow(10,-lod_thld)
        
        # define KEY exp
        exp_list = []
        exp_ = LOD.objects.filter(locus_identifier = search_gene,gxe=False,experiment_name=experiment)
        for exp in exp_:
            exp_dic = {}
            exp_dic['gene'] = exp.locus_identifier.locus_identifier
            exp_dic['marker'] = exp.marker_name.marker_name#.encode('ascii','ignore')
            exp_dic['lod'] = float(exp.LOD_score)
            exp_list.append(exp_dic)
        output_dic['exp'] = exp_list 
        toc6 = time.time()
        print 'LOD curve in %f seconds' % (toc6-tic)   
        output_dic_js = json.dumps(output_dic)
        search_gene = json.dumps(search_gene.upper())
        
        ########## output json file #########
        #with open('test1212.json','wb') as output:
        #    json.dump(output_dic,output,indent = 4)
        return render_to_response('qtl/eQTL.html',{'output_dic_js':output_dic_js,
                                                   'search_gene': search_gene,
                                                   'lod_thld': lod_thld,
                                                   'nr_eQTL':nr_eQTL,
                                                   'nr_gene':nr_gene,
                                                   'p':pv,
                                                   #'peaks_list':peaks_list,
                                                   'peaks_gene_list':peaks_gene_list})
    
    else:
        return render_to_response('qtl/gene.html',{})

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
    

    
    
######################################
# main #
######################################
#if __name__ == "__main__":
    

