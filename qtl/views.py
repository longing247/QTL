import time,os
import csv
import json
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
from ViewsFunc import getMarkersList,getMarkerNamesList,getMarkerNamesList,getChrMarkers,geneUpload,markerUpload,lodUpload,parentUpload,rilUpload
from ViewsFunc import getAll,getAllMetabolite,getData,metaboliteUpload,metaboliteParentUpload,metaboliteRILUpload,metaboliteLODUpload,envLODUpload,envMLODUpload
from ViewsFunc import findPeak,findMPeak,yield_gene,marker_plot,m_marker_plot,geneUpdate,decimalFormat,genotypeUpload,getExperimentName

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
                    js_peak_lod_env = json.dumps(float('{0:.2f}'.format(peak_lod_env)),cls=DjangoJSONEncoder)
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
                          
            exps = Experiment.objects.exclude(experiment_name = experiment)
            if exps.exists():
                exps_name_list = list (exps.values_list('experiment_name',flat = True))
                exps_lod_list = []
                exp_mul_list = []
                markers_list = list(Marker.objects.all().order_by('marker_chromosome','marker_cm').values_list('marker_name', flat = True))
                for exp_name in exps_name_list:
                    one_lod_list = []
                    if LOD.objects.filter(locus_identifier = search_gene,experiment_name = exp_name,gxe = 0).exists():
                        for marker in markers_list:
                            if LOD.objects.filter(locus_identifier = search_gene, marker_name = marker,experiment_name = exp_name,gxe =0).exists():
                                gxe_lod = LOD.objects.get(locus_identifier = search_gene, marker_name = marker,experiment_name = exp_name,gxe =0).LOD_score
                                one_lod_list.append(float('{0:.2f}'.format(gxe_lod)))
                            else:
                                one_lod_list.append(None) 
                        exps_lod_list.append(one_lod_list)
                        one_lod_list = []
                        exp_mul_list.append(exp_name)
                    if LOD.objects.filter(locus_identifier = search_gene,experiment_name = exp_name,gxe = 1).exists():
                        for marker in markers_list:
                            if LOD.objects.filter(locus_identifier = search_gene, marker_name = marker,experiment_name = exp_name,gxe =1).exists():
                                gxe_lod = LOD.objects.get(locus_identifier = search_gene, marker_name = marker,experiment_name = exp_name,gxe =1).LOD_score
                                one_lod_list.append(float('{0:.2f}'.format(gxe_lod)))
                            else:
                                one_lod_list.append(None)
                        exps_lod_list.append(one_lod_list)
                        one_lod_list = []
                        exp_mul_list.append(exp_name+'GxE')       
                if is_lod_exp:
                    if is_gxp_exp:
                        one_lod_list = []
                        for marker in markers_list:
                            if LOD.objects.filter(locus_identifier = search_gene, marker_name = marker,experiment_name = experiment,gxe =0).exists():
                                gxe_lod = LOD.objects.get(locus_identifier = search_gene, marker_name = marker,experiment_name = experiment,gxe =0).LOD_score
                                one_lod_list.append(float('{0:.2f}'.format(gxe_lod)))
                            else:
                                one_lod_list.append(None)
                        exps_lod_list.append(one_lod_list)
                        exp_mul_list.append(experiment)
                    if is_gxe_exp:
                        one_lod_list = []
                        for marker in markers_list:
                            if LOD.objects.filter(locus_identifier = search_gene, marker_name = marker,experiment_name = experiment,gxe =1).exists():
                                gxe_lod = LOD.objects.get(locus_identifier = search_gene, marker_name = marker,experiment_name = experiment,gxe =1).LOD_score
                                one_lod_list.append(float('{0:.2f}'.format(gxe_lod)))
                                one_lod_list.append(None)
                        exps_lod_list.append(one_lod_list)
                        exp_mul_list.append(experiment+'GxE')
                js_marker_mul_list= json.dumps(markers_list)
                js_lod_mul_list = json.dumps(exps_lod_list,cls=DjangoJSONEncoder)
                js_exp_mul_list = json.dumps(exp_mul_list)
                render_dic['js_marker_mul_list']=js_marker_mul_list#markers along the chromosome 
                render_dic['js_lod_mul_list']=js_lod_mul_list# he corresponding LOD expression value of the searched gene/trait with environmental interaction against the marker list.
                render_dic['js_exp_mul_list'] = js_exp_mul_list
                render_dic['is_multi_exp'] = True
                
            
            
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

def upload_success(request):
    return render_to_response('qtl/success.html')
    
    
######################################
# main #
######################################
#if __name__ == "__main__":
    

