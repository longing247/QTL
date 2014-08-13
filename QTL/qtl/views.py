import time,os
#import numpy as np
import csv

from django.shortcuts import render,get_object_or_404,render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.http import HttpRequest as request
from django.template import RequestContext
from django.views.generic import FormView
from django.conf import settings 
from django.core.context_processors import csrf
from django.contrib import messages

from .models import Experiment,Gene,Marker,LOD
from .forms import GeneUploadFileForm,MarkerUploadFileForm, EXPUploadFileForm

'''
WUR Seed Lab eQTL expression dataset was used for test purpose. 
The derived files: gene_test.txt, marker_test.txt, and lod_test.txt are all from the original datadataset.
All the files were normalized use tab delimited text file format with extension (.txt) and were converted from Unicode to UTF8 by Notepad in order to be compatible with MySQL.
Python built-in module csv was invoked to parse the Django upoloaded files.
'''


def geneuploadView(request):
    '''
    handling gene file upload event through HTTP request method
    URL: ~/qtl/gene/
    '''
    
    args = {}
    args.update(csrf(request))
    if request.method=='POST':
        form = GeneUploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            geneUpload(request.FILES['genefile'])# request.FILES['genefile'] is from <input name='*' /> and need to be identical to form attribute(genefile).
            return HttpResponseRedirect('success/')

        else:
            print 'form is not valid'
            messages.error(request,'Error')
    else:
        form = GeneUploadFileForm()      
    return render_to_response('qtl/geneupload.html',args)

def markeruploadView(request):
    '''
    handling marker file upload event through HTTP request method
    URL: ~/qtl/marker/
    '''
    
    args = {}
    args.update(csrf(request))
    if request.method=='POST':
        form = MarkerUploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            markerUpload(request.FILES['markerfile'])
            return HttpResponseRedirect('success/')
        else:
            print 'form is not valid'
            messages.error(request,'Error')
    else:
        form = MarkerUploadFileForm()
    return render_to_response('qtl/markerupload.html',args)

def experimentuploadView(request):
    '''
    handling expriment file upload event through Http request method
    URL: ~/qtl/experiment/ 
    '''
    
    args = {}
    args.update(csrf(request))
    if request.method=='POST':
        #file_instance = UploadFile()
        form = EXPUploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            expUpload(request.FILES['expfile'])
            return HttpResponseRedirect('success/')

        else:
            messages.error(request,'Error')
    else:
        form = EXPUploadFileForm()       
    return render_to_response('qtl/experimentupload.html',args)

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
    populating makrer data into MySQL under table qtl_marker table.  
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
def expUpload(f):   
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
    Solution:
    delete from qtl_lod where locus_identifier_id = 'AT2G36270';
    split the data set from that row into another one.
    Anyhow this part of function needs to be optimized by either using mmap or generator.
    '''
    tic = time.clock()   
    csv_reader = csv.reader(f,delimiter='\t')
    print f.name
    add_exp = Experiment()
    add_exp.experiment_name = f.name
    add_exp.save()
    marker_list = []  
    i=0
    for line in csv_reader:
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
    
def upload_success(request):
    return render_to_response('qtl/success.html')


    
######################################
# main #
######################################
#if __name__ == "__main__":
    

