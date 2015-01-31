from django.conf.urls import patterns, include, url

from django.contrib import admin
from qtl.views import *
admin.autodiscover()


urlpatterns = patterns('',


    url(r'^admin/', include(admin.site.urls)),
    #url(r'^qtl/', include('qtl.urls')),
    url(r'^$', searchGeneView,name='index'),
    url(r'^about/$',aboutView,name='aboutView'),
    url(r'^upload/$',uploadView,name='uploadView'),   
    url(r'^upload/success/$',upload_success,name='uploadSuccessView'),                 
    url(r'^gene/success/$',upload_success,name='uploadSuccessView'),  
    url(r'^marker/success/$',upload_success,name='uploadSuccessView'), 
    url(r'^experiment/success/$',upload_success,name='uploadSuccessView'), 
    url(r'^parent/success/$',upload_success,name='uploadSuccessView'),
    url(r'^gene/overlap/$',searchOverlapQTLView,name='search_overlap'), #/qtl/gene/overlap?lod_threshold=&trait=
    url(r'^gene/$',searchGeneView,name='search_gene'),
    url(r'^marker/$',searchMarkerView,name='search_marker'),
    url(r'^metabolite/$',searchMetaboliteView,name='search_metabolite'),
    url(r'^gene/chromosome/$',chromosomeView,name='chromosome'),
    url(r'^gene/eQTLPlot/$',eQTLPlotView,name='eQTL'),
    url(r'^chromosome/success/$',upload_success,name='success'),
    url(r'^documentation/$',documentationView,name='documentation'),
    
)
