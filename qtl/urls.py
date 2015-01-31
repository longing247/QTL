'''
Created on Jul 30, 2014

@author: jiao
'''
from django.conf.urls import patterns,url
from qtl import views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = patterns('',
    url(r'^$',views.searchGeneView,name='search_gene'),
    url(r'^about/$',views.aboutView,name='aboutView'),
    url(r'^upload/$',views.uploadView,name='uploadView'),   
    url(r'^upload/success/$',views.upload_success,name='uploadSuccessView'),                 
    url(r'^gene/success/$',views.upload_success,name='uploadSuccessView'),  
    url(r'^marker/success/$',views.upload_success,name='uploadSuccessView'), 
    url(r'^experiment/success/$',views.upload_success,name='uploadSuccessView'), 
    url(r'^parent/success/$',views.upload_success,name='uploadSuccessView'),
    url(r'^gene/overlap/$',views.searchOverlapQTLView,name='search_overlap'), #/qtl/gene/overlap?lod_threshold=&trait=
    url(r'^gene/$',views.searchGeneView,name='search_gene'),
    url(r'^marker/$',views.searchMarkerView,name='search_marker'),
    url(r'^metabolite/$',views.searchMetaboliteView,name='search_metabolite'),
    url(r'^gene/chromosome/$',views.chromosomeView,name='chromosome'),
    url(r'^gene/eQTLPlot/$',views.eQTLPlotView,name='eQTL'),
    url(r'^chromosome/success/$',views.upload_success,name='success'),
    url(r'^documentation/$',views.documentationView,name='documentation'),
    
    ) 
#+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)