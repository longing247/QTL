'''
Created on Jul 30, 2014

@author: jiao
'''
from django.conf.urls import patterns,url
from qtl import views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = patterns('',
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
    url(r'^chromosome/$',views.chromosomeView,name='index'),
    url(r'^chromosome/success/$',views.upload_success,name='success'),
    ) 
#+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)