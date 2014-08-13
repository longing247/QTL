'''
Created on Jul 30, 2014

@author: jiao
'''
from django.conf.urls import patterns,url
from qtl import views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = patterns('',
    url(r'^gene/success/$',views.upload_success,name='uploadSuccessView'),  
    url(r'^gene/$',views.geneuploadView,name='gene_upload'),  
    url(r'^marker/$',views.markeruploadView,name='marker_upload'), 
    url(r'^marker/success/$',views.upload_success,name='uploadSuccessView'), 
    url(r'^experiment/$',views.experimentuploadView,name='experiment_upload'), 
    url(r'^experiment/success/$',views.upload_success,name='uploadSuccessView'), 
    url(r'^parent/$',views.parentuploadView,name='parent_upload'), 
    url(r'^parent/success/$',views.upload_success,name='uploadSuccessView'),
    ) 
#+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)