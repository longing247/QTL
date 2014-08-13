'''
Created on Jul 30, 2014

@author: jiao
'''
from django.conf.urls import patterns,url
from qtl import views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = patterns('',

                       ) 
#+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)