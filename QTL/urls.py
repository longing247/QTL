from django.conf.urls import patterns, include, url

from django.contrib import admin
from qtl import views
admin.autodiscover()


urlpatterns = patterns('',


    url(r'^admin/', include(admin.site.urls)),
    url(r'^qtl/', include('qtl.urls')),
    #url(r'^$', views.indexView,name='index'),

    
)
