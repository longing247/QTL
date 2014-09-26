from django.conf.urls import patterns, include, url

from django.contrib import admin
from qtl import views
admin.autodiscover()


urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'QTL.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),

    url(r'^admin/', include(admin.site.urls)),
    url(r'^qtl/', include('qtl.urls')),
    url(r'^$', views.indexView,name='index'),
#    url(r'^qtl/',views.qtlView,name='qtlView'),
    url(r'^search/', include('haystack.urls')), 
    
)