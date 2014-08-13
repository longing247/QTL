from django.conf.urls import patterns, include, url

from django.contrib import admin
from qtl import views
admin.autodiscover()


urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'QTL.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),

    url(r'^admin/', include(admin.site.urls)),
#    url(r'^qtl/',views.qtlView,name='qtlView'),
    url(r'^qtl/gene/success/$',views.upload_success,name='uploadSuccessView'),  
    url(r'^qtl/gene/$',views.geneuploadView,name='gene_upload'),  
    url(r'^qtl/marker/$',views.markeruploadView,name='marker_upload'), 
    url(r'^qtl/marker/success/$',views.upload_success,name='uploadSuccessView'), 
    url(r'^qtl/experiment/$',views.experimentuploadView,name='experiment_upload'),    
    
)
