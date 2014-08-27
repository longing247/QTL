from django.contrib import admin

from .models import Gene,Marker,LOD

'''
class GeneInLine(admin.StackedInline):
    model = Gene
    extra = 1
    #fields = [ 'experiment_name_id','locus_identifier_id', 'marker_name_id']

class MarkerInLine(admin.StackedInline):
    model = Marker
    extra = 1
    #fields = [ 'experiment_name_id','locus_identifier_id', 'marker_name_id']    
'''    
class LODAdmin(admin.ModelAdmin):

    list_display = ( 'experiment_name','locus_identifier', 'marker_name','LOD_score')
    search_fields = ['locus_identifier__locus_identifier']

admin.site.register(LOD,LODAdmin)   