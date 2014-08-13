

from django.contrib import admin

from .models import Gene

class GeneAdmin(admin.ModelAdmin):
    fields = [ 'id','experiment_name']
    
admin.site.register(Gene,GeneAdmin)    