'''
Created on Jul 30, 2014

@author: jiao
'''
from django import forms

class GeneUploadFileForm(forms.Form):
    genefile = forms.FileField()

class MarkerUploadFileForm(forms.Form):
    markerfile = forms.FileField()
    
class LODUploadFileForm(forms.Form):
    expfile = forms.FileField()

class ParentUploadFileForm(forms.Form):
    parentfile = forms.FileField()
    
class RILUploadFileForm(forms.Form):
    rilfile = forms.FileField()

class MetaboliteUploadFileForm(forms.Form):
    metabolitefile = forms.FileField()
    
class MParentUploadFileForm(forms.Form):
    MParentFile = forms.FileField()

class MRILUploadFileForm(forms.Form):
    MRILFile = forms.FileField()

class MLODUploadFileForm(forms.Form):
    MLODFile = forms.FileField()

class ENVLODUploadFileForm(forms.Form):
    envLODFile = forms.FileField()
    
class ENVMLODUploadFileForm(forms.Form):
    envMLODFile = forms.FileField()
    
class GeneUpdateFileForm(forms.Form):    
    geneUpDateFile = forms.FileField()