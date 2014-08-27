'''
Created on Jul 30, 2014

@author: jiao
'''
from django import forms


class GeneUploadFileForm(forms.Form):
    '''
    classdocs
    '''
    #title = forms.CharField(max_length=50)
    #file = forms.FileField()
    genefile = forms.FileField()

class MarkerUploadFileForm(forms.Form):
    '''
    classdocs
    '''
    #title = forms.CharField(max_length=50)
    #file = forms.FileField()
    markerfile = forms.FileField()
    
class EXPUploadFileForm(forms.Form):
    '''
    classdocs
    '''
    #title = forms.CharField(max_length=50)
    #file = forms.FileField()
    expfile = forms.FileField()

class ParentUploadFileForm(forms.Form):
    '''
    classdocs
    '''
    #title = forms.CharField(max_length=50)
    #file = forms.FileField()
    parentfile = forms.FileField()
