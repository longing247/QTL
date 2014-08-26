'''
Created on Aug 13, 2014

@author: jiao
'''
import datetime
from haystack import indexes
from .models import Parent
    
    
class ParentIndex(indexes.SearchIndex, indexes.Indexable):
    '''
    parent_type = models.CharField(max_length=20)
    expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    locus_identifier = models.ForeignKey(Gene)
    '''
    text = indexes.CharField(document = True, use_template = True)
    parent = indexes.CharField(model_attr='parent_type')
    exp = indexes.DecimalField(model_attr='expression')
  

    def get_model(self):
        return Parent
    
    def index_queryset(self,using=None):
        '''
        Used when the entire index for model is updated.
        Further it can be customized to acts objects.filter()
        '''
        return self.get_model().objects.all()
        