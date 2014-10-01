'''
Created on Oct 1, 2014

@author: jiao
'''

from django import template

register = template.Library()

@register.filter
def get_key(dict_name,arg):
    return dict_name.get(arg.encode('ascii','ignore'),'')