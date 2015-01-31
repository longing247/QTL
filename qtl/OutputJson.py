'''
Created on Oct 15, 2014

@author: jiao
'''

import json
from django.core.serializers.json import DjangoJSONEncoder

def outputJson(data,file_name):
    with open(file_name,'w') as outfile:
        json.dump(data,outfile)


if __name__ == '__main__':
    test = {"gene":{"gene1":[1,2,3,4,5]}}
    outputJson(test,'test.json')