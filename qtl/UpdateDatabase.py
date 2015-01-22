'''
Created on Jan 13, 2015

@author: jiao
'''

import sys
import os

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,Parent,Experiment

def updateSE(f):
    with open(f,'r') as fi:
        i = 0
        for line in fi:
            if i ==0:
                linele = line.strip().split('\t') 
                for j in range(1,9):
                    p = Parent.objects.get(locus_identifier = linele[0][3:].strip().encode('ascii','ignore'),parent_type = linele[j].strip().encode('ascii','ignore'))    
                    p.se = float(linele[j+8].strip())
                    p.save()
                i+=1
            else:
                linele = line.strip().split('\t') 
                for j in range(1,9):
                    p = Parent.objects.get(locus_identifier = linele[0].strip().encode('ascii','ignore'),parent_type = linele[j].strip().encode('ascii','ignore'))    
                    p.se = float(linele[j+8].strip())
                    p.save()
if __name__ == '__main__':
    updateSE('UpdateSEUTF8.txt')