'''
Created on Jan 13, 2015

@author: jiao
'''

import sys
import os
import csv

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,Parent,Experiment

def updateSE(f):
    csv_reader = csv.reader(f,delimiter='\t')
    i = 0
    
    for line in csv_reader:
        i+=1
        linele = line.strip().split('\t') 
        print linele[0]
        for j in range(1,9):
            g = Gene.objects.get(locus_identifier = linele[0].strip(),parent_type = linele[j].strip())    
            g.se = linele[j+8].strip()
            g.save()


if __name__ == '__main__':
    updateSE('')