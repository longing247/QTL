'''
Created on Sep 23, 2014

@author: jiao
'''
import time
import sys
import os
import math

sys.path.append('/home/jiao/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,ExperimentMarker,RIL
from django.db.models import Q


def mysqlCorrelationAll(gene_name):
    '''
    calculate pearson correlation coefficient from MySQL server side for a certain query gene.
    execution time was reduced to 87 sec, but still too long.
    '''
    tic = time.time()
    #the primary key field can not be leaved out.
    query_script = '''SELECT id, name,
                ((psum - (sum1 * sum2 / n)) / sqrt((sum1sq - pow(sum1, 2.0) / n) * (sum2sq - pow(sum2, 2.0) / n))) AS r
                FROM 
                    (SELECT
                        target_gene.id AS id,
                        target_gene.locus_identifier_id AS name,
                        SUM(query_gene.ril_exp) AS sum1,
                        SUM(target_gene.ril_exp) AS sum2,
                        SUM(query_gene.ril_exp * query_gene.ril_exp) AS sum1sq,
                        SUM(target_gene.ril_exp * target_gene.ril_exp) AS sum2sq,
                        SUM(query_gene.ril_exp * target_gene.ril_exp) AS psum,
                        COUNT(*) AS n  
                    FROM
                        qtl_ril AS query_gene
                    LEFT JOIN
                        qtl_ril AS target_gene
                    ON
                        query_gene.ril_name = target_gene.ril_name
                    WHERE
                        query_gene.locus_identifier_id = %s AND query_gene.locus_identifier_id <> target_gene.locus_identifier_id
                    GROUP BY
                        query_gene.locus_identifier_id, target_gene.locus_identifier_id) AS CORR
                ORDER BY r DESC
                '''
    gene_corr = RIL.objects.raw(query_script,[gene_name])  
    
    toc = time.time()
    print 'in %f seconds' % (toc-tic)
    #return gene_list,corr_list
    return gene_corr

def mysqlCorrelationSingle(gene_name,target_gene_name):
    '''
    MySQL backed script to calculate Pearson correlation of expression value of two genes
    '''
    tic = time.time()
    #the primary key field can not be leaved out.
    query_script = '''SELECT id, name,
                ((psum - (sum1 * sum2 / n)) / sqrt((sum1sq - pow(sum1, 2.0) / n) * (sum2sq - pow(sum2, 2.0) / n))) AS r
                FROM 
                    (SELECT
                        target_gene.id AS id,
                        target_gene.locus_identifier_id AS name,
                        SUM(query_gene.ril_exp) AS sum1,
                        SUM(target_gene.ril_exp) AS sum2,
                        SUM(query_gene.ril_exp * query_gene.ril_exp) AS sum1sq,
                        SUM(target_gene.ril_exp * target_gene.ril_exp) AS sum2sq,
                        SUM(query_gene.ril_exp * target_gene.ril_exp) AS psum,
                        COUNT(*) AS n  
                    FROM
                        qtl_ril AS query_gene
                    LEFT JOIN
                        qtl_ril AS target_gene
                    ON
                        query_gene.ril_name = target_gene.ril_name
                    WHERE
                        query_gene.locus_identifier_id = %s AND target_gene.locus_identifier_id = %s
                    GROUP BY
                        query_gene.locus_identifier_id, target_gene.locus_identifier_id) AS CORR
                ORDER BY r DESC
                '''
    gene_corr = RIL.objects.raw(query_script,[gene_name,target_gene_name])  
    #gene_list = []
    #corr_list = []
    #for gene in gene_corr:
    #    gene_list.append(gene.name)
    #    corr_list.append(gene.r)
    toc = time.time()
    print 'in %f seconds' % (toc-tic)
    #return gene_list,corr_list
    return gene_corr     
   
if __name__=="__main__":   
    query = mysqlCorrelationSingle('AT1G09950','AT5G53700')
    print query[0].r