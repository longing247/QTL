from django.db import models


class Gene(models.Model):
    
    
    locus_identifier = models.CharField(max_length=30,primary_key=True) #AT1G01480
    gene_model_name = models.CharField(max_length=30,blank = True)#AT1G01480.1
    gene_model_description = models.TextField(blank = True)
    gene_model_type = models.CharField(max_length = 40,blank = True)
    primary_gene_symbol = models.TextField(blank = True)
    all_gene_symbols = models.TextField(blank = True)
    #tair_accession = models.CharField(max_length=30) #Locus:2025361
    #gene__type = models.CharField(max_length=15) # protein_coding
    #associated_loci = models.TextField()
    #TEST FOR SYNCHRONIZATION
    
    def __unicode__(self):
        return self.locus_identifier
    
    
class Marker(models.Model):
    
    
    marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
    marker_chromosome = models.IntegerField()#1
    marker_cm = models.DecimalField(max_digits = 3, decimal_places =1)#64.6
    marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639
    #associated_locus = models.ForeignKey(Locus) #AT1G01480
    #marker_aliases = models.CharField(max_length=15) #PVV4
    #tair_accession = models.CharField(max_length=30) #GeneticMarker:1945638
    #marker_type = models.CharField(max_length=10) #CAPS
    #marker_length = models.CharField() #The data format of marker_length is like 1.000 (bp) 
    #is_PCR_marker = models.BooleanField() #true
    #special_condition = models.CharField()
    #chromosome = models.IntegerField() #1
    
    #map = models.ManyToManyField(Map) or add an instance in between to avoid of many-to-many relationship
        
    def __unicode__(self):
        return self.marker_name

class Experiment(models.Model):
    
    
    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    def __unicode__(self):
        return self.experiment_name
    
class LOD(models.Model):
    
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    locus_identifier = models.ForeignKey(Gene)
    marker_name = models.ForeignKey(Marker)
    
    def __unicode__(self):
        return self.LOD_score
    

    
#class GO(models.Model):
# is it necessary to split it into 3 categories


#    GO_accession_number = models.IntegerField()
#    GO_term = models.TextField()
    
#    def __unicode__(self):
#        return self.marker_name

    
    
#class Polymorphism(models.Model):
#Polymorphism specie variant

    
    #marker = models.ForeignKey(Marker)
    #variant = models.TextField()

#class Map(models.Model):    
#Types of maps in TAIR. (null,genetic, nuc_sequence,physical_framework, etc..)
    #map = models.ForeignKey(Map)
    #map_name = models.CharField(max_length=15)
    #map_type = models.charField(max_length=20)
    #map_link = models.charField()

#    def __unicode__(self):
#        return self.title



