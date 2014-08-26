from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome

class Arabidopsis:
    
    
    def draw(self,features_list):
    #The features can either be SeqFeature objects, or tuples of values:
    #start (int), end (int), strand (+1, -1, O or None), label (string),
    #ReportLab color (string or object), and optional ReportLab fill color.
    #for example: feature = {'Chr1':[(8600,8700,None, 'MSAT100008','blue'),()],...}
        entries = [("Chr I", 30427671),#SeqIO.read('NC_003070.gbk',"genbank"): 30427671
                   ("Chr II", 19698289),#19698289
                   ("Chr III", 23459830),#23459830
                   ("Chr IV", 18585056),#18585056
                   ("Chr V", 26985502 )]#26975502 
        
        max_len = 30432563 #Could compute this
        telomere_length = 1000000 #For illustration
        chr_diagram = BasicChromosome.Organism(output_format='jpg')
        #chr_diagram = BasicChromosome.Organism()
        #chr_diagram.page_size = (29.7*cm, 21*cm) #A4 landscape
        #chr_diagram.page_size = (42*cm, 29.7*cm) #A5 landscape
        chr_diagram.page_size = (59.4*cm, 42*cm) #A5 landscape
        for index, (name, length) in enumerate(entries):
            features = features_list[name] 
            
            #features = [f for f in record.features if f.type==name]
            
            #Record an Artemis style integer color in the feature's qualifiers,
            #1 = Black, 2 = Red, 3 = Green, 4 = blue, 5 =cyan, 6 = purple 
            #for f in features: f.qualifiers["color"] = [index+2]
        
            cur_chromosome = BasicChromosome.Chromosome(name)
            #Set the scale to the MAXIMUM length plus the two telomeres in bp,
            #want the same scale used on all five chromosomes so they can be
            #compared to each other
            cur_chromosome.scale_num = max_len + 2 * telomere_length
        
            #Add an opening telomere
            start = BasicChromosome.TelomereSegment()
            start.scale = telomere_length
            cur_chromosome.add(start)
        
            #Add a body - again using bp as the scale length here.
        
            
            body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
            body.scale = length
            cur_chromosome.add(body)
        
            #Add a closing telomere
            end = BasicChromosome.TelomereSegment(inverted=True)
            end.scale = telomere_length
            cur_chromosome.add(end)
        
            #This chromosome is done
            chr_diagram.add(cur_chromosome)
        
        chr_diagram.draw("arabidopsis_chrom_marker.jpg", "Arabidopsis thaliana")
        #Image("arabidopsis_chrom_marker.png")
            