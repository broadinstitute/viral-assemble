"""Extract genes from a sequin file, and output in quast's four-column format"""

import fileinput

seqName = None
geneBeg = None
geneEnd = None
genesSeen = set()
for line in fileinput.input():
    L = line.strip().split()
    #print(L)
    if len(L)==2 and L[0]=='>Feature' and L[1].startswith('gb|'):
        seqName = L[1].split('|')[1]
    elif geneBeg is None and geneEnd is None and len(L)==3 and L[0].isdigit() and L[1].isdigit() and L[2]=='gene':
        geneBeg, geneEnd, gene = line.strip().split()
    elif all((seqName,geneBeg,geneEnd)) and len(L)==2 and L[0]=='gene':
        geneName = L[1]
        assert geneName not in genesSeen
        genesSeen.add( geneName )
        #print("seqName=",seqName,'geneName=',geneName,'geneBeg=',geneBeg,'geneEnd=',geneEnd)
        print('\t'.join((seqName,geneName,geneBeg,geneEnd)))
        geneBeg = None
        geneEnd = None
        

        
    
        
        
