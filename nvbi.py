from Bio import Entrez
import xml.etree.ElementTree as ET

refgenome=''
search_id="GSE145919"
print("Input GEOID is",search_id)
def idselector(tree):
    ''' Returns top 1 id after searching'''
    x=0
    req_id=''
    for node in tree.iter("IdList"):
        for elem in node.iter():
            if(x==1):
                print(elem.tag, elem.text) #1st Result found
                id=elem.text
                break
            x+=1
    return id

Entrez.email = "bhavay18384@iiitd.ac.in"
handle = Entrez.esearch(db="gds", term=search_id)
pp=handle.read()
#print(pp)
tree = ET.fromstring(pp)
print("Selecting the following ID")
id=idselector(tree)
handle = Entrez.esummary(db="gds", id=id)
pp=handle.read()
#print(pp)
tree = ET.fromstring(pp)
targetsra=''
for item in tree.findall('DocSum'):
    for item2 in item.findall("Item"):
        if item2.attrib['Name']=='taxon':
            refgenome=item2.text
        if item2.attrib['Name']=='ExtRelations':
            for item3 in item2.findall("Item"):
                if item3.attrib['Name']=='ExtRelation':
                    for item4 in item3.findall("Item"):
                        if item4.attrib['Name']=='TargetObject':
                            targetsra=item4.text

print("SRA in relation is ",targetsra) #SRA to  be searched again
handle = Entrez.esearch(db="sra", term=targetsra)
pp=handle.read()
#print(pp)
tree = ET.fromstring(pp)
print("Selecting the following ID")
id=idselector(tree)
variable=Entrez.efetch(db="sra",id=id,rettype="full")
#print(variable.read())

#import os
#os.system("pysradb metadata SRP250724")
print("Organism is",refgenome)
from pysradb import sraweb,download 
db = sraweb.SRAweb()
db.download(targetsra,skip_confirmation=True)