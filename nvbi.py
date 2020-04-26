from Bio import Entrez
import xml.etree.ElementTree as ET
import sys
from pysradb import sraweb, download
import os

# This function have been moved up
def idselector(tree):
    ''' Returns top 1 id after searching'''
    x = 0
    req_id = ''
    for node in tree.iter("IdList"):
        for elem in node.iter():
            if (x == 1):
                print(elem.tag, elem.text)  # 1st Result found
                id = elem.text
                break
            x += 1
    return id

def downloadRefGenome():
    pass;

def createIndex(filename:str):
    os.system("hisat2-build "+filename);

def startAlignment(indexLocation, firstName, secondName):
    os.system("hisat2 -x "+indexLocation+" -m1 "+firstName+" -m2 "+secondName+" -S output.sam")

def convertSamToBam(samFilename):
    os.system("samtools view -S -b "+samFilename+" > output.bam")

def sortBamFiles(filename):
    os.system("samtools sort "+filename+" -o sample.sorted.bam")

def createCountMatrix():
    # http://genomespot.blogspot.com/2015/01/generate-rna-seq-count-matrix-with.html
    pass

if __name__ == '__main__':
    refgenome = ''

    # Example ID
    search_id = "GSE145919"

    # Over-writting with user Entered ID, only if present
    if(len(sys.argv)>1):
        search_id=sys.argv[1]

    print("Input GEOID is", search_id)


    Entrez.email = "bhavay18384@iiitd.ac.in"
    handle = Entrez.esearch(db="gds", term=search_id)
    pp = handle.read()
    # print(pp)
    tree = ET.fromstring(pp)
    print("Selecting the following ID")
    id = idselector(tree)
    handle = Entrez.esummary(db="gds", id=id)
    pp = handle.read()
    # print(pp)
    tree = ET.fromstring(pp)
    targetsra = ''
    for item in tree.findall('DocSum'):
        for item2 in item.findall("Item"):
            if item2.attrib['Name'] == 'taxon':
                refgenome = item2.text
            if item2.attrib['Name'] == 'ExtRelations':
                for item3 in item2.findall("Item"):
                    if item3.attrib['Name'] == 'ExtRelation':
                        for item4 in item3.findall("Item"):
                            if item4.attrib['Name'] == 'TargetObject':
                                targetsra = item4.text


    # Fetching Target SRA
    print("SRA in relation is ", targetsra)  # SRA to  be searched again
    handle = Entrez.esearch(db="sra", term=targetsra)
    pp = handle.read()
    # print(pp)
    tree = ET.fromstring(pp)
    print("Selecting the following ID")
    id = idselector(tree)
    variable = Entrez.efetch(db="sra", id=id, rettype="full")
    # print(variable.read())

    # import os
    # os.system("pysradb metadata SRP250724")
    print("Organism is", refgenome)

    targetsra="app*"
    # Downloading .SRA file
    # db = sraweb.SRAweb()
    # db.download(targetsra, skip_confirmation=True);
    # SRR11550936
    # Searching for Downloaded SRA File in File System
    os.system("find . -name "+targetsra +"> downloadPath.txt")
    downloadFileLocation = open("downloadPath.txt","r").readline();

    # Creating fastq file from Downloaded SRA File
    os.system("fastq-dump --split-3 "+downloadFileLocation)

    # Quality Control
    os.system("fastqc *fastq") # Can be Modified

    # To Download Refernce Genome
    filename = downloadRefGenome()

    # Creating Index for Hisat2
    createIndex(filename)

    # Starting Alignment
    startAlignment("indexName","firstName","secondName")

    # Converting Sam to Bam
    convertSamToBam("samFileName")

    # Sort Bam Files
    sortBamFiles("BAMFilename")

    # Create Count Matrix (User Option)
    createCountMatrix()