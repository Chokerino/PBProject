from Bio import Entrez
import xml.etree.ElementTree as ET
import sys
from pysradb import sraweb, download
import os
import multiprocessing


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

def createFastqFiles(sraFileNames):
    oneNames=[];secondNames=[]
    for i in range(len(sraFileNames)):
        command = "fastq-dump --split-3 {}".format(sraFileNames[i])
        os.system(command)
        if(os.path.exists(sraFileNames[i]+"_1.fastq")):
            oneNames.append(sraFileNames[i]+"_1.fastq")
            secondNames.append(sraFileNames[i]+"_2.fastq")
        else:
            oneNames.append(sraFileNames[i]+".fastq")
    return oneNames,secondNames

def createIndex(refGenome:str,annotations:str):
    command = "hisat-build -p {} {} {}".format(multiprocessing.cpu_count(),refGenome,annotations)
    os.system(command)
    return annotations

def startAlignment(indexLocation:str, firstFileNamesList:list, secondFilesNameList=[]):
    outputFilenames = [];
    if(len(secondFilesNameList)==0):
        for i in range(len(firstFileNamesList)):
            command = "hisat2 -x {} -p {} -U {} -S output_{}.sam".format(indexLocation, multiprocessing.cpu_count(), firstFileNamesList[i],i+1)
            os.system(command)
            outputName = "output_{}.sam".format(i+1);
            outputFilenames.append(outputFilenames)
    else:
        for i in range(len(firstFileNamesList)):
            command = "hisat2 -x {} -p {} -1 {} -2 {} -S output_{}.sam".format(indexLocation, multiprocessing.cpu_count(), firstFileNamesList[i], secondFilesNameList[i],i+1)
            os.system(command)
            outputName = "output_{}.sam".format(i + 1);
            outputFilenames.append(outputFilenames)
    return outputFilenames


def convertSamToBam(samFilenames):
    bamFileNames = []
    for i in range(len(samFilenames)):
        command = "samtools view -S -b {} > {}".format(samFilenames[i],samFilenames[i][0:-3]+"bam")
        os.system(command)
        bamFileNames.append(samFilenames[i][0:-3]+"bam")
    return bamFileNames


def sortBamFiles(BamfileNames: list):
    sortedBamFileNames = []
    for i in range(len(BamfileNames)):
        command = "samtools sort {} -o {}".format(BamfileNames[i],BamfileNames[i][0:-3]+"sorted.bam")
        os.system(command)
        sortedBamFileNames.append(BamfileNames[i][0:-3]+"sorted.bam")
    return sortedBamFileNames

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

    sraFileNames = ["Put List of All SRA File Names"]

    # Creating fastq file from Downloaded SRA File
    firstList,secondList = createFastqFiles(sraFileNames)

    # Quality Control
    os.system("fastqc *fastq") # Can be Modified

    # To Download Refernce Genome
    fnaName,gtfName = downloadRefGenome()

    # Creating Index for Hisat2
    index_name = createIndex(fnaName,gtfName)

    # Starting Alignment //Done
    listOfOutputFiles = startAlignment(index_name,firstList,secondList)

    # Converting Sam to Bam // Done
    BamFileNames = convertSamToBam(listOfOutputFiles)

    # Sort Bam Files // Done
    sortedBamFileNames = sortBamFiles(BamFileNames)

    # Create Count Matrix (User Option)
    createCountMatrix()