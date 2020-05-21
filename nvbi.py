from Bio import Entrez
import xml.etree.ElementTree as ET
import sys
from pysradb import sraweb, download
import os
import multiprocessing
from Bio import SeqIO
import re


# esearch -db assembly -query "hg19" | esummary
# wget

# This function have been moved up
def idselector(tree):
    """ Returns top 1 id after searching"""
    x = 0
    req_id = ""
    for node in tree.iter("IdList"):
        for elem in node.iter():
            if x == 1:
                print(elem.tag, elem.text)  # 1st Result found
                id = elem.text
                break
            x += 1
    return id


def writeMessage(geoid,message):
    #pythonoutputfile = open("~/{}/python_script_realtime_log.txt".format(geoid),"a")
    #pythonoutputfile.write(message+"\n")
    #pythonoutputfile.close()
    return 0

"""
    given GEOID (GSE)
    This Function esearch it's corresponding GSM ID
    and then downloads it's html file and returns that
"""


def getHTML(GEOID):
    command = "(esearch -db gds --query {} | esummary) > sampleScraped.xml".format(
        GEOID
    )
    os.system(command)
    sampleScrapedFile = open("sampleScraped.xml", "r")
    root = ET.parse(sampleScrapedFile)
    SampleID = None
    for id in root.iter("Accession"):
        SampleID = id.text
        if SampleID[0:3] == "GSM":
            break
    command = "wget https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={} -O experimentSample.html".format(
        SampleID
    )
    os.system(command)
    sampleScrapedFile.close()
    htmlFile = open("experimentSample.html", "r")
    return htmlFile


"""
    Still Experimental (Under Build)
    For Given HTML File this function scrapes that HTML file
    for genome_build and if it exists returns genome_build
    else None
"""


def findGenomeBuild(htmlFile):
    genomeBuild = None
    listContainingBuild = None
    for line in htmlFile:
        if "genome_build" in line.lower() or "genome build" in line.lower():
            listContainingBuild = line.lower().split("<")
            break
    words = []
    if listContainingBuild == None:
        return None
    else:
        for sentences in listContainingBuild:
            if "genome_build" in sentences or "genome build" in sentences:
                words = sentences.split(" ")
    for i in range(len(words)):
        if "genome_build" in words[i]:
            genomeBuild = words[i + 1]
        if "genome" in words[i] and "build" in words[i + 1]:
            genomeBuild = words[i + 2]
    return genomeBuild


"""
    This Function Intakes orgainsm name and then returns
    refseq_id for it's latest genome build 
"""


def getLatestRefGenomeName(org_name: str):
    command = "./datasets assembly_descriptors tax_name '{}' --refseq > orgRefSeqInfo.txt".format(
        org_name
    )
    refSeqAccessionID = None
    file = open("orgRefSeqInfo.txt", "r")
    textOfAccessionFile = []
    for line in file:
        textOfAccessionFile.extend(line.split('"'))

    for words in range(len(textOfAccessionFile)):
        if "assembly_accession" in textOfAccessionFile[words]:
            refSeqAccessionID = textOfAccessionFile[words + 2]
            break
    return refSeqAccessionID


"""
    For Given Refseq id this function downloads reference Genome
    and it's annotations and returns the name of the downloaded
    return annotation File name, reference Genome Filename 
"""


def downloadGenomeUsingWget(refSeqId):
    command = "(esearch -db assembly --query '{}' | esummary) > refGenomeScrape.xml".format(
        refSeqId
    )
    os.system(command)
    ftpPathToGenome = None
    refGenomeScrape = open("refGenomeScrape.xml", "r")
    root = ET.parse(refGenomeScrape)
    for i in root.iter("FtpPath_RefSeq"):
        if len(str(i) != 0):
            ftpPathToGenome = i
    refGenomeScrape.close()
    name = ftpPathToGenome.split("/")
    command = "wget -cNrv -t 45 '{}/{}_genomic.gtf.gz' -O {}_annotations.gz".format(
        ftpPathToGenome, name[-1], name[-1]
    )
    os.system(command)
    annotationFileName = "{}_annotations.gz".format(name[-1])
    command = "wget -cNrv -t 45 '{}/{}_genomic.fna.gz' -O {}_genome.gz".format(
        ftpPathToGenome, name[-1], name[-1]
    )
    os.system(command)
    refgenomeFileName = "{}_genome.gz".format(name[-1])
    return annotationFileName, refgenomeFileName


"""
    Put GEOID->eSearch->  get GSM_ID from Here
    Put GSM_ID -> wget -> get HTML Page
    Scrape HTML_PAGE -> Genome_Build
    If Genome_Build not Exist
    do 
        ./dataset with organism name -> Refseq Accession ID  (GCF_XXXX )
        use wget to download Reference genome using Accession ID
    else if Genome Build Exists
    do
        esearch -db assembly --query build_name -> Refseq Accession ID 
        use wget to download Reference genome using Accession ID
"""


def downloadRefGenome(GEOID: str, OrganismName):
    annotationFileName, refGenomeFileName = None, None
    htmlFile = getHTML(GEOID)
    genomeBuild = findGenomeBuild(htmlFile)
    if genomeBuild == None:
        refGenomName = getLatestRefGenomeName(OrganismName)
        annotationFileName, refGenomeFileName = downloadGenomeUsingWget(refGenomName)
        return annotationFileName, refGenomeFileName
    else:
        print(genomeBuild, "<- Genome Build")
        # Search for Refseq id in esearch assmbly Results
        command = "(esearch -db assembly --query '{}' | esummary)> refGenomeScrape.xml".format(
            genomeBuild
        )  
        os.system(command)
        refGenomeScrapeFile = open("refGenomeScrape.xml", "r")
        root = ET.parse(refGenomeScrapeFile)
        refGenomeName = None
        for i in root.iter("FtpPath_RefSeq"):
            if len(str(i) != 0):
                refGenomeName = i
        # If RefSeq id doesn't Exists in esearch, then download latest genome
        if refGenomeName == None:
            refGenomName = getLatestRefGenomeName(OrganismName)
            annotationFileName, refGenomeFileName = downloadGenomeUsingWget(
                refGenomName
            )
            return annotationFileName, refGenomeFileName
        else:  # If Refseq Id Exists in Esearch scrape then download that build
            annotationFileName, refGenomeFileName = downloadGenomeUsingWget(
                refGenomeName
            )
            return annotationFileName, refGenomeFileName


def copyFilesFromSubdirs(sraFileNames):
    listOfDirs = os.listdir()
    for i in range(len(listOfDirs)):
        os.system("mv {}/* {}".format(listOfDirs[i], os.getcwd()))
    for dirs in listOfDirs:
        os.system("rm -rf {}".format(dirs))
    sraFileNames = os.listdir()
    return sraFileNames

def createFastqFiles(allSRAFiles):
    allSRAFiles = copyFilesFromSubdirs(allSRAFiles)
    oneNames = []
    secondNames = []
    for i in range(len(allSRAFiles)):
        command = "fastq-dump --split-3 {}".format(allSRAFiles[i])
        os.system(command)
        if os.path.exists(allSRAFiles[i] + "_1.fastq"):
            oneNames.append(allSRAFiles[i] + "_1.fastq")
            secondNames.append(allSRAFiles[i] + "_2.fastq")
        else:
            oneNames.append(allSRAFiles[i] + ".fastq")
    return oneNames, secondNames


def preprocess(firstList, secondList, corr_rRNA, corr_trim, sortmernaDbDir):
    """
    corr_rRNA trims out rRNA reads
    corr_trim trims out bases with QS less than threshold
    """
    if corr_rRNA == True:
        firstList, secondList = removerRnaContamination(
            firstList, secondList, sortmernaDbDir, [], []
        )
    if corr_trim == True:
        firstName, secondName = trimmingBadReads(
            firstName=firstList, secondName=secondList
        )
    return firstList, secondList


def removerRnaContamination(
    firstName: list,
    secondName: list,
    sortmernaDbDir,
    outputFirstName: list,
    outputSecondName: list,
):
    """
    Uses sortmerna library to filter out rRNA contamination
    """
    os.system("cp {}/* {}".format(sortmernaDbDir, os.getcwd()))
    referenceStr = "--ref "  # Complete this string
    if secondName == []:
        for i, read in enumerate(firstName):
            os.system(
                "sortmerna {} --reads {} --fastx --other {}/rRNAfiltered{}".format(
                    referenceStr, read, os.getcwd(), i
                )
            )
            outputFirstName.append("rRNAfiltered{}.fastq".format(i))
    else:
        for i in range(len(firstName)):
            os.system(
                "sortmerna {} --reads {} --reads {} --fastx --other {}/rRNAfiltered{} --out2".format(
                    referenceStr, firstName[i], secondName[i], os.getcwd(), i
                )
            )
            outputFirstName.append("rRNAfiltered{}_fwd.fastq".format(i))
            outputSecondName.append("rRNAfiltered{}_rev.fastq".format(i))
    os.system("rm -f {}").format(referenceStr)  # Complete this
    return (outputFirstName, outputSecondName)


def trimmingBadReads(phred=33, firstName=[], secondName=[]):
    """
    uses trimmomatic to trim out bases with low scores,
    uses default parameters i.e. removes bases with quality score < 20 in windows of 5
    """
    illuminaclip_adapters = (
        "ILLUMINACLIP:/usr/local/TRIMHOME/adapters/TruSeq3-SE.fa:2:30:10"
    )
    illuminaclip_Attribute = "SLIDINGWINDOW:5:20 MINLEN:30"
    attribute = ""
    if secondName == []:
        attribute = "SE -threads {} -phred{}".format(multiprocessing.cpu_count(), phred)
    else:
        attribute = "PE -threads {} -phred{}".format(multiprocessing.cpu_count(), phred)
    if secondName == []:
        firstNameTrimmed = []
        secondNameTrimmed = []
        for i in range(len(firstName)):
            os.system(
                "java -jar /usr/local/TRIMHOME/trimmomatic.jar {} {}/{}.fastq {}/trimmedRead{}.fastq {} {}".format(
                    attribute,
                    os.getcwd,
                    firstName[i],
                    os.getcwd(),
                    firstName,
                    illuminaclip_adapters,
                    illuminaclip_Attribute,
                )
            )
            firstNameTrimmed.append("trimmedRead{}".format(firstName[i]))
        return firstNameTrimmed, secondNameTrimmed
    else:
        firstNameTrimmed = []
        secondNameTrimmed = []
        for i in range(len(firstName)):
            fwd_input = os.getcwd() + "/" + firstName[i]
            rev_input = os.getcwd() + "/" + secondName[i]
            output_forward_paired = os.getcwd() + "/trimmedRead" + firstName[i]
            output_forward_unpaired = os.getcwd() + "/" + "removethis_1.fastq"
            output_reverse_paired = os.getcwd() + "/trimmedRead" + secondName[i]
            output_reverse_unpaired = os.getcwd() + "/" + "removethis_2.fastq"
            os.system(
                "java -jar /usr/local/TRIMHOME/trimmomatic.jar {} {}.fastq {}.fastq {}.fastq {}.fastq {}.fastq {}.fastq {} {}".format(
                    attribute,
                    fwd_input,
                    rev_input,
                    output_forward_paired,
                    output_forward_unpaired,
                    output_reverse_paired,
                    output_reverse_unpaired,
                    illuminaclip_adapters,
                    illuminaclip_Attribute,
                )
            )
            os.system("rm -f *removethis_1.fastq *removethis_2.fastq")
            firstNameTrimmed.append("trimmedRead{}".format(firstName[i]))
            secondNameTrimmed.append("trimmedRead{}".format(secondName[i]))
        return firstNameTrimmed, secondNameTrimmed


def createIndex(refGenome: str, annotations: str):
    command = "hisat-build -p {} {} {}".format(
        multiprocessing.cpu_count(), refGenome, annotations
    )
    os.system(command)
    return annotations


def startAlignment(
    indexLocation: str, firstFileNamesList: list, secondFilesNameList=[]
):
    outputFilenames = []
    if len(secondFilesNameList) == 0:
        for i in range(len(firstFileNamesList)):
            command = "hisat2 -x {} -p {} -U {} -S output_{}.sam".format(
                indexLocation, multiprocessing.cpu_count(), firstFileNamesList[i], i + 1
            )
            os.system(command)
            outputName = "output_{}.sam".format(i + 1)
            outputFilenames.append(outputFilenames)
    else:
        for i in range(len(firstFileNamesList)):
            command = "hisat2 -x {} -p {} -1 {} -2 {} -S output_{}.sam".format(
                indexLocation,
                multiprocessing.cpu_count(),
                firstFileNamesList[i],
                secondFilesNameList[i],
                i + 1,
            )
            os.system(command)
            outputName = "output_{}.sam".format(i + 1)
            outputFilenames.append(outputFilenames)
    return outputFilenames


def convertSamToBam(samFilenames):
    bamFileNames = []
    for i in range(len(samFilenames)):
        command = "samtools view -S -b {} > {}".format(
            samFilenames[i], samFilenames[i][0:-3] + "bam"
        )
        os.system(command)
        bamFileNames.append(samFilenames[i][0:-3] + "bam")
    return bamFileNames


def sortBamFiles(BamfileNames: list):
    sortedBamFileNames = []
    for i in range(len(BamfileNames)):
        command = "samtools sort {} -o {}".format(
            BamfileNames[i], BamfileNames[i][0:-3] + "sorted.bam"
        )
        os.system(command)
        sortedBamFileNames.append(BamfileNames[i][0:-3] + "sorted.bam")
    return sortedBamFileNames


def qualityControl():
    cwd = os.getcwd()
    os.system("mkdir qcReports")
    os.system("fastqc -t {} -o qcReports/ *.fastq".format(multiprocessing.cpu_count()))
    os.chdir("qcReports")
    os.system("multiqc .")
    # send(multiqc_report.html)     IMPLEMENT THIS


def createCountMatrix(pathToGTF, bamInputFile):
    """
    Exon specific expression and alternative splicing not implemented
    Uses featureCounts

    Reference - http://genomespot.blogspot.com/2015/01/generate-rna-seq-count-matrix-with.html
    """
    os.system(
        "featureCounts -Q 10 -T {} -a {} -o countMatrix.txt {}.bam".format(
            multiprocessing.cpu_count(), pathToGTF, bamInputFile
        )
    )


if __name__ == "__main__":
    print("Process ID: {}".format(os.getpid()))
    refgenome = ""

    # Example IDGSE138181
    search_id = "GSE29968"
    breakpoint = 2
    try:
        # Over-writting with user Entered ID, only if present
        if len(sys.argv) > 1:
            breakpoint = int(sys.argv[1])
            search_id = sys.argv[2]
    except:
        print("Exception Occured during argument Parsing")
        writeMessage(search_id,"Process ID: {}".format(os.getpid()))
        writeMessage(search_id,"Exception Occured during argument Parsing")


    print("Input GEOID is", search_id)
    print("Breaking step is {}".format(breakpoint))
    writeMessage(search_id,"Input GEOID is {}".format(search_id))


    Entrez.email = "bhavay18384@iiitd.ac.in"
    handle = Entrez.esearch(db="gds", term=search_id)
    pp = handle.read()
    # print(pp)
    tree = ET.fromstring(pp)
    print("Selecting the following ID")
    writeMessage(search_id,"Selecting the following ID")
    id = idselector(tree)
    writeMessage(search_id,str(id))
    handle = Entrez.esummary(db="gds", id=id)
    pp = handle.read()
    # print(pp)
    tree = ET.fromstring(pp)
    targetsra = ""
    for item in tree.findall("DocSum"):
        for item2 in item.findall("Item"):
            if item2.attrib["Name"] == "taxon":
                refgenome = item2.text
            if item2.attrib["Name"] == "ExtRelations":
                for item3 in item2.findall("Item"):
                    if item3.attrib["Name"] == "ExtRelation":
                        for item4 in item3.findall("Item"):
                            if item4.attrib["Name"] == "TargetObject":
                                targetsra = item4.text

    # Fetching Target SRA
    print("SRA in relation is ", targetsra)  # SRA to  be searched again
    writeMessage(search_id,"SRA in relation is {}".format(targetsra))
    #writeMessage()
    handle = Entrez.esearch(db="sra", term=targetsra)
    pp = handle.read()
    # print(pp)
    tree = ET.fromstring(pp)
    print("Selecting the following ID")
    writeMessage(search_id,"Selecting the following ID")
    id = idselector(tree)
    writeMessage(search_id,str(id))
    variable = Entrez.efetch(db="sra", id=id, rettype="full")
    # print(variable.read())

    # import os
    # os.system("pysradb metadata SRP250724")
    print("Organism is", refgenome)
    writeMessage(search_id,"Organism is {}".format(refgenome))

    # Downloading .SRA file
    db = sraweb.SRAweb()
    db.download(targetsra, skip_confirmation=True)
    # SRR11550936
    # Searching for Downloaded SRA File in File System
    os.system("find . -name " + targetsra + "> downloadPath.txt")
    downloadFileLocation = open("downloadPath.txt", "r").readline()
    #sraFileNames = ["Put List of All SRA File Names"]
    
    # Creating fastq file from Downloaded SRA File
    #firstList, secondList = createFastqFiles([])

    # To Download Refernce Genome and it's Annotations
    gtfName, fnaName = downloadRefGenome("GSE146443", "Caenorhabditis elegans")

    # To Preprocess the data (rRNA contamination removal, trimming)
    #preprocess(firstList, secondList)

    # Creating Index for Hisat2
    #index_name = createIndex(fnaName, gtfName)

    # Starting Alignment //Done
    #listOfOutputFiles = startAlignment(index_name, firstList, secondList)

    # Converting Sam to Bam // Done
    #BamFileNames = convertSamToBam(listOfOutputFiles)

    # Sort Bam Files // Done
    #sortedBamFileNames = sortBamFiles(BamFileNames)

    # Quality Control
    #qc = qualityControl()

    # Create Count Matrix (User Option)
    #createCountMatrix()
