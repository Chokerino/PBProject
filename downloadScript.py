import os
import pandas
from pysradb.sraweb import SRAweb
import sys, getopt

class Download_fq_from_sra:
    
    def __init__(self, sra_id):
        self.db = SRAweb()
        self.sra_id = sra_id
    
    def get_metadata(self):
        df = self.db.sra_metadata(self.sra_id, detailed=True)
        return df
    
    def download_fq_file(self):
        print(os.getcwd())
        os.system('mkdir {}'.format(self.sra_id))
        metadata = self.get_metadata()
        print(metadata)
        os.chdir(self.sra_id)
        for run_acc in metadata.loc[:,"run_accession"]:
            print(run_acc)
            return_value = os.system("fasterq-dump {}".format(str(run_acc)))
            print(return_value)
        self.build_collections(metadata)

    def build_collections(self, df):
        all_types = set(list(df.loc[:,"source_name"]))
        for collec in all_types:
            os.system("mkdir {}_{}".format(collec, self.sra_id))
            for i, curr_collec in enumerate(df.loc[:,"source_name"]):
                if collec == curr_collec:
                    curr_run = str(df.loc[i,"run_accession"])
                    os.system('mv {}* {}_{}/'.format(curr_run, collec, self.sra_id))

 
if __name__ == "__main__":
    download_object = Download_fq_from_sra("SRP007169")
    download_object.get_metadata()
    download_object.download_fq_file()

