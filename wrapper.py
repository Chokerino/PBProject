import subprocess


process = subprocess.Popen(("fasterq-dump {} -p -t $HOME/temp_files".format("SRR278178")).split(), stdout=subprocess.PIPE)
process.wait()