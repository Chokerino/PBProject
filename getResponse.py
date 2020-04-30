import urllib3

req = urllib3.Request('ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL24nnn/GPL24676/soft/GPL24676_family.soft.gz')
response = urllib3.urlopen(req)
the_page = response.read()
