#!/usr/bin/python
# python thisfile.py SSR123456

import os,sys

ID = sys.argv[1]
H1 = ID[0:3]
H2 = ID[0:6]
os.system('/home/k821209/.aspera/connect/bin/ascp  -i /home/k821209/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l80m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra ./'%(H1,H2,ID,ID))
