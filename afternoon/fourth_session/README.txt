#To work on this tutorial, you will need ms and seq-gen. These can be downloaded from:
# ms:
# http://home.uchicago.edu/rhudson1/source/mksamples.html
# seq-gen:
# http://tree.bio.ed.ac.uk/software/seqgen/

#creating simulated dataset to test migrate-n

#for a single pop and seven loci
#run ms
ms 20 7 -T | tail +4 | grep -v // >treefile.singlepop

#run seq-gen
seq-gen -mF84 -l 100 -s .2 -t 2 <treefile.singlepop >seqfile.singlepop

#run seqgen2mig.py
./migrate_sim.py seqfile.singlepop infile.singlepop "20"

#run migrate-n
migrate-n parmfile.singlepop

#for two pops, and seven loci
#run ms
ms 40 7 -T -I 2 20 20 4 | tail +4 | grep -v // >treefile.twopop

#run seq-gen
seq-gen -mF84 -l 100 -s .2 -t 2 <treefile.twopop >seqfile.twopop

#run seqgen2mig.py
./seqgen2mig.py seqfile.twopop infile.twopop "20 20"

#run migrate-n
migrate-n parmfile.twopop

#for two pops and 20 loci
#run ms
ms 40 20 -T -I 2 20 20 4 | tail +4 | grep -v // >treefile.twopop.20loci

#run seq-gen
seq-gen -mF84 -l 100 -s .2 -t 2 <treefile.twopop.20loci >seqfile.twopop.20loci

#run seqgen2mig.py
./seqgen2mig.py seqfile.twopop.20loci infile.twopop.20loci "20 20"

#run migrate-n
migrate-n parmfile.twopop.20loci