Emily Zhu | BIO312 | Stony Brook University | Dr. Joshua Rest
# Phylogenetic Analysis of the Proteasome Gene 
This repository is a tutorial to conduct phylogeneic analysis of the Proteasome Gene through methods such as gene alignment, gene tree reconciliation, and predicting protein domains. The bioinformatic software to conduct the methods to which they were installed into the AWS system by Associate Professor, Dr. Josuha Rest.

In this guide, you will be working with bioinformatic formats such as: GFF, FASTA, GenBank, Newick, Thirdkind, etc.

## Contents

a: Introduction + Getting Started

1: Determining Homologous Sequences 

2: Gene Alignment

3: Gene Tree Reconciliation

4: Protein Domain Predictions


# Introduction + Getting Started
In this analysis, we will be using various bioinformatics software to conduct analysis on gene families. Bioinformatic packages should be downloaded to the system prior the conduction of the lab.

Here's a guide on [installing packages](https://docs.github.com/en/packages/learn-github-packages/installing-a-package)

In this analaysis, we will work with 9 animals, 5 protostomes and 4 duterostomes: *Homo sapiens, Branchiostoma belcheri, Mizuhopecten yessoensis, Acanthaster planci, Ciona intestinalis, Drosophila melanogaster, Adineta vaga, Echinococcus granulosus, Lingula anatina*

Between every section, push the commands to the repository on github to view the files and save your work. Learn how to do so [here](https://docs.github.com/en/get-started/importing-your-projects-to-github/importing-source-code-to-github/adding-locally-hosted-code-to-github). 
## Getting started
In this lab, $MYGIT is used in lieu of of the github username. To create this shortcut, add the variable $MYGIT to your .bash_profile. 
Here's a [guide](https://stackoverflow.com/questions/14524590/how-to-add-export-statement-in-a-bash-profile-file) on how to do that. 

If you do *not* want to use this shortcut, simply replace $MYGIT with your github username. 

REMEMBER: Between every section, **push** the commands to the repository on github to view the files and save your work. Learn how to do so [here](https://docs.github.com/en/get-started/importing-your-projects-to-github/importing-source-code-to-github/adding-locally-hosted-code-to-github). 

# 01 Determining Homologous Sequences
Use BLAST to find orthologs and paralogs to the Proteasome Subunit Alpha Type-3-like gene in bilaterians. This is done through working with proteoms. 
The folder containing all the data from this section is named **01-psma-homologs**.

## Begin by creating a BLAST database

We will be searching for homologs in proteomes from several species. 

Begin in the repository.
```bash
cd ~/labs/projectrepo-$MYGIT
```
Uncompress the proteomes using the following command:
```bash
gunzip proteomes/*.gz
```
Using the cat command, move all protein sequences into a single file
```bash
cat  proteomes/* > allprotein.fas
```
Now with these sequences, build a blast database.
```bash
makeblastdb -in allprotein.fas -dbtype prot
```
## Conduct the BLAST search: BLAST a Proteasome Subunit Alpha Type-3-like against the database to identify homologs. 
Create a folder for the Proteasome BLAST search using the`mkdir` command:
```bash
mkdir /home/ec2-user/labs/projectrepo-$MYGIT/01-psma-homologs
```
Go to the folder: 
```bash
cd ~/labs/projectrepo-$MYGIT/01-psma-homologs
```
Download the query protein. 
```bash
ncbi-acc-download -F fasta -m protein XP_021370576.1
```
Perform the BLAST search. 
```bash
blastp -db ../allprotein.fas -query XP_021370576.1.fa -outfmt 0 -max_hsps 1 -out psma.blastp.typical.out
```
Analyze the output using the `less` command.

## Perform a BLAST search with a tabular output
Create a more detailed output of the same analysis using the following command:
```bash
blastp -db ../allprotein.fas -query XP_021370576.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out psma.blastp.detail.out
```
Look at the output in ``psma.blastp.detail.out`` using the ``less -S`` command.
To determine the number of hits for each species, use the `grep -c` command folllowed by the abbreviated species name followed by psma.blastp.detail.out 
For example: to determine the number of total human hits in the file use the command 
```bash
grep -c Hsapiens psma.blastp.detail.out
```
## Filtering the BLAST output.
This is done to identify only high-scoring putative homologs.
Using the e-value cutoff of <1e-35 the resulting output will contain homologs that have a relatively high sequence identity and a high alignment score as well as minimize the possibility that we include false homologs. 

The e-value cut off used is 1e-35.

Use this command to filter out the output file with our e-value cutoff:
```bash
awk '{if ($6< 1e-35 )print $1 }' psma.blastp.detail.out > psma.blastp.detail.filtered.out
```

Now, determine the total number of hits in the BLAST result using the `wc` command:
```bash
wc -l psma.blastp.detail.filtered.out
```
To determine the number of paralogs in each species
Use the`grep` command.
```bash
grep -o -E "^[A-Z][a-z]+\." psma.blastp.detail.filtered.out  | sort | uniq -c
```
Push the repository. 

# 02 Gene Alignment
Align the gene family using MUSCLE and Seqkit. 
The folder containing all the data from this section is named **02-psma-alignment**.
## Install software
The following utility software was installed in order to run the gene alignment:
[[aha](https://github.com/theZiz/aha)] convert text-based sequence alignment into html.
[[conda](https://docs.conda.io/en/latest/)] an open source package and enviornment management system to install aha.
[[a2ps](https://www.gnu.org/software/a2ps/)] to convert the html file to a poscript file.
[[yum](https://blog.packagecloud.io/what-is-yum-package-manager/)] a package manager installed to Amazon Linux 2 OS to download a2ps.
[[alignbuddy](https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy)] from [[buddysuite](https://github.com/biologyguy/BuddySuite)] to work with alignments
[[pip](https://pypi.org/project/pip/)] a Python package manager will be used to install alignbuddy. 

This is the command to install the software, only run this once.
```bash
conda install -y -n base -c conda-forge aha
sudo yum install -y a2ps
pip install buddysuite
```
## Begin the gene family alignment
Begin in the repository. 
```bash
cd ~/labs/projectrepo-$MYGIT
```
Create a folder for the Proteasome alignment using the `mkdir` command:
```bash
mkdir /home/ec2-user/labs/projectrepo-$MYGIT/02-psma-alignment
```
Go to the folder: 
```bash
cd ~/labs/projectrepo-$MYGIT/02-psma-alignment
```
Use the `seqkit` command to obtain the sequences in the BLAST output file created in the previous section:
```bash
seqkit grep --pattern-file ~/labs/projectrepo-$MYGIT/01-psma-homologs/psma.blastp.detail.filtered.out ~/labs/projectrepo-$MYGIT/01-psma-homologs/allprotein.fas > ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.fas
```
## Using the Muscle program, perform a global multiple sequence alignment
This will align all the entire length of the sequences with one another. 

Make a multiple sequence alignments using [MUSCLE](https://drive5.com/muscle/manual/):
```bash
muscle -in ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.fas -out ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.fas
```
View the alignment in alv, use arrow keys to move around the alignment:
```bash
alv -kli  ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.fas | less -RS
```
Use the "majority"  in alv to hilight where the most common amino acids were found in 50% of sequences.
```bash
alv -kli --majority ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.fas | less -RS
```
Convert the alignment into a pdf to be viewed in github.
```bash
alv -ki -w 100 ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.fas | aha > ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.html
a2ps -r --columns=1 ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.html -o ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.ps
ps2pdf ~/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.ps ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.pdf
```
## Alignment informatinon
Determine the length and average precent identity.

Calculate the length of the alignment:
```bash
alignbuddy  -al  ~/labs/projectrepo-$MYGIT/psma-alignment/psma.homologs.al.fas
```
Remove columns with gaps and calculate the length of the alignment:
```bash
alignbuddy -trm all  ~/labs/projectrepo-$MYGIT/psma-alignment/psma.homologs.al.fas | alignbuddy  -al
```
Remove the completely conserved portions and calculate the length of the alignment:
```bash
alignbuddy -dinv 'ambig' ~/labs/projectrepo-$MYGIT/psma-alignment/psma.homologs.al.fas | alignbuddy  -al
```
Use [[t_coffee](https://www.tcoffee.org/Projects/tcoffee/)] to calculate the average percent identity 
```bash
t_coffee -other_pg seq_reformat -in ~/labs/projectrepo-$MYGIT/psma-alignment/psma.homologs.al.fas -output sim
```
Use [[alignbuddy](https://github.com/biologyguy/BuddySuite/blob/master/buddysuite/AlignBuddy.py)] to calculate the average percent identity 
```bash
 alignbuddy -pi ~/labs/projectrepo-$MYGIT/psma-alignment/psma.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
# 03 Gene Tree Reconciliation 
The folders containing all the data from this section are named **03-psma-tree** and **04-psma-tree-reconciled**. 

Create a folder for the Proteasome gene tree using the `mkdir` command:
```bash
mkdir /home/ec2-user/labs/projectrepo-$MYGIT/03-psma-tree
```
Go to the folder: 
```bash
cd ~/labs/projectrepo-$MYGIT/03-psma-tree
```
Download psma.homologs.al.fas from the previous section.
```bash
cp ~/labs/projectrepo-$MYGIT/02-psma-alignment/psma.homologs.al.fas ~/labs/projectrepo-$MYGIT/03-psma-tree/psma.homologs.al.fas  
```
## Construct a Phylogenetic Tree for Opsin Homologs from Sequence Data
Use the software [IQ-TREE] to infer the optimal phylogenetic tree based on a sequence alignment using the data in the previous section. 
Copy the alignment file from the previous section; IQTree will place the output files in the same directories as the input file. This will run for ~10 minutes.
    If it gets stopped, you can restart it by runing the same command again from the same place.
```bash
iqtree -s ~/labs/projectrepo-$MYGIT/03-psma-tree/psma/psma.homologs.al.fas -bb 1000 -nt 2
```
Look at the Substitution model, use up down arrows to scroll. 
```bash
nano ~/labs/projectrepo-$MYGIT/03-psma-tree/psma/psma.homologs.al.fas.iqtree
```
## Look at the unrooted tree
Display the newick formatted .treefile:
```bash
nw_display ~/labs/03-psma-tree/psma/psma.homologs.al.fas.treefile
```
Look at the unrooted files using an R script.
```bash
Rscript --vanilla ~/labs/03-psma-tree/plotUnrooted.R  ~/labs/03-psma-tree/psma/psma.homologs.al.fas.treefile ~/labs/03-psma-tree/psma/psma.homologs.al.fas.treefile.pdf 0.4
```
## Midpoint rooting
Use [gotree](https://github.com/evolbioinfo/gotree) to reroot the tree.
```bash
gotree reroot midpoint -i ~/labs/03-psma-tree/psma/psma.homologs.al.fas.treefile -o ~/labs/03-psma-tree/psma/psma.homologs.al.mid.treefile
```
Look at the midpoint rooted tree by displaying the entire ASCII Image:
```bash
nw_order -c n ~/labs/03-psma-tree/psma/psma.homologs.al.mid.treefile  | nw_display -
```
Create an output as a graphic: 
nw_order -c n ~/labs/03-psma-tree/psma/psma.homologs.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/labs/03-psma-tree/psma/psma.homologs.al.mid.treefile.svg -

#### Branch lengths
Switch the view to a **cladogram**: 
```bash
nw_order -c n ~/labs/03-psma-tree/psma/psma.homologs.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/labs/03-psma-tree/psma/psma.homologs.al.midCl.treefile.svg -
```
Push files and examine the proteasome gene family!

## Reconciling the gene tree within the species tree.

Begin in the repository 
cd ~/labs/04-psma-tree-reconciled/psma

Download the software package [Notung](http://www.cs.cmu.edu/~durand/Notung/)
```bash
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar --help  
```
Make a copy of the midpoint rooted gene tree from the previous section into this folder
```bash
cp ~/labs/03-psma-tree/psma/psma.homologs.al.mid.treefile ~/labs/04-psma-tree-reconciled/psma/psma.homologs.al.mid.treefile
```
Perform the reconcillation
```bash
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ~/labs/03-psma-tree/species.tre -g ~/labs/04-psma-tree-reconciled/psma/psma.homologs.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/labs/04-psma-tree-reconciled/psma/
```
Examine the file psma.homologs.al.mid.treefile.reconciled.events.txt
using`less -S`.

Look at the species tree
nw_display ~/labs/03-psma-tree/species.tre

To determine the lineages, use:
grep NOTUNG-SPECIES-TREE ~/labs/04-psma-tree-reconciled/psma/psma.homologs.al.mid.treefile.reconciled | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -

Push the files to look at the graphics.

View the gene-within species tree using thirdkind.
Begin by creating a RecPhyloXML object.
```bash
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/labs/04-psma-tree-reconciled/psma/psma.homologs.al.mid.treefile.reconciled --include.species
```
Create a gene reconcillation within species tree reconcillation using thirdkind.
```bash
thirdkind -Iie -D 40 -f ~/labs/04-psma-tree-reconciled/psma/psma.homologs.al.mid.treefile.reconciled.xml -o  ~/labs/04-psma-tree-reconciled/psma/psma.homologs.al.mid.treefile.reconciled.svg
```
Push to view the reconcillations

# 04 Protein Domain Predictions
The folder containing all the data from this section is named **05-psma-domain-prediction**.

Create a folder for the Proteasome gene protein domain using the `mkdir` command:
```bash
mkdir /home/ec2-user/labs/projectrepo-$MYGIT/05-psma-domain-prediction/psma
```
Go to the folder: 
```bash
cd ~/labs/projectrepo-$MYGIT/05-psma-domain-prediction/psma
```
Begin with making a copy of the raw unaligned seqeunce from the previous section.
```bash
sed 's/*//' ~/labs/02-psma-alignment/psma/psma.homologs.fas > ~/labs/05-psma-domain-prediction/psma/psma.homologs.fas
```
Download the [Pfam](https://www.ebi.ac.uk/interpro/) database.
```bash
wget -O ~/data/Pfam_LE.tar.gz ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz && tar xfvz ~/data/Pfam_LE.tar.gz  -C ~/data
```
Run RPS-BLAST
```bash
rpsblast -query ~/labs/05-psma-domain-prediction/psma/psma.homologs.fas -db ~/data/Pfam -out ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
Plot the predicted protein domains onto the phylogeny.
```bash
sudo /usr/local/bin/Rscript  --vanilla ~/labs/05-psma-domain-prediction/plotTreeAndDomains.r ~/labs/05-psma-domain-prediction/psma/psma.homologs.al.mid.treefile ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out ~/labs/05-psma-domain-prediction/psma/psma.tree.rps.pdf
```
Look at the annotations in psma.rps-blast.out
```bash 
mlr --inidx --ifs "\t" --opprint  cat ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out | tail -n +2 | less -S
```
Determine the most common Pfam domain using this command:
```bash
cut -f 1 ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out | sort | uniq -c
```
Determine the longest annotated Pfam domain using this command:
```bash
cut -f 6 ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out | sort | uniq -c
```
Dwtermine the shortest annotated protein domain?
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out |  sort  -k2nr
```
Determine the protein domain with the best e-value using this command:
```bash 
sort  -k5rg ~/labs/05-psma-domain-prediction/psma/psma.rps-blast.out | less -S
```
[Push]([here](https://docs.github.com/en/get-started/importing-your-projects-to-github/importing-source-code-to-github/adding-locally-hosted-code-to-github)) everything to the repository. 
