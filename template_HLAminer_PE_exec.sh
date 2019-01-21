#! /bin/bash

##from HugoMananet-HLAminer-master-hlaminer.1.4.def.simg /opt/bin
cd EXECDIR
{
##HPRArnaseq_classI-II.sh
echo -e"$(date)\nRunning HPRArnaseq_classI-II.sh..."
### Run bwa or your favorite short read aligner
echo "Running bwa..."
bwa aln -e 0 -o 0 /opt/database/HLA-I_II_CDS.fasta FASTQ1 > aln_1.sai
bwa aln -e 0 -o 0 /opt/database/HLA-I_II_CDS.fasta FASTQ2 > aln_2.sai
bwa sampe -o 1000 /opt/database/HLA-I_II_CDS.fasta aln_1.sai aln_2.sai FASTQ1 FASTQ2 > SAMPLEID.aln.sam
### Predict HLA
echo "Predicting HLA..."
/opt/bin/HLAminer.pl -p /opt/database/hla_nom_p.txt -a SAMPLEID.aln.sam -h /opt/database/HLA-I_II_CDS.fasta -s 500
##clean up
rm *sai
samtools view -hC -T /opt/database/HLA-I_II_CDS.fasta SAMPLEID.aln.sam

##HPTASRrnaseq_classI-II.sh
echo -e "$(date)\nRunning HPTASRrnaseq_classI-II.sh..."
###Run TASR
echo "Running TASR..."
#TASR Default is -k 15 for recruiting reads. You may increase k, as long as k < L/2 where L is the minimum shotgun read length
echo FASTQ1 > fastq.files
echo FASTQ2 >> fastq.files

/opt/bin/TASR -f fastq.files -m 20 -k 20 -s /opt/database/HLA-I_II_CDS.fasta -i 1 -b TASRhla -w 1
###Restrict 200nt+ contigs
cat TASRhla.contigs |perl -ne 'if(/size(\d+)/){if($1>=200){$flag=1;print;}else{$flag=0;}}else{print if($flag);}' > TASRhla200.contigs
###Create a [NCBI] blastable database
echo "Formatting blastable database..."
/opt/bin/formatdb -p F -i TASRhla200.contigs
###Align HLA contigs to references
echo "Aligning TASR contigs to HLA references..."
/opt/bin/parseXMLblast.pl -c /opt/bin/ncbiBlastConfig.txt -d /opt/database/HLA-I_II_CDS.fasta -i TASRhla200.contigs -o 0 -a 1 > tig_vs_hla-ncbi.coord
###Align HLA references to contigs
echo "Aligning HLA references to TASR contigs (go have a coffee, it may take a while)..."
/opt/bin/parseXMLblast.pl -c /opt/bin/ncbiBlastConfig.txt -i /opt/database/HLA-I_II_CDS.fasta -d TASRhla200.contigs -o 0 > hla_vs_tig-ncbi.coord
###Predict HLA alleles
echo "Predicting HLA alleles..."
/opt/bin/HLAminer.pl -p /opt/database/hla_nom_p.txt -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla200.contigs -h /opt/database/HLA-I_II_CDS.fasta
} 2>&1 | tee SAMPLEID.HPRA-HPTA.$(date +"%Y%M%d").log.txt
