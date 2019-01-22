#!/usr/bin/env nextflow

/* Take input of sampleMap file (csv)
* run HLAminer inside Singularity, output HLAtypes tsv per sample
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '------------------------------------'
  log.info 'NEXTFLOW: RUN HLAminer ON EXOME DATA'
  log.info '------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run exome_HLA.simg.nf \
    --sampleCsv "send/location/sampleMap.csv" \
    -with-timeline exome_HLA.timeline \
    -with-report exome_HLA.report" \
    -c exome_HLA.simg.nextflow.config'

  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    sampleCsv     FILE      CSV format, header to include "sampleID,read1,read2" in case of paired data (NB single-end not implemented yet)'
  log.info ''
  exit 1
}

/* 1.0: HLAminer
*/
params.outDir = "analysis"

Channel.fromPath("$params.sampleCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.sampleID, file(row.read1), file(row.read2)] }
       .set { hlaminer }

process hlamin {

  publishDir path: "$params.outDir/$sampleID", mode: 'copy', pattern: "*[.csv,.log]"

  input:
  set val(sampleID), file(read1), file(read2) from hlaminer

  output:
  file("$sampleID*csv") into parsing

  """
  bwa aln -e 0 -o 0 /opt/database/HLA-I_II_CDS.fasta $read1 > aln_1.sai
  bwa aln -e 0 -o 0 /opt/database/HLA-I_II_CDS.fasta $read2 > aln_2.sai
  bwa sampe -o 1000 /opt/database/HLA-I_II_CDS.fasta aln_1.sai aln_2.sai $read1 $read2 > $sampleID".aln.sam"
  /opt/bin/HLAminer.pl -p /opt/database/hla_nom_p.txt -a $sampleID".aln.sam" -h /opt/database/HLA-I_II_CDS.fasta -s 500
  rm *sai
  samtools view -hC -T /opt/database/HLA-I_II_CDS.fasta $sampleID".aln.sam" > $sampleID".aln.cram"
  rm $sampleID".aln.sam"

  echo $read1 > fastq.files
  echo $read2 >> fastq.files

  /opt/bin/TASR -f fastq.files -m 20 -k 20 -s /opt/database/HLA-I_II_CDS.fasta -i 1 -b TASRhla -w 1
  cat TASRhla.contigs |perl -ne 'if(/size(\\d+)/){if(\$1>=200){\$flag=1;print;}else{\$flag=0;}}else{print if(\$flag);}' > TASRhla200.contigs
  /opt/bin/formatdb -p F -i TASRhla200.contigs
  echo "Aligning TASR contigs to HLA references..."
  /opt/bin/parseXMLblast.pl -c /opt/bin/ncbiBlastConfig.txt -d /opt/database/HLA-I_II_CDS.fasta -i TASRhla200.contigs -o 0 -a 1 > tig_vs_hla-ncbi.coord
  /opt/bin/parseXMLblast.pl -c /opt/bin/ncbiBlastConfig.txt -i /opt/database/HLA-I_II_CDS.fasta -d TASRhla200.contigs -o 0 > hla_vs_tig-ncbi.coord
  /opt/bin/HLAminer.pl -p /opt/database/hla_nom_p.txt -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla200.contigs -h /opt/database/HLA-I_II_CDS.fasta

  for x in *_H*csv; do
    mv \$x $sampleID"."\$x
  done
  """
}

/* 2.: Parse outputs
*/
process pars {

  publishDir "$params.outDir/$sampleID", mode:"copy", pattern: "*tsv"

  input:
  file(csv) from parsing

  output:
  file('*tsv') into completedPars

  shell:
  '''
  perl parse_HLAminer_output.pl !{csv}
  '''
}
