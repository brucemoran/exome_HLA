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
       .into { hpraminer; hptasrminer }

process hpra {

  publishDir path: "$params.outDir/HLAminer/$sampleID", mode: 'copy', pattern: "*[.csv,.log]"

  input:
  set val(sampleID), file(read1), file(read2) from hpraminer

  output:
  set val(sampleID), file("$sampleID*csv") into parsehpra

  """
  bwa aln -e 0 -o 0 /opt/database/HLA-I_II_CDS.fasta $read1 > aln_1.sai
  bwa aln -e 0 -o 0 /opt/database/HLA-I_II_CDS.fasta $read2 > aln_2.sai
  bwa sampe -o 1000 /opt/database/HLA-I_II_CDS.fasta aln_1.sai aln_2.sai $read1 $read2 > $sampleID".aln.sam"
  /opt/bin/HLAminer.pl -p /opt/database/hla_nom_p.txt -a $sampleID".aln.sam" -h /opt/database/HLA-I_II_CDS.fasta -s 500
  rm *sai
  samtools view -hC -T /opt/database/HLA-I_II_CDS.fasta $sampleID".aln.sam" > $sampleID".aln.cram"
  rm $sampleID".aln.sam"

  mv HLAminer_HPRA.csv $sampleID".HLAminer_HPRA.csv"
  mv HLAminer_HPRA.log $sampleID".HLAminer_HPRA.log"
  """
}

process hptasr {

  publishDir path: "$params.outDir/HLAminer/$sampleID", mode: 'copy', pattern: "*[.csv,.log]"

  input:
  set val(sampleID), file(read1), file(read2) from hptasrminer

  output:
  set val(sampleID), file("$sampleID*csv") into parsehptasr

  """
  echo $read1 > fastq.files
  echo $read2 >> fastq.files

  /opt/bin/TASR -f fastq.files -m 20 -k 20 -s /opt/database/HLA-I_II_CDS.fasta -i 1 -b TASRhla -w 1
  cat TASRhla.contigs |perl -ne 'if(/size(\\d+)/){if(\$1>=200){\$flag=1;print;}else{\$flag=0;}}else{print if(\$flag);}' > TASRhla200.contigs
  /opt/bin/formatdb -p F -i TASRhla200.contigs
  echo "Aligning TASR contigs to HLA references..."
  /opt/bin/parseXMLblast.pl -c /opt/bin/ncbiBlastConfig.txt -d /opt/database/HLA-I_II_CDS.fasta -i TASRhla200.contigs -o 0 -a 1 > tig_vs_hla-ncbi.coord
  /opt/bin/parseXMLblast.pl -c /opt/bin/ncbiBlastConfig.txt -i /opt/database/HLA-I_II_CDS.fasta -d TASRhla200.contigs -o 0 > hla_vs_tig-ncbi.coord
  /opt/bin/HLAminer.pl -p /opt/database/hla_nom_p.txt -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla200.contigs -h /opt/database/HLA-I_II_CDS.fasta

  mv HLAminer_HPTASR.csv $sampleID".HLAminer_HPTASR.csv"
  mv HLAminer_HPTASR.log $sampleID".HLAminer_HPTASR.log"
  """
}

/* 2.: Parse outputs
*/
parsehpra.join(parsehptasr).set { parseboth }
process pars {

  publishDir "$params.outDir/HLAminer/$sampleID", mode:"copy", pattern: "*tsv"

  input:
  set val(sampleID), file(csv1), file(csv2) from parseboth

  output:
  file('*tsv') into rprocess

  shell:
  '''
  #! /usr/bin/perl

  use strict;
  use warnings;

  ##parser for HLAminer to write to standardised format
  ##parses file name for sampleID (splits on first "." keeping [0])
  ##format outputs:
  #HLA-XYZ XYZ*01:001 for all precitions
  ##strategy: if previous line started ^HLA, \
  #and this line contains "Prediction", \
  #then until next line !~ "\\*", \
  #split on "," and return [0]=~s/\\s+//;

  ##file naming
  my @inputs = ("!{csv1}", "!{csv2}");
  foreach my $s (@inputs){
    my $fileinp = $s;
    my @outname = split(/\\./, $fileinp);
    my $fileout = "!{sampleID}" . "." . $outname[1] . ".tsv";

    ##iterate over file
    open(IN, $fileinp);

    ##to write output
    open(OUT, ">$fileout");
    print OUT "SAMPLEID\\tHLAminer_Call\\n";

    ##flags to determine course of action
    my $HLAline=0;
    my $PRDline=0;
    my $KEPline=0;

    while(<IN>){
      if($_=~m/^HLA/){
        $HLAline=1;
        next;
      }
      if($_=~m/Prediction/){
        $PRDline=1;
        next;
      }
      ##we have a line of interest; split and print OUT
      if(($PRDline==1) && ($_ ne "")){
        if($_=~m/:/){
          my @sp=split(/\\,/,$_);
          my $HLAcall=$sp[0];
          $HLAcall=~s/\\s+//g;
          print OUT "!{sampleID}" . "\\t" . $HLAcall . "\\n";
          next;
        }
      }
    }

    close IN;
    close OUT;
  }
  '''
}

/* 3.: Concatenate into R object; indicate frequency of each set
*/

process rfreq {

  publishDir "$params.outDir/R", mode:"copy"

  input:
  file(mixr) from rprocess.collect()

  output:
  file('*RData') into completeR

  shell:
  '''
  #! /usr/bin/Rscript

  ##script to read in all tsv for HPRA, HPTASR tsv files
  ##combines, then creates table of frequency of all HLAminer calls

  library(tidyverse)

  ##iterate over HPRA, HPTASR inputs
  lapply(c("HPRA", "HPTASR"), function(HPTYPE){

    ##create list to hold output in single object
    all_list <- as.list(1,2,3)

    ##read in all matches to HPTYPE.tsv
    all_list[[1]] <- map_df(list.files(pattern = paste0(HPTYPE,".tsv"),
                                full.names = TRUE),
                      read_tsv)
    ##tabulate
    all_table <- table(all_list[[1]]$HLAminer_Call) %>% sort()
    ##frequency and table in tibble
    all_list[[2]] <- tibble(HLAminer_Call = names(all_table),
                            count = c(all_table),
                            frequency = c(all_table/length(unique(all_list[[1]]$SAMPLEID)))) %>%
                     arrange(., desc(frequency))

    ##rename list
    names(all_list) <- paste(HPTYPE, c("tsv", "table"), sep="_")

    ##assign name to object for HPTYPE
    assignName <- paste0(HPTYPE, "_list")
    assign(assignName, value = all_list)

    ##save output
    save(list = assignName, file = paste0(HPTYPE, ".table.RData"))
  })
  '''
}
completeR.subscribe { println "Completed R processing: " + it }
