#! /bin/bash

############
## SETUP ##
##########

##inputs
BASEDIR=$1
SAMPMAP=$2

##ensure specified
if [[ $# < 2 ]];then
  echo -e "Two inputs required: \n\t(i) full path to base directory (created if not extant)\n\t(ii) full path to sample map file (tab-delimited, sampleID -> fastq location (comma-sep if paired-end))"
fi

##check dirs, samplemap
if [[ ! -d $BASEDIR ]];then
  echo "$BASEDIR is not extant, creating..."
  mkdir -p $BASEDIR
fi

cd $BASEDIR
singularity pull shub://HugoMananet/HLAminer:hlaminer.1.4.def

if [[ ! -e $SAMPMAP ]];then
  echo "No sample map file found at $SAMPMAP"
  exit 127;
fi

#################
## RUN FASTQS ##
###############

##check inputs exist; create script and output directory for sample
DATE=$(date +"%Y%M%d")
cat $SAMPMAP | while read LINE; do
  SAMPLEID=$(echo $LINE | cut -f 1)
  FASTQ=$(echo $LINE | cut -f 2)

  ##dir to hold output
  EXECDIR="$BASEDIR/analysis/$SAMPLEID"
  if [[ ! -d "$EXECDIR" ]];then
    mkdir -p "$EXECDIR"
  fi

  cd $EXECDIR

  ##test paired-end
  if [[ "$FASTQ" =~ "," ]];then
    echo "Found paired-end data..."
    FASTQ1=$(echo $FASTQ | cut -d "," -f 1)
    FASTQ2=$(echo $FASTQ | cut -d "," -f 2)
    if [[ -e $FASTQ1 && -e $FASTQ2 ]]; then
      sed "s/FASTQ1/$FASTQ1/g" template_HLAminer_PE_exec.sh | \
      sed "s/FASTQ2/$FASTQ2/g" | \
      sed "s/SAMPLEID/$SAMPLEID/g" | \
      sed "s/EXECDIR/$EXECDIR/g" > "$SAMPLEID.HLAminer_PE_exec.$DATE.sh"
      singularity exec $BASEDIR/HugoMananet-HLAminer-master-hlaminer.1.4.def.simg sh "$SAMPLEID.HLAminer_PE_exec.$DATE.sh"
    else
      echo "Could not find $FASTQ1 and/or $FASTQ2, ensure paths are correct in $SAMPMAP"
    fi
  else
    if [[ -e $FASTQ ]]; then
      echo "Found single-end data..."
      sed "s/FASTQ/$FASTQ/g" template_HLAminer_PE_exec.sh | \
      sed "s/SAMPLEID/$SAMPLEID/g" | \
      sed "s/EXECDIR/$EXECDIR/g" > "$SAMPLEID.HLAminer_SE_exec.$DATE.sh"
      $BASEDIR/HugoMananet-HLAminer-master-hlaminer.1.4.def.simg sh "$SAMPLEID.HLAminer_PE_exec.$DATE.sh"
    else
      echo "Could not find $FASTQ..."
    fi
  fi

  ####################
  ## PARSE OUTPUTS ##
  ##################

  ##in the spirit of HLAminer (Perl)...
  for x in *_HP*.csv;do
    mv $x $SAMPLEID.$x
    perl $BASEDIR/parse_HLminer_output.pl $SAMPLEID.$x
  done
done
