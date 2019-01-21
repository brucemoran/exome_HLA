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
if [[ ! -e $BASEDIR/HugoMananet-HLAminer-master-hlaminer.1.4.def.simg ]];then
  singularity pull shub://HugoMananet/HLAminer:hlaminer.1.4.def
fi

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
  SAMPLEID=$(echo $LINE | perl -ane 'print $F[0];')
  FASTQ=$(echo $LINE | perl -ane 'print $F[1];')

  ##dir to hold output
  EXECDIR="$BASEDIR/analysis/$SAMPLEID"
  echo "Creating $EXECDIR..."
  if [[ ! -d "$EXECDIR" ]];then
    mkdir -p "$EXECDIR"
  fi

  cd $EXECDIR

  ##test paired-end
  if [[ "$FASTQ" =~ "," ]];then
    echo "Found paired-end data..."
    FASTQ1=$(echo $FASTQ | perl -ane '@s=split(/\,/,$_);print $s[0];')
    FASTQ2=$(echo $FASTQ | perl -ane '@s=split(/\,/,$_);print $s[1];')
    if [[ -e $FASTQ1 && -e $FASTQ2 ]]; then
      cat template_HLAminer_PE_exec.sh | sed "s/FASTQ1/$FASTQ1/g" | sed "s/FASTQ2/$FASTQ2/g" | sed "s/SAMPLEID/$SAMPLEID/g" | sed "s/EXECDIR/$EXECDIR/g" > "$EXECDIR/$SAMPLEID.HLAminer_PE_exec.$DATE.sh"
      #singularity exec $BASEDIR/HugoMananet-HLAminer-master-hlaminer.1.4.def.simg sh "$SAMPLEID.HLAminer_PE_exec.$DATE.sh"
    else
      echo "Could not find $FASTQ1 and/or $FASTQ2, ensure paths are correct in $SAMPMAP"
    fi
  else
    if [[ -e $FASTQ ]]; then
      echo "Found single-end data..."
      cat template_HLAminer_PE_exec.sh | sed "s/FASTQ/$FASTQ/g" |    sed "s/SAMPLEID/$SAMPLEID/g" | sed "s/EXECDIR/$EXECDIR/g" > "$EXECDIR/$SAMPLEID.HLAminer_SE_exec.$DATE.sh"
      #singularity exec $BASEDIR/HugoMananet-HLAminer-master-hlaminer.1.4.def.simg sh "$SAMPLEID.HLAminer_PE_exec.$DATE.sh"
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
