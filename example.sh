# Everything plugged together.  Best uses as a template.
# Assumes one sample per input file, assumed components have not been
# installed properly...

HG19=/path/to/hg19
TMP=/var/tmp

# this starts a number of workers assuming you run SGE 
qsub -N workers -t 1-50 -pe smp 1- \
    network-aware-bwa/bwa worker -t \$NSLOTS -h `hostname` -p 6789

for f in "$@" ; do
    pipeline/mergeTrimReadsBAM --keepOrig -u -o /dev/stdout "${f}" | \
        network-aware-bwa/bwa bam2bam -g ${HG19} -p 6789 -t 0 - \
        -f ${f%bam}u.bam__
done

qdel workers

for f in "$@" ; do
    samtools sort -o ${f%bam}u.bam ${TMP}/${f} | \
        biohazard/dist/build/bam-rmdup/bam-rmdup --keep \
        -o ${f%bam}hg19.bam && rm ${f%bam}u.bam 
done

