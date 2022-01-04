#loads RefSeq nucleotide and protein sequences into RGD database
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
ELIST=mtutaj@mcw.edu

echo  "starting EntrezGene pipeline in RefSeqLoad mode"
cd $HOMEDIR
java -Dspring.config=../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$HOMEDIR/properties/log4j2.xml \
    -jar EntrezGeneLoading.jar -refseq_load -species rat | tee refseq_load.log
mailx -s "[$SERVER] REFSEQ Load EntrezGene pipeline OK" mtutaj@mcw.edu < refseq_load.log
