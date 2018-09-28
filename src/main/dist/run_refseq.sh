#checks genes in rgd againt refseq files on NCBI FTP site
#
echo  "starting EntrezGene pipeline in RefSeq mode"
cd /home/rgddata/pipelines/EntrezGeneLoading
java -Dspring.config=../properties/default_db.xml -jar EntrezGeneLoading.jar -refseq -species rat > refseq.log
mailx -s "[KIRWAN] REFSEQ EntrezGene pipeline finished running" mtutaj@mcw.edu < refseq.log
