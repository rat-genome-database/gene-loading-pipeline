#on demand load of all mitochondrial genes
#
echo  "starting EntrezGene pipeline in mitochondrial mode"
cd /home/rgddata/pipelines/EntrezGeneLoading
java -Dspring.config=../properties/default_db.xml -jar EntrezGeneLoading.jar -mitochondrial -species rat > mtrat.log
mailx -s "[KIRWAN] Rat EntrezGene pipeline for mitochondrial genes finished running" rgd.developers@mcw.edu < mtrat.log
#
java -Dspring.config=../properties/default_db.xml -jar EntrezGeneLoading.jar -mitochondrial -species mouse > mtmouse.log
mailx -s "[KIRWAN] Mouse EntrezGene pipeline for mitochondrial genes finished running" rgd.developers@mcw.edu < mtmouse.log
#
java -Dspring.config=../properties/default_db.xml -jar EntrezGeneLoading.jar -mitochondrial -species human > mthuman.log
mailx -s "[KIRWAN] Human EntrezGene pipeline for mitochondrial genes finished running" rgd.developers@mcw.edu < mthuman.log

