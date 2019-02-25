#!/usr/bin/env bash

#delete aliases that are the same as gene symbol or name
#
HOMEDIR=/home/rgddata/pipelines/EntrezGeneLoading
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

ELIST=mtutaj@mcw.edu

echo  "starting EntrezGene pipeline"
cd $HOMEDIR
java -Dspring.config=../properties/default_db2.xml \
    -Dlog4j.configuration=file://$HOMEDIR/properties/log4j.properties \
    -jar lib/EntrezGeneLoading.jar \
    -delete_redundant_aliases > deleted_redundant_aliases.log
mailx -s "[$SERVER] EntrezGene pipeline finished running" $ELIST < deleted_redundant_aliases.log

