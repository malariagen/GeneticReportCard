#!/bin/sh

SUFFIX=crt

#5'
#Flank seq:     	TATTATTTATTTAAGTGTA
#Rev. Flank seq:	TACACTTAAATAAATAATA
#
#3'
#Flank seq:     	ATTTTTGCTAAAAGAAC
#Rev. Flank seq:	GTTCTTTTAGCAAAAAT
#
REGEX="TATTATTTATTTAAGTGTA\|TACACTTAAATAAATAATA\|ATTTTTGCTAAAAGAAC\|GTTCTTTTAGCAAAAAT"

#REGEX1=TATTATTTATTTAAGTGTA
#REGEX2=TACACTTAAATAAATAATA
#REGEX3=ATTTTTGCTAAAAGAAC
#REGEX4=GTTCTTTTAGCAAAAAT

samtools="/nfs/users/nfs_o/om1/samtools-1.2/samtools"
SAMPLE=$1
ALN=$2

echo $SAMPLE
# Get reads from BAM file for de novo alignment
$samtools view $ALN | grep $REGEX >  ./flankReads/$SAMPLE-$SUFFIX-flankReads.sam
#$samtools view $ALN | grep $REGEX1 >  ./flankReads/$SAMPLE-$SUFFIX-flankReads.sam
#$samtools view $ALN | grep $REGEX2 >> ./flankReads/$SAMPLE-$SUFFIX-flankReads.sam
#$samtools view $ALN | grep $REGEX3 >> ./flankReads/$SAMPLE-$SUFFIX-flankReads.sam
#$samtools view $ALN | grep $REGEX4 >> ./flankReads/$SAMPLE-$SUFFIX-flankReads.sam
