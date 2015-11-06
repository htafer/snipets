#!/bin/bash                                                                                                                                                                                                                                                                      

function usage
{
    echo "usage: simplify the preparation of the tbl ";
    echo "files used to generate the sqn file necessary for the ncbi submission\n";
}

while [ "$1" != "" ]; do
    case $1 in
        -s | --swissProtBLAST ) shift
                                swissProtBLAST=$1
                                ;;   
        -f | --functional)      shift
                                functional=$1
                                ;;   
        -S | --swissProtFASTA)  shift
                                swissProtFASTA=$1
                                ;;   
        -N | --names)           shift
                                names=$1
                                ;;   
        -F | --genomeSeq)       shift
                                genomeSeq=$1
                                ;;   
        -G | --genomeAnn)       shift
                                genomeAnn=$1
                                ;;   
        -a | --annPrefix)       shift
                                annPrefix=$1
                                ;;   
        -h | --help )           usage
                                exit
                                ;;   
        * )                     usage
                                exit 1
    esac
    shift
done


parseNameFromIPRSP.pl -i $functional -s $swissProtFASTA -b $swissProtBLAST  > $names
tblFile=`echo $genomeSeq | sed -r 's/\.f.+/.tbl/'`;
echo "cat $genomeAnn | enumerateCDS.pl | testFeatures2.pl -F $genomeSeq -g /dev/stdin -a $annPrefix -f $functional -n $names > $tblFile";
cat $genomeAnn | enumerateCDS.pl | testFeatures2.pl -F $genomeSeq -g /dev/stdin -a $annPrefix -f $functional -n $names > $tblFile

printf "Now you are ready to pass the .tbl and .fsa file to linux64.tbl2asn\n"
printf "linux64.tbl2asn -t template.sbt -p . -a s -V v -j \"[organism=Debaryomyces fabryi][strain=CBS789]\" -i $genomeSeq"







