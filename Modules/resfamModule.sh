#!/bin/bash 
#SBATCH -J resfamModule
#SBATCH -o resfamModule.out
#SBATCH -e resfamModule.err
#SBATCH -t RESFAMTIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=RESFAMPARTITION
#SBATCH -N 1

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

#part1 #is necesary create identity value from ANI output(?)

source CHAINREACTION

CORES=$(grep -c ^processor /proc/cpuinfo)
THEFAA=FAAPATH
	

if ! [ -d $THEFAA ];then
	echo "* $THEFAA folder doesn't exist"
	exit
fi


gnum=$(ls -1 $THEFAA/* |wc -l |awk '{print $1}')
ls -1 $THEFAA/* > faasamples.txt

echo "#!/bin/bash
#SBATCH -J resfam
#SBATCH -o resfam_%a.out
#SBATCH -e resfam_%a.err
#SBATCH -t RESFAMTIME
#SBATCH --partition=RESFAMPARTITION
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --array=1-$gnum

i=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p faasamples.txt)
name=\$(echo \"\$i\" |rev |cut -d \"/\" -f 1 |rev)

HMMERHOME/bin/hmmscan -o \$name.hmmoutput --tblout \$name.hmmtable --cut_ga --cpu $CORES RESFAMHMMINDEX \$i

" > resfam.sh

callSLURM resfam.sh "Resfam"

rm -f faasamples.txt

if mkdir AR > /dev/null 2>&1; then
	cd AR
else
	cd AR
	rm -rf *
fi

mv ../*.hmm[ot]* .
mv ../resfam*.* .

#get unique Resfam ID
for i in $(ls -1 *.hmmtable)
do
	grep -v "#" $i |awk '{print $2}' |sort |uniq -c |awk '{print $2, $1}' > $i.unique.txt
done

ls -1 *.unique.txt > thetables.txt
cat *unique.txt |awk '{print $1}' |sort |uniq > common.txt

#Make ID Matrix
echo "RFID" > header
cat header common.txt > permanente
rm header
while read thetable
do
    newname=$(echo "$thetable" | rev |cut -d "/" -f 1 |rev |awk -F"." '{print $1}')
    touch permanente
    awk -v name="$newname" 'BEGIN{print name;n["default_var"]=0}{if(FNR==NR){n[$1]=$2}else{if($1 in n){print n[$1]}else{print "0"}}}' $thetable common.txt > tmpactual
    paste permanente tmpactual > tmp
    rm permanente
    mv tmp permanente

done < thetables.txt
mv permanente arMatrix.txt
rm -f tmpactual common.txt thetables.txt *unique.txt

#conver ID to Description

grep -v "#" *.hmmtable > allhmmtables

awk '{printf "%s\t",$2;for(i=19;i<NF;i++){printf "%s ",$i}print $NF}' allhmmtables > id2desc.tsv

sed -i "s/\//:/g" id2desc.tsv 
cat id2desc.tsv |while read line
do 
	id=$(echo $line |awk '{print $1}')
	desc=$(echo "$line" |awk -F"\t" '{print $2}')
	sed -i "s/$id/$desc/g" arMatrix.txt
done

rm -f allhmmtables id2desc.tsv

#plot arMatrix

echo '
library(pheatmap)

df<-read.table("arMatrix.txt",header = T,sep = "\t",row.names = 1,check.names = F)

wd<- 1+ncol(df)
ht<- 2+round(nrow(df)/2)

pdf(file="arMatrix.pdf", width = wd, height = ht, onefile=FALSE)
pheatmap(df,cluster_rows = F, cluster_cols = F, 
         cellwidth = 30,cellheight = 30,
         main = "Antibiotic Resistance")
dev.off()
' > pheatmap.R

RHOME/bin/Rscript pheatmap.R && rm -f pheatmap.R
cd ..

chainReaction "RESFAM" $1 $DEBUG

