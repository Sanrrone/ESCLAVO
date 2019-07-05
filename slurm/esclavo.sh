#!/usr/bin/env bash

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

# $1 is the config file

if [ "$1" == "" ];then
	echo "* config file is mandatory"
	echo "usage: bash esclavo.sh config.conf"
	exit
fi

grep -v "^#" $1 |sed "s/#.*//g" |awk '{if($0!="")print}' > tmpconf

while read parameter
do
	Pname=$(echo "$parameter" |awk -F"=" '{print $1}')	
	case $Pname in
		"ESCLAVOPATH")
			ESCLAVOPATH=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"MODULE")
			MODULE=$(echo "$parameter" | awk -F"=" '{print $2}' | sed "s/,/ /g")
		;;
		"PROKKAHOME")
			PROKKAHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"GENOMESPATH")
			GENOMESPATH=$(echo "$parameter" | awk -F"=" '{print $2}')
			GENOMESPATH=$(readlink -f $GENOMESPATH)
		;;
		"PROKKAEXTRACODE")
			PROKKAEXTRACODE=$(echo "$parameter" | awk -F"=" '{gsub("\"","");print $2}')	
		;;
		"GETHOMOLOGUESHOME")
			GETHOMOLOGUESHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"GTEXTRACODE")
			GTEXTRACODE=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"ANIHOME")
			ANIHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"MUMMERHOME")
			MUMMERHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"AMPHORA2HOME")
			AMPHORA2HOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PYTHONBIN")
			PYTHONBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PERLBIN")
			PERLBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"USERMAIL")
			USERMAIL=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PROKKATIME")
			PROKKATIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PROKKAPARTITIONS")
			PROKKAPARTITIONS=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"FAAPATH")
			FAAPATH=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"GETHTIME")
			GETHTIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"GETHPARTITIONS")
			GETHPARTITIONS=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"GETHIDENTITY")
			GETHIDENTITY=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"ANITIME")
			ANITIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"ANIPARTITIONS")
			ANIPARTITIONS=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"AMPHORA2TIME")
			AMPHORA2TIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"AMPHORA2PARTITIONS")
			AMPHORA2PARTITIONS=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"TRANSXBIN")
			TRANSXBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"EXTRACTSEQ")
			EXTRACTSEQ=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"RESFAMTIME")
			RESFAMTIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"RESFAMPARTITION")
			RESFAMPARTITION=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"HMMERHOME")
			HMMERHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"RESFAMHMMINDEX")
			RESFAMHMMINDEX=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"RHOME")
			RHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PRINSEQHOME")
			PRINSEQHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PRINSEQEXTRACODE")
			PRINSEQEXTRACODE=$(echo "$parameter" |awk -F"=" '{print $2}')
		;;
		"INPUTREADSLIST")
			INPUTREADSLIST=$(echo "$parameter" |awk -F"=" '{print $2}')
		;;
		"QCTIME")
			QCTIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"QCPARTITION")
			QCPARTITION=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PEARBIN")
			PEARBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"MERGERTIME")
			MERGERTIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"MERGERPARTITION")
			MERGERPARTITION=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PEARQ")
			PEARQ=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PEARMEMORY")
			PEARMEMORY=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"PEAREXTRACODE")
			PEAREXTRACODE=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"SPADESBIN")
			SPADESBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"SPADESMEMORY")
			SPADESMEMORY=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"SPADESEXTRACODE")
			SPADESEXTRACODE=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"ASSEMBLYTIME")
			ASSEMBLYTIME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"ASSEMBLYPARTITION")
			ASSEMBLYPARTITION=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
		"KRONAHOME")
			KRONAHOME=$(echo "$parameter" | awk -F"=" '{print $2}')
		;;
	esac
done < <(grep "" tmpconf)
rm -f tmpconf


#functions
function qualityControlFunction {

	if [ -d QC ];then
		echo "* Another instance of esclavo (with QC module) was performed"
		echo "* Check the config file or delete QC folder"
		exit
	fi

	GNUMS=$(wc -l $INPUTREADSLIST |awk '{print $1}')
	seq 1 1 $GNUMS > jobids.txt
	paste $INPUTREADSLIST jobids.txt > $INPUTREADSLIST.esclavocopy && rm jobids.txt
	INPUTREADSLIST="$INPUTREADSLIST.esclavocopy"

	sed "s/QCTIME/$QCTIME/g" ${ESCLAVOPATH}/Modules/qcModule.sh > qcModule.sh
	sed -i "s/QCPARTITION/$QCPARTITION/g" qcModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" qcModule.sh
	sed -i "s:PERLBIN:$PERLBIN:g" qcModule.sh
	sed -i "s:PRINSEQHOME:$PRINSEQHOME:g" qcModule.sh
	sed -i "s/PRINSEQEXTRACODE/$PRINSEQEXTRACODE/g" qcModule.sh
	sed -i "s:INPUTREADSLIST:$INPUTREADSLIST:g" qcModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" qcModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" qcModule.sh


	sed -i "s:GNUMS:$GNUMS:g" qcModule.sh

	rm -f QC_JOBS.txt
	sbatch qcModule.sh $1 $DEBUG
}

function mergerReadsFunction {

	if [ -d MERGE ];then
		echo "* Another instance of esclavo (with MERGE module) was performed"
		echo "* Check the config file or delete MERGED folder"
		exit
	fi

	sed "s/MERGERTIME/$MERGERTIME/g" ${ESCLAVOPATH}/Modules/mergerModule.sh > mergerModule.sh
	sed -i "s/MERGERPARTITION/$MERGERPARTITION/g" mergerModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" mergerModule.sh
	sed -i "s:PEARBIN:$PEARBIN:g" mergerModule.sh
	sed -i "s/PEARQ/$PEARQ/g" mergerModule.sh
	sed -i "s/PEARMEMORY/$PEARMEMORY/g" mergerModule.sh
	sed -i "s/PEAREXTRACODE/$PEAREXTRACODE/g" mergerModule.sh
	sed -i "s:INPUTREADSLIST:$INPUTREADSLIST:g" mergerModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" mergerModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" mergerModule.sh

	GNUMS=$(wc -l $INPUTREADSLIST |awk '{print $1}')
	seq 1 1 $GNUMS > jobids.txt
	paste $INPUTREADSLIST jobids.txt > $INPUTREADSLIST.esclavocopy && rm jobids.txt
	INPUTREADSLIST="$INPUTREADSLIST.esclavocopy"
	sed -i "s:GNUMS:$GNUMS:g" qcModule.sh

	rm -f MERGE_JOBS.txt

	sbatch mergerModule.sh $1 $DEBUG

}
function assemblyFunction {

	if [ -d ASSEMBLY ];then
		echo "* Another instance of esclavo (with ASSEMBLY module) was performed"
		echo "* Check the config file or delete ASSEMBLY folder"
		exit
	fi

	sed "s/ASSEMBLYTIME/$ASSEMBLYTIME/g" ${ESCLAVOPATH}/Modules/assemblyModule.sh > assemblyModule.sh
	sed -i "s/ASSEMBLYPARTITION/$ASSEMBLYPARTITION/g" assemblyModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" assemblyModule.sh
	sed -i "s:SPADESBIN:$SPADESBIN:g" assemblyModule.sh
	sed -i "s/SPADESMEMORY/$SPADESMEMORY/g" assemblyModule.sh
	sed -i "s/SPADESEXTRACODE/$SPADESEXTRACODE/g" assemblyModule.sh
	sed -i "s:INPUTREADSLIST:$INPUTREADSLIST:g" assemblyModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" assemblyModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" assemblyModule.sh

	sbatch assemblyModule.sh $1 $DEBUG

}
function prokkaFunction {

	sed "s/PROKKATIME/$PROKKATIME/g" ${ESCLAVOPATH}/Modules/prokkaModule.sh > prokkaModule.sh
	sed -i "s/PROKKAPARTITIONS/$PROKKAPARTITIONS/g" prokkaModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" prokkaModule.sh
	sed -i "s:PERLBIN:$PERLBIN:g" prokkaModule.sh
	sed -i "s:GENOMESPATH:$GENOMESPATH:g" prokkaModule.sh
	sed -i "s:PROKKAHOME:$PROKKAHOME:g" prokkaModule.sh
	sed -i "s/PROKKAEXTRACODE/$PROKKAEXTRACODE/g" prokkaModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" prokkaModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" prokkaModule.sh


	sbatch prokkaModule.sh $1 $DEBUG

}
function gethomologuesFunction {

	sed "s/GETHTIME/$GETHTIME/g" ${ESCLAVOPATH}/Modules/getHModule.sh > getHModule.sh
	sed -i "s/GETHPARTITIONS/$GETHPARTITIONS/g" getHModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" getHModule.sh
	sed -i "s:PERLBIN:$PERLBIN:g" getHModule.sh
	sed -i "s:FAAPATH:$FAAPATH:g" getHModule.sh
	sed -i "s:GETHOMOLOGUESHOME:$GETHOMOLOGUESHOME:g" getHModule.sh
	sed -i "s:GETHIDENTITY:$GETHIDENTITY:g" getHModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" getHModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" getHModule.sh
	
	sbatch getHModule.sh $1 $DEBUG

}
function aniFunction {

	sed "s/ANITIME/$ANITIME/g" ${ESCLAVOPATH}/Modules/aniModule.sh > aniModule.sh
	sed -i "s/ANIPARTITIONS/$ANIPARTITIONS/g" aniModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" aniModule.sh
	sed -i "s:ANIHOME:$ANIHOME:g" aniModule.sh
	cd $GENOMESPATH
	GENOMESPATH=$(pwd)
	cd $OLDPWD
	sed -i "s:GENOMESPATH:$GENOMESPATH:g" aniModule.sh
	sed -i "s:MUMMERHOME:$MUMMERHOME:g" aniModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" aniModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" aniModule.sh

	sbatch aniModule.sh $1 $DEBUG

}
function amphora2Function {

	sed "s/AMPHORA2TIME/$AMPHORA2TIME/g" ${ESCLAVOPATH}/Modules/amphora2Module.sh > amphora2Module.sh
	sed -i "s/AMPHORA2PARTITIONS/$AMPHORA2PARTITIONS/g" amphora2Module.sh
	sed -i "s/USERMAIL/$USERMAIL/g" amphora2Module.sh
	sed -i "s:PERLBIN:$PERLBIN:g" amphora2Module.sh
	sed -i "s:AMPHORA2HOME:$AMPHORA2HOME:g" amphora2Module.sh
	sed -i "s:GENOMESPATH:$GENOMESPATH:g" amphora2Module.sh
	sed -i "s:TRANSXBIN:$TRANSXBIN:g" amphora2Module.sh
	sed -i "s:EXTRACTSEQ:$EXTRACTSEQ:g" amphora2Module.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" amphora2Module.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" amphora2Module.sh

	#ignoring all fastas with *.orf extension
	rm -f $GENOMESPATH/*.orf
	ls -1 $GENOMESPATH/* > samples.txt
	GNUMS=$(wc -l samples.txt |awk '{print $1}')

	seq 1 1 $GNUMS > jobids.txt
	paste samples.txt jobids.txt > tmpsamples && rm samples.txt jobids.txt
	mv tmpsamples samples.txt

	sed -i "s:GNUMS:$GNUMS:g" amphora2Module.sh

	#hours=$(wc -l samples.txt |awk '{print int($1*10/60)}')
	#minutes=$(wc -l samples.txt |awk -v hrs=$hours '{if(hours==0){print int($1*10)+2}else{print (hrs*60)-int($1*10)}}' |awk '{if($1<10){print "10"}else{print $1}}')

	#if [ $((hours)) -le 4 ];then
	#	partition="debug,defq,short"
	#else
	#	if [ $((hours)) -gt 48 ];then
	#		partition="defq"
	#	else
	#		partition="defq,short"
	#	fi
	#fi

	sbatch amphora2Module.sh $1 $DEBUG

}
function kroningFunction {

	sed "s/KRONINGTIME/$KRONINGTIME/g" ${ESCLAVOPATH}/Modules/kroningModule.sh > kroningModule.sh
	sed -i "s/KRONINGPARTITIONS/$KRONINGPARTITIONS/g" kroningModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" kroningModule.sh
	sed -i "s/KRONAHOME/$KRONAHOME/g" kroningModule.sh
	sed -i "s/GFFFILE/$GFFFILE/g" kroningModule.sh
	sed -i "s/SAMFILE/$SAMFILE/g" kroningModule.sh
	#sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" kroningModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" kroningModule.sh



	if mkdir Kroning > /dev/null 2>&1; then
		cd KRONING
	else
		cd KRONING
		rm -rf *
	fi

	mv ../kroningModule.sh .
	sbatch kroningModule.sh ../$1 $DEBUG

}
function resfamFunction {
	
	sed "s/RESFAMTIME/$RESFAMTIME/g" ${ESCLAVOPATH}/Modules/resfamModule.sh > resfamModule.sh
	sed -i "s/RESFAMPARTITION/$RESFAMPARTITION/g" resfamModule.sh
	sed -i "s/USERMAIL/$USERMAIL/g" resfamModule.sh
	sed -i "s:FAAPATH:$FAAPATH:g" resfamModule.sh	
	sed -i "s:HMMERHOME:$HMMERHOME:g" resfamModule.sh
	sed -i "s:RESFAMHMMINDEX:$RESFAMHMMINDEX:g" resfamModule.sh
	sed -i "s:RHOME:$RHOME:g" resfamModule.sh
	sed -i "s:CALLSLURM:$ESCLAVOPATH/Modules/callSLURM.sh:g" resfamModule.sh
	sed -i "s:CHAINREACTION:$ESCLAVOPATH/Modules/chainReaction.sh:g" resfamModule.sh

	sbatch resfamModule.sh $1 $DEBUG

}

function reorderModulesFunction {
	
	IFS=' ' read -a array <<< "$MODULE"
	MODULE=""
	#make order module
	if echo ${array[@]} | grep -q -w "QC"; then 
		MODULE=$(echo "QC $MODULE")
	fi
	if echo ${array[@]} | grep -q -w "MERGE"; then 
		MODULE=$(echo "$MODULE MERGE")
	fi	
	if echo ${array[@]} | grep -q -w "ASSEMBLY"; then 
		MODULE=$(echo "$MODULE ASSEMBLY")
	fi	
	if echo ${array[@]} | grep -q -w "ANI"; then 
		MODULE=$(echo "$MODULE ANI")
	fi	
	if echo ${array[@]} | grep -q -w "AMPHORA2"; then 
		MODULE=$(echo "$MODULE AMPHORA2")
	fi	
	if echo ${array[@]} | grep -q -w "KRONING"; then 
		MODULE=$(echo "$MODULE KRONING")
	fi	
	if echo ${array[@]} | grep -q -w "PROKKA"; then 
		MODULE=$(echo "$MODULE PROKKA")
	fi	
	if echo ${array[@]} | grep -q -w "GETHOMOLOGUES"; then 
		MODULE=$(echo "$MODULE GETHOMOLOGUES")
	fi	
	if echo ${array[@]} | grep -q -w "RESFAM"; then 
		MODULE=$(echo "$MODULE RESFAM")
	fi	

}

function criticalvariablesFunction {

	pass=0
	chain=0
	echo "Checking critical variables:"

	#to ensure the correcto order of modules
	reorderModulesFunction


	if [ "$ESCLAVOPATH" == "" ];then
		echo "* You must provide the ESCLAVOPATH in config file"
		pass=$((pass+1))
	fi

	if [ ! -f "$PYTHONBIN" ];then
		echo "* You must provide the python binary e.g. /usr/bin/python"
		pass=$((pass+1))
	fi

	if [ ! -f "$PERLBIN" ];then
		echo "* You must provide the perl binary e.g. /usr/bin/perl"
		pass=$((pass+1))
	fi

	if [[ "$MODULE" =~ "QC" ]];then
		if [ "$PRINSEQHOME" == "" ];then
			echo "* no PRINSEQHOME was specified in the config file"
			pass=$((pass+1))
		fi
		if ! [ -f ${PRINSEQHOME}/prinseq-lite.pl ];then
			echo "* prinseq-lite.pl no exist in $PRINSEQHOME"
			pass=$((pass+1))
		fi
		if ! [ -f $INPUTREADSLIST ] || [ "$INPUTREADSLIST" == "" ];then
			echo "* invalid INPUTREADSLIST, check if the file exist or if the field is not empty"
			pass=$((pass+1))
		fi
		chain=$((chain+1))
		
	fi

	if [[ "$MODULE" =~ "MERGE" ]];then
		if ! [ -f $PEARBIN ];then
			echo "* PEARBIN no exist"
			pass=$((pass+1))
		fi
		if [ "$INPUTREADSLIST" == "" ];then
			echo "* provide an INPUTREADSLIST is mandatory"
			pass=$((pass+1))
		fi
		if [ "$PEARMEMORY" == "" ];then
			PEARMEMORY="64G"
		fi
		chain=$((chain+1))
		
	fi

	if [[ "$MODULE" =~ "ASSEMBLY" ]];then
		if ! [ -f $SPADESBIN ];then
			echo "* SPADESBIN no exist"
			pass=$((pass+1))
		fi
		if [ "$INPUTREADSLIST" == "" ];then
			echo "* provide an INPUTREADSLIST is mandatory"
			pass=$((pass+1))
		fi
		if [ "$SPADESMEMORY" == "" ];then
			SPADESMEMORY="64G"
		fi
		chain=$((chain+1))
		
	fi

	if [[ "$MODULE" =~ "PROKKA" ]] || [[ "$MODULE" =~ "ANI" ]] || [[ "$MODULE" =~ "AMPHORA2" ]];then
		if ! [[ "$MODULE" =~ "ASSEMBLY" ]];then
			if [ ! -d "$GENOMESPATH" ];then
				echo "* Prokka, Amphora2 and Pyani needs GENOMESPATH variable in the config file (no exist)"
				pass=$((pass+1))
			fi
		fi
	fi

	if [[ "$MODULE" =~ "PROKKA" ]]; then

		if [ "$PROKKAHOME" == "" ];then
			echo "* No PROKKAHOME is specified"
			pass=$((pass+1))
		fi
		if ! [ -f ${PROKKAHOME}/bin/prokka ];then
			echo "* prokka no exist in $PROKKAHOME/bin"
			pass=$((pass+1))
		fi
		chain=$((chain+1))

	fi

	if [[ "$MODULE" =~ "GETHOMOLOGUES" ]]; then

		if [ "$GETHOMOLOGUESHOME" == "" ];then
			echo "* No GETHOMOLOGUESHOME is specified"
			pass=$((pass+1))
		fi
		if ! [ -f $GETHOMOLOGUESHOME/get_homologues.pl ] || ! [ -f $GETHOMOLOGUESHOME/compare_clusters.pl ] || ! [ -f $GETHOMOLOGUESHOME/parse_pangenome_matrix.pl ]
		then
			echo "* get_homologues.pl or compare_clusters.pl or parse_pangenome_matrix.pl no exist in $GETHOMOLOGUESHOME"
			pass=$((pass+1))
		fi
	fi

	if [[ "$MODULE" =~ "ANI" ]]; then

		if [ "$ANIHOME" == "" ];then
			echo "* No ANIHOME is specified"
			pass=$((pass+1))
		fi
		if ! [ -f $ANIHOME/average_nucleotide_identity.py ];then
			echo "* average_nucleotide_identity.py no exist in $ANIHOME"
			pass=$((pass+1))
		fi
	fi

	if [[ "$MODULE" =~ "AMPHORA2" ]]; then

		if [ "$AMPHORA2HOME" == "" ];then
			echo "* No AMPHORA2HOME is specified"
			pass=$((pass+1))
		fi

		if ! [ -f $AMPHORA2HOME/Scripts/MarkerScanner.pl ] || ! [ -f $AMPHORA2HOME/Scripts/MarkerAlignTrim.pl ] || ! [ -f $AMPHORA2HOME/Scripts/Phylotyping.pl ];then
			echo "* MarkerScanner.pl or MarkerAlignTrim or Phylotyping.pl no exist in $AMPHORA2HOME/Scripts"
			pass=$((pass+1))
		fi
	fi


	if [[ "$MODULE" =~ "KRONING" ]]; then
		if [ -f "$SAMFILE" ];then
			SAMFILE=$(echo "$i" |rev |cut -d "/" -f 1 |rev)
			SAMFILEDIR=$(echo "$i" |rev |cut -d "/" -f 2- |rev)
			cd $SAMFILEDIR
			dbpath=$(pwd)
			SAMFILE=$(echo "$dbpath/$SAMFILE")
			cd $OLDPWD
		else
			echo "* $SAMFILE no exist"
			exit
		fi

		if [ -f "$GFFFILE" ];then
			GFFFILE=`echo "$i" |rev |cut -d "/" -f 1 |rev`
			GFFFILEDIR=`echo "$i" |rev |cut -d "/" -f 2- |rev`
			cd $GFFFILEDIR
			dbpath=`pwd`
			GFFFILE=`echo "$dbpath/$GFFFILE"`
			cd $OLDPWD
		else
			echo "* $GFFFILE no exist"
			exit
		fi
	fi

	if [[ "$MODULE" =~ "RESFAM" ]]; then
		if ! [ -f "$HMMERHOME"/bin/hmmscan ] || [ "$RESFAMHMMINDEX" == "" ];then
			echo "* $HMMERHOME/bin/hmmscan not found or RESFAMHMMINDEX empty in the config file"
			pass=$((pass+1))
		fi

		if ! [ -f "$RESFAMHMMINDEX" ] ;then
			echo "* $RESFAMHMMINDEX not index found"
			pass=$((pass+1))
		else
			ok=$(ls ${RESFAMHMMINDEX}.* |wc -l |awk '{print $1}')
			if ! [ $((ok)) -ge 3 ];then
				echo "* $RESFAMHMMINDEX is not prepared, use hmmpress in your index"
				pass=$((pass+1))
			fi
		fi

		if ! [ -f "$RHOME"/bin/Rscript ];then
			echo "* $RHOME/bin/Rscript not found"
			pass=$((pass+1))
		fi
	fi

	if [[ "$USERMAIL" != ?*@?*.?* ]];then
		echo "* Valid mail is mandatory"
		pass=$((pass+1))
	fi

	if [ $((pass)) -eq 0 ];then
		echo "* All parameters ok"
	else
		exit
	fi

}


#begin the code

criticalvariablesFunction

echo "* Launch Zone"
for d in $MODULE
do
	case $d in
		"QC")
			qualityControlFunction $1
			exit
		;;
		"MERGE")
			mergerReadsFunction $1
			exit
		;;
		"ASSEMBLY")
			assemblyFunction $1
			exit
		;;
		"PROKKA")
			prokkaFunction $1
			exit
	   	;;
	   	"GETHOMOLOGUES")
			gethomologuesFunction $1
	   	;;
		"ANI")
			aniFunction $1
	   	;;
	   	"AMPHORA2")
			amphora2Function $1
	   	;;
	   	"KRONING")
			kroningFunction $1
		;;
		"RESFAM")
			resfamFunction $1
		;;
	esac
done

echo "* Done [ go for a coffee :D ]"
