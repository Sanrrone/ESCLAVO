if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

function chainReaction {
	echo "chain Reaction from $1"
	MODULESOURCE=$1
	CFILE=$2
	grep -v "^#" $CFILE |sed "s/#.*//g" |awk '{if($0!="")print}' > tmpconf

	while read parameter
	do
		Pname=$(echo "$parameter" |awk -F"=" '{print $1}')	
		case $Pname in
			"ESCLAVOPATH")
				ESCLAVOPATH=$(echo "$parameter" | awk -F"=" '{print $2}' )			
			;;
			"MODULE")
				MODULE=$(echo "$parameter" | awk -F"=" '{print $2}' | sed "s/,/ /g")
			;;
			"GENOMESPATH")
				GENOMESPATH=$(echo "$parameter" | awk -F"=" '{print $2}')
				GENOMESPATH=$(readlink -f $GENOMESPATH)
			;;
			"EXTRACTSEQ")
				EXTRACTSEQ=$(echo "$parameter" | awk -F"=" '{print $2}')
			;;
			"TRANSXBIN")
				TRANSXBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
			;;
			"PERLBIN")
				PERLBIN=$(echo "$parameter" | awk -F"=" '{print $2}')
			;;
		esac
	done < <(grep "" tmpconf)
	rm -f tmpconf

	case $MODULESOURCE in
		"qcModule")
			mkdir QC

			mv *.*pass*.fastq QC/.
			mv qc*.* QC/.

			ls -1 QC/*.pass_1.fastq > pass1
			ls -1 QC/*.pass_2.fastq > pass2
			paste pass1 pass2 >QCoutpass.txt
			rm pass[12]
			
			CHAIN=$(echo "$MODULE" |awk '{if($0~"MERGE" || $1~"ASSEMBLY"){print "CALL"}else{print "EXIT"}}')
			if [ "$CHAIN" == "CALL" ];then

				sed "s/INPUTREADSLIST=.*/INPUTREADSLIST=QCoutpass.txt/g" $CFILE > $MODULESOURCE.config
				sed -i "s/QC,//g" $MODULESOURCE.config
				CFILE=$MODULESOURCE.config
			fi
		;;
		"mergeMod")
			mkdir MERGE
			mv *.fastq MERGE/.
			mv merge* MERGE/.

			ls -1 MERGED/*assembled.fastq > MERGEDassembled.txt
			CHAIN=$(echo "$MODULE" |awk '{if($0~"ASSEMBLY"){print "CALL"}else{print "EXIT"}}')
			if [ "$CHAIN" == "CALL" ];then

				sed "s/INPUTREADSLIST=.*/INPUTREADSLIST=MERGEDassembled.txt/g" $CFILE > $MODULESOURCE.config
				sed -i "s/MERGE,\?//g" $MODULESOURCE.config


				CFILE=$MODULESOURCE.config
			fi
		;;
		"assembly")
			
			mkdir ASSEMBLY
			mv *.assembly ASSEMBLY/.
			mv assembly*.* ASSEMBLY/.

			CHAIN=$(echo "$MODULE" |awk '{if($0~"PROKKA" || $0~"ANI" || $0~"AMPHORA2" || $0~"KRONING"){print "CALL"}else{print "EXIT"}}')
			if [ "$CHAIN" == "CALL" ];then
				echo $(pwd)
				mkdir genomes
				cd ASSEMBLY
				for i in $(ls -1d *.assembly)
				do
					name=$(echo "$i" |sed "s/\.assembly//g" |sed "s/\.assembled//g" |sed "s/\.pass_[12]//g" |sed "s/R[12]//g" |sed "s/\.\.//g" |sed "s/__//g" |sed "s/\.fasta//g" |sed "s/\.fastq//g")
					cp $i/scaffolds_filter500.fasta ../genomes/$name.fasta
				done
				cd ..

				sed "s:GENOMESPATH=.*:GENOMESPATH=genomes:g" $CFILE > $MODULESOURCE.config
				sed -i "s/ASSEMBLY,\?//g" $MODULESOURCE.config

				CFILE=$MODULESOURCE.config
			fi
		;;
		"prokMod")
			CHAIN=$(echo "$MODULE" |awk '{if($0~"GETHOMOLOGUES" || $0~"RESFAM"){print "CALL"}else{print "EXIT"}}')
			if [ "$CHAIN" == "CALL" ];then
				mkdir FAA
				cd ANNOT
				for i in $(ls -1d *_annot)
				do
					cp $i/*.faa ../FAA/.
				done
				cd ..

				sed "s:FAAPATH=.*:FAAPATH=FAA:g" $CFILE > $MODULESOURCE.config
				sed -i "s/ASSEMBLY,\?//g" $MODULESOURCE.config
				CFILE=$MODULESOURCE.config
			fi
		;;
		"AMPH2P1")
			rm AMPHORA2P1_JOBS.txt
			mkdir AMPHORA2
			mv amphora* AMPHORA2/.
			cd AMPHORA2
			rm -f GENOMESPATH/*.orf
			########################################################################################################################################
			
			echo "part2 nucl from pep"
			#part2 get nucl from prot (7gen x min)
			
			for i in $(ls -1d amphora_*);   
			do  
				cd $i;   
				
				for j in $(ls *.pep);   
				do   
					grep ">" $j |head -n1 |sed "s/_[0-9]\+ / /g" |awk '{gsub(">|\\[|- |\\]|\\(|\\)","");print}' |awk -v fna="$i" -v gen="$j" '{if($2>$3){mayor=$2;menor=$3}else{mayor=$3;menor=$2};if($0 ~ "REVERSE SENSE"){print $1, menor, mayor, "-", gen, fna}else{print $1, menor, mayor, "+", gen, fna}}' ; 
				done ;  
				
				cd .. ;   
			done > coords.pep
			
			
			#USE IN SLURM
			#~5 seconds per pep
			hours=$(ls -1 amphora_*/*.pep |wc -l |awk '{print int($1/60/60)}')
			minutes=$(ls -1 amphora_*/*.pep |wc -l |awk -v hrs=$hours '{if(hours==0){print int($1/60/60)+2}else{print (hrs*60)-int($1/60/60)}}'  |awk '{if($1<10){print "5"}else{print $1}}')
			pepnum=$(wc -l coords.pep |awk '{print $1}')
			
			if [ $((pepnum)) -eq 0 ];then
				echo "* no .pep was found, exiting."
				exit
			fi

			seq 1 1 $pepnum > jobids.txt
			paste coords.pep jobids.txt > tmp && rm coords.pep  jobids.txt
			mv tmp coords.pep
	
	
			if [ $((hours)) -le 4 ];then
				partition="debug,defq,short"
			else
				if [ $((hours)) -gt 48 ];then
					partition="defq"
				else
					partition="defq,short"
				fi
			fi

			if [ $((pepnum)) -ge 1000 ];then
				#~5 seconds per pep (to make sure: 30)
				hours=$(ls -1 amphora_*/*.pep |wc -l |awk '{print int($1*30/60/60)}')
				minutes=$(ls -1 amphora_*/*.pep |wc -l |awk -v hrs=$hours '{if(hours==0){print int($1*30/60/60)+2}else{print (hrs*60)-int($1*30/60/60)}}'  |awk '{if($1<10){print "5"}else{print $1}}')

				if [ $((hours)) -le 4 ];then
					partition="debug,defq,short"
				else
					if [ $((hours)) -gt 48 ];then
						partition="defq"
					else
						partition="defq,short"
					fi
				fi
				
				echo "#!/bin/bash
#SBATCH -J AMPH2P2
#SBATCH -o pepxnuc.out
#SBATCH -e pepxnuc.err
#SBATCH -t $hours:$minutes:00
#SBATCH --partition=$partition
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL

source $ESCLAVOPATH/Modules/chainReaction.sh
CFILE=\$1

cat coords.pep |while read c1 c2 c3 c4 c5 c6 jobid
do

    cd \$c6
	fna=\$(echo \$c6 |sed \"s/amphora_//g\")
	bash $EXTRACTSEQ \"\$c1\" \$c2 \$c3 \$c4 \$c5 $GENOMESPATH/\$fna
    cd ..

echo \"\$jobid\" >> AMPHORA2P2_JOBS.txt

done

sbatch $ESCLAVOPATH/Modules/callSLURM.sh AMPH2P2 AMPHORA2P2_JOBS.txt $pepnum ../\$CFILE $DEBUG

" > amphora2P2.sh
			sbatch amphora2P2.sh $CFILE				

			else
			
				echo "#!/bin/bash
#SBATCH -J AMPH2P2
#SBATCH -o pepxnuc.out
#SBATCH -e pepxnuc.err
#SBATCH -t $hours:$minutes:00
#SBATCH --partition=$partition
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --array=1-$pepnum

source $ESCLAVOPATH/Modules/chainReaction.sh
CFILE=\$1

coord=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p coords.pep)
c1=\$(echo \$coord |awk '{print \$1}')
c2=\$(echo \$coord |awk '{print \$2}')
c3=\$(echo \$coord |awk '{print \$3}')
c4=\$(echo \$coord |awk '{print \$4}')
c5=\$(echo \$coord |awk '{print \$5}')
c6=\$(echo \$coord |awk '{print \$6}')
jobid=\$(echo \$coord |awk '{print \$7}')

    cd \$c6

	fna=\$(echo \$c6 |sed \"s/amphora_//g\")
	bash $EXTRACTSEQ \"\$c1\" \$c2 \$c3 \$c4 \$c5 $GENOMESPATH/\$fna

    cd ..

echo \"\$jobid\" >> AMPHORA2P2_JOBS.txt

if [ \$((jobid)) -eq $pepnum ];then

	completedjobs=\$(wc -l AMPHORA2P2_JOBS.txt |awk '{print \$1}')
	if [ \$((completedjobs)) -eq $pepnum ];then
		chainReaction \"AMPH2P2\" ../\$CFILE $DEBUG
	else
		sbatch $ESCLAVOPATH/Modules/callSLURM.sh AMPH2P2 AMPHORA2P2_JOBS.txt $pepnum ../\$CFILE $DEBUG
	fi
fi

" > amphora2P2.sh
			sbatch amphora2P2.sh $CFILE
		fi
		exit

		;;
		"AMPH2P2")
			for i in $(ls -1d amphora_*)
			do  
				cd $i
				ls -1 *.pep.fasta > thepeps.txt
				cd ..
			done   
		
			first=0
			for i in $(ls -1d amphora_*);   
			do  
		
				cd $i
				if [ $((first)) -eq 0 ];then
					first=1
					back="${i}/thepeps.txt"
				else
					grep -w -F -x -f ../$back thepeps.txt > tmp.txt
					cp tmp.txt ../common.txt && rm tmp.txt
					back="common.txt"
				fi
				cd ..
			done  
		
			mkdir groupedGenes
		
			cat common.txt |while read gene
			do
				for i in $(ls -1d amphora_*);
				do  
					cd $i
					cat $gene >> ../groupedGenes/$gene
					echo "" >> ../groupedGenes/$gene
					cd ..
				done
			done
		
			cd groupedGenes
			#remove space from genes
			for i in $(ls -1)
			do
				sed "s/ /_/g" $i > tmp
				rm $i
				mv tmp $i
			done
		
			#align 1seg x gen x genome
			hours=$(ls -1 ../amphora_*/*.pep |wc -l |awk '{print int($1/60/60)}')
			minutes=$(ls -1 ../amphora_*/*.pep |wc -l |awk -v hrs=$hours '{if(hours==0){print $1/60}else{print (hrs*60)- int($1/60)}}'  |awk '{if($1<10){print "10"}else{print $1}}')
		
			if [ $((hours)) -le 4 ];then
				partition="debug,defq,short"
			else
				if [ $((hours)) -gt 48 ];then
					partition="defq"
				else
					partition="defq,short"
				fi
			fi
		
			echo "#!/bin/bash
#SBATCH -J AMPH2P3
#SBATCH -o align.out
#SBATCH -e align.err
#SBATCH -t $hours:$minutes:00
#SBATCH --partition=$partition
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAI
#SBATCH --mail-user=USERMAIL

source $ESCLAVOPATH/Modules/chainReaction.sh

for i in \$(ls -1 *.pep.fasta)
do
	name=\$(echo \$i |sed \"s/\.pep\.fasta//g\")
	$PERLBIN $TRANSXBIN -i \$i -o \$name -p F
done

mkdir aligned
mv *.nt_ali.fasta aligned/.

cd aligned
sed -i '/^>/ s/_c\?[0-9]*:[0-9]*$//g' *.nt_ali.fasta
paste -d'' *.nt_ali.fasta > wholealign.fasta
cd ..

mv aligned ../.
rm -f ../../samples.txt
$ESCLAVOPATH/Modules/chainReaction \"AMPH2P3\" ../\$1 $DEBUG

 " > amphora2P3.sh

 			sbatch amphora2P3.sh $CFILE
			exit
		;;
		"GETH")
			#########FOR THE FUTURE
			#########ALIGN COREGENOME FOR PHYLOSEQ 1 hora x 1000 genomas (~ 16 seg/archivo)
			#corelist=pangenome_matrix_t0__core_list.txt
			#ffnpath=/home/ecastron/lustre/Sandro/gingi/GETH/FFN
			#
			#rm -rf theffn *.ffn
			#mkdir theffn
			#
			#cat $corelist |while read fasta
			#do
			#        grep ">" $fasta |while read line
			#        do
			#                idtokeep=$(echo $line |awk -F"_" '{gsub(">","");print $1"_"$2}')
			#                ffn=$(echo $line |awk -F"_" 'BEGIN{toprint=""}{for(i=3;i<NF;i++){toprint=toprint$i"_"}toprint=toprint$NF}END{print toprint}')
			#
			#                awk -v elid=$idtokeep -v elfasta=$ffn 'BEGIN{RS=">"}{if($0~elid){print ">"elfasta"_"$0;exit}}' $ffnpath/$ffn.fna.ffn
			#        done > $fasta.ffn
			#done
			#
			#mv *.ffn theffn/.
			#cd theffn

			#TRANSLATORX 1200 ALIGN IN 30 MIN
			#for i in $(ls -1 *.ffn)
			#do
			#~/programs/perl/bin/perl ~/programs/translatorx_vLocal.pl -i $i -p F -o $i.msa
			#done
			#cd ..
		;;
		*)
			CHAIN="EXIT"
		;;
	esac

	if [ "$CHAIN" == "CALL" ];then
		bash $ESCLAVOPATH/esclavo.sh $CFILE.config $DEBUG
	else
		LASTMODULE=$(echo "$MODULE" |awk '{print $NF}')
		if [ "$LASTMODULE" == "$MODULESOURCE" ];then
			echo "No more modules to launch :D"
			rm -f samples.txt ASSEMBLY.config QCoutpass.txt QC.config MERGEDassembled.txt MERGE.config
		fi
	fi
}