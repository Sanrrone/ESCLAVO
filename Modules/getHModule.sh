#!/bin/bash 
#SBATCH -J getHModule
#SBATCH -o getHModule.out
#SBATCH -e getHModule.err
#SBATCH -t GETHTIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=GETHPARTITIONS
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


ok=$(ls -1 ${THEFAA}/ |wc -l)
mkdir GETH
if [ $((ok)) -gt 2 ]; then

	#-S identity, -C coverage (75 default)
	PERLBIN GETHOMOLOGUESHOME/get_homologues.pl -d ${THEFAA} -n $CORES -S GETHIDENTITY -C 75 -t 0 -c -z -G
	PERLBIN GETHOMOLOGUESHOME/get_homologues.pl -d ${THEFAA} -n $CORES -S GETHIDENTITY -C 75 -t 0 -c -z -M

else
	echo "* Not sufficient faa files (minimal: 3), check $THEFAA outputs"
	mv getH*.* GETH/.
	chainReaction "GETHOMOLOGUES" $1 $DEBUG
	exit
fi

#part2

cd ${THEFAA}_homologues

	cpath=$(ls -d1 *_/ |awk '{printf "%s,",$1}END{printf "\n"}' |sed 's/.$//')

	PERLBIN GETHOMOLOGUESHOME/compare_clusters.pl -d $cpath -o intersection -m -T

	cd intersection

		PERLBIN GETHOMOLOGUESHOME/parse_pangenome_matrix.pl -m pangenome_matrix_t0.tab -s
		cp *.pdf ../../GETH/.
		cp pangenome_matrix_t0.tab ../../GETH/.
	cd ..
cd ..

mv ${THEFAA}_homologues GETH/.
mv getH*.* GETH/.

chainReaction "GETHOMOLOGUES" $1 $DEBUG
