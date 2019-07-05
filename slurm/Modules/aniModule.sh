#!/bin/bash 
#SBATCH -J aniModule
#SBATCH -o aniModule.out
#SBATCH -e aniModule.err
#SBATCH -t ANITIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=ANIPARTITIONS
#SBATCH -N 1

	#part1 #is necesary create identity value from ANI output(?)
if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

source CHAINREACTION


			echo "#!/bin/bash
#SBATCH -J ani
#SBATCH -o ani.out
#SBATCH -e ani.err
#SBATCH -t ANITIME
#SBATCH --partition=ANIPARTITIONS
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL

ANIHOME/average_nucleotide_identity.py -i GENOMESPATH -o ANIm --nucmer_exe MUMMERHOME/nucmer -m ANIm -g --gformat pdf

" > ani.sh
	callSLURM ani.sh "ANIm (pyani)"

	mv ani*.* ANIm/.

	chainReaction "ANI" $1 $DEBUG