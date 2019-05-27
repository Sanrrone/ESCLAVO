#!/bin/bash 
#SBATCH -J prokMod
#SBATCH -o prokkaModule.out
#SBATCH -e prokkaModule.err
#SBATCH -t PROKKATIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=PROKKAPARTITIONS
#SBATCH -N 1

if [[ "$@" =~ "--debug" ]]; then
		DEBUG="--debug"
        set -ex
else
        set -e
        DEBUG=""
fi

source CHAINREACTION

CORES=$(grep -c ^processor /proc/cpuinfo)
ls -1 GENOMESPATH/* > samples.txt
gnum=$(wc -l samples.txt |awk '{print $1}')

			echo "#!/bin/bash
#SBATCH -J prokka
#SBATCH -o prokka_%a.out
#SBATCH -e prokka_%a.err
#SBATCH -t PROKKATIME
#SBATCH --partition=PROKKAPARTITIONS
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --array=1-$gnum

module load parallel

i=\$(sed -n \"\$SLURM_ARRAY_TASK_ID\"p samples.txt)
name=\$(echo \"\$i\" |rev |cut -d \"/\" -f 1 |rev)

PERLBIN PROKKAHOME/bin/prokka \$i --outdir \${name}_annot --prefix \$name --addgenes --centre c --locustag l --evalue 1e-6 --cpus $CORES PROKKAEXTRACODE

" > prokka.sh

	callSLURM prokka.sh "Prokka"

	mkdir ANNOT
	mv *_annot ANNOT/.
	mv prokka*.* ANNOT/.

	chainReaction "prokMod" $1 $DEBUG

