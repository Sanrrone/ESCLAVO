#!/bin/bash
set -e
shopt -s direxpand
shopt -s expand_aliases
if [ -f ~/.bash_profile ]; then source ~/.bash_profile; fi
if [ -f ~/.bashrc ]; then source ~/.bashrc; fi
if [ -f ~/.bash_alias ]; then source ~/.bash_alias; fi

export ESCLAVOHOME=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
source $ESCLAVOHOME/modules/checkVariables.sh
source $ESCLAVOHOME/modules/statusb.sh
source $ESCLAVOHOME/modules/humanDecont.sh
source $ESCLAVOHOME/modules/qc.sh
source $ESCLAVOHOME/modules/statusa.sh
source $ESCLAVOHOME/modules/assignTaxonomy.sh

#usage: bash 16s18sits.sh --force -p /home/sandro/Programas/ESCLAVO/projects -f /home/sandro/Programas/ESCLAVO/0-raw -pt .fastq.gz
#software requirements
# FASTQC (better download and create alias), MULTIQC, DADA2 (r), silvaDB https://zenodo.org/record/1172783#.XUMZwfx7k5k

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -p|--projectfolder)
    PROJECTFOLDER="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--fastqfolder)
    FASTQFOLDER="$2"
    shift # past argument
    shift # past value
    ;;
    -pt|--fastqpattern)
    PATTERN="$2"
    shift # past argument
    shift # past value
    ;;
    -to|--tolerance)
    TOLERANCE="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--module)
    MODULE="$2"
    shift # past argument
    shift # past value
    ;;
    --force)
    FORCE=true
    shift # past argument
    ;;
    --debug)
    DEBUG=true
    set -ex
    shift # past argument
    ;;
    -h|--help)
    echo "usage: bash 16s18sits.sh -p [project folder] -f [fastq foder] -pt [fastq pattern]"
    echo "example: bash 16s18sits.sh -p mynewproject -f 0-raw -pt .fastq.gz"
    echo -e "\noptions available:\n"
    echo "-p Project folder, if no exist, the Pipeline will assume is the actual folder"
    exit
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
if [ "$TOLERANCE" == "" ];then
	TOLERANCE=0.1
fi
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ -n $1 ]]; then
    echo "the following argument lack of parameter: $1"
    exit
fi

checkVariables


cd $PROJECTFOLDER
PNAME=$(pwd | awk -F"/" '{print $NF"_eConf.tsv"}')
export PCONF=$(pwd | awk -v pname=$PNAME '{print $0"/"pname}')

if [ ! -f "$PNAME" ] || [ $FORCE ]; then
	echo "	value" > $PNAME
	echo "pfolder	$(pwd)" >> $PNAME
	echo "ffolder	$FASTQFOLDER" >> $PNAME
	echo "fqpattern	$PATTERN" >> $PNAME
	echo "analysis	16s18sits" >> $PNAME
	echo "aversion	1.0" >> $PNAME
	echo "created	$(date --iso-8601)" >> $PNAME
	echo "status	open" >> $PNAME
	echo "pPercent	0" >> $PNAME
	echo "lastStep	statusb" >> $PNAME
fi

for mod in $MODULE
do
    case $mod in
        "statusb")
            statusb $FORCE
        ;;
        "humanDecont")
            humanDecont
        ;;
        "qc")
            qc
        ;;
        "statusa")
            statusa $FORCE
        ;;
        "assignTaxonomy")
            assignTaxonomy 2-taxInsight
        ;;
        *)
            echo "Module $mod not recognized"
        ;;
    esac
    cd $PROJECTFOLDER
done

echo "ESCLAVO: Pipeline done :)"

    
#         ___     _,.--.,_                ___     _,.--.,_                 ___     _,.--.,_               ___     _,.--.,_    
#      .-~   ~--"~-.   ._ "-.          .-~   ~--"~-.   ._ "-.          .-~   ~--"~-.   ._ "-.          .-~   ~--"~-.   ._ "-.    
#     /      ./_    Y    "-. \        /      ./_    Y    "-. \        /      ./_    Y    "-. \        /      ./_    Y    "-. \    
#    Y       :~     !         Y      Y       :~     !         Y      Y       :~     !         Y      Y       :~     !         Y    
#    lq p    |     /         .|      lq p    |     /         .|      lq p    |     /         .|      lq p    |     /         .|    
# _   \. .-, l    /          |j   _   \. .-, l    /          |j   _   \. .-, l    /          |j   _   \. .-, l    /          |j    
#()\___) |/   \_/";          !   ()\___) |/   \_/";          !   ()\___) |/   \_/";          !   ()\___) |/   \_/";          !    
# \._____.-~\  .  ~\.      ./     \._____.-~\  .  ~\.      ./     \._____.-~\  .  ~\.      ./     \._____.-~\  .  ~\.      ./    
#            Y_ Y_. "vr"~  T                 Y_ Y_. "vr"~  T                 Y_ Y_. "vr"~  T                 Y_ Y_. "vr"~  T    
#            (  (    |L    j                 (  (    |L    j                 (  (    |L    j                 (  (    |L    j    
#            [nn[nn..][nn..]                 [nn[nn..][nn..]                 [nn[nn..][nn..]                 [nn[nn..][nn..]    
#        ~~~~~~~~~~~~~~~~~~~~~~~         ~~~~~~~~~~~~~~~~~~~~~~~         ~~~~~~~~~~~~~~~~~~~~~~~         ~~~~~~~~~~~~~~~~~~~~~~~   