function makeReport {
	set -e

	rawfolder=$1
	qcfolder=$2
	taxfolder=$3
	proyectname=$(echo $PNAME |awk -F"_" '{print $(NF-1)}')

	cp $ESCLAVOHOME/modules/reportTemplate.Rmd report.Rmd
	sed "s:FASTQFOLDER:$rawfolder:g" $ESCLAVOHOME/modules/reportTemplate.R > report.R
	sed -i "s:PFOLDER:$PROJECTFOLDER:g" report.R
	sed -i "s:PNAME:$proyectname:g" report.R
	sed -i "s:MULTIQCHTML:$qcfolder/multiqc_report.html:g" report.R
	sed -i "s:ABUNDANCEFILE:$taxfolder/abundance.tsv:g" report.R
	sed -i "s:ABUNDANCEFILE:$taxfolder/abundance.tsv:g" report.Rmd
	sed -i "s:ESCLAVOLOGO:$ESCLAVOHOME/DB/esclavologo.jpg:g" report.Rmd

	Rscript report.R

	rm -f report.R report.Rmd report.aux report.dvi report.toc report.out

}