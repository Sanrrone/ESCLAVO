function checkVariables {
	set -e

if [ "$PROJECTFOLDER" == "" ];then
	export PROJECTFOLDER="."
fi

if [ "$MODULE" == "" ] || [ "$MODULE" == "all" ] || [ "$MODULE" == "ALL" ]; then
	MODULE="statusb humanDecont qc statusa assignTaxonomy report"
fi

}