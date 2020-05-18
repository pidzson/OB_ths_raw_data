#!/bin/bash
#defaults
fields=2
atest=logistic

if [ $# -gt 0 ]; then
    fields=$1
    if [ $# -gt 1 ]; then
        atest=$2
    fi
fi

let digs=$fields*2

case "$atest" in
    logistic)
        mdl=additive
        pval=P_Logit
        ;;
    fisher)
        mdl=allelic
        pval=P_FET
        ;;
    chisq)
        mdl=allelic
        pval=P_Chisq
        ;;
    *)
        echo "cannot use $atest as test type"
        exit 1
esac

mkdir -p results_$fields

echo "preproc/${fields}field_pyhla_input.tsv"

python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage --digit $digs --assoc --test $atest --model $mdl --adjust Bonferroni --print --out results_$fields/assoc.txt 
#--combinations DRB1,DQB1 --freq 0.01

if [ ! -f results_$fields/assoc.txt ]; then
    echo "PyHLA was unsuccessful"
else
    echo -e "\nsorting...\n\n"

    scol=$(head -n1 results_$fields/assoc.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
    head -n1 results_$fields/assoc.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/assoc_sorted.txt
    tail -n +2 results_$fields/assoc.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/assoc_sorted.txt
fi

