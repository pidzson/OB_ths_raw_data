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
        #pval=P_Logit
        pval=P_adj
        ;;
    fisher)
        mdl=allelic
        #pval=P_FET
        pval=P_adj
        ;;
    chisq)
        mdl=allelic
        #pval=P_Chisq
        pval=P_adj
        ;;
    *)
        echo "cannot use $atest as test type"
        exit 1
esac

mkdir -p results_$fields

python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test $atest --model $mdl --adjust Bonferroni --print --out results_$fields/assoc_${atest}.txt

if [ ! -f results_$fields/assoc_${atest}.txt ]; then
    echo "PyHLA was unsuccessful"
else
    echo -e "\nsorting...\n\n"

    scol=$(head -n1 results_$fields/assoc_${atest}.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
    head -n1 results_$fields/assoc_${atest}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/assoc_${atest}_sorted.txt
    tail -n +2 results_$fields/assoc_${atest}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/assoc_${atest}_sorted.txt

    echo -e "set term png\nset out \"results_$fields/assoc_${atest}.png\"\nplot [0:10] \"results_$fields/assoc_${atest}_sorted.txt\" u 0:(-log(\$$scol)) with lp title \"\", \"results_$fields/assoc_${atest}_sorted.txt\" u 0:(-log(\$$scol)+.05):1 with labels title \"\"\nset term svg\nset out \"results_$fields/assoc_${atest}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
    gnuplot -p -c vis.plt
    rm vis.plt
fi

