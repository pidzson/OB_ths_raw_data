#!/bin/bash
#defaults
fields=2
atest=logistic
freq=0.075
covar=
ucovar=

if [ $# -gt 0 ]; then
    fields=$1
    if [ $# -gt 1 ]; then
        atest=$2
        if [ $# -gt 2 ]; then
            covar="_covar"
            ucovar=" using covariance"
        fi
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

if [ "$covar" == "_covar" ]; then
    python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test $atest --model $mdl --adjust Bonferroni --print --out results_$fields/assoc_${atest}${covar}.txt --freq $freq
else
    python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --digit $digs --assoc --test $atest --model $mdl --adjust Bonferroni --print --out results_$fields/assoc_${atest}${covar}.txt --freq $freq
fi

if [ ! -f results_$fields/assoc_${atest}${covar}.txt ]; then
    echo "PyHLA was unsuccessful"
else
    echo -e "\nsorting...\n\n"

    scol=$(head -n1 results_$fields/assoc_${atest}${covar}.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
    head -n1 results_$fields/assoc_${atest}${covar}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/assoc_${atest}${covar}_sorted.txt
    tail -n +2 results_$fields/assoc_${atest}${covar}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/assoc_${atest}${covar}_sorted.txt
    labels=$(cat results_$fields/assoc_${atest}${covar}_sorted.txt | awk '{s+=1; if(s>1){ if(s!=2) printf ","; printf "\"%s\" %d",$1,s-2 }}')
    echo -e "set term png\nset out \"results_$fields/assoc_${atest}${covar}.png\"\nset xla 'Loci'\nset title 'Results of association test with $atest${ucovar} on $fields fields'\nset yla 'Adjusted P-value (log)'\nset key autotitle columnhead\nunset xtics\nset xtics rotate by -90 ($labels)\nset grid xtics\nplot \"results_$fields/assoc_${atest}${covar}_sorted.txt\" u 0:(-log(\$$scol)):(log(\$8/$freq)+.5) with lp lc 2 pt 4 ps variable title \"\", -log(0.05) with l lc 7 title \"Confidence threshold\"\nset term svg\nset out \"results_$fields/assoc_${atest}${covar}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
    gnuplot -p -c vis.plt
    rm vis.plt
fi

