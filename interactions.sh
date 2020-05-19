#!/bin/bash
#defaults
fields=2
itest=fisher

if [ $# -gt 0 ]; then
    fields=$1
    if [ $# -gt 1 ]; then
        itest=$2
    fi
fi

let digs=$fields*2

case "$itest" in
    fisher)
        #pval=P_FET
        pval=P_adj
        ;;
    chisq)
        #pval=P_Chisq
        pval=P_adj
        ;;
    *)
        echo "cannot use $itest as test type"
        exit 1
esac

mkdir -p results_$fields

python ../PyHLA/PyHLA.py --interaction --level allele --input preproc/${fields}field_pyhla_input.tsv --digit $digs --test $itest --model allelic --adjust Bonferroni --print --out results_$fields/interaction_${itest}.txt 

if [ ! -f results_$fields/interaction_${itest}.txt ]; then
    echo "PyHLA was unsuccessful"
else
    #echo -e "\nfiltering for significant interactions\n\n"
    #cat results_$fields/interaction_${itest}.txt | tr -s ' ' | sed 's/ /\t/g' | cut -f1-10 | awk '{ if($1=="ID1") print $0; else{ printf "%s\t%s" ,$1,$2; for(i=3; i<=10; i++) if($i < 0.05) printf "\t%d",1; else printf "\t%d",0; printf "\n"}}' | tr -s ' ' | sed 's/ /\t/g' | tee results_$fields/interaction_${itest}_graph.txt
    python gen_heatmap.py results_$fields/interaction_${itest}.txt | tee results_$fields/interaction_${itest}_graph.txt
    echo -e "set term png\nset out \"results_$fields/interaction_${itest}.png\"\nset xla 'Loci'\nset title 'Results of interaction test with $itest on $fields fields'\nset yla 'Loci'\nset xtics rotate by -90\nset grid xtics ytics\nset palette defined ( 0 \"red\", 1 \"white\", 2 \"green\" )\nset cbrange [-10:10]\nplot \"results_$fields/interaction_${itest}_graph.txt\" matrix rowheaders columnheaders u 1:2:3 w image title \"\"\nset term svg\nset out \"results_$fields/interaction_${itest}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
    gnuplot -p -c vis.plt
    rm vis.plt
fi

