#!/bin/bash
#defaults
fields=2
#pval=P_Logit
pval=P_adj

if [ $# -gt 0 ]; then
    fields=$1
fi

let digs=$fields*2

mkdir -p results_$fields

loci=("DQA1" "DQB1" "DRB1")

#pairs
for (( li=0; li<${#loci[@]}; li++ )); do
    for (( lj=$li+1; lj<${#loci[@]}; lj++ )); do
        ptest=${loci[$li]}_${loci[$lj]}
        python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out results_$fields/pairwise_${ptest}.txt --combinations ${loci[$li]},${loci[$lj]} --freq 0.03

        if [ ! -f results_$fields/pairwise_${ptest}.txt ]; then
            echo "PyHLA was unsuccessful"
        else
            echo -e "\nsorting...\n\n"

            scol=$(head -n1 results_$fields/pairwise_${ptest}.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
            head -n1 results_$fields/pairwise_${ptest}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/pairwise_${ptest}_sorted.txt
            tail -n +2 results_$fields/pairwise_${ptest}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/pairwise_${ptest}_sorted.txt

            echo -e "set term png\nset out \"results_$fields/pairwise_${ptest}.png\"\nplot [0:10] \"results_$fields/pairwise_${ptest}_sorted.txt\" u 0:(-log(\$$scol)) with lp title \"\", \"results_$fields/pairwise_${ptest}_sorted.txt\" u 0:(-log(\$$scol)+.05):1 with labels title \"\"\nset term svg\nset out \"results_$fields/pairwise_${ptest}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
            gnuplot -p -c vis.plt
        fi
    done
done

#triplet
ptest=${loci[0]}_${loci[1]}_${loci[2]}
python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out results_$fields/pairwise_${ptest}.txt --combinations ${loci[0]},${loci[1]},${loci[2]} --freq 0.01

if [ ! -f results_$fields/pairwise_${ptest}.txt ]; then
    echo "PyHLA was unsuccessful"
else
    echo -e "\nsorting...\n\n"

    scol=$(head -n1 results_$fields/pairwise_${ptest}.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
    head -n1 results_$fields/pairwise_${ptest}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/pairwise_${ptest}_sorted.txt
    tail -n +2 results_$fields/pairwise_${ptest}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/pairwise_${ptest}_sorted.txt

    echo -e "set term png\nset out \"results_$fields/pairwise_${ptest}.png\"\nplot [0:10] \"results_$fields/pairwise_${ptest}_sorted.txt\" u 0:(-log(\$$scol)) with lp title \"\", \"results_$fields/pairwise_${ptest}_sorted.txt\" u 0:(-log(\$$scol)+.05):1 with labels title \"\"\nset term svg\nset out \"results_$fields/pairwise_${ptest}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
    gnuplot -p -c vis.plt
    rm vis.plt
fi

