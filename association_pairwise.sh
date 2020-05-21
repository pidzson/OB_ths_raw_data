#!/bin/bash
#defaults
fields=2
#pval=P_Logit
pval=P_adj
let freq="5-$fields"
freq="0.0${freq}5"
covar=
ucovar=
tcovar=", -log(0.05) with l lc 7 title \"Confidence threshold\""

if [ $# -gt 0 ]; then
    fields=$1
    if [ $# -gt 1 ]; then
        covar="_covar"
        ucovar=" using covariance"
        tcovar=
    fi
fi

let digs=$fields*2

mkdir -p results_$fields

loci=("DQA1" "DQB1" "DRB1")

#pairs
for (( li=0; li<${#loci[@]}; li++ )); do
    for (( lj=$li+1; lj<${#loci[@]}; lj++ )); do
        ptest=${loci[$li]}_${loci[$lj]}
        if [ "$covar" == "_covar" ]; then
            python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out results_$fields/pairwise_${ptest}${covar}.txt --combinations ${loci[$li]},${loci[$lj]} --freq $freq
        else
            python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out results_$fields/pairwise_${ptest}${covar}.txt --combinations ${loci[$li]},${loci[$lj]} --freq $freq
        fi

        if [ ! -f results_$fields/pairwise_${ptest}${covar}.txt ]; then
            echo "PyHLA was unsuccessful"
        else
            echo -e "\nsorting...\n\n"

            scol=$(head -n1 results_$fields/pairwise_${ptest}${covar}.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
            head -n1 results_$fields/pairwise_${ptest}${covar}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/pairwise_${ptest}${covar}_sorted.txt
            tail -n +2 results_$fields/pairwise_${ptest}${covar}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/pairwise_${ptest}${covar}_sorted.txt
            labels=$(tail -n +2 results_$fields/pairwise_${ptest}${covar}_sorted.txt | awk '{s+=1; if(s>1) printf ","; printf "\"%s\" %d",$1,s-1 }')
            cat results_$fields/pairwise_${ptest}${covar}_sorted.txt | awk -v scol=$scol '{print $1","$scol","$8}' > results_$fields/pairwise_${ptest}${covar}_plot_data.csv
            echo -e "set term png\nset out \"results_$fields/pairwise_${ptest}${covar}.png\"\nset xla 'Loci'\nset title 'Results of pairwise association test with logistic${ucovar} on $fields fields'\nset yla 'Adjusted P-value (log)'\nset key autotitle columnhead\nunset xtics\nset xtics rotate by -90 ($labels)\nset grid xtics\nplot \"results_$fields/pairwise_${ptest}${covar}_sorted.txt\" u 0:(-log(\$$scol)):(log(\$8/$freq)+.5) with lp lc 2 pt 4 ps variable title \"\"${tcovar}\nset term svg\nset out \"results_$fields/pairwise_${ptest}${covar}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
            gnuplot -p -c vis.plt
            rm vis.plt
        fi
    done
done

#triplet
let freq="30 - 5 * $fields"
freq="0.0${freq}"
ptest=${loci[0]}_${loci[1]}_${loci[2]}
if [ "$covar" == "_covar" ]; then
    python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --covar preproc/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out results_$fields/pairwise_${ptest}${covar}.txt --combinations ${loci[0]},${loci[1]},${loci[2]} --freq $freq
else
    python ../PyHLA/PyHLA.py --input preproc/${fields}field_pyhla_input.tsv --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out results_$fields/pairwise_${ptest}${covar}.txt --combinations ${loci[0]},${loci[1]},${loci[2]} --freq $freq
fi

if [ ! -f results_$fields/pairwise_${ptest}${covar}.txt ]; then
    echo "PyHLA was unsuccessful"
else
    echo -e "\nsorting...\n\n"

    scol=$(head -n1 results_$fields/pairwise_${ptest}${covar}.txt | awk -v pcol=$pval '{for(i=1; i<= NF; ++i) if($i == pcol) {print i; break;}}')
    head -n1 results_$fields/pairwise_${ptest}${covar}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee results_$fields/pairwise_${ptest}${covar}_sorted.txt
    tail -n +2 results_$fields/pairwise_${ptest}${covar}.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a results_$fields/pairwise_${ptest}${covar}_sorted.txt
    labels=$(tail -n +2 results_$fields/pairwise_${ptest}${covar}_sorted.txt | awk '{s+=1; if(s>1) printf ","; printf "\"%s\" %d",$1,s-1 }')
    cat results_$fields/pairwise_${ptest}${covar}_sorted.txt | awk -v scol=$scol '{print $1","$scol","$8}' > results_$fields/pairwise_${ptest}${covar}_plot_data.csv
    echo -e "set term png\nset out \"results_$fields/pairwise_${ptest}${covar}.png\"\nset xla 'Loci'\nset title 'Results of pairwise association test with logistic${ucovar} on $fields fields'\nset yla 'Adjusted P-value (log)'\nset key autotitle columnhead\nunset xtics\nset xtics rotate by -90 ($labels)\nset grid xtics\nplot \"results_$fields/pairwise_${ptest}${covar}_sorted.txt\" u 0:(-log(\$$scol)):(log(\$8/$freq)+.5) with lp lc 2 pt 4 ps variable title \"\"${tcovar}\nset term svg\nset out \"results_$fields/pairwise_${ptest}${covar}.svg\"\nrep\nset term qt\nset out\nrep" > vis.plt
    gnuplot -p -c vis.plt
    rm vis.plt
fi

