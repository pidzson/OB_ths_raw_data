#!/bin/bash
if [ $# -gt 0 ]; then
    digs=$1
else
    digs=4
fi

rm -f test_result_$digs.txt test_result_sorted_$digs.txt

#python PyHLA/PyHLA.py --input data/formatted_reftyping.tsv --covar data/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage --digit $digs --assoc --test logistic --model additive --adjust FDR --print --out test_result_$digs.txt 
python PyHLA/PyHLA.py --input OB_ths_raw_data/preproc/3field_pyhla_input.tsv --covar data/formatted_pheno.tsv --covar-name risk_group,gender,ALL_lineage,age_10 --digit $digs --assoc --test logistic --model additive --adjust Bonferroni --print --out test_result_$digs.txt --freq 0.015 --combinations DQB1,DQA1,DRB1
#python PyHLA/PyHLA.py --input data/formatted_reftyping.tsv --digit 8 --assoc --test logistic --model additive --adjust FDR --print --out test_result.txt 

if [ ! -f test_result_$digs.txt ]; then
    echo "PyHLA was unsuccessful"
else

    echo -e "\nsorting...\n\n"

    scol=$(head -n1 test_result_$digs.txt | awk '{for(i=1; i<= NF; ++i) if($i == "P_Logit") {print i; break;}}')
    head -n1 test_result_$digs.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | tee test_result_sorted_$digs.txt
    tail -n +2 test_result_$digs.txt | tr -s " " | sed 's/^ //' | sed 's/ /\t/g' | awk -v scol=$scol '{print $scol "\t" $0}' | sort -g | sed 's/^[^\t]*\t//' | tee -a test_result_sorted_$digs.txt

    echo "plot [0:10] \"test_result_sorted_$digs.txt\" u 0:(-log(\$$scol)) with lp title \"\", \"test_result_sorted_$digs.txt\" u 0:(-log(\$$scol)+.05):1 with labels title \"\"" > vis.plt
    #gnuplot -p -c vis.plt
fi

