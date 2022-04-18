#!/bin/sh
source activate tf

case=$1
control=$2
dir=~/Desktop/Halogen_Health_work/AUC_based_estimator_of_predictive_GWAS_model
python SNP_simulation_model.py --PRS_file ~/Desktop/Halogen_Health_work/AUC_based_estimator_of_predictive_GWAS_model/PGS000147_formatted.txt --outpath ~/Desktop/Halogen_Health_work/AUC_based_estimator_of_predictive_GWAS_model --Ncases $case --Ncontrol $control

for i in $(seq 1 $(($(awk < ${dir}/Simulated_genotype_data.txt '{print NF}'| head -n1)-2)));do
ncase=$((2 * case))
ncontrl=$((2*control))

tail -n+2 Simulated_genotype_data.txt | awk -v i=$((i+2)) -v ncase=${ncase}  -v ncontr=${ncontrl} '{if($2==1){case+=$i} else {contr+=$i}} END {print ncase/2,ncontr/2,(case*(ncontr-contr))/(contr*(ncase-case)) }'; done

#tail -n+2 Simulated_genotype_data.txt | awk -v i=$((i+2)) -v ncase=${ncase}  -v ncontr=${ncontrl} '{if($2==1){case+=$i} else {contr+=$i}} END {print ncase/2,ncontr/2,contr/ncontr,(case+contr)/(ncontr+ncase) }';done

#tail -n+2 Simulated_genotype_data.txt | awk -v i=$((i+2)) -v ncase=${ncase}  -v ncontr=${ncontrl} '{if($2==1){case+=$i} else {contr+=$i}} END {print case,"\t",ncase-case,"\n",contr,"\t",ncontr-contr,"\n","OR=","\t", (case*(ncontr-contr))/(contr*(ncase-case)) }'; done
