#！/bin/bash

#ref 1000G download from https://github.com/mancusolab/ma-focus

conda activate ma-focus

# Import PrediXcan db as Focus db
focus import IBD_Expr_EUR_driven.db predixcan --tissue Colon --name EUR --assay rnaseq --output IBD_Expr_EUR_driven_focus --use-ens-id --from-gencode --verbose

# Run Focus
focus finemap /Example_Data/EUR_IBD_perchr/chr${i}.focus.sumstats.gz /Focus_1000G/1000GP3_multiPop_allelesAligned/EUR/1000G.EUR.QC.allelesAligned.${i} /Example_Data/IBD_Expr_EUR_driven_focus.db --tissue Colon --chr ${i} --out ./Expr/EUR_IBD_${i} --locations 38:EUR