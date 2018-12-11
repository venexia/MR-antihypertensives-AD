# Repurposing antihypertensive drugs for the prevention of Alzheimer’s disease: a Mendelian Randomization study

This respository contains the code to reproduce the analysis from the paper 'Repurposing antihypertensive drugs for the prevention of Alzheimer’s disease: a Mendelian Randomization study', which will be avaliable shortly from biorxiv.

## Abstract

Background: Evidence concerning the potential repurposing of antihypertensives for Alzheimer’s disease prevention is inconclusive. We used Mendelian randomization, which can be more robust to confounding by indication and patient characteristics, to investigate the effects of lowering systolic blood pressure (SBP), via different antihypertensive drug classes, on Alzheimer’s disease.
 
Methods: We used summary statistics from genome wide association studies of SBP (from UK Biobank) and Alzheimer’s disease (from the International Genomics of Alzheimer's Project) in a two-sample Mendelian randomization analysis. We identified single nucleotide polymorphisms (SNPs) that mimic the action of antihypertensive targets and estimated the effect of lowering SBP, via antihypertensive drug classes, on Alzheimer’s disease. We also report the effect of lowering SBP on Alzheimer’s disease by combining all drug targets and without consideration of the associated drugs.
 
Results: There was limited evidence that lowering SBP, via antihypertensive drug classes, affected Alzheimer’s disease risk. For example, calcium channel blockers had an odds ratio (OR) per 10mmHg lower SBP of 1.53 (95% confidence interval (CI): 0.94 to 2.49; p=0.09; SNPs=17). We also found limited evidence for an effect of lowering SBP on Alzheimer’s disease when combining all drug targets (OR per 10mmHg lower SBP: 1.14; 95%CI: 0.83 to 1.56; p=0.41; SNPs=59) and without consideration of the associated drug targets (OR per 10mmHg lower SBP: 1.04; 95%CI: 0.95 to 1.13; p=0.45; SNPs=153).
 
Conclusions: Lowering SBP itself is unlikely to affect risk of developing Alzheimer’s disease. Consequently, if specific antihypertensive drug classes do affect risk of Alzheimer’s disease, they are unlikely to do so via SBP.


## Using this code

To run this code, set your working directy in the file 'config.R' and run this file. All other files are called when required from this file. If you would like more information concerning the data setup, please contact venexia.walker@bristol.ac.uk. 

## Availability of data

The data used in this project are publicly avaliable and have been obtained from the following sources:

MR-Base - http://mrbase.org

DrugBank - https://www.drugbank.ca/releases/latest

GTEx - https://gtexportal.org/home/datasets

Neale lab systolic blood pressure GWAS - http://www.nealelab.is/uk-biobank/

Larsson et al instruments - https://www.bmj.com/content/359/bmj.j5375

Ostergaard et al instruments - https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001841

Gill et al instruments - https://www.biorxiv.org/content/early/2018/11/05/460543

## Funding statement

This work was supported by the Perros Trust and the Integrative Epidemiology Unit. The Integrative Epidemiology Unit is supported by the Medical Research Council and the University of Bristol [grant number MC_UU_00011/1, MC_UU_00011/3]. 
