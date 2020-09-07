# Drug-Target-Interaction
- Generally, to formulate a cure for a disease requires lot of wet-lab experiments
- This can cost 1-2 Billion Dollars per disease and upto 1-2 years
- I trained my model on 50,000 such interactions
  - Positive Interactions
    - When a drug binds to a protein
  - Negative Interaction
    - When a drug does not bind to a protein
- Given a Protein and a Drug, My model predicts if a drug is a binder to a Protein or not.
- Hence saving a lot of time and money

# Web interface
## Single Ligand Prediction
- Input:
  - PDB File
  - Fasta File
  - Ligand File
- Output
  - Is the ligand a binder to the pocket or not

## Entire DrugBank Prediction
- Input:
  - PDB File
  - Fasta File
- Output
  - DTI on entrire DrugBank Database

## Approved Drugs Prediction
- Input:
  - PDB File
  - Fasta File
- Output
  - DTI On Approved Drugs from DrugBank

## To execute:
- change pwd variable inside flaskblog/drug_target_interaction/routes.py to your current directory
- python run.py


<br>
<br>
<br>
- Web Interface template taken from https://github.com/CoreyMSchafer/code_snippets/tree/master/Python/Flask_Blog
