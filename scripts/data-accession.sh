#!bin/bash/

curl -o dbCAN-PUL_v5 https://aca.unl.edu/static/DBCAN-PUL/dbCAN-PUL_v5_05-10-2023.xlsx
curl -o dbCAN-PUL.substrate.mapping https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.substrate.mapping.xls
# converted to CSV

for PUL in /env/cns/proj/agc/scratch_microscope/Data/DBCAN/899/dbCAN-PUL/*;do echo -n "$(basename $PUL)"$','; cut -f3,9 $PUL/cgc.gff | cut -d ';' -f1 | tr '\n' ',';echo;done > pul_data


