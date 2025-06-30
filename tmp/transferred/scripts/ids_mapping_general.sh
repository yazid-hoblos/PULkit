#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 10

grep -P "\tCDS\t" $1/*/genomic.gff |
awk '{
  match($0, /locus_tag=([^;]+)/, a);
  match($0, /protein_id=([^;]+)/, b);
  if (a[1] && b[1]) print a[1] "\t" b[1];
}' > ids_mapping

