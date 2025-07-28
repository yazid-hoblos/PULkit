from Bio import SeqIO
import sys 

gff_lines = ["##gff-version 3"]
start = 1
gap = 1000  # fake gene spacing

file = sys.argv[1] if len(sys.argv) > 1 else "xylan_pul_proteins.faa"
output = sys.argv[2] if len(sys.argv) > 2 else "synthetic.gff"

for record in SeqIO.parse(file, "fasta"):
    contig = record.id.split("|")[0]
    gene_id = record.id
    length = len(record.seq) * 3  # approximate CDS length in bp
    end = start + length - 1
    strand = "+"  # or random.choice(["+", "-"]) 

    gff_line = f"{contig}\tdbcan\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={gene_id}"
    gff_lines.append(gff_line)
    start = end + gap  # space next gene

with open(output, "w") as f:
    f.write("\n".join(gff_lines))

