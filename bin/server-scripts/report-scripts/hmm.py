import argparse
import os
import logging
import psutil
import pyhmmer

# Constants for your domain-specific filtering, update or remove if not needed
GT2_PREFIX = "GT2"
GT2_FAMILY_NAME = "GT2_Family"

def run_hmmsearch(hmm_file, faa_file, output_file, e_value_threshold=1e-4, coverage_threshold=0.35, cpu=4):
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    available_memory = psutil.virtual_memory().available
    target_size = os.stat(faa_file).st_size
    results = []

    with pyhmmer.plan7.HMMFile(hmm_file) as hmm_files:
        with pyhmmer.easel.SequenceFile(faa_file, digital=True) as seqs:
            # Decide whether to load all sequences at once or stream
            targets = seqs.read_block() if target_size < available_memory * 0.1 else seqs

            logging.info(f"Running hmmsearch on {faa_file} with {hmm_file}...")
            for hits in pyhmmer.hmmsearch(hmm_files, targets, cpus=cpu, domE=e_value_threshold):
                for hit in hits:
                    for domain in hit.domains.included:
                        alignment = domain.alignment
                        coverage = (alignment.hmm_to - alignment.hmm_from + 1) / alignment.hmm_length
                        hmm_name = alignment.hmm_name.decode('utf-8')
                        if GT2_PREFIX in hmm_name:
                            hmm_name = GT2_FAMILY_NAME
                        i_evalue = domain.i_evalue

                        if i_evalue < e_value_threshold and coverage > coverage_threshold:
                            results.append([
                                hmm_name,
                                alignment.hmm_length,
                                alignment.target_name.decode('utf-8'),
                                alignment.target_length,
                                i_evalue,
                                alignment.hmm_from,
                                alignment.hmm_to,
                                alignment.target_from,
                                alignment.target_to,
                                coverage,
                                os.path.basename(hmm_file).split(".")[0]
                            ])

    logging.info(f"Search completed. Found {len(results)} hits passing filters.")

    # Write results to TSV
    header = [
        "hmm_name", "hmm_length", "target_name", "target_length", "i_evalue",
        "hmm_from", "hmm_to", "target_from", "target_to", "coverage", "hmm_file"
    ]
    with open(output_file, 'w') as out:
        out.write("\t".join(header) + "\n")
        for row in results:
            out.write("\t".join(map(str, row)) + "\n")

    logging.info(f"Results written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Run hmmsearch with PyHMMER, filter by coverage and e-value.")
    parser.add_argument("hmm_file", help="Path to HMM database file (e.g. STP.hmm)")
    parser.add_argument("faa_file", help="Path to protein FASTA file")
    parser.add_argument("-o", "--output", default="results.tsv", help="Output TSV file name (default: results.tsv)")
    parser.add_argument("--evalue", type=float, default=1e-4, help="E-value threshold (default: 1e-4)")
    parser.add_argument("--coverage", type=float, default=0.35, help="Coverage threshold (default: 0.35)")
    parser.add_argument("--cpu", type=int, default=4, help="Number of CPUs to use (default: 4)")
    args = parser.parse_args()

    run_hmmsearch(args.hmm_file, args.faa_file, args.output, args.evalue, args.coverage, args.cpu)

if __name__ == "__main__":
    main()

