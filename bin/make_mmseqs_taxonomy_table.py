#!/usr/bin/env python
from utils import save_excel
from Bio import SeqIO
import pandas as pd
import argparse
import logging
import errno
import glob
import gzip
import sys
import re
import os

logger = logging.getLogger()

def parser_args(args=None):
    Description = "Create table containing mapped taxonomy information."
    Epilog = """Example usage: python make_variants_long_table.py --mmseqs2_dir ./mmseqs2/ --contigs_dir ./contigs/"""
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-p",
        "--pattern_sample",
        type=str,
        default=None,
        help="Pattern to sort samples by their numeric features (default: None).",
    )
    parser.add_argument(
        "-mm",
        "--mmseqs2_dir",
        type=str,
        default="./mmseqs2",
        help="Directory containing output of mmseqs2 for each sample (default: './mmseqs2').",
    )
    parser.add_argument(
        "-cd",
        "--contigs_dir",
        type=str,
        default="./contigs",
        help="Directory containing assembled contigs for each sample (default: './contigs').",
    )
    parser.add_argument(
        "-as",
        "--aln_file_suffix",
        type=str,
        default="_tophit_aln",
        help="Suffix to trim off tophit alignment file name to obtain sample name (default: '_tophit_aln').",
    )
    parser.add_argument(
        "-rs",
        "--report_file_suffix",
        type=str,
        default="_tophit_report",
        help="Suffix to trim off tophit report file name to obtain sample name (default: '_tophit_report').",
    )
    parser.add_argument(
        "-mmp",
        "--mmseqs_file_prefix",
        type=str,
        default="",
        help="Prefix to trim off mmseqs result file names to obtain sample name (default: '').",
    )
    parser.add_argument(
        "-cs",
        "--contigs_file_suffix",
        type=str,
        default=".contigs.fa.gz",
        help="Suffix to trim off contigs file name to obtain sample name (default: '.contigs.fa.gz').",
    )
    parser.add_argument(
        "-cp",
        "--contigs_file_prefix",
        type=str,
        default="",
        help="Prefix to trim off contigs file name to obtain sample name. Similar to 'assembler-' (default: '').",
    )
    parser.add_argument(
        "-of",
        "--output_file",
        type=str,
        default="contigs.taxonomy.csv",
        help="Full path to output file (default: 'contigs.taxonomy.csv').",
    )
    parser.add_argument(
        "--excel",
        action="store_true",
        default=False,
        help="Create corresponding excel table."
    )
    parser.add_argument(
        "-db", "--database", type=str, help="Database used to map the taxonomy."
    )
    parser.add_argument(
        "-st", "--search_type", type=int, default = 2, help="Search type of mmseqs used to map the taxonomy (default: 2)."
    )
    return parser.parse_args(args)

def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

def get_file_dict(file_dir, file_suffix, file_prefix="", pattern=None):
    files = glob.glob(os.path.join(file_dir, f"*{file_suffix}"))
    if pattern:
        files.sort(key=lambda x: eval_sample_num(pattern=pattern, string=x))
    samples = [os.path.basename(x).removeprefix(f"{file_prefix}").removesuffix(f"{file_suffix}") for x in files]
    return dict(zip(samples, files))

# ".*[a-zA-Z]-[A-Z]?([0-9]+)[A-Z]?_S[0-9]+"
def eval_sample_num(pattern: str, string: str) -> int:
    match = re.search(pattern, string)
    num = None
    if match:
        num = match.group(1)
    if num and num.isdigit():
        return int(num)
    logger.error(f"Found {num} in {string} which is not an integer!")
    sys.exit(1)

def search_type_int_to_str(stype: int) -> str:
    stype_dict = {
        0: "auto",
        1: "amino acid",
        2: "translated",
        3: "nucleotide",
        4: "translated nucleotide alignment"
    }

    return stype_dict[stype]

def process_fasta(fasta) -> pd.DataFrame:
    fasta_df = pd.concat( [ pd.Series({'qseqid' : x.id,
                'kmer_cov' : re.search('multi=(.+?) ' if x.id.startswith("k") else 'cov_(.+?)$', x.description).group(1), # megahit or spades kmer coverage
                'seq' : x.seq.__str__()}) for x in SeqIO.parse(fasta, 'fasta') ], axis=1).T
    fasta_df["kmer_cov"] = fasta_df["kmer_cov"].astype(float)
    fasta_df = fasta_df.loc[fasta_df.groupby("qseqid")["kmer_cov"].idxmax()]
    return fasta_df

def contig_fasta_to_table(fasta_file) -> pd.DataFrame:
    try:
        if fasta_file.endswith(".gz"):
            with gzip.open(fasta_file, "rt") as fasta_unziped:
                fasta_df = process_fasta(fasta_unziped)
        else:
            fasta_df = process_fasta(fasta_file)
        return fasta_df.dropna()
    except Exception as e:
        logger.error(f'Failed to parse fasta from {fasta_file}. Error: {e}')
        print(e)
        sys.exit(1)

def tophit_aln_to_table(aln_file):
    aln_header = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    try:
        aln_df = pd.read_table(aln_file, header=None, names=aln_header)
        return aln_df.dropna()
    except Exception as e:
        logger.warning(f'{aln_file} is empty.Error: {e}')
        pass

def tophit_report_to_table(report_file):
    rep_header = ['sseqid', 'numbers', 'uniq_cov', 'uniq_cov_aln', 'average_ident', 'taxid', 'rank', 'name', 'taxonomy']
    try:
        rep_df = pd.read_table(report_file, header=None, names=rep_header)
        return rep_df.dropna()
    except Exception as e:
        print(f'{report_file} is empty. Error: {e}')
        pass

def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    allowed_search_types = [2, 3, 4]
    if args.search_type not in allowed_search_types:
        logger.error(
            f"Invalid option '--search_type {args.search_type}'. Valid options: " | ", ".join(allowed_search_types)
        )
        sys.exit(1)

    aln_files = get_file_dict(args.mmseqs2_dir, args.aln_file_suffix, args.mmseqs_file_prefix, args.pattern_sample)
    report_files = get_file_dict(args.mmseqs2_dir, args.report_file_suffix, args.mmseqs_file_prefix, args.pattern_sample)
    contigs_files = get_file_dict(args.contigs_dir, args.contigs_file_suffix, args.contigs_file_prefix, args.pattern_sample)

    if set(aln_files) != set(report_files):
        logger.error(
            f"Number of tophit alignment ({len(aln_files)}) and tophit report ({len(report_files)}) files do not match!"
        )
        sys.exit(1)
    else: 
        if set(aln_files) != set(contigs_files):
                logger.error(
                    f"Number of tophit alignment ({len(aln_files)}) and contigs ({len(contigs_files)}) files do not match!"
                )
                sys.exit(1)

    sample_tables = []

    for sample in aln_files:
        tophit_aln_table = tophit_aln_to_table(aln_files[sample])

        if not tophit_aln_table.empty:
            tophit_report_table = tophit_report_to_table(report_files[sample])
            contigs_table = contig_fasta_to_table(contigs_files[sample])

            merged_table = pd.DataFrame(data=tophit_aln_table)
            merged_table.insert(0, "sample", sample)
            merged_table = pd.merge(merged_table, tophit_report_table, on="sseqid")
            merged_table = pd.merge(merged_table, contigs_table, on="qseqid")
            merged_table["search_type"] = search_type_int_to_str(args.search_type)

            merged_table["database"] = args.database

            sample_tables.append(merged_table)
        
    if sample_tables:
        
        cols_order: list[str] = [
            'sample', 'qseqid', 'contig_length', 'kmer_cov', 'seq', 'sseqid', 'pident',
            'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
            'evalue', 'bitscore', 'numbers', 'uniq_cov', 'uniq_cov_aln',
            'average_ident', 'taxid', 'rank', 'name', 'taxonomy', 'search_type', 'database'
        ]

        taxonomic_order: list[str] = [
            'life', 'domain', 'kingdom', 'phylum', 'class', 'order',
            'family', 'genus', 'species'
        ]

        samples: list[str] = [sample for sample in aln_files]

        merged_tables = pd.concat(sample_tables)

        merged_tables['contig_length'] = merged_tables['seq'].str.len().values
        ## Added tax order for future sort
        merged_tables['rank'] = pd.Categorical(merged_tables['rank'], categories=taxonomic_order, ordered=True)
        merged_tables['sample'] = pd.Categorical(merged_tables['sample'], categories=samples, ordered=True)
        merged_tables = merged_tables.sort_values(["sample", "kmer_cov", "length", "evalue"], ascending=[True, False, False, True])

        merged_tables = merged_tables[cols_order]
        merged_tables.to_csv(args.output_file, index=False, encoding="utf-8-sig")
        if args.excel:
            excel_name = args.output_file.replace("csv", "xlsx")
            save_excel(merged_tables, outname=excel_name, skip_adjust=["seq"])


if __name__ == "__main__":
    sys.exit(main())

