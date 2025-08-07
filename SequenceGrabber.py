
from collections import defaultdict
from typing import OrderedDict
import pandas as pd
from os import path
import numpy as np
import time
import logging
import argparse
import re


class SequenceGrabber():

    name = "grabbersequences"

    def __init__(self, chrom_file, name:str ="None"):
        self.seqs = defaultdict(lambda: str)
        if name != "None":
            self.name = name
        self.chrom_dict = self.__parse_chromosome_file(chrom_file)

    def grab_sequences(self, loci_file, cutoff=300):
        logging.info(f"Grabbing sequences with {cutoff}bp on either side")
        print(type(cutoff))
        df = pd.read_excel(loci_file, header=0)
        # Remove Unnamed columns sometimes imported from excel file
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        # Remove where the whole row is NA when spacers are added in
        df = df.dropna(how='all')
        try:
            for i, row in df.iterrows():
                logging.debug(f"Working on {row['SNP']}")
                chrom_sequence = self.chrom_dict.get(str(row["Chromosome"]))
                position = int(row["Position"])
                low_cutoff = position-cutoff
                logging.debug(f"Low cutoff is: {low_cutoff}")
                high_cutoff = position+cutoff
                logging.debug(f"High cutoff is: {high_cutoff}")
                ref = row["REF"]
                alt = row["ALT"]
                self.array = np.fromiter(chrom_sequence.keys(), dtype=int)
                keys = self.__find_nearest(
                    low_cutoff=low_cutoff, high_cutoff=high_cutoff)
                seq = ''
                for key in keys:
                    value = chrom_sequence[key]
                    if key == keys[0]:
                        value = value[low_cutoff-key:]
                    elif key == keys[-1]:
                        value = value[:1+(high_cutoff-key)]
                    seq += value
                seq_list = list(seq)
                logging.debug(
                    f"The list of bases is {len(seq_list)} bases long")
                seq_list[(len(seq_list)-1)//2:(len(seq_list)+2) //
                         2] = f"[{ref}/{alt}]"
                final_seq = "".join(seq_list)
                df.loc[i, 'Sequence'] = final_seq
                self.seqs[row["SNP"]] = final_seq
            self.df = df
        except KeyError as e:
            print(
                f"I'm looking for a header in your loci file called '{e.args[0]}' and I didn't find it. Please add it!")
            exit(1)
        return df

    def __parse_chromosome_file(self, chrom_file):
        chrom_dict = defaultdict(OrderedDict)
        try:
            with open(path.join(chrom_file), 'r') as infile:
                lines = infile.readlines()
                n = 1
                chrom = None
                for line in lines:
                    if line.startswith(">"):
                        try:
                            chrom = int(line.strip(">chr").rstrip())
                            chrom = line.strip(">chr").rstrip()
                            n = 1
                        except ValueError:
                            word = 'chromosome'
                            pattern = rf"{word}\s+(\S+)"
                            match = re.search(pattern, line)
                            if match:
                                chrom = match.group(1)
                                logging.debug(chrom)
                            else:
                                chrom = line.strip(">chr").rstrip()
                            n = 1
                    if not line.startswith(">"):
                        l = line.strip()
                        chrom_dict[chrom][n] = l
                        n += len(l)
            try:
                scaffold_phrases = ("scaffold", "Scaf")
                scaffold_count = sum(1 for key in chrom_dict if any(
                    phrase in key for phrase in scaffold_phrases))
                logging.debug(f"ID'ed Scaffolds: {scaffold_count}")
            except TypeError:
                logging.debug(f'I ran into a type error but idk why')
            chromosome_phrases = ("chromosome", 'chr', "Ca")
            chromosome_count = [key for key in chrom_dict if any(
                phrase in key for phrase in chromosome_phrases)]
            if chromosome_count == 0:
                logging.warning(f"I didn't find any chromosomes.. The tags should start with chromosome or chr")
            logging.debug(f"ID'ed Chromosome Count: {len(chromosome_count)}")
            logging.debug(f"ID'ed Chromosomes:")
            for i in chromosome_count:
                last_item = list(chrom_dict.get(i).items())[-1]
                first_item = list(chrom_dict.get(i).items())[0]
                logging.debug(
                    f"Chromosome {i} starts at: {first_item[0]+len(first_item[1])} and goes to: {last_item[0]+len(last_item[1])} bases long")
            logging.debug(f"ID'ed Chromosomes:")
            return chrom_dict
        except FileNotFoundError as e:
            print(
                f"I can't find the sequence file in {path.join(self.sequence_directory,chrom_file)}")
            exit(1)

    def __find_nearest(self, high_cutoff, low_cutoff):
        # To find the high index, we subtract each element in the array by the high_cutoff
        # and get the absolute value of that. Then we get the minimum of that. That gets
        # us really close to where we need to be and gives us the index of that value
        high_idx = np.abs(self.array - high_cutoff).argmin()
        logging.debug(f"High idx: {high_idx}")
        # Check to see if the high index is greater than the high cutoff. If it is
        # set the index -1
        if self.array[high_idx] > high_cutoff:
            logging.debug(
                f"High_idx: {self.array[high_idx]} is greater than high_cutoff")
            high_idx -= 1
        # Same function here, just in the reverse direction. Find the low cutoff
        # value
        low_idx = np.abs(self.array - low_cutoff).argmin()
        logging.debug(f"Low idx: {low_idx}")
        if self.array[low_idx] > low_cutoff:
            logging.debug(f"low_idx is greater than low_cutoff")
            low_idx -= 1
        return self.array[low_idx:high_idx+1]

    def __insert_snp(self, seq, snp):
        row = self.vcf_df.loc[self.vcf_df['ID'] == snp]
        if row.empty:
            logging.debug(f"SNP for marker {snp} not found in VCF file!")
            ref = "unk"
            alt = "unk"
        else:
            ref = row.iloc[0]['REF']
            alt = row.iloc[0]['ALT']
        seq_list = list(seq)
        seq_list[(len(seq_list)-1)//2:(len(seq_list)+2)//2] = f"[{ref}/{alt}]"
        seq = "".join(seq_list)
        return seq


    def write_out(self, outdir, write_excel=False):
        if not write_excel:
            with open(path.join(outdir, f"{time.strftime('%Y%m%d')}_{self.name}_sequences.txt"), 'w') as outfile:
                for k, v in self.seqs.items():
                    outfile.write(f">{k}\n{v}\n\n")
        else:
            writer = pd.ExcelWriter(
                path.join(outdir, f"{time.strftime('%Y%m%d')}_{self.name}_sequences.xlsx"))
            self.df.to_excel(writer)
            writer.close()


def readArgs():
    """Reads the arguments provided by the user."""
    desc = """
    This program parses through an excel file of SNPs and grabs the position and surrounding sequence (with the REF/ALT) inserted into the sequence
    """
    p = argparse.ArgumentParser(description=desc, prog="SequenceGrabber.py")
    # Required arguments
    p.add_argument("-o", "--outdir", required=True,
                   help="Specifies the directory where the data is written.")
    p.add_argument("-x", "--write_excel", default=False,
                   help="Output the excel dataframe", action="store_true")
    p.add_argument("-i", "--snpfile", required=True,
                   help="The Input SNP file.. Should be xlsx with required colums: SNP, REF, ALT, Chromosome, Position")
    p.add_argument("-g", "--genomefile", required=True,
                   help="The Input Genome file.. Should be a fasta file")
    
    # Optional arguments
    p.add_argument("-n", "--name", required=False, default="None",
                   help="The Name ")
    p.add_argument("-b", "--bp_cutoff", required=False, default=300, type=int,
                   help="Specify how many bp on either side of the SNP you want -- defaults to 300")    
    p.add_argument("-v", "--verbose", default=False,
                   help="Output verbose", action="store_true")
    args = p.parse_args()
    return args


def main():
    args = readArgs()
    if args.verbose:
        logging.basicConfig(format='%(asctime)s [%(funcName)s] %(levelname)s - %(message)s',
                            datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG)
        logging.debug("Outputting verbose")
    else:
        logging.basicConfig(format='%(asctime)s [%(funcName)s] %(levelname)s - %(message)s',
                            datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
    OUTDIR = args.outdir
    SNPFILE = args.snpfile
    GENOMEFILE = args.genomefile
    NAME = args.name
    WRITE_EXCEL = args.write_excel
    cutoff = args.bp_cutoff

    s = SequenceGrabber(chrom_file=GENOMEFILE, name=NAME)
    df = s.new_grab_seqs(SNPFILE, cutoff=cutoff)
    s.write_out(OUTDIR, write_excel=WRITE_EXCEL)


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    main()
