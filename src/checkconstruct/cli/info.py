"""
Info
"""
import primer3
import pandas as pd
import logging
import itertools
from Bio.Seq import Seq
import sys


logger = logging.getLogger(__name__)

import checkconstruct.util as checkconstruct


def add_info(primers, info_function, keyword_args=dict(), column_name="Unknown", decimals=1):
    """
    Flexible function to add information to pandas dataframe.
    :param primers: Dataframe containing sequences for primers/oligos under column 'Sequence".
    :param info_function: Method/function to be applied on sequence to extract information
    :param keyword_args: Dict containing keywords and values relating to info_function
    :param column_name: Name of column to bee added.
    :param decimals: Number of decimals to round to.
    """
    info_list = []
    for index, row in primers.iterrows():
        info = info_function(row["Sequence"], **keyword_args)
        info_list.append(round(info, decimals))

    primers[column_name] = info_list


def analyze_homostructures(primers, homodimer_threshold=-5000, hairpin_threshold=-1000, keyword_args=dict()):

    structures = []

    for name, primer in primers.iterrows():
        if len(primer["Sequence"]) > 60:
            logger.warning(f'Primer {name} is too long ({len(primer["Sequence"])}bp > 60 bp) for analysis of homostructures. ')
            continue

        hairpin = primer3.calcHairpin(primer["Sequence"], **keyword_args)

        if hairpin.structure_found and hairpin.dg < hairpin_threshold:
            structures.append({
                "Type": "Hairpin",
                "dG": round(hairpin.dg,1),
                "Tm": round(hairpin.tm,1),
                "Name": name
            })

        homodimer = primer3.calcHomodimer(primer["Sequence"], **keyword_args)

        if homodimer.structure_found and homodimer.dg < homodimer_threshold:
            structures.append({
                "Type": "Self-dimer",
                "dG (cal/mol)": round(homodimer.dg,1),
                "Tm (C)": round(homodimer.tm,1),
                "Name": name
            })

    return structures


def analyze_heterostructures(primers, heterodimer_threshold=-5000, keyword_args=dict()):

    structures = []

    names = primers.index.tolist()
    for name, name2 in itertools.combinations(names, 2):

        seq = primers.at[name, "Sequence"]
        seq2 = primers.at[name2, "Sequence"]
        # At least one primer needs to be shorter than 60 bp.
        if len(seq) > 60 and len(seq2) > 60:
            logging.warning(f'Primers {name} and {name2} are too long (> 60 bp) for analysis of heterostructures. ')
            continue

        heterodimer = primer3.calcHeterodimer(seq, seq2, **keyword_args)

        if heterodimer.structure_found and heterodimer.dg < heterodimer_threshold:
            structures.append({
                "Type": "Heterodimer",
                "dG (cal/mol)": round(heterodimer.dg,1),
                "Tm (C)": round(heterodimer.tm,1),
                "Name": f'{name} and {name2}'
            })

        seq2 = Seq(seq2)
        seq2_revcomp = seq2.reverse_complement()

        heterodimer = primer3.calcHeterodimer(seq, str(seq2_revcomp), **keyword_args)

        if heterodimer.structure_found and heterodimer.dg < heterodimer_threshold:
            structures.append({
                "Type": "Heterodimer",
                "dG (cal/mol)": round(heterodimer.dg,1),
                "Tm (C)": round(heterodimer.tm,1),
                "Name": f'{name} and revcomp({name2})'
            })

    return structures


def max_dg(sequence, mv_conc=None, dv_conc=None, dna_conc=None):
    seq = Seq(sequence)
    revcomp_seq = seq.reverse_complement()

    if len(sequence) > 60:
        return 0

    heterodimer = primer3.calcHeterodimer(str(seq), str(revcomp_seq), mv_conc=mv_conc, dv_conc=dv_conc,
                                          dna_conc=dna_conc)

    return heterodimer.dg


def main(args):
    print("SETTINGS USED:")
    print("-" * 30)
    for option, value in sorted(vars(args).items()):
        print(f"{option}: {value}")
    print("-"*30)

    primers = checkconstruct.import_primers(args.input_file)

    # Add size information
    add_info(primers, len, column_name="Size (bp)")

    # Add melting temperature information
    add_info(primers, primer3.calcTm, column_name="Tm (C)", keyword_args={'mv_conc': args.mv_conc,
                                                                          'dv_conc': args.dv_conc,
                                                                          'dntp_conc': args.dntp_conc,
                                                                          'dna_conc': args.oligo_conc,
                                                                          'salt_corrections_method': args.salt_corr_method})

    # Add info about maximum delta G for sequence (binding to complement).
    add_info(primers, max_dg, column_name="Max dG (cal/mol)", keyword_args={'mv_conc': args.mv_conc,
                                                                            'dv_conc': args.dv_conc,
                                                                            'dna_conc': args.oligo_conc})

    homostructures = analyze_homostructures(primers, homodimer_threshold=args.homodimer_threshold,
                                            hairpin_threshold=args.hairpin_threshold,
                                            keyword_args={'mv_conc': args.mv_conc,
                                                          'dv_conc': args.dv_conc,
                                                          'dna_conc': args.oligo_conc})

    heterostructures = analyze_heterostructures(primers, heterodimer_threshold=args.heterodimer_threshold,
                                                keyword_args={'mv_conc': args.mv_conc,
                                                              'dv_conc': args.dv_conc,
                                                              'dna_conc': args.oligo_conc})
    if homostructures or heterostructures:
        struct = homostructures + heterostructures
        struct_df = pd.DataFrame(data=struct, columns=["Name", "Type", "dG (cal/mol)", "Tm (C)"]).set_index("Name")

        with pd.option_context('display.max_rows', 200, 'display.max_columns', None):
            print(struct_df.to_string())
    else:
        logging.info("No homostructures or heterostructures detected.")

    print("-" * 10)
    print('RESULT')
    print("-" * 10)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', 30):
        print(primers.to_string())


def add_arguments(parser):
    parser.add_argument("input_file",
                        help="Input tsv file with primers. Frist column contain primer name and second primer sequence.")

    parser.add_argument('-c', '--oligo_conc',
                        help="Oligo conc (nM)",
                        default=200,
                        type=float)

    parser.add_argument('-m', '--mv_conc',
                        help="Monovalent cation conc (mM)",
                        default=50,
                        type=float)

    parser.add_argument('-d', '--dv_conc',
                        help="Divalent cation conc (mM)",
                        default=1.5,
                        type=float)

    parser.add_argument('-n', '--dntp_conc',
                        help="dNTPs conc (mM)",
                        default=0.2,
                        type=float)

    parser.add_argument('-s', '--salt_corr_method',
                        help="Salt correction method",
                        default='owczarzy',
                        choices=['schildkraut', 'owczarzy', 'santalucia'],
                        type=str)

    parser.add_argument('--hairpin_threshold',
                        help="deltaG (cal/mol) threshold for alerting about hairpin.",
                        default=-2000,
                        type=int)

    parser.add_argument('--homodimer_threshold',
                        help="deltaG (cal/mol) threshold for alerting about homodimer.",
                        default=-5000,
                        type=int)

    parser.add_argument('--heterodimer_threshold',
                        help="deltaG (cal/mol) threshold for alerting about heterodimer.",
                        default=-5000,
                        type=int)