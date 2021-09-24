"""
Filter FASTA with primer sequences
"""
from collections import Counter

from tqdm import tqdm
from Bio import SeqIO
import primer3

from checkconstruct.util import smart_open, Primer, print_stats


def add_arguments(parser):
    parser.add_argument("input", help="FASTA with primer sequences to filter.")
    parser.add_argument(
        "-o",
        "--output",
        help="Write FASTA with filtered sequences to file instead of stdout.",
        default="-",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Upper limit for repeats. Default: %(default)s",
    )
    parser.add_argument(
        "--gc-upper",
        type=int,
        default=60,
        help="Upper limit for GC. Default: %(default)s",
    )
    parser.add_argument(
        "--gc-lower",
        type=int,
        default=40,
        help="Lower limit for GC. Default: %(default)s",
    )
    parser.add_argument(
        "--gc-clamp",
        type=bool,
        default=True,
        help="Presence of GC-clamp (one of the two 3′ terminal bases must be G or C). Default: %(default)s",
    )
    parser.add_argument(
        "--gc-uniformity",
        type=float,
        default=0.4,
        help="Threshold for variation (CV) of GC along the sequence. Default: %(default)s",
    )
    parser.add_argument(
        "--hairpin_threshold",
        default=-3000,
        type=int,
        help="deltaG (cal/mol) threshold for hairpin. Default: %(default)s.",
    )
    parser.add_argument(
        "--homodimer_threshold",
        default=-6000,
        type=int,
        help="deltaG (cal/mol) threshold for homodimer. Default: %(default)s.",
    )
    parser.add_argument(
        "-c",
        "--dna_conc",
        default=200,
        type=float,
        help="DNA conc (nM) for hairpin/homodimer calculation. Default: %(default)s",
    )
    parser.add_argument(
        "-m",
        "--mv_conc",
        default=50,
        type=float,
        help="Monovalent cation conc (mM) for hairpin/homodimer calculation. Default: %(default)s",
    )
    parser.add_argument(
        "-d",
        "--dv_conc",
        default=1.5,
        type=float,
        help="Divalent cation conc (mM) for hairpin/homodimer calculation. Default: %(default)s",
    )


# TODO Add option to blast seqeunces


def main(args):
    run_filter(
        input=args.input,
        output=args.output,
        repeats_threshold=args.repeats,
        gc_upper_threshold=args.gc_upper,
        gc_lower_threshold=args.gc_lower,
        gc_clamp=args.gc_clamp,
        gc_uniformity_threshold=args.gc_uniformity,
        hairpin_threshold=args.hairpin_threshold,
        homodimer_threshold=args.homodimer_threshold,
        dna_conc=args.dna_conc,
        mv_conc=args.mv_conc,
        dv_conc=args.dv_conc,
    )


def run_filter(
    input: str,
    output: str,
    repeats_threshold: int,
    gc_upper_threshold: float,
    gc_lower_threshold: float,
    gc_clamp: bool,
    gc_uniformity_threshold: float,
    hairpin_threshold: int,
    homodimer_threshold: int,
    dna_conc: float,
    mv_conc: float,
    dv_conc: float,
):

    summary = Counter()

    struct_configs = {
        "mv_conc": mv_conc,
        "dv_conc": dv_conc,
        "dna_conc": dna_conc,
    }

    with smart_open(input, mode="r") as infile, smart_open(output, mode="w") as outfile:
        for entry in tqdm(SeqIO.parse(infile, "fasta"), desc="Parsing primers"):
            summary["Primers reads"] += 1
            primer = Primer(name=entry.name, seq=str(entry.seq))

            if not (gc_lower_threshold <= primer.gc <= gc_upper_threshold):
                summary["Out of GC-range"] += 1
                continue

            if primer.has_repeats(threshold=repeats_threshold):
                summary["Repeats"] += 1
                continue

            if gc_clamp and not primer.has_gc_clamp():
                summary["No GC-clamp"] += 1
                continue

            if primer.has_uniform_gc(window=5, threshold=gc_uniformity_threshold):
                summary["Non-uniform GC rate"] += 1
                continue

            hairpin = primer3.calcHairpin(primer.seq, **struct_configs)

            if hairpin.structure_found and hairpin.dg < hairpin_threshold:
                summary["Hairpin"] += 1
                continue

            homodimer = primer3.calcHomodimer(primer.seq, **struct_configs)

            if homodimer.structure_found and homodimer.dg < homodimer_threshold:
                summary["Homodimer"] += 1
                continue

            summary["Primers written"] += 1
            SeqIO.write(entry, outfile, format="fasta")

    print_stats(summary)
