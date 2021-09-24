"""
Analyce set of primers
"""
import primer3
import pandas as pd
import logging
import itertools

from checkconstruct.util import Primer

logger = logging.getLogger(__name__)


def analyze_homostructures(
    primers,
    structures,
    homodimer_threshold=-5000,
    hairpin_threshold=-1000,
    mv_conc=None,
    dv_conc=None,
):
    for primer in primers:
        checked_seq = primer.seq
        if len(primer) > 60:
            logger.warning(
                f"Primer {primer.name} is too long ({len(primer)}bp > 60 bp) for analysis of "
                f"homostructures. "
            )
            logger.info(f"Trimming primer to 60 bp by removing bases from 5'-end")
            checked_seq = primer.seq[-60:]

        keyword_args = {"mv_conc": mv_conc, "dv_conc": dv_conc, "dna_conc": primer.conc}

        hairpin = primer3.calcHairpin(checked_seq, **keyword_args)

        if hairpin.structure_found and hairpin.dg < hairpin_threshold:
            structures.add_struture(
                primer.name, "Hairpin", round(hairpin.dg, 1), round(hairpin.tm, 1)
            )

        homodimer = primer3.calcHomodimer(checked_seq, **keyword_args)

        if homodimer.structure_found and homodimer.dg < homodimer_threshold:
            structures.add_struture(
                primer.name,
                "Self-dimer",
                round(homodimer.dg, 1),
                round(homodimer.tm, 1),
            )


def analyze_heterostructures(
    primers, structures, heterodimer_threshold=-5000, mv_conc=None, dv_conc=None
):
    for primer1, primer2 in itertools.combinations(primers, 2):
        # At least one primer needs to be shorter than 60 bp.
        if len(primer1) > 60 and len(primer2) > 60:
            logging.warning(
                f"Primers {primer1.name} and {primer2.name} are too long (> 60 bp) for analysis of "
                f"heterostructures. "
            )
            continue

        keyword_args = {
            "mv_conc": mv_conc,
            "dv_conc": dv_conc,
            "dna_conc": max(primer1.conc, primer2.conc),
        }

        heterodimer = primer3.calcHeterodimer(primer1.seq, primer2.seq, **keyword_args)

        if heterodimer.structure_found and heterodimer.dg < heterodimer_threshold:
            structures.add_struture(
                f"{primer1.name} and {primer2.name}",
                "Heterodimer",
                round(heterodimer.dg, 1),
                round(heterodimer.tm, 1),
            )

        heterodimer2 = primer3.calcHeterodimer(
            primer1.seq, primer2.revcomp, **keyword_args
        )

        if heterodimer2.structure_found and heterodimer2.dg < heterodimer_threshold:
            structures.add_struture(
                f"{primer1.name} and {primer2.name}'",
                "Heterodimer",
                round(heterodimer2.dg, 1),
                round(heterodimer2.tm, 1),
            )


def import_primers(filename, dna_conc):
    primers = []
    with open(filename, "r") as file:
        for line in file:
            values = line.strip().split("\t")
            name, sequence = (None, None)
            conc = dna_conc
            if len(values) == 3:
                name, sequence, conc = values
                conc = float(conc)
            elif len(values) == 2:
                name, sequence = values
            else:
                logging.error("Wrong file format")
                continue

            primers.append(Primer(name=name, seq=sequence, conc=conc))

    return primers


class Structures:
    def __init__(self):
        self.list = list()
        self.columns = ["Name", "Type", "dG (cal/mol)", "Tm (°C)"]

    def add_struture(self, name, structure_type, delta_g, melt_temp):
        self.list.append(
            {
                "Name": name,
                "Type": structure_type,
                "dG (cal/mol)": delta_g,
                "Tm (°C)": melt_temp,
            }
        )

    def to_df(self):
        return pd.DataFrame(data=self.list, columns=self.columns).set_index("Name")

    def print(self, sort_by="Tm (°C)", ascending=False):
        with pd.option_context("display.max_rows", 200, "display.max_columns", None):
            print(self.to_df().sort_values(by=sort_by, ascending=ascending).to_string())


def main(args):
    Primer.TM_CONFIGS = {
        "mv_conc": args.mv_conc,
        "dv_conc": args.dv_conc,
        "dntp_conc": args.dntp_conc,
        "dna_conc": args.oligo_conc,
        "salt_corrections_method": args.salt_corr_method,
    }

    Primer.STRUCT_CONFIGS = {
        "mv_conc": args.mv_conc,
        "dv_conc": args.dv_conc,
        "dna_conc": args.oligo_conc,
    }

    primers = import_primers(args.input_file, dna_conc=args.oligo_conc)

    structures = Structures()

    analyze_homostructures(
        primers,
        structures,
        homodimer_threshold=args.homodimer_threshold,
        hairpin_threshold=args.hairpin_threshold,
        mv_conc=args.mv_conc,
        dv_conc=args.dv_conc,
    )
    if not args.skip_hetero:
        analyze_heterostructures(
            primers,
            structures,
            heterodimer_threshold=args.heterodimer_threshold,
            mv_conc=args.mv_conc,
            dv_conc=args.dv_conc,
        )

    if structures.list:
        structures.print()
    elif args.skip_hetero:
        logging.info("No homostructures detected.")
    else:
        logging.info("No homostructures or heterostructures detected.")

    print("-" * 10)
    print("RESULT")
    print("-" * 10)
    primers_df = pd.DataFrame.from_records([p.to_dict() for p in primers])
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.max_colwidth",
        30,
    ):
        print(primers_df.to_string())


def add_arguments(parser):
    parser.add_argument(
        "input_file",
        help="Input tsv file with primers. Frist column contain primer name and second primer "
        "sequence. There is also the option to add a third column with the concentration "
        "(in nM).",
    )

    parser.add_argument(
        "-c",
        "--oligo_conc",
        help="Oligo conc (nM). Can also be set for each primer individually, see input_file.",
        default=200,
        type=float,
    )

    parser.add_argument(
        "-m", "--mv_conc", help="Monovalent cation conc (mM)", default=50, type=float
    )

    parser.add_argument(
        "-d", "--dv_conc", help="Divalent cation conc (mM)", default=1.5, type=float
    )

    parser.add_argument(
        "-n", "--dntp_conc", help="dNTPs conc (mM)", default=0.2, type=float
    )

    parser.add_argument(
        "-s",
        "--salt_corr_method",
        help="Salt correction method",
        default="owczarzy",
        choices=["schildkraut", "owczarzy", "santalucia"],
        type=str,
    )

    parser.add_argument(
        "--hairpin_threshold",
        help="deltaG (cal/mol) threshold for alerting about hairpin.",
        default=-3000,
        type=int,
    )

    parser.add_argument(
        "--homodimer_threshold",
        help="deltaG (cal/mol) threshold for alerting about homodimer.",
        default=-6000,
        type=int,
    )

    parser.add_argument(
        "--heterodimer_threshold",
        help="deltaG (cal/mol) threshold for alerting about heterodimer.",
        default=-6000,
        type=int,
    )

    parser.add_argument(
        "--skip-hetero",
        default=False,
        action="store_true",
        help="Skip heterodimer analysis.",
    )
