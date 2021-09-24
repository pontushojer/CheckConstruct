"""
Genereate sets of handles that don't from dimers with eachother.
"""
import logging
import itertools
from collections import defaultdict, namedtuple, Counter
from multiprocessing import Pool
from functools import partial

import primer3
from tqdm import tqdm

from checkconstruct.util import Primer, print_stats


logger = logging.getLogger(__name__)

summary = Counter()


def parse_tsv(iterator):
    Sequnce = namedtuple("Sequence", ["name", "seq"])
    for line in iterator:
        name, sequence = line.strip().split("\t")
        yield Sequnce(name, sequence)


def parse_fasta(iterator):
    Sequnce = namedtuple("Sequence", ["name", "seq"])
    while True:
        try:
            name = next(iterator)
            seq = next(iterator)
        except StopIteration:
            break

        yield Sequnce(name.replace(">", "").strip(), seq.strip())


def parse_primers(
    filename, struct_configs, hairpin_threshold=-3000, homodimer_threshold=-6000
):
    with open(filename, "r") as file:
        if filename.endswith("tsv"):
            parser = parse_tsv
        elif filename.endswith(".fa") or filename.endswith(".fasta"):
            parser = parse_fasta
        else:
            print("Error: wrong file format")
            exit(1)

        for primer in tqdm(parser(file), desc="Parsing file"):
            summary["Handles"] += 1
            yield Primer(name=primer.name, seq=primer.seq)


def main(args):
    struct_configs = {
        "mv_conc": args.mv_conc,
        "dv_conc": args.dv_conc,
        "dna_conc": args.oligo_conc,
    }

    primers = list(parse_primers(
        args.input_file,
        struct_configs,
        args.hairpin_threshold,
        args.homodimer_threshold,
    ))

    primer_combos = defaultdict(set)
    pairs = itertools.combinations(primers, 2)
    if args.threads == 1:
        primer_combos = get_combos(pairs, struct_configs, args.heterodimer_threshold)
    else:
        with Pool(processes=args.threads) as pool:
            func = partial(get_combo, struct_configs=struct_configs, heterodimer_threshold=args.heterodimer_threshold)
            for p1, p2 in tqdm(pool.imap_unordered(func, pairs, chunksize=10)):
                if p1 is not None and p2 is not None:
                    primer_combos[p1].add(p2)
                    primer_combos[p2].add(p1)

    # TODO Fix this as it is both slow and not working properly.
    # Find sets of primers in which all pairs do not form dimers.
    sets = set()
    for primer, neighbours in tqdm(primer_combos.items(), desc="Building sets"):
        for neighbour in neighbours:
            shared = {primer, neighbour}
            neighbour_neighbours = primer_combos.get(neighbour, set())
            for nn in neighbours & neighbour_neighbours:
                nn_neigbours = primer_combos.get(nn, set())
                if shared.issubset(nn_neigbours):
                    shared.add(nn)

            sets.add(frozenset(shared))

    summary["Unique sets"] = len(sets)
    for all_primers in sets:
        sorted_primers = sorted(all_primers)
        print(",".join([p.name for p in sorted_primers]), end="\t")
        print(",".join([p.seq for p in sorted_primers]))

    print_stats(summary)


def get_combo(pair, struct_configs, heterodimer_threshold):
    if forms_dimer(*pair, struct_configs, heterodimer_threshold):
        return None, None
    return pair


def get_combos(pairs, struct_configs, heterodimer_threshold):
    combos = defaultdict(set)
    for p1, p2 in tqdm(pairs):
        if forms_dimer(p1, p2, struct_configs, heterodimer_threshold):
            continue
        combos[p1].add(p2)
        combos[p2].add(p1)
    return combos


def forms_dimer(primer1, primer2, struct_configs, heterodimer_threshold):
    heterodimer = primer3.calcHeterodimer(primer1.seq, primer2.seq, **struct_configs)

    if heterodimer.structure_found and heterodimer.dg < heterodimer_threshold:
        return True

    heterodimer2 = primer3.calcHeterodimer(
        primer1.seq, primer2.revcomp, **struct_configs
    )

    if heterodimer2.structure_found and heterodimer2.dg < heterodimer_threshold:
        return True
    return False


def add_arguments(parser):
    parser.add_argument(
        "input_file",
        help="Input TSV file with primers. First column contain primer name and second primer "
        "sequence.",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Nr of threads.")
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
