import logging
import contextlib
import sys
import statistics
import dataclasses

from Bio.Seq import Seq
import primer3

logger = logging.getLogger(__name__)


@dataclasses.dataclass(order=False, frozen=True)
class Primer:
    name: str
    seq: str
    conc: float = None

    TM_CONFIGS = None
    STRUCT_CONFIGS = None

    def __lt__(self, other):
        return self.seq > other.seq

    def __eq__(self, other):
        return self.name == other.name and self.seq == other.seq

    def __len__(self):
        return len(self.seq)

    @property
    def length(self):
        return len(self)

    @property
    def tm(self, **kwargs) -> float:
        configs = Primer.TM_CONFIGS.copy()
        if self.conc is not None:
            configs["dna_conc"] = self.conc
        return round(primer3.calcTm(self.seq, **configs), 1)

    @property
    def revcomp(self) -> str:
        return str(Seq(self.seq).reverse_complement())

    @property
    def gc(self) -> float:
        return round(self._gc(self.seq), 1)

    @property
    def max_dg(self, **kwargs) -> float:
        if self.length > 60:
            return 0

        configs = Primer.STRUCT_CONFIGS.copy()
        if self.conc is not None:
            configs["dna_conc"] = self.conc

        return round(primer3.calcHeterodimer(self.seq, self.revcomp, **configs).dg, 1)

    def has_repeats(self, threshold: int = 3) -> bool:
        current, count, max_count = (None, 1, 1)
        for c in self.seq:
            if c == current:
                count += 1
                max_count = max(count, max_count)
                if max_count >= threshold:
                    return True
            else:
                current = c
                count = 1
        return False

    def has_gc_clamp(self) -> bool:
        return self.seq[-1] in "GC" or self.seq[-2] in "GC"

    def has_uniform_gc(self, window: int = 5, threshold=0.4) -> bool:
        gc_range = [
            self._gc(self.seq[i : i + window]) for i in range(0, self.length - window)
        ]
        # Coeff. of variation < threadshols
        return statistics.stdev(gc_range) / statistics.mean(gc_range) < threshold

    @classmethod
    def _gc(cls, sequence) -> float:
        return 100 * sum(1 for s in sequence if s in ["G", "C"]) / len(sequence)

    def to_dict(self):
        return {
            "Name": self.name,
            "Sequence": self.seq,
            "Length (bp)": self.length,
            "Tm (°C)": self.tm,
            "GC (%)": self.gc,
            "Max dG (cal/mol)": self.max_dg,
        }


def simple_tm(dna_string):
    """
    This formula is valid for oligos <14 bases and assumes that the reaction is carried out
    in the presence of 50mM monovalent cations.
    :param dna_string:
    :return:
    """
    if len(dna_string) > 14:
        logger.warning("Sequence longer than the recommended 14 bp")

    count_at = 0
    count_gc = 0
    for nt in dna_string.upper():
        if nt == "A" or nt == "T":
            count_at += 1
        elif nt == "G" or nt == "C":
            count_gc += 1
        else:
            logger.warning("Sequence contains non-standard letter {nt}")
            return None

    return 2 * count_at + 4 * count_gc


@contextlib.contextmanager
def smart_open(filename=None, mode="w"):
    # From https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    if filename and filename is not None and filename != "-":
        fh = open(filename, mode)
    else:
        fh = sys.stdout if mode == "w" else sys.stdin

    try:
        yield fh
    finally:
        if fh is not sys.stdout or fh is not sys.stdin:
            fh.close()


def print_stats(stats):
    k_width = max(map(len, stats.keys()))
    v_width = max(map(len, map(str, stats.values())))
    tot_width = k_width + v_width + 1

    print("=" * tot_width, file=sys.stderr)
    print("Summary stats", file=sys.stderr)
    print("-" * tot_width, file=sys.stderr)
    for k, v in stats.items():
        print(f"{k:<{k_width}} {v:>{v_width}}", file=sys.stderr)
    print("=" * tot_width, file=sys.stderr)
