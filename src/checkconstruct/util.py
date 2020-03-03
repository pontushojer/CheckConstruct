import logging

logger = logging.getLogger(__name__)


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
