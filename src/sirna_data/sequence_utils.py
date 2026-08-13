"""Nucleotide sequence utilities: alphabet conversion (`to_rna`/`to_dna`) and
real DNA-template transcription (`transcribe_template_to_mrna`).

`to_rna`/`to_dna` are a plain notation swap, NOT transcription: they change
'T' <-> 'U' in place, on the same strand, in the same order. That's the
right operation when a sequence is already in its final, correct
orientation -- e.g. an siRNA guide, or an mRNA window -- and some downstream
tool just expects the other alphabet's spelling of that same sequence.

Real transcription is a different operation. RNA polymerase reads a DNA
*template* strand 3'->5' and synthesizes a complementary RNA 5'->3' -- so
the resulting mRNA is the template strand's REVERSE COMPLEMENT, with U in
place of T, not just a letter swap. (If you instead already have the
*coding/sense* strand -- the one written the same direction as, and with
the same sequence as, the mRNA except for T vs. U -- then a plain `to_rna`
swap with no complementing IS correct, because the coding strand and the
mRNA are the same sequence by convention.) Use `transcribe_template_to_mrna`
when you specifically have a DNA template strand and need the mRNA it
produces; use `to_rna`/`to_dna` when you already have the right strand and
just need the other alphabet's spelling of it.
"""
from __future__ import annotations

_DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def to_rna(seq: str) -> str:
    """Uppercase `seq` and replace every 'T' with 'U' (DNA -> RNA alphabet).
    A no-op (aside from uppercasing) on a sequence that's already RNA.

    This is a notation swap only -- NOT transcription. It assumes `seq` is
    already the strand you want (e.g. a coding/sense strand, whose sequence
    already equals the mRNA's except for T vs. U), and just re-spells it in
    the RNA alphabet in the same order. If you actually have a DNA
    *template* strand and want the mRNA it transcribes to, use
    `transcribe_template_to_mrna` instead -- that reverse-complements,
    which this function deliberately does not do.
    """
    return seq.upper().replace("T", "U")


def to_dna(seq: str) -> str:
    """Uppercase `seq` and replace every 'U' with 'T' (RNA -> DNA alphabet).
    A no-op (aside from uppercasing) on a sequence that's already DNA.

    Notation swap only, same caveat as `to_rna`: no complementing, so `seq`
    should already be the strand/orientation you want.
    """
    return seq.upper().replace("U", "T")


def transcribe_template_to_mrna(template_strand: str) -> str:
    """Transcribe a DNA template (antisense) strand into the mRNA sequence
    it produces: the reverse complement of `template_strand`, in the RNA
    alphabet.

    `template_strand` is expected written 5'->3' (the standard convention
    for writing out any nucleotide sequence). RNA polymerase reads the
    template strand 3'->5' and synthesizes RNA 5'->3', so the resulting
    mRNA -- also returned 5'->3' -- is `template_strand`'s reverse
    complement with U in place of T, not merely a T->U letter swap.

    Unlike `to_rna`/`to_dna`, this performs a real biological operation on
    strictly A/C/G/T input, so it validates rather than silently passing
    through anything unrecognized: lowercase is accepted and normalized,
    but ambiguity codes (e.g. 'N'), RNA characters ('U'), or any other
    character raise `ValueError` -- there's no well-defined complement for
    them here, and guessing would risk silently producing a wrong mRNA
    sequence.
    """
    template = template_strand.upper()
    try:
        complement = "".join(_DNA_COMPLEMENT[base] for base in template)
    except KeyError as e:
        raise ValueError(
            f"template_strand must be a DNA sequence (A/C/G/T only); "
            f"found unexpected character {e.args[0]!r}"
        ) from e
    return to_rna(complement[::-1])
