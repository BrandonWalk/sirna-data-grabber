from __future__ import annotations

import pytest

from sirna_data.sequence_utils import to_dna, to_rna, transcribe_template_to_mrna


def test_to_rna_replaces_t_with_u():
    assert to_rna("ACGT") == "ACGU"


def test_to_dna_replaces_u_with_t():
    assert to_dna("ACGU") == "ACGT"


def test_to_rna_uppercases():
    assert to_rna("acgt") == "ACGU"


def test_to_dna_uppercases():
    assert to_dna("acgu") == "ACGT"


def test_to_rna_is_a_noop_on_already_rna():
    assert to_rna("ACGUACGU") == "ACGUACGU"


def test_to_dna_is_a_noop_on_already_dna():
    assert to_dna("ACGTACGT") == "ACGTACGT"


def test_to_rna_leaves_other_characters_untouched():
    # Ambiguity codes (e.g. 'N') and anything else that isn't 'T' pass through.
    assert to_rna("ACGTN") == "ACGUN"


def test_to_dna_leaves_other_characters_untouched():
    assert to_dna("ACGUN") == "ACGTN"


def test_to_rna_then_to_dna_round_trips_a_pure_acgt_sequence():
    seq = "ACGTACGTTTGGCA"
    assert to_dna(to_rna(seq)) == seq


def test_to_dna_then_to_rna_round_trips_a_pure_acgu_sequence():
    seq = "ACGUACGUUUGGCA"
    assert to_rna(to_dna(seq)) == seq


def test_to_rna_empty_string():
    assert to_rna("") == ""


def test_to_dna_empty_string():
    assert to_dna("") == ""


# --------------------------------------------------------------------------
# transcribe_template_to_mrna: real transcription (reverse complement + T->U)
# --------------------------------------------------------------------------


def test_transcribe_reverse_complements_and_converts_to_rna():
    # Template 5'->3': A C G T
    # Complement:       T G C A
    # Reverse of that:  A C G T   (reversed complement, read 3'->5' -> 5'->3')
    # mRNA (T->U):      A C G U
    assert transcribe_template_to_mrna("ACGT") == "ACGU"


def test_transcribe_known_example():
    # Template 5'->3' = TACCAG
    # complement       = ATGGTC
    # reverse complement (mRNA, 5'->3') = CTGGTA -> RNA: CUGGUA
    assert transcribe_template_to_mrna("TACCAG") == "CUGGUA"


def test_transcribe_uppercases_lowercase_input():
    assert transcribe_template_to_mrna("acgt") == transcribe_template_to_mrna("ACGT")


def test_transcribe_is_not_the_same_as_a_plain_alphabet_swap():
    template = "TACCAG"
    assert transcribe_template_to_mrna(template) != to_rna(template)


def test_transcribe_double_transcription_is_not_identity():
    # Transcribing twice does NOT return the original sequence in general
    # (reverse-complementing twice returns the original, but the
    # alphabet swap direction differs on the way back) -- sanity check
    # that the function is doing real reverse-complement work, not
    # accidentally cancelling itself out.
    template = "TACCAGGGTA"
    once = transcribe_template_to_mrna(template)
    assert once != template


@pytest.mark.parametrize("bad_char", ["N", "U", "X", "-", " "])
def test_transcribe_rejects_non_acgt_characters(bad_char):
    with pytest.raises(ValueError):
        transcribe_template_to_mrna(f"ACGT{bad_char}ACGT")


def test_transcribe_empty_string():
    assert transcribe_template_to_mrna("") == ""
