import sys

import bakta.constants as bc
import bakta.config as cfg
import bakta.utils as bu
from bakta.features.r_rna import R_RNA_PROFILES
from bakta.features.ori_arch import _parse_attributes, _parse_orics


def test_archaeal_argument_selects_archaeal_profile(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['bakta', '--organism', 'archaea', 'genome.fna'])
    assert bu.parse_arguments().organism == bc.ORGANISM_ARCHAEA


def test_archaeal_rrna_models_are_domain_specific():
    models = R_RNA_PROFILES[bc.ORGANISM_ARCHAEA]['models']
    assert 'RF01959' in models
    assert 'RF02540' in models
    assert 'RF00177' not in models


def test_archaeal_recoding_profile_includes_pyrrolysine():
    assert bc.ORGANISM_PROFILES[bc.ORGANISM_ARCHAEA]['recoding_codons']['TAG'] == 'pyrrolysine'


def test_archaeal_profile_ignores_gram_option():
    assert cfg.normalize_gram(bc.ORGANISM_ARCHAEA, bc.GRAM_POSITIVE) == bc.GRAM_UNKNOWN
    assert cfg.normalize_gram(bc.ORGANISM_BACTERIA, bc.GRAM_POSITIVE) == bc.GRAM_POSITIVE


def test_ori_finder_arch_gff_attributes_are_preserved():
    assert _parse_attributes("type='Type I (RIP-adjacent)';RIP=CDC6") == {'type': 'Type I (RIP-adjacent)', 'RIP': 'CDC6'}


def test_ori_finder_arch_gff_is_mapped_to_bakta_feature(tmp_path):
    result = tmp_path / 'ori.gff3'
    result.write_text("seq1\tOri-Finder-Arch\trep_origin\t10\t20\t42.0\t.\t0\ttype='Type I';RIP=CDC6\n")
    orics = _parse_orics(result, {'seq1': {'id': 'seq1', 'nt': 'A' * 100, 'length': 100}})
    assert orics[0]['type'] == bc.FEATURE_ORIC
    assert orics[0]['score'] == 42.0
    assert orics[0]['attributes']['RIP'] == 'CDC6'
