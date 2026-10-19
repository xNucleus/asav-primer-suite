"""
Unit-Tests fuer die wissenschaftliche Kernlogik (asav_core).
Ausfuehren mit:  python -m pytest -q   (oder ./venv/Scripts/python.exe -m pytest -q)
"""

import pytest

from asav_core import (
    aggregate_carrier_probability,
    assess_primer_snp,
    calculate_metrics,
    cross_dimer_warning,
    exon_target_rel,
    genomic_to_rel,
    is_snp_relevant,
    maf_from_populations,
    parse_hgvsg_to_genomic,
    primer_pair_snp_load,
    primer_warnings,
    recommend_traffic_light,
    rel_to_genomic,
    secondary_structure_warnings,
    snp_biology,
    snp_carrier_probability,
    snp_is_pathogenic,
    snp_penalty,
    snp_presence,
    snps_in_footprint,
    target_well_placed,
)

PLUS = {"start": 1000, "end": 2000, "strand": 1}
MINUS = {"start": 1000, "end": 2000, "strand": -1}


# ----------------------------- Koordinaten -----------------------------
def test_rel_to_genomic_plus_strand():
    assert rel_to_genomic(PLUS, 0, 20) == (1000, 1019)
    assert rel_to_genomic(PLUS, 10, 20) == (1010, 1029)


def test_rel_to_genomic_minus_strand():
    # 5'-Ende des Minus-Strang-Gens liegt am genomischen end.
    assert rel_to_genomic(MINUS, 0, 20) == (1981, 2000)
    assert rel_to_genomic(MINUS, 10, 20) == (1971, 1990)


@pytest.mark.parametrize("meta", [PLUS, MINUS])
def test_coordinate_roundtrip(meta):
    # rel -> genomic (erste Base) -> rel muss identisch sein, fuer beide Straenge.
    for rel in (0, 5, 123, 999):
        g_start, g_end = rel_to_genomic(meta, rel, 20)
        first_base = g_start if meta["strand"] == 1 else g_end
        assert genomic_to_rel(meta, first_base) == rel


# ----------------------------- Metriken -----------------------------
def test_gc_content_and_case_insensitivity():
    tm_upper, gc_upper = calculate_metrics("ACGTACGTACGTACGTACGT")
    tm_lower, gc_lower = calculate_metrics("acgtacgtacgtacgtacgt")
    assert gc_upper == 50.0
    assert (tm_upper, gc_upper) == (tm_lower, gc_lower)  # Gross-/Kleinschreibung egal


def test_tm_is_nearest_neighbor_not_wallace():
    # Bei langen Primern (hier 30 nt) ueberschaetzt die Wallace-Faustregel die
    # Tm massiv; das NN-Modell liefert einen deutlich niedrigeren, realistischen Wert.
    seq = "ACGTACGTACGTACGTACGTACGTACGTAC"  # 30 nt, 50 % GC
    tm, _ = calculate_metrics(seq)
    wallace = 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))
    assert wallace == 90          # Faustregel ueberschaetzt grob
    assert 55 < tm < 80           # NN liefert plausiblen PCR-Bereich
    assert wallace - tm > 5       # NN deutlich unter der Faustregel


def test_empty_sequence():
    assert calculate_metrics("") == (0.0, 0.0)


# ----------------------------- Primer-Qualitaet -----------------------------
def test_missing_gc_clamp_flagged():
    assert any("GC-Clamp" in w for w in primer_warnings("ACGTACGTACGTACGTACGT"))


def test_polyx_run_flagged():
    assert any("Mononukleotid" in w for w in primer_warnings("AAAAACGTACGTACGTACGC"))


def test_clean_primer_has_no_gc_or_clamp_warning():
    # 50 % GC, endet auf G, kein langer Run.
    warns = primer_warnings("ATCGATCGATGCATCGATCG")
    assert not any("GC-Clamp" in w for w in warns)
    assert not any("außerhalb" in w for w in warns)


# ----------------------------- SNP-Biologie -----------------------------
def test_snp_biology_pathogenic_wins():
    cls, _, _ = snp_biology({"clinical_significance": ["pathogenic"],
                             "consequence_type": "missense_variant"})
    assert cls.startswith("PATHOGEN")


def test_snp_biology_high_impact():
    cls, _, _ = snp_biology({"consequence_type": "stop_gained"})
    assert cls == "HIGH IMPACT"


def test_snp_biology_moderate_and_benign_and_unknown():
    assert snp_biology({"consequence_type": "missense_variant"})[0] == "MODERATE"
    assert snp_biology({"clinical_significance": ["benign"],
                        "consequence_type": "intron_variant"})[0].startswith("BENIGN")
    assert snp_biology({"consequence_type": "intron_variant"})[0].startswith("MODIFIER")


# ----------------------------- Primer-SNP-Bewertung -----------------------------
def test_snp_at_3prime_end_is_high_risk_forward():
    # Forward-Primer rel 100..119, 3'-Ende genomisch 1119.
    snp = {"start": 1119, "consequence_type": "missense_variant"}
    res = assess_primer_snp(snp, PLUS, primer_rel_start=100, p_len=20,
                            crit_zone=5, is_reverse=False)
    assert res["dist"] == 0
    assert res["risk"] == "HIGH RISK"


def test_snp_at_5prime_end_is_low_risk_forward():
    snp = {"start": 1100, "consequence_type": "intron_variant"}  # rel 100, weit vom 3'-Ende
    res = assess_primer_snp(snp, PLUS, primer_rel_start=100, p_len=20,
                            crit_zone=5, is_reverse=False)
    assert res["dist"] == 19
    assert res["risk"] == "LOW RISK"


def test_reverse_primer_3prime_distance():
    # Reverse-Primer: 3'-Ende an der niedrigsten rel-Position (rel 100).
    snp = {"start": 1100, "consequence_type": "missense_variant"}  # rel 100
    res = assess_primer_snp(snp, PLUS, primer_rel_start=100, p_len=20,
                            crit_zone=5, is_reverse=True)
    assert res["dist"] == 0
    assert res["risk"] == "HIGH RISK"


# ----------------------------- Frequenz (MAF) -----------------------------
def test_maf_from_populations_picks_minor_allele():
    pops = [
        {"population": "1000GENOMES:phase_3:ALL", "frequency": 0.543},
        {"population": "1000GENOMES:phase_3:ALL", "frequency": 0.457},
        {"population": "1000GENOMES:phase_3:EUR", "frequency": 0.30},  # andere Pop -> ignoriert
    ]
    assert maf_from_populations(pops, "1000GENOMES:phase_3:ALL") == 0.457


def test_maf_from_populations_none_without_data():
    assert maf_from_populations([]) is None
    assert maf_from_populations(None) is None
    assert maf_from_populations([{"population": "1000GENOMES:phase_3:ALL",
                                  "frequency": 1.0}]) is None  # nur ein Allel


def test_snp_presence_carrier_probability():
    # Haeufig: P(>=1 Minor-Allel) = 1 - (1-0.4)^2 = 0.64.
    carrier, desc = snp_presence({"maf": 0.4})
    assert abs(carrier - 0.64) < 1e-9
    assert "häufig" in desc
    # Selten: kleine Wahrscheinlichkeit, als "1 von N" formuliert.
    carrier_r, desc_r = snp_presence({"maf": 0.0005})
    assert 0 < carrier_r < 0.01
    assert "1 von" in desc_r


def test_snp_presence_no_data_never_claims_absence():
    carrier, desc = snp_presence({})
    assert carrier is None
    assert "nicht ausgeschlossen" in desc   # nie "kommt nicht vor"


def test_is_snp_relevant_only_common():
    # Ranking-Relevanz = HAEUFIG (MAF >= Schwelle). Seltene – auch pathogene –
    # zaehlen NICHT (werden separat geflaggt, dominieren aber das Ranking nicht).
    assert is_snp_relevant({"maf": 0.2}, 0.01) is True
    assert is_snp_relevant({"maf": 0.001}, 0.01) is False
    assert is_snp_relevant({"maf": 0.001, "clinical_significance": ["pathogenic"]},
                           0.01) is False
    assert is_snp_relevant({}, 0.01) is False


def test_snp_is_pathogenic_flag():
    assert snp_is_pathogenic({"clinical_significance": ["pathogenic"]}) is True
    assert snp_is_pathogenic({"clinical_significance": ["likely_pathogenic"]}) is True
    assert snp_is_pathogenic({"clinical_significance": ["benign"]}) is False
    assert snp_is_pathogenic({}) is False


# ----------------------------- SNP-Last (Ranking) -----------------------------
def test_snp_penalty_3prime_worse_than_5prime():
    # Beide haeufig (MAF 0.2) -> relevant; 3'-Ende muss teurer sein als 5'-Bereich.
    snp_3 = {"start": 1119, "consequence_type": "missense_variant", "maf": 0.2}
    snp_5 = {"start": 1100, "consequence_type": "intron_variant", "maf": 0.2}
    pen_3 = snp_penalty(snp_3, PLUS, 100, 20, 5, 0.01)
    pen_5 = snp_penalty(snp_5, PLUS, 100, 20, 5, 0.01)
    assert pen_3 == 6.0   # 5 (3'-Zone) + 1 (missense)
    assert pen_5 == 1.0   # 1 (5'-Bereich) + 0 (intron)


def test_rare_nonclinical_snp_has_zero_penalty():
    # Genau der Fall des dichten dbSNP-Rauschens: selten + harmlos -> 0.
    rare = {"start": 1119, "consequence_type": "missense_variant", "maf": 0.0005}
    assert snp_penalty(rare, PLUS, 100, 20, 5, 0.01) == 0.0


def test_rare_pathogenic_not_in_ranking():
    # Selten + pathogen -> NICHT ranking-relevant -> Penalty 0 (wird separat geflaggt).
    rare_patho = {"start": 1119, "clinical_significance": ["pathogenic"],
                  "consequence_type": "missense_variant", "maf": 0.0005}
    assert snp_penalty(rare_patho, PLUS, 100, 20, 5, 0.01) == 0.0


def test_common_pathogenic_gets_surcharge():
    # Haeufig UND pathogen am 3'-Ende -> relevant + Aufschlag: 5 (3') + 4 (pathogen).
    common_patho = {"start": 1119, "clinical_significance": ["pathogenic"],
                    "consequence_type": "missense_variant", "maf": 0.2}
    assert snp_penalty(common_patho, PLUS, 100, 20, 5, 0.01) == 9.0


def test_snps_in_footprint_filters_by_position():
    # Forward-Primer rel 100..119 -> genomisch 1100..1119 (PLUS-Strang).
    all_snps = [{"start": 1110}, {"start": 1200}, {"start": 1105}, {"start": 1099}]
    inside = snps_in_footprint(all_snps, PLUS, 100, 20)
    starts = sorted(s["start"] for s in inside)
    assert starts == [1105, 1110]   # 1099 und 1200 liegen ausserhalb


def test_primer_pair_load_counts_rare_but_does_not_penalize():
    # 1 haeufiger (zaehlt) + 1 seltener (nur Info) -> Penalty nur vom haeufigen,
    # aber total = 2 (nichts verschwiegen).
    common = {"start": 1119, "consequence_type": "missense_variant", "maf": 0.2}  # 6.0
    rare = {"start": 1100, "consequence_type": "intron_variant", "maf": 0.0001}   # 0.0
    pen, relevant, total = primer_pair_snp_load(
        [common, rare], [], PLUS, 100, 100, 20, 5, 0.01
    )
    assert pen == 6.0
    assert relevant == 1
    assert total == 2


def test_snp_free_pair_has_zero_load():
    assert primer_pair_snp_load([], [], PLUS, 100, 200, 20, 5, 0.01) == (0, 0, 0)


# ============================ NEU: primer3-Sekundaerstruktur ============================
def test_secondary_structure_flags_palindrome_but_not_clean():
    bad = secondary_structure_warnings("ATCGATCGATGCATCGATCG", annealing_tm=55.0)
    good = secondary_structure_warnings("CCTTCAGCTCAAGACTTAAC", annealing_tm=55.0)
    assert bad != []      # echte Hairpin/Dimer-Struktur erkannt
    assert good == []     # sauberer Primer wird NICHT ueberflaggt


def test_cross_dimer_detects_complementary_pair():
    w = cross_dimer_warning("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT")
    assert any("Primer-Dimer" in x for x in w)


# ============================ NEU: Traegerwahrscheinlichkeit / Aggregat (J) ============================
def test_snp_carrier_probability():
    assert abs(snp_carrier_probability({"maf": 0.4}) - 0.64) < 1e-9
    assert snp_carrier_probability({}) is None


def test_aggregate_carrier_probability_combines_and_counts_missing():
    p, n_with, n_without = aggregate_carrier_probability(
        [{"maf": 0.1}, {"maf": 0.1}, {"maf": None}])
    assert abs(p - 0.3439) < 1e-6     # 1 - (0.81 * 0.81)
    assert (n_with, n_without) == (2, 1)


def test_aggregate_none_when_no_frequency_data():
    p, n_with, n_without = aggregate_carrier_probability([{"maf": None}, {}])
    assert p is None and n_with == 0 and n_without == 2


# ============================ NEU: Ampel (H) ============================
def test_traffic_light_green_yellow_red():
    base = {"tm_f": 60, "tm_r": 60, "delta": 0.5, "penalty": 0,
            "snp_relevant": 0, "snp_3prime": 0}
    assert recommend_traffic_light(base, 60, 5)[0] == "grün"
    assert recommend_traffic_light({**base, "snp_relevant": 1}, 60, 5)[0] == "gelb"
    assert recommend_traffic_light({**base, "snp_patho": 1}, 60, 5)[0] == "gelb"
    assert recommend_traffic_light({**base, "notable": 1}, 60, 5)[0] == "gelb"  # nennenswerter Hinweis
    assert recommend_traffic_light({**base, "penalty": 3}, 60, 5)[0] == "grün"  # nur kosmetisch (GC-Clamp) -> bleibt grün
    assert recommend_traffic_light({**base, "snp_3prime": 1}, 60, 5)[0] == "rot"
    assert recommend_traffic_light({**base, "serious": 1}, 60, 5)[0] == "rot"   # Hairpin/starker Dimer
    assert recommend_traffic_light({**base, "spec": "offtarget"}, 60, 5)[0] == "rot"


# ============================ NEU: Ziel-Aufloesung (A) ============================
def test_exon_target_rel_plus_and_minus_strand():
    plus = {"start": 1000, "end": 2000, "strand": 1,
            "exons": [{"start": 1000, "end": 1100}, {"start": 1500, "end": 1600}]}
    assert exon_target_rel(plus, 1) == 50      # Exon1 +Strang = niedrigste Koord., Mitte 1050
    minus = {"start": 1000, "end": 2000, "strand": -1,
             "exons": [{"start": 1000, "end": 1100}, {"start": 1500, "end": 1600}]}
    assert exon_target_rel(minus, 1) == 450    # Exon1 -Strang = hoechste Koord., Mitte 1550
    assert exon_target_rel(plus, 9) is None     # Exon ausserhalb


def test_parse_hgvsg_to_genomic():
    assert parse_hgvsg_to_genomic("NC_000017.11:g.7676154G>C") == ("17", 7676154)
    assert parse_hgvsg_to_genomic("NC_000023.11:g.100A>T") == ("X", 100)
    assert parse_hgvsg_to_genomic("kaputt") == (None, None)


# ============================ NEU: Amplikon-Platzierung (B) ============================
def test_target_well_placed():
    assert target_well_placed(110, 40, 160, 20, 50) is True   # Ziel mittig, >=50bp Abstand
    assert target_well_placed(65, 40, 160, 20, 50) is False   # zu nah am Forward-3'-Ende
    assert target_well_placed(None, 40, 160, 20, 50) is True  # kein Ziel -> immer ok
