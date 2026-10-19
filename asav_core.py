"""
asav_core – streamlit-freie Kernlogik der ASAV-Suite.

Hier liegt die wissenschaftlich relevante, testbare Logik:
Koordinaten-Umrechnung (strangkorrekt), Schmelztemperatur (Nearest-Neighbor),
Primer-Qualitaetskriterien und SNP-Bewertung. Bewusst ohne UI-Abhaengigkeiten,
damit sie per pytest geprueft werden kann (siehe test_core.py).
"""

import re

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

try:                       # primer3 ist die bevorzugte Engine; Fallback bleibt Tm_NN
    import primer3
    _HAS_PRIMER3 = True
except Exception:          # pragma: no cover
    primer3 = None
    _HAS_PRIMER3 = False

# ---- PCR-/Thermodynamik-Parameter (dokumentierte Annahmen) ----
NA_CONC = 50.0     # mM einwertige Kationen (Na+/K+)
MG_CONC = 1.5      # mM Mg2+
DNTP_CONC = 0.6    # mM dNTPs
PRIMER_CONC = 250  # nM Primer-Konzentration

GC_MIN, GC_MAX = 40.0, 60.0   # akzeptabler GC-Bereich (%)

# VEP-Konsequenzen nach Ensembl-IMPACT-Rating (Auszug der relevanten Klassen)
CONS_HIGH = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification",
}
CONS_MODERATE = {
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "regulatory_region_ablation",
}

# Gewichte fuer die SNP-Last im Kandidaten-Ranking (Position x Biologie).
# Position: ein SNP nahe dem 3'-Ende verhindert die Amplifikation am staerksten.
SNP_W_3PRIME = 5.0      # SNP in der kritischen 3'-Zone
SNP_W_MID = 2.0         # SNP im mittleren Primer-Bereich
SNP_W_5PRIME = 1.0      # SNP im 5'-Bereich (am wenigsten kritisch)
# Biologie: klinisch relevante Varianten zusaetzlich bestrafen (Aufschlag).
SNP_W_PATHOGENIC = 4.0  # pathogen / likely_pathogenic (ClinVar)
SNP_W_HIGH = 2.0        # High-Impact-Konsequenz (z. B. stop_gained)
SNP_W_MODERATE = 1.0    # Moderate-Impact-Konsequenz (z. B. missense)

# Schwelle, ab der eine Variante als "haeufig" gilt (Standard-Definition: 1 %).
MAF_COMMON_DEFAULT = 0.01
# Default-Population fuer die Frequenz. gnomAD GENOMES (gnomADg) deckt auch
# intronische Regionen ab und ist daher fuer genomweites Primer-Design besser
# geeignet als die Exome-Kohorte. Im UI umschaltbar (POPULATION_OPTIONS).
MAF_POPULATION = "gnomADg:ALL"
POPULATION_OPTIONS = {
    "Global (gnomAD)": "gnomADg:ALL",
    "Europäisch (NFE)": "gnomADg:nfe",
    "Afrikanisch": "gnomADg:afr",
    "Ostasiatisch": "gnomADg:eas",
    "Südasiatisch": "gnomADg:sas",
    "Amerikanisch (AMR)": "gnomADg:amr",
}

# Sekundaerstruktur-Schwellen (primer3), klinik-konservativ.
# Hairpin nur kritisch, wenn seine Tm die ANNEALING-Temperatur erreicht (sonst
# schmilzt er vorher und konkurriert nicht mit der Template-Bindung).
DEFAULT_ANNEALING_TM = 55.0  # °C (entspricht Ziel-Tm 60 - 5)
DIMER_DG_WARN = -7000.0      # cal/mol – ab hier moeglicher Dimer
DIMER_DG_STRONG = -9000.0    # cal/mol – starker Dimer (IDT-Richtwert ~-9 kcal/mol)


# ============================ Koordinaten ============================
def rel_to_genomic(meta, rel_start, length):
    """Relative Position (0-based, im 5'->3'-Genstrang) -> genomische 1-based
    Koordinaten (g_start <= g_end). Strangkorrekt.

    Ensembl liefert die Sequenz eines Minus-Strang-Gens bereits in 5'->3'-
    Orientierung des Gens; `start`/`end` aus dem Lookup sind aber stets
    genomisch aufsteigend. Diese Funktion ueberbrueckt beides.
    """
    if meta.get("strand", 1) == 1:
        g_start = meta["start"] + rel_start
        g_end = g_start + length - 1
    else:
        g_end = meta["end"] - rel_start
        g_start = g_end - length + 1
    return g_start, g_end


def genomic_to_rel(meta, g_pos):
    """Genomische 1-based Position -> relative Position im 5'->3'-Genstrang."""
    if g_pos is None:
        return 0
    if meta.get("strand", 1) == 1:
        return g_pos - meta["start"]
    return meta["end"] - g_pos


# ============================ Metriken ============================
def _p3_kwargs():
    """primer3-Bedingungen aus unseren PCR-Konstanten (mM bzw. nM)."""
    return dict(mv_conc=NA_CONC, dv_conc=MG_CONC, dntp_conc=DNTP_CONC, dna_conc=PRIMER_CONC)


def calculate_metrics(seq_str):
    """(Tm, GC%). Tm via primer3 (Industriestandard, SantaLucia-NN) mit unseren
    PCR-Bedingungen; Fallback auf Biopython Tm_NN, dann Tm_GC bei nicht-ACGT-Basen."""
    if not seq_str:
        return 0.0, 0.0
    s = seq_str.upper()
    gc = 100.0 * (s.count("G") + s.count("C")) / len(s)
    tm = None
    if _HAS_PRIMER3:
        try:
            tm = primer3.calc_tm(s, **_p3_kwargs())
        except Exception:
            tm = None
    if tm is None:
        try:
            tm = mt.Tm_NN(Seq(s), Na=NA_CONC, Mg=MG_CONC, dNTPs=DNTP_CONC,
                          dnac1=PRIMER_CONC, dnac2=PRIMER_CONC, saltcorr=7)
        except (ValueError, IndexError):
            tm = mt.Tm_GC(Seq(s), Na=NA_CONC)
    return round(float(tm), 1), round(gc, 1)


def primer_warnings(seq):
    """Standard-Primer-Qualitaetskriterien (dependency-frei). Liste von Hinweisen."""
    s = seq.upper()
    w = []
    if not s:
        return w
    gc = 100.0 * (s.count("G") + s.count("C")) / len(s)
    if not (GC_MIN <= gc <= GC_MAX):
        w.append(f"GC {gc:.0f}% außerhalb {GC_MIN:.0f}–{GC_MAX:.0f}%")
    if s[-1] not in "GC":
        w.append("kein GC-Clamp (3'-Base A/T)")
    if sum(1 for b in s[-5:] if b in "GC") >= 4:
        w.append("3'-Ende GC-überladen (>3/5)")
    if re.search(r"(A{5,}|T{5,}|G{5,}|C{5,})", s):
        w.append("Mononukleotid-Run ≥5")
    if re.search(r"(AT){4,}|(TA){4,}|(GC){4,}|(CG){4,}", s):
        w.append("Dinukleotid-Repeat")
    return w


def secondary_structure_warnings(seq, annealing_tm=DEFAULT_ANNEALING_TM):
    """primer3-Warnungen fuer Hairpin + Self-Dimer eines einzelnen Primers.
    Ein Hairpin wird nur gemeldet, wenn seine Tm die Annealing-Temperatur erreicht.
    Leere Liste, wenn primer3 fehlt (dann greifen nur die Heuristiken oben)."""
    if not _HAS_PRIMER3 or not seq:
        return []
    s = seq.upper()
    w = []
    try:
        h = primer3.calc_hairpin(s, **_p3_kwargs())
        if h.structure_found and h.tm >= annealing_tm:
            w.append(f"Hairpin (Tm {h.tm:.0f}°C ≥ Annealing)")
        d = primer3.calc_homodimer(s, **_p3_kwargs())
        if d.structure_found and d.dg <= DIMER_DG_WARN:
            tag = "stark" if d.dg <= DIMER_DG_STRONG else "moderat"
            w.append(f"Self-Dimer {tag} (ΔG {d.dg / 1000:.1f} kcal/mol)")
    except Exception:      # pragma: no cover
        pass
    return w


def cross_dimer_warning(fwd, rev):
    """primer3-Warnung fuer Heterodimer zwischen Forward- und Reverse-Primer."""
    if not _HAS_PRIMER3 or not fwd or not rev:
        return []
    try:
        d = primer3.calc_heterodimer(fwd.upper(), rev.upper(), **_p3_kwargs())
        if d.structure_found and d.dg <= DIMER_DG_WARN:
            tag = "stark" if d.dg <= DIMER_DG_STRONG else "moderat"
            return [f"Primer-Dimer F/R {tag} (ΔG {d.dg / 1000:.1f} kcal/mol)"]
    except Exception:      # pragma: no cover
        pass
    return []


# ============================ SNP-Bewertung ============================
def snp_biology(snp):
    """(Klasse, Farbe, Konsequenz) aus clinical_significance + consequence_type."""
    clin = [c.lower() for c in (snp.get("clinical_significance") or [])]
    cons = snp.get("consequence_type", "") or ""
    if any("pathogenic" in c for c in clin):           # pathogenic / likely_pathogenic
        return "PATHOGEN (ClinVar)", "#B91C1C", cons
    if cons in CONS_HIGH:
        return "HIGH IMPACT", "#EF4444", cons
    if cons in CONS_MODERATE:
        return "MODERATE", "#F59E0B", cons
    if any("benign" in c for c in clin):
        return "BENIGN (ClinVar)", "#10B981", cons
    return "MODIFIER/unbekannt", "#6B7280", cons


def assess_primer_snp(snp, meta, primer_rel_start, p_len, crit_zone, is_reverse=False):
    """Bewertet einen SNP relativ zum Primer.
    Die Position (Abstand vom 3'-Ende) bestimmt das PCR-Risiko;
    die Biologie (ClinVar/Konsequenz) ergaenzt den Kontext."""
    snp_rel = genomic_to_rel(meta, snp.get("start"))
    if not is_reverse:                                  # Forward: 3' an hoechster rel-Pos.
        dist = (primer_rel_start + p_len - 1) - snp_rel
    else:                                               # Reverse: 3' an niedrigster rel-Pos.
        dist = snp_rel - primer_rel_start
    dist = abs(dist)

    bio_class, bio_color, cons = snp_biology(snp)

    if dist <= crit_zone:
        risk, color = "HIGH RISK", "#FF4B4B"
        reason = f"Kritisch: {dist} bp vom 3'-Ende"
    elif dist <= crit_zone * 2:
        risk, color = "MED RISK", "#FFA500"
        reason = f"Warnung: {dist} bp vom 3'-Ende"
    else:
        risk, color = "LOW RISK", "#3FA34D"
        reason = f"Stabil: {dist} bp vom 3'-Ende (5'-Bereich)"

    style = (
        f"background-color:{color}; color:white; padding:4px 10px; "
        f"border-radius:8px; font-weight:bold; display:inline-block; margin-bottom:8px;"
    )
    return {
        "risk": risk, "reason": reason, "style": style,
        "bio_class": bio_class, "bio_color": bio_color, "consequence": cons,
        "snp_rel": snp_rel, "dist": dist,
    }


# ============================ Frequenz (MAF) ============================
def maf_from_populations(populations, pop=MAF_POPULATION):
    """Minor-Allele-Frequency (Frequenz des zweithaeufigsten Allels) der
    angegebenen Population, oder None falls keine Frequenzdaten vorliegen."""
    freqs = [p.get("frequency") for p in (populations or [])
             if p.get("population") == pop and p.get("frequency") is not None]
    if len(freqs) < 2:
        return None
    return sorted(freqs, reverse=True)[1]


def snp_frequency_label(snp):
    """Anzeige-Text fuer die Frequenz eines SNP (fuer den Bericht/das Audit)."""
    maf = snp.get("maf")
    if maf is None:
        return "keine Frequenzdaten (selten)"
    return f"MAF {maf * 100:.2f}%"


def snp_carrier_probability(snp):
    """P(zufaelliger Mensch traegt >=1 Minor-Allel) = 1 - (1 - MAF)^2
    (Hardy-Weinberg), oder None wenn keine Frequenzdaten vorliegen.
    Dieser Fall ist fuer die PCR relevant, weil schon ein betroffenes Allel das
    Priming stoeren kann (Allel-Dropout auch bei Heterozygoten)."""
    maf = snp.get("maf")
    if maf is None:
        return None
    return 1 - (1 - maf) ** 2


def snp_presence(snp):
    """Wie wahrscheinlich traegt ein zufaelliger Mensch diese Variante?
    Liefert (traegerwahrscheinlichkeit_oder_None, beschreibung)."""
    carrier = snp_carrier_probability(snp)
    if carrier is None:
        return None, "in der gewählten Population nicht beobachtet – sehr selten, " \
                     "Vorkommen aber nicht ausgeschlossen"

    maf = snp.get("maf")
    if maf >= 0.05:
        lvl = "häufig"
    elif maf >= 0.01:
        lvl = "gelegentlich"
    elif maf >= 0.001:
        lvl = "selten"
    else:
        lvl = "sehr selten"

    pct = carrier * 100
    if pct >= 1:
        freq = f"~{pct:.0f}% tragen sie"
    else:
        freq = f"~1 von {round(1 / carrier):,} tragen sie".replace(",", ".")
    return carrier, f"{lvl} · {freq}"


def aggregate_carrier_probability(snps):
    """(P_gesamt, n_mit_daten, n_ohne_daten) fuer eine Menge Varianten:
    Wahrscheinlichkeit, dass ein zufaelliger Patient MINDESTENS EINE davon traegt
    = 1 - prod(1 - carrier_i). Varianten ohne Frequenzdaten koennen nicht in die
    Zahl eingehen, werden aber SEPARAT gezaehlt (nicht als 0 verschwiegen)."""
    p_none = 1.0
    n_with = n_without = 0
    for s in snps:
        c = snp_carrier_probability(s)
        if c is None:
            n_without += 1
        else:
            p_none *= (1 - c)
            n_with += 1
    p_any = (1 - p_none) if n_with else None
    return p_any, n_with, n_without


def is_snp_relevant(snp, maf_threshold):
    """Relevant fuers RANKING = HAEUFIG (MAF >= Schwelle). Nur haeufige Varianten
    betreffen einen messbaren Anteil der Bevoelkerung und sollen das Ranking
    steuern. Seltene Varianten – auch pathogene – werden ANGEZEIGT und (wenn
    pathogen) separat geflaggt, dominieren aber das Ranking NICHT. Sonst waeren an
    dicht annotierten Genen wie TP53 (sehr viele ultra-seltene pathogene Eintraege)
    alle Primer 'rot' und das Werkzeug unbrauchbar."""
    maf = snp.get("maf")
    return maf is not None and maf >= maf_threshold


def snp_is_pathogenic(snp):
    """True, wenn ClinVar die Variante als (likely) pathogen fuehrt – unabhaengig
    von der Frequenz. Solche Stellen werden immer sichtbar gemacht (Sicherheit)."""
    clin = [c.lower() for c in (snp.get("clinical_significance") or [])]
    return any("pathogenic" in c for c in clin)


# ============================ SNP-Last (Ranking) ============================
def snp_penalty(snp, meta, primer_rel_start, p_len, crit_zone, maf_threshold,
                is_reverse=False):
    """Numerische Strafe fuer EINEN SNP unter einem Primer. Seltene, klinisch
    unauffaellige Varianten bekommen 0 (werden nur angezeigt). Relevante:
    Positionsgewicht (Naehe zum 3'-Ende) plus biologischer Aufschlag."""
    if not is_snp_relevant(snp, maf_threshold):
        return 0.0

    a = assess_primer_snp(snp, meta, primer_rel_start, p_len, crit_zone, is_reverse)
    if a["dist"] <= crit_zone:
        pen = SNP_W_3PRIME
    elif a["dist"] <= 2 * crit_zone:
        pen = SNP_W_MID
    else:
        pen = SNP_W_5PRIME

    bio = a["bio_class"]
    if bio.startswith("PATHOGEN"):
        pen += SNP_W_PATHOGENIC
    elif bio == "HIGH IMPACT":
        pen += SNP_W_HIGH
    elif bio == "MODERATE":
        pen += SNP_W_MODERATE
    return pen


def snps_in_footprint(all_snps, meta, rel_start, p_len):
    """Filtert vorab geholte SNPs auf das genomische Fenster EINES Primers.
    Ueberlappung ueber den GESAMTEN Variantenbereich (start..end), damit auch
    Indels/groessere Varianten erfasst werden, nicht nur Punkt-SNPs."""
    g0, g1 = rel_to_genomic(meta, rel_start, p_len)
    lo, hi = (g0, g1) if g0 <= g1 else (g1, g0)
    out = []
    for s in all_snps:
        s0 = s.get("start")
        if s0 is None:
            continue
        s1 = s.get("end", s0)
        if s1 is None:
            s1 = s0
        a, b = (s0, s1) if s0 <= s1 else (s1, s0)
        if a <= hi and b >= lo:        # Bereichs-Ueberlappung Variante x Primer
            out.append(s)
    return out


def primer_pair_snp_load(snps_f, snps_r, meta, f_rel, r_rel, p_len, crit_zone,
                         maf_threshold):
    """(Penalty, Anzahl_relevant, Anzahl_gesamt) fuer ein Primer-Paar.
    Penalty zaehlt nur relevante (haeufige/pathogene) Varianten; die
    Gesamtzahl wird zusaetzlich gemeldet, damit nichts verschwiegen wird."""
    pen = sum(snp_penalty(s, meta, f_rel, p_len, crit_zone, maf_threshold, False)
              for s in snps_f)
    pen += sum(snp_penalty(s, meta, r_rel, p_len, crit_zone, maf_threshold, True)
               for s in snps_r)
    relevant = (sum(1 for s in snps_f if is_snp_relevant(s, maf_threshold))
                + sum(1 for s in snps_r if is_snp_relevant(s, maf_threshold)))
    total = len(snps_f) + len(snps_r)
    return round(pen, 1), relevant, total


# ============================ Ampel-Empfehlung (H) ============================
def recommend_traffic_light(cand, target_tm, tm_tolerance):
    """Ampel-Empfehlung ('grün'/'gelb'/'rot') + Begruendung fuer einen Kandidaten.
    Klinik-konservativ: im Zweifel nicht gruen.
    Erwartete cand-Felder: tm_f, tm_r, delta, penalty (Qualitaets-Warnungen),
    snp_relevant, snp_3prime (relevante Varianten in der 3'-Zone),
    spec (optional: 'ok' | 'offtarget' | None)."""
    red, yellow, reasons = False, False, []

    if cand.get("snp_3prime", 0) > 0:
        red = True
        reasons.append("häufige Variante in 3'-Zone")
    if cand.get("spec") == "offtarget":
        red = True
        reasons.append("unspezifisch (Off-Target)")
    if cand.get("snp_relevant", 0) > 0:
        yellow = True
        reasons.append(f"{cand['snp_relevant']} häufige Variante(n) im Primer")
    if cand.get("snp_patho", 0) > 0:
        yellow = True
        reasons.append(f"{cand['snp_patho']} pathogene Stelle(n) im Primer (selten – beachten)")

    for tm in (cand.get("tm_f"), cand.get("tm_r")):
        if tm is None:
            continue
        d = abs(tm - target_tm)
        if d > tm_tolerance:
            red = True
            reasons.append("Tm außerhalb Fenster")
        elif d > tm_tolerance * 0.7:
            yellow = True

    delta = cand.get("delta", 0)
    if delta > 3:
        red = True
        reasons.append(f"ΔTm {delta}°C")
    elif delta > 1.5:
        yellow = True

    # Ernste Sekundaerstruktur (Hairpin am Annealing / starker Dimer) -> rot.
    if cand.get("serious", 0) > 0:
        red = True
        reasons.append("Sekundärstruktur (Hairpin/starker Dimer)")
    # Nennenswerte Qualitaets-Hinweise (PolyX, GC-Bereich, moderate Dimere) -> gelb.
    # Rein kosmetische (fehlender GC-Clamp) zaehlen NICHT, sonst ist 'grün' nie
    # erreichbar (ein fehlender GC-Clamp ist nur eine Praeferenz, kein Defekt).
    elif cand.get("notable", 0) > 0:
        yellow = True
        reasons.append(f"{cand['notable']} Qualitäts-Hinweis(e)")

    if red:
        return "rot", "; ".join(reasons) or "kritisch"
    if yellow:
        return "gelb", "; ".join(reasons) or "mit Einschränkung"
    return "grün", "alle Kriterien erfüllt"


# ============================ Ziel-Auflösung (A) & Platzierung (B) ============================
def exon_target_rel(meta, exon_number):
    """Mittlere rel-Position (5'->3') des angegebenen Exons (1-based) oder None.
    Exon-Nummerierung in Transkriptionsrichtung: bei +Strang aufsteigende
    Genom-Koordinate = Exon 1, bei -Strang absteigend (5'-Ende an hoher Koordinate)."""
    exons = meta.get("exons", [])
    if not exons or exon_number < 1 or exon_number > len(exons):
        return None
    rev = meta.get("strand", 1) != 1
    ordered = sorted(exons, key=lambda e: e["start"], reverse=rev)
    ex = ordered[exon_number - 1]
    mid_genomic = (ex["start"] + ex["end"]) // 2
    return genomic_to_rel(meta, mid_genomic)


def parse_hgvsg_to_genomic(hgvsg):
    """Parst eine genomische HGVS wie 'NC_000017.11:g.7676154G>C' zu (chrom, pos).
    Gibt (None, None) zurueck, wenn das Format nicht passt."""
    m = re.match(r"NC_0*(\d+)\.\d+:g\.(\d+)", hgvsg or "")
    if not m:
        return None, None
    num = int(m.group(1))
    chrom = {23: "X", 24: "Y", 12920: "MT"}.get(num, str(num))
    return chrom, int(m.group(2))


def target_well_placed(target_rel, f_rel, r_rel, p_len, margin):
    """True, wenn die Zielposition zwischen den Primern liegt und mindestens
    `margin` bp von beiden 3'-Enden entfernt ist (saubere Sanger-Lesbarkeit).
    target_rel=None -> keine Zielvorgabe, immer True."""
    if target_rel is None:
        return True
    f_three = f_rel + p_len - 1     # 3'-Ende Forward (hohe rel-Position)
    r_three = r_rel                 # 3'-Ende Reverse (niedrige rel-Position)
    if not (f_three < target_rel < r_three):
        return False
    return (target_rel - f_three) >= margin and (r_three - target_rel) >= margin
