"""
ASAV Pro – Automated SNP-Aware Validation Suite
================================================
Primer-Design & -Validierung auf Basis echter Ensembl-Daten (GRCh38).

Wissenschaftliche Kernpunkte:
  * Schmelztemperatur + Sekundaerstruktur via primer3 (SantaLucia-NN),
    Fallback Biopython Tm_NN – nicht die Wallace-Faustregel.
  * Strangkorrekte Umrechnung relativer Sequenz- in genomische Koordinaten,
    damit alles auch fuer Minus-Strang-Gene (z. B. TP53) stimmt.
  * Ziel waehlbar per rsID / Exon-Nummer / genomischer Koordinate / HGVS
    (HGVS via Ensembl Variant Recoder); Amplikon umschliesst das Ziel mit
    Mindestabstand zu beiden Primern (saubere Sanger-Lesbarkeit).
  * SNP-Ranking nach Bevoelkerungs-Haeufigkeit (gnomAD, Population waehlbar);
    seltene und pathogene Varianten werden gezeigt + geflaggt, dominieren das
    Ranking aber nicht. Aggregat-Risiko = P(Patient traegt >=1 Variante).
  * Ampel-Empfehlung (gruen/gelb/rot) pro Kandidat; CSV-Export der Kandidaten.
  * On-Demand-Spezifitaet via NCBI BLAST gegen das humane Genom.
  * Offline-/Demo-Cache + Ensembl-Release-Anzeige (Reproduzierbarkeit).

Reine, testbare Logik liegt in asav_core.py (siehe test_core.py).
"""

import hashlib
import json
import os
import re

import pandas as pd
import plotly.graph_objects as go
import requests
import streamlit as st
from Bio.Seq import Seq
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from asav_core import (
    MAF_COMMON_DEFAULT, MG_CONC, NA_CONC, POPULATION_OPTIONS, PRIMER_CONC,
    aggregate_carrier_probability, assess_primer_snp, calculate_metrics,
    cross_dimer_warning, exon_target_rel, genomic_to_rel, is_snp_relevant,
    maf_from_populations, parse_hgvsg_to_genomic, primer_pair_snp_load,
    primer_warnings, recommend_traffic_light, rel_to_genomic,
    secondary_structure_warnings, snp_frequency_label, snp_is_pathogenic,
    snp_presence, snps_in_footprint, target_well_placed,
)

# UI-/App-spezifische Konstanten
MAX_PINS = 8                  # max. markierbare SNPs fuer den Bericht
ENSEMBL_SERVER = "https://rest.ensembl.org"
SANGER_MARGIN = 50            # bp Mindestabstand Ziel <-> Primer-3'-Ende (Standard)

# Beispiel-Presets (klinisch bekannte Loci) – für schnellen Einstieg & Demo.
PRESETS = {
    "— eigenes —": None,
    "BRCA2 (Brustkrebs)": {"ids": "ENSG00000139618", "mode": "Exon-Nr", "value": "11"},
    "TP53 (Li-Fraumeni)": {"ids": "ENSG00000141510", "mode": "rsID", "value": "rs1042522"},
    "CFTR (Mukoviszidose)": {"ids": "ENSG00000001626", "mode": "Exon-Nr", "value": "11"},
    "HBB (Sichelzell)": {"ids": "ENSG00000244734", "mode": "rsID", "value": "rs334"},
}
TRAFFIC_COLORS = {"grün": "#16A34A", "gelb": "#F59E0B", "rot": "#DC2626"}

# =================================================================
# 1. SYSTEM-KONFIGURATION & STYLING
# =================================================================
st.set_page_config(
    page_title="ASAV Pro – SNP-Aware Primer Suite",
    page_icon="🧬",
    layout="wide",
)

st.markdown(
    """
    <style>
    [data-testid="stMetricValue"] { font-size: 1.8rem !important; white-space: nowrap; }
    [data-testid="stMetricLabel"] { font-size: 1rem !important; }
    .stExpander { border: 1px solid #f0f2f6; border-radius: 8px; margin-bottom: 10px; }

    .seq-context {
        font-family: 'Courier New', monospace;
        font-size: 1.1rem;
        line-height: 1.6;
        word-wrap: break-word;
        background-color: #121212;
        color: #e0e0e0;
        padding: 20px;
        border-radius: 8px;
        border: 2px solid #333;
        letter-spacing: 1px;
    }
    .highlight-fwd {
        background-color: rgba(16, 185, 129, 0.2);
        border-bottom: 4px solid #10B981 !important;
        padding-bottom: 2px;
    }
    .highlight-rev {
        background-color: rgba(59, 130, 246, 0.2);
        border-bottom: 4px solid #3B82F6 !important;
        padding-bottom: 2px;
    }
    .highlight-snp-inline {
        color: #ff4d4d !important;
        font-weight: 900 !important;
        text-decoration: underline wavy #ff4d4d;
        cursor: help;
    }
    .stPopover { width: 100%; }
    </style>
    """,
    unsafe_allow_html=True,
)

# =================================================================
# 2. API-SCHNITTSTELLE (robust, mit Retry & Rate-Limit-Handling)
# =================================================================
SESSION = requests.Session()
SESSION.headers.update({
    "Content-Type": "application/json",
    "User-Agent": "ASAV-Primer-Suite/2.0 (Bachelorarbeit)",
})

# --- Demo-/Offline-Cache: erfolgreiche Antworten auf Platte spiegeln, bei
# Netzwerk-Ausfall daraus bedienen (z. B. fuer die Praesentation ohne WLAN). ---
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".asav_cache")


def _cache_path(tag, payload):
    h = hashlib.sha1(f"{tag}|{payload}".encode("utf-8")).hexdigest()[:16]
    return os.path.join(CACHE_DIR, f"{tag}_{h}.json")


def _cache_write(path, data):
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f)
    except OSError:
        pass


def _cache_read(path):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


@retry(
    reraise=True,
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=8),
    retry=retry_if_exception_type(requests.RequestException),
)
def _ensembl_request(ext):
    """Einzelner GET gegen die Ensembl-REST-API. Wirft bei Fehlern (fuer Retry)."""
    r = SESSION.get(ENSEMBL_SERVER + ext, timeout=20)
    if r.status_code == 429:  # Rate-Limit -> erneut versuchen
        raise requests.RequestException("429 Too Many Requests")
    r.raise_for_status()
    return r.json()


def _ensembl_get(ext):
    """Liefert JSON oder None. Erfolg -> Disk-Cache aktualisieren; bei Netz-Fehler
    -> zuletzt gecachte Antwort (Offline-/Demo-Modus)."""
    path = _cache_path("get", ext)
    try:
        data = _ensembl_request(ext)
        _cache_write(path, data)
        return data
    except requests.RequestException:
        return _cache_read(path)


@retry(
    reraise=True,
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=8),
    retry=retry_if_exception_type(requests.RequestException),
)
def _ensembl_post_request(ext, payload):
    """POST gegen die Ensembl-REST-API (fuer Bulk-Abfragen)."""
    r = SESSION.post(ENSEMBL_SERVER + ext, data=json.dumps(payload), timeout=30)
    if r.status_code == 429:
        raise requests.RequestException("429 Too Many Requests")
    r.raise_for_status()
    return r.json()


def _ensembl_post(ext, payload):
    path = _cache_path("post", ext + json.dumps(payload, sort_keys=True))
    try:
        data = _ensembl_post_request(ext, payload)
        _cache_write(path, data)
        return data
    except requests.RequestException:
        return _cache_read(path)


@st.cache_data(ttl=3600)
def get_gene_metadata(ensembl_id):
    data = _ensembl_get(f"/lookup/id/{ensembl_id}?expand=1")
    if not data:
        return None
    transcripts = data.get("Transcript", [])
    # Kanonisches Transkript bevorzugen, sonst das laengste.
    canonical = next((t for t in transcripts if t.get("is_canonical") == 1), None)
    if canonical:
        chosen = canonical
    elif transcripts:
        chosen = max(transcripts, key=lambda x: x.get("end", 0) - x.get("start", 0))
    else:
        chosen = {}
    return {
        "symbol": data.get("display_name", "N/A"),
        "desc": data.get("description", "Keine Beschreibung verfügbar"),
        "exons": chosen.get("Exon", []),
        "start": data.get("start"),
        "end": data.get("end"),
        "strand": data.get("strand", 1),
        "chrom": data.get("seq_region_name"),
        "transcript": chosen.get("id", "N/A"),
        "is_canonical": bool(canonical),
    }


@st.cache_data(ttl=3600)
def get_sequence_data(ensembl_id):
    """Gen-Sequenz in 5'->3'-Orientierung des Gens (strangkorrekt von Ensembl)."""
    data = _ensembl_get(f"/sequence/id/{ensembl_id}")
    if data and data.get("seq"):
        return data["seq"].upper()
    return None


@st.cache_data(ttl=3600)
def get_real_snps(chrom, g_start, g_end):
    """Variationen in einer genomischen Region (1-based, g_start<=g_end)."""
    if g_start is None or g_end is None:
        return []
    if g_start > g_end:
        g_start, g_end = g_end, g_start
    data = _ensembl_get(
        f"/overlap/region/human/{chrom}:{int(g_start)}-{int(g_end)}?feature=variation"
    )
    return data if isinstance(data, list) else []


@st.cache_data(ttl=3600)
def get_snp_frequencies(rsids, pop):
    """Bulk-Abruf der MAF (gewaehlte Population) fuer viele rsIDs (ein POST je 200).
    rsids: hashbares Tupel. Liefert {rsid: maf_oder_None}."""
    ids = [r for r in rsids if r and r.startswith("rs")]
    out = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        data = _ensembl_post("/variation/human?pops=1", {"ids": chunk})
        if isinstance(data, dict):
            for rsid, v in data.items():
                out[rsid] = maf_from_populations(v.get("populations"), pop)
    return out


@st.cache_data(ttl=3600)
def get_ensembl_release():
    """Ensembl-Release-Nummer (fuer Reproduzierbarkeit/Zitierbarkeit)."""
    data = _ensembl_get("/info/data")
    if isinstance(data, dict):
        rel = data.get("releases", [])
        return rel[0] if rel else None
    return None


@st.cache_data(ttl=3600)
def resolve_rsid(rsid):
    """rsID -> (chrom, pos) ueber das erste genomische Mapping; sonst (None, None)."""
    data = _ensembl_get(f"/variation/human/{rsid}")
    if isinstance(data, dict):
        for m in data.get("mappings", []):
            if m.get("seq_region_name") and m.get("start"):
                return str(m["seq_region_name"]), int(m["start"])
    return None, None


@st.cache_data(ttl=3600)
def resolve_hgvs(hgvs):
    """HGVS (genomisch oder c./p. mit Transkript) -> (chrom, pos) via Variant Recoder."""
    data = _ensembl_get(f"/variant_recoder/human/{hgvs}")
    if isinstance(data, list) and data:
        for allele in data[0].values():
            for g in (allele.get("hgvsg") or []):
                chrom, pos = parse_hgvsg_to_genomic(g)
                if chrom:
                    return chrom, pos
    return None, None


@st.cache_data(ttl=86400, show_spinner=False)
def check_specificity_blast(fwd, rev):
    """On-Demand-Spezifitaet via NCBI BLAST (blastn-short) gegen das humane Genom.
    UCSC In-Silico-PCR ist inzwischen CAPTCHA-gesperrt, daher dieser Weg.
    Langsam (1-3 min) und auf NCBI-Verfuegbarkeit angewiesen -> nur on demand.
    Liefert dict {status, fwd_hits, rev_hits, shared_loci, note}."""
    try:
        from io import StringIO
        from Bio.Blast import NCBIWWW, NCBIXML
    except Exception as e:  # pragma: no cover
        return {"status": "unavailable", "note": f"BLAST-Modul fehlt: {e}"}

    def _loci(seq):
        # blastn-short gegen humanes Genom (GRCh38, top-level), nur beste Treffer
        handle = NCBIWWW.qblast(
            "blastn", "GPIPE/9606/current/ref_top_level", seq,
            expect=1000, word_size=7, hitlist_size=50,
            megablast=False,
        )
        rec = NCBIXML.read(handle)
        handle.close()
        loci = []
        for aln in rec.alignments:
            for hsp in aln.hsps:
                # nahezu vollstaendige, hoch-identische Treffer = potenzielle Bindestelle
                ident = hsp.identities / hsp.align_length
                if hsp.align_length >= len(seq) - 2 and ident >= 0.90:
                    chrom = aln.hit_def
                    pos = min(hsp.sbjct_start, hsp.sbjct_end)
                    strand = "+" if hsp.sbjct_start <= hsp.sbjct_end else "-"
                    loci.append((chrom, pos, strand))
        return loci

    try:
        f_loci = _loci(fwd)
        r_loci = _loci(rev)
    except Exception as e:
        return {"status": "error", "note": f"NCBI nicht erreichbar/Timeout: {e}"}

    # gemeinsame Chromosomen, auf denen beide Primer binden + plausibler Amplikon-Abstand
    shared = 0
    for fc, fp, fs in f_loci:
        for rc, rp, rs in r_loci:
            if fc == rc and fs != rs and 0 < abs(rp - fp) <= 5000:
                shared += 1
    status = "ok" if shared <= 1 else "offtarget"
    return {"status": status, "fwd_hits": len(f_loci), "rev_hits": len(r_loci),
            "shared_loci": shared,
            "note": ("spezifisch (genau ein Amplikon erwartet)" if status == "ok"
                     else f"{shared} potenzielle Amplikon-Loci im Genom")}


# =================================================================
# 3. DARSTELLUNGS-HELFER (Genkarte & Sequenz-Viewer)
#    Reine Logik (Koordinaten, Tm, Bewertung) liegt in asav_core.py.
# =================================================================
def draw_clinical_gene_map(gene_id, metadata, marked_data, zoom_mode=False):
    """Genkarte in genomischer Forward-Koordinate (bp ab Gen-Start)."""
    fig = go.Figure()
    gene_len = metadata["end"] - metadata["start"]

    fig.add_trace(go.Scatter(
        x=[0, gene_len], y=[0, 0], mode="lines",
        line=dict(color="#D1D5DB", width=2), hoverinfo="skip",
    ))
    for exon in metadata.get("exons", []):
        ex_s = exon["start"] - metadata["start"]
        ex_e = exon["end"] - metadata["start"]
        fig.add_trace(go.Scatter(
            x=[ex_s, ex_e], y=[0, 0], mode="lines",
            line=dict(color="#3B82F6", width=14), name="Exon",
        ))

    relevant = [info for info in marked_data.values() if info["gene"] == gene_id]
    snp_x = []
    if relevant:
        snp_x = [m["map_x"] for m in relevant]
        snp_labels = [m["id"] for m in relevant]
        fig.add_trace(go.Scatter(
            x=snp_x, y=[0.2] * len(snp_x), mode="markers+text",
            marker=dict(symbol="diamond-tall", size=15, color="#EF4444",
                        line=dict(width=1, color="white")),
            text=snp_labels, textposition="top center",
        ))

    x_range = [-gene_len * 0.02, gene_len * 1.02]
    if zoom_mode and snp_x:
        pad = max(2000, gene_len * 0.05)
        x_range = [max(0, min(snp_x) - pad), min(gene_len, max(snp_x) + pad)]

    strand_lbl = "+" if metadata.get("strand", 1) == 1 else "−"
    fig.update_layout(
        height=240, margin=dict(l=10, r=10, t=40, b=10), template="plotly_white",
        showlegend=False,
        xaxis=dict(title=f"bp ab Gen-Start (Strang {strand_lbl})", range=x_range),
        yaxis=dict(showticklabels=False, range=[-0.7, 0.7]),
    )
    return fig


def generate_highlighted_sequence_html(window_seq, win_start_rel, fwd_pos, rev_pos,
                                        p_len, snps_list, meta):
    """Sequenz-Fenster mit Primer- und SNP-Highlights (alles in rel-Koordinaten)."""
    snp_rel_positions = {genomic_to_rel(meta, s.get("start")): s.get("id") for s in snps_list}

    out = []
    for i, char in enumerate(window_seq):
        rel_pos = win_start_rel + i
        res = char
        if rel_pos in snp_rel_positions:
            res = f'<span class="highlight-snp-inline" title="{snp_rel_positions[rel_pos]}">{res}</span>'
        if fwd_pos <= rel_pos < fwd_pos + p_len:
            res = f'<span class="highlight-fwd">{res}</span>'
        elif rev_pos <= rel_pos < rev_pos + p_len:
            res = f'<span class="highlight-rev">{res}</span>'
        out.append(res)
    return f'<div class="seq-context">{"".join(out)}</div>'


# =================================================================
# 4. SESSION STATE & SIDEBAR
# =================================================================
st.session_state.setdefault("batch_data", {})
st.session_state.setdefault("marked_snps_data", {})
st.session_state.setdefault("run_analysis", False)
st.session_state.setdefault("zoom_state", {})
st.session_state.setdefault("blast_results", {})
st.session_state.setdefault("gene_ids", "ENSG00000139618\nENSG00000141510")
st.session_state.setdefault("target_mode", "Position (bp)")
st.session_state.setdefault("target_value", "1000")
st.session_state.setdefault("preset_choice", "— eigenes —")


def _apply_preset():
    p = PRESETS.get(st.session_state.get("preset_choice"))
    if p:
        st.session_state["gene_ids"] = p["ids"]
        st.session_state["target_mode"] = p["mode"]
        st.session_state["target_value"] = p["value"]


with st.sidebar:
    st.header("🎯 Ziel & Population")
    st.selectbox("Beispiel-Preset (Schnellstart)", list(PRESETS),
                 key="preset_choice", on_change=_apply_preset)
    pop_label = st.selectbox("Population (Häufigkeit)", list(POPULATION_OPTIONS), index=0)
    pop_string = POPULATION_OPTIONS[pop_label]

    st.selectbox("Ziel angeben per", ["Position (bp)", "rsID", "Exon-Nr", "HGVS"],
                 key="target_mode")
    st.text_input(
        "Zielwert", key="target_value",
        help="Position: bp im Gen · rsID: z.B. rs1042522 · Exon-Nr: z.B. 11 · "
             "HGVS: z.B. NC_000017.11:g.7676154G>C (oder eine rsID).",
    )
    target_mode = st.session_state["target_mode"]
    target_value = st.session_state["target_value"].strip()
    num_candidates = st.number_input("Anzahl Kandidaten pro Gen", 1, 15, 5)

    with st.expander("🔧 Erweiterte Parameter"):
        p_len = st.slider("Primer-Länge (nt)", 18, 30, 20)
        target_tm = st.slider("Ziel-Schmelztemp. (°C)", 50, 70, 60)
        tm_tolerance = st.slider("Tm-Toleranz zum Ziel (°C)", 1.0, 10.0, 5.0, step=0.5)
        max_tm_diff = st.slider("Max. Tm-Differenz F/R (°C)", 0.5, 5.0, 2.0)
        crit_zone = st.slider("Kritische 3'-Zone (bp)", 3, 15, 5,
                              help="SNPs in dieser Distanz zum 3'-Ende gelten als HIGH RISK.")
        maf_pct = st.slider(
            "Häufigkeits-Schwelle MAF (%)", 0.1, 5.0, MAF_COMMON_DEFAULT * 100, step=0.1,
            help="Ab dieser Frequenz gilt eine Variante als 'häufig' und geht ins "
                 "Ranking ein. Seltenere werden angezeigt, aber nicht gewertet.",
        )
        min_amplicon = st.number_input("Min. Amplikon (bp)", 80, 2000, 140, step=10)
        max_amplicon = st.number_input("Max. Amplikon (bp)", 100, 4000, 500, step=10)
        sanger_margin = st.slider(
            "Mindestabstand Ziel↔Primer (bp)", 0, 150, SANGER_MARGIN, step=5,
            help="Die Zielvariante soll so weit von beiden Primern liegen "
                 "(saubere Sanger-Lesbarkeit des Amplikons).",
        )
    maf_threshold = maf_pct / 100.0

    if st.button("🗑️ Komplett-Reset", use_container_width=True):
        st.session_state.batch_data = {}
        st.session_state.marked_snps_data = {}
        st.session_state.run_analysis = False
        st.session_state.zoom_state = {}
        st.session_state.blast_results = {}
        st.rerun()

    st.divider()
    _release = get_ensembl_release()
    st.caption(
        f"Tm: primer3 (SantaLucia-NN) · Na⁺ {NA_CONC:.0f} mM · Mg²⁺ {MG_CONC} mM · "
        f"Primer {PRIMER_CONC} nM · Population {pop_label}"
    )
    st.caption(
        "Daten: Ensembl REST"
        + (f" (Release {_release})" if _release else "")
        + " · GRCh38 · gnomAD/ClinVar"
    )

# =================================================================
# 5. PHASE 1: VALIDIERUNG & GENKARTEN
# =================================================================
st.title("🧬 ASAV Pro – SNP-Aware Primer Suite")

st.header("1️⃣ Gen-Validierung & Sektoren-Mapping")
with st.container(border=True):
    st.text_area("Ensembl-Gen-IDs", key="gene_ids")
    if st.button("Gen-Strukturen laden", use_container_width=True):
        raw_ids = re.split(r"[,\n\s]+", st.session_state["gene_ids"])
        clean_ids = [i.strip() for i in raw_ids if i.strip()]
        new_batch = {}
        errors = []
        for eid in clean_ids:
            meta = get_gene_metadata(eid)
            seq = get_sequence_data(eid)
            if meta and seq:
                new_batch[eid] = {**meta, "seq": seq}
                st.session_state.zoom_state[eid] = False
            else:
                errors.append(eid)
        st.session_state.batch_data = new_batch
        st.session_state.run_analysis = False
        if errors:
            st.error(
                "Konnte folgende IDs nicht laden (ungültig oder Ensembl nicht "
                f"erreichbar): {', '.join(errors)}"
            )

if st.session_state.batch_data:
    for eid, info in st.session_state.batch_data.items():
        with st.container(border=True):
            col_metrics, col_map = st.columns([1.5, 2])
            with col_metrics:
                strand_lbl = "+" if info.get("strand", 1) == 1 else "−"
                st.subheader(f"{info['symbol']} ({eid})")
                st.caption(
                    f"Transkript {info['transcript']} "
                    f"{'(kanonisch)' if info['is_canonical'] else '(längstes)'} · "
                    f"Strang {strand_lbl}"
                )
                m1, m2 = st.columns(2)
                m1.metric("Chr", f"Chr {info['chrom']}")
                m2.metric("Sektoren", f"{len(info['exons'])} Exons")
                st.metric("Start-Position", f"{info['start']:,} bp")

                nav1, nav2 = st.columns(2)
                if nav1.button("🔍 SNP-Fokus", key=f"z_{eid}"):
                    st.session_state.zoom_state[eid] = not st.session_state.zoom_state.get(eid, False)
                    st.rerun()
                if nav2.button("🗑️ Pins löschen", key=f"c_{eid}"):
                    st.session_state.marked_snps_data = {
                        k: v for k, v in st.session_state.marked_snps_data.items()
                        if v["gene"] != eid
                    }
                    st.rerun()
            with col_map:
                st.plotly_chart(
                    draw_clinical_gene_map(
                        eid, info, st.session_state.marked_snps_data,
                        st.session_state.zoom_state.get(eid, False),
                    ),
                    use_container_width=True,
                )

    # =================================================================
    # 6. PHASE 2: PRIMER-ANALYSE & AMPLIKON-AUDIT
    # =================================================================
    st.divider()
    st.header("2️⃣ Primer-Analyse & Amplikon-Audit")
    if st.button("🚀 Primer-Berechnung & Sequenz-Audit starten", use_container_width=True):
        st.session_state.run_analysis = True

    if st.session_state.run_analysis:
        # Ziel ggf. EINMAL genomisch aufloesen (rsID/HGVS gelten genomweit).
        tgt_chrom, tgt_pos = None, None
        if target_mode == "rsID" and target_value:
            tgt_chrom, tgt_pos = resolve_rsid(target_value)
        elif target_mode == "HGVS" and target_value:
            tgt_chrom, tgt_pos = resolve_hgvs(target_value)

        def target_rel_for(info):
            """Liefert (rel_position_oder_None, hinweis) fuer dieses Gen."""
            if target_mode == "Position (bp)":
                try:
                    return int(float(target_value)), None
                except ValueError:
                    return None, "ungültige Position"
            if target_mode == "Exon-Nr":
                try:
                    n = int(float(target_value))
                except ValueError:
                    return None, "ungültige Exon-Nummer"
                r = exon_target_rel(info, n)
                return (r, None) if r is not None else (None, f"Exon {n} nicht vorhanden")
            # rsID / HGVS
            if tgt_chrom is None:
                return None, "Ziel nicht auflösbar"
            if str(tgt_chrom) != str(info["chrom"]) or not (info["start"] <= tgt_pos <= info["end"]):
                return None, "Ziel liegt nicht in diesem Gen"
            return genomic_to_rel(info, tgt_pos), None

        annealing_tm = target_tm - 5      # Referenz fuer Hairpin-Bewertung
        csv_rows = []

        for eid, info in st.session_state.batch_data.items():
            st.subheader(f"Kandidaten für {info['symbol']}")
            full_seq = info["seq"]

            target_rel, note = target_rel_for(info)
            if target_rel is None:
                st.warning(f"Kein Ziel für {info['symbol']}: {note}.")
                continue
            st.caption(f"Ziel: {target_mode} = {target_value} → rel-Position {target_rel} bp")

            # SNPs der gesamten Suchregion einmal holen; lokal je Primer filtern.
            margin = sanger_margin
            search_span = num_candidates * 8 * 4
            reach = margin + p_len + search_span
            lo_rel = max(0, target_rel - reach)
            hi_rel = min(len(full_seq), target_rel + reach + p_len)
            gA, gB = rel_to_genomic(info, lo_rel, hi_rel - lo_rel)
            region_snps = get_real_snps(info["chrom"], gA, gB)
            maf_map = get_snp_frequencies(
                tuple(sorted({s.get("id") for s in region_snps if s.get("id")})), pop_string
            )
            region_snps = [{**s, "maf": maf_map.get(s.get("id"))} for s in region_snps]

            candidate_pairs = []
            # Forward upstream, Reverse downstream; 3'-Enden mind. `margin` bp vom
            # Ziel entfernt, nach aussen wandernd -> Amplikon waechst kontrolliert,
            # Ziel liegt garantiert lesbar im Amplikon.
            for i in range(num_candidates * 8):
                offset = i * 4
                f_s = target_rel - margin - p_len - offset
                r_o = target_rel + margin + offset
                if f_s < 0 or (r_o + p_len) > len(full_seq):
                    continue
                amplicon_len = (r_o + p_len) - f_s
                if not (min_amplicon <= amplicon_len <= max_amplicon):
                    continue
                if not target_well_placed(target_rel, f_s, r_o, p_len, margin):
                    continue

                f_seq = full_seq[f_s:f_s + p_len].upper()
                r_seq = str(Seq(full_seq[r_o:r_o + p_len].upper()).reverse_complement())
                tm_f, gc_f = calculate_metrics(f_seq)
                tm_r, gc_r = calculate_metrics(r_seq)

                if not (abs(tm_f - tm_r) <= max_tm_diff
                        and abs(tm_f - target_tm) <= tm_tolerance
                        and abs(tm_r - target_tm) <= tm_tolerance):
                    continue

                # Qualitaet: Heuristik + primer3-Sekundaerstruktur + Cross-Dimer
                warns_f = primer_warnings(f_seq) + secondary_structure_warnings(f_seq, annealing_tm)
                warns_r = primer_warnings(r_seq) + secondary_structure_warnings(r_seq, annealing_tm)
                cross = cross_dimer_warning(f_seq, r_seq)

                snps_f = snps_in_footprint(region_snps, info, f_s, p_len)
                snps_r = snps_in_footprint(region_snps, info, r_o, p_len)
                snp_pen, snp_rel, snp_tot = primer_pair_snp_load(
                    snps_f, snps_r, info, f_s, r_o, p_len, crit_zone, maf_threshold
                )
                # relevante SNPs in der kritischen 3'-Zone (fuer Ampel)
                snp_3prime = 0
                for s in snps_f:
                    if is_snp_relevant(s, maf_threshold) and \
                       assess_primer_snp(s, info, f_s, p_len, crit_zone, False)["dist"] <= crit_zone:
                        snp_3prime += 1
                for s in snps_r:
                    if is_snp_relevant(s, maf_threshold) and \
                       assess_primer_snp(s, info, r_o, p_len, crit_zone, True)["dist"] <= crit_zone:
                        snp_3prime += 1
                snp_patho = sum(1 for s in (snps_f + snps_r) if snp_is_pathogenic(s))
                p_any, n_with, n_without = aggregate_carrier_probability(snps_f + snps_r)
                # Warnungen klassifizieren -> Ampel-Logik:
                #   serious = Hairpin/starker Dimer -> rot
                #   notable = nennenswert (PolyX, GC-Bereich, moderate Dimere) -> gelb
                #   kosmetisch (fehlender GC-Clamp) -> nur Info, beeinflusst Ampel nicht
                all_w = warns_f + warns_r + cross
                serious = sum(1 for w in all_w if "Hairpin" in w or "stark" in w)
                notable = sum(1 for w in all_w
                              if "Hairpin" not in w and "stark" not in w and "GC-Clamp" not in w)

                cand = {
                    "f_seq": f_seq, "r_seq": r_seq,
                    "tm_f": tm_f, "tm_r": tm_r, "gc_f": gc_f, "gc_r": gc_r,
                    "f_pos": f_s, "r_pos": r_o, "amplicon": amplicon_len,
                    "delta": round(abs(tm_f - tm_r), 1),
                    "warns_f": warns_f, "warns_r": warns_r, "cross": cross,
                    "penalty": len(warns_f) + len(warns_r) + len(cross),
                    "serious": serious, "notable": notable,
                    "snp_penalty": snp_pen, "snp_relevant": snp_rel,
                    "snp_total": snp_tot, "snp_3prime": snp_3prime,
                    "snp_patho": snp_patho,
                    "snps_f": snps_f, "snps_r": snps_r,
                    "risk_any": p_any, "risk_with": n_with, "risk_without": n_without,
                }
                light, light_reason = recommend_traffic_light(cand, target_tm, tm_tolerance)
                cand["light"], cand["light_reason"] = light, light_reason
                candidate_pairs.append(cand)

            if not candidate_pairs:
                st.warning(
                    "Keine Kandidaten im Filter. Tipp: Tm-Toleranz erhöhen, "
                    "Amplikon-Grenzen lockern, Mindestabstand senken oder Ziel prüfen."
                )
                continue

            # Ranking: Ampel (grün<gelb<rot), dann SNP-Last, Qualitaet, Tm-Differenz.
            light_order = {"grün": 0, "gelb": 1, "rot": 2}
            ranked = sorted(
                candidate_pairs,
                key=lambda x: (light_order[x["light"]], x["snp_penalty"],
                               x["penalty"], x["delta"]),
            )

            emoji_map = {"grün": "🟢", "gelb": "🟡", "rot": "🔴"}
            for idx, c in enumerate(ranked[:num_candidates]):
                rare = c["snp_total"] - c["snp_relevant"]
                rare_txt = f" ({rare} selten)" if rare else ""
                if c["snp_relevant"] == 0:
                    snp_flag = f"🧬 keine häufigen SNP{rare_txt}"
                else:
                    snp_flag = (f"🧬 {c['snp_relevant']} häufig · "
                                f"Last {c['snp_penalty']}{rare_txt}")
                if c["snp_patho"]:
                    snp_flag += f" · ⚠{c['snp_patho']} patho"
                # Aggregat-Risiko-Text (Vorhandensein mind. einer Variante)
                if c["risk_any"] is not None:
                    risk_txt = (f"{c['risk_any'] * 100:.1f}% Wahrscheinlichkeit, dass ein "
                                f"Patient ≥1 Variante unter den Primern trägt")
                    if c["risk_without"]:
                        risk_txt += f" (+{c['risk_without']} ohne Frequenzdaten)"
                else:
                    n = c["risk_without"]
                    risk_txt = ("keine Frequenzdaten für die Varianten unter den Primern"
                                + (f" ({n} Varianten)" if n else " – keine Varianten"))

                csv_rows.append({
                    "Gen": info["symbol"], "Ensembl_ID": eid, "Ampel": c["light"],
                    "Forward_5_3": c["f_seq"], "Reverse_5_3": c["r_seq"],
                    "Tm_F": c["tm_f"], "Tm_R": c["tm_r"], "GC_F": c["gc_f"], "GC_R": c["gc_r"],
                    "Amplikon_bp": c["amplicon"], "dTm": c["delta"],
                    "SNP_haeufig": c["snp_relevant"], "SNP_gesamt": c["snp_total"],
                    "Pathogene_Stellen": c["snp_patho"],
                    "Risiko_Prozent": (round(c["risk_any"] * 100, 1)
                                       if c["risk_any"] is not None else ""),
                    "Warnungen": "; ".join(c["warns_f"] + c["warns_r"] + c["cross"]) or "-",
                    "Empfehlung": c["light_reason"],
                })

                emoji = emoji_map[c["light"]]
                title = (f"{emoji} Kandidat #{idx + 1} | {snp_flag} | "
                         f"Amplikon {c['amplicon']} bp | ΔTm {c['delta']}°C")
                with st.expander(title, expanded=(idx == 0)):
                    color = TRAFFIC_COLORS[c["light"]]
                    st.markdown(
                        f"<div style='background:{color}1a; border-left:6px solid {color}; "
                        f"padding:8px 12px; border-radius:6px; margin-bottom:8px;'>"
                        f"<b>{emoji} Empfehlung: {c['light'].upper()}</b> — {c['light_reason']}</div>",
                        unsafe_allow_html=True,
                    )
                    st.caption(f"🧬 Vorhandensein: {risk_txt}")
                    if c["snp_patho"]:
                        st.warning(
                            f"⚠️ {c['snp_patho']} bekannte **pathogene Stelle(n)** unter den "
                            f"Primern (meist selten – im SNP-Audit prüfen, Primer ggf. "
                            f"verschieben). Wird nicht versteckt, fließt aber nicht ins Ranking."
                        )
                    for w in c["cross"]:
                        st.caption(f"⚠️ {w}")

                    # On-Demand-Spezifitaet (NCBI BLAST gegen humanes Genom)
                    bkey = f"{eid}_{idx}"
                    bres = st.session_state.blast_results.get(bkey)
                    bc1, bc2 = st.columns([1.4, 3])
                    if bc1.button("🎯 Spezifität prüfen (BLAST)", key=f"blast_{bkey}",
                                  help="NCBI BLAST gegen das humane Genom – dauert 1–3 min"):
                        with st.spinner("BLAST gegen humanes Genom läuft (1–3 min)…"):
                            st.session_state.blast_results[bkey] = check_specificity_blast(
                                c["f_seq"], c["r_seq"])
                        st.rerun()
                    if bres:
                        msg = bres.get("note", "")
                        if bres["status"] == "ok":
                            bc2.success(f"Spezifität: {msg}")
                        elif bres["status"] == "offtarget":
                            bc2.error(f"Spezifität: {msg}")
                        else:
                            bc2.info(f"Spezifität: {msg}")

                    col1, col2, col3 = st.columns([2, 2, 1.2])
                    with col1:
                        st.markdown("**Forward (5'-3')**")
                        st.code(c["f_seq"])
                        st.write(f"Tm: **{c['tm_f']}°C** | GC: **{c['gc_f']}%**")
                        for w in c["warns_f"]:
                            st.caption(f"⚠️ {w}")
                    with col2:
                        st.markdown("**Reverse (5'-3')**")
                        st.code(c["r_seq"])
                        st.write(f"Tm: **{c['tm_r']}°C** | GC: **{c['gc_r']}%**")
                        for w in c["warns_r"]:
                            st.caption(f"⚠️ {w}")
                    with col3:
                        st.write(" ")
                        gene_marked = [
                            sid for sid, m in st.session_state.marked_snps_data.items()
                            if m["gene"] == eid
                        ]
                        with st.popover(f"⚙️ SNP-Audit ({len(gene_marked)})",
                                        use_container_width=True):
                            # bereits bei der Berechnung ermittelt (kein erneuter API-Call)
                            snps_f = c["snps_f"]
                            snps_r = c["snps_r"]

                            t1, t2 = st.tabs(["Forward SNPs", "Reverse SNPs"])
                            audit_iter = zip(
                                [t1, t2], [snps_f, snps_r], ["F", "R"],
                                [c["f_pos"], c["r_pos"]],
                            )
                            for tab, s_list, p_type, p_rel in audit_iter:
                                with tab:
                                    if not s_list:
                                        st.success("Keine bekannten Varianten unter diesem Primer.")
                                    for s in s_list:
                                        sid = s.get("id")
                                        a = assess_primer_snp(
                                            s, info, p_rel, p_len, crit_zone,
                                            is_reverse=(p_type == "R"),
                                        )
                                        relevant = is_snp_relevant(s, maf_threshold)
                                        freq_txt = snp_frequency_label(s)
                                        _, presence_txt = snp_presence(s)
                                        freq_color = "#1D4ED8" if relevant else "#9CA3AF"

                                        cl, cr = st.columns([3, 1])
                                        cl.write(f"**{sid}**" + ("" if relevant else " · _nur Info_"))
                                        is_p = sid in st.session_state.marked_snps_data
                                        marks_full = len(st.session_state.marked_snps_data) >= MAX_PINS

                                        btn_key = f"b_{eid}_{idx}_{p_type}_{sid}_{p_rel}"
                                        if cr.button("📍" if is_p else "⚪", key=btn_key):
                                            if is_p:
                                                del st.session_state.marked_snps_data[sid]
                                            elif not marks_full:
                                                st.session_state.marked_snps_data[sid] = {
                                                    "id": sid, "gene": eid,
                                                    "rel_pos": a["snp_rel"],
                                                    "map_x": (s.get("start") or info["start"]) - info["start"],
                                                    "alleles": s.get("alleles"),
                                                    "status": a["risk"],
                                                    "bio": a["bio_class"],
                                                    "consequence": a["consequence"],
                                                    "frequenz": freq_txt,
                                                    "vorkommen": presence_txt,
                                                }
                                            else:
                                                st.warning(f"Max. {MAX_PINS} Markierungen erreicht.")
                                            st.rerun()

                                        st.markdown(
                                            f"<span style='{a['style']}'>{a['risk']}</span> "
                                            f"<span style='background-color:{a['bio_color']}; color:white; "
                                            f"padding:4px 10px; border-radius:8px; font-weight:bold;'>"
                                            f"{a['bio_class']}</span> "
                                            f"<span style='background-color:{freq_color}; color:white; "
                                            f"padding:4px 10px; border-radius:8px; font-weight:bold;'>"
                                            f"{freq_txt}</span>",
                                            unsafe_allow_html=True,
                                        )
                                        ranking_note = ("häufig → fließt ins Ranking ein" if relevant
                                                        else "selten → nur Info, nicht im Ranking")
                                        if snp_is_pathogenic(s):
                                            ranking_note += " · ⚠ ClinVar-pathogen (trotz Seltenheit beachten)"
                                        st.caption(
                                            f"Grund: {a['reason']} · "
                                            f"Konsequenz: {a['consequence'] or 'n/a'}"
                                        )
                                        st.caption(f"Vorkommen: {presence_txt} · {ranking_note}")
                                        st.write("")

                    # --- Amplikon-Kontext / dynamischer Zoom ---
                    st.divider()
                    btn_text = (
                        f"🔬 Fokus-Ansicht: {len(gene_marked)} SNP(s) markiert"
                        if gene_marked else "🔬 Standard-Amplikom anzeigen"
                    )
                    with st.popover(btn_text, use_container_width=True):
                        if not gene_marked:
                            st.info("💡 Markiere SNPs im SNP-Audit für dynamische Detail-Zooms.")
                            st.markdown("**Komplettes Amplikom (Standardansicht)**")

                            win_start_rel = max(0, min(c["f_pos"], c["r_pos"]) - 50)
                            win_end_rel = min(
                                len(full_seq),
                                max(c["f_pos"] + p_len, c["r_pos"] + p_len) + 50,
                            )
                            window_seq_raw = full_seq[win_start_rel:win_end_rel]
                            _, win_gc = calculate_metrics(window_seq_raw)
                            cw1, cw2 = st.columns(2)
                            cw1.metric("Fenster-Länge", f"{len(window_seq_raw)} bp")
                            cw2.metric("Lokaler GC-Gehalt", f"{win_gc}%")

                            g0, g1 = rel_to_genomic(info, win_start_rel,
                                                    win_end_rel - win_start_rel)
                            snps_in_window = get_real_snps(info["chrom"], g0, g1)
                            st.markdown(
                                generate_highlighted_sequence_html(
                                    window_seq_raw, win_start_rel, c["f_pos"],
                                    c["r_pos"], p_len, snps_in_window, info,
                                ),
                                unsafe_allow_html=True,
                            )
                        else:
                            st.success(f"✅ {len(gene_marked)} dynamische Zoom-Fenster generiert.")
                            for m_sid in gene_marked[:MAX_PINS]:
                                m_snp = st.session_state.marked_snps_data[m_sid]
                                st.markdown(f"#### 🎯 Zoom-Fokus: {m_snp['id']}")
                                snp_rel = m_snp["rel_pos"]

                                win_start_rel = max(0, snp_rel - 50)
                                win_end_rel = min(len(full_seq), snp_rel + 50)
                                window_seq_raw = full_seq[win_start_rel:win_end_rel]
                                _, win_gc = calculate_metrics(window_seq_raw)
                                cw1, cw2 = st.columns(2)
                                cw1.metric("Fenster-Länge", f"{len(window_seq_raw)} bp")
                                cw2.metric("Lokaler GC-Gehalt", f"{win_gc}%")

                                g0, g1 = rel_to_genomic(info, win_start_rel,
                                                        win_end_rel - win_start_rel)
                                snps_in_window = get_real_snps(info["chrom"], g0, g1)
                                st.markdown(
                                    generate_highlighted_sequence_html(
                                        window_seq_raw, win_start_rel, c["f_pos"],
                                        c["r_pos"], p_len, snps_in_window, info,
                                    ),
                                    unsafe_allow_html=True,
                                )
                                st.write("---")

                        st.caption(
                            "Legende: <span style='border-bottom:3px solid #10B981'>Forward</span> | "
                            "<span style='border-bottom:3px solid #3B82F6'>Reverse</span> | "
                            "<span style='color:#ff4d4d; font-weight:bold;'>SNP</span> (Hover für rsID)",
                            unsafe_allow_html=True,
                        )

        if csv_rows:
            st.divider()
            cdf = pd.DataFrame(csv_rows)
            st.download_button(
                "📥 Primer-Kandidaten als CSV (alle Gene)",
                cdf.to_csv(index=False).encode("utf-8"),
                "ASAV_Primer_Kandidaten.csv", "text/csv", use_container_width=True,
            )

# =================================================================
# 7. PHASE 3: LABOR-BERICHT
# =================================================================
if st.session_state.marked_snps_data:
    st.divider()
    st.header("3️⃣ Validierter klinischer Labor-Bericht")
    report_df = pd.DataFrame.from_dict(st.session_state.marked_snps_data, orient="index")

    if "alleles" in report_df:
        def _fmt_alleles(a):
            text = ",".join(a) if isinstance(a, list) else str(a)
            return text[:40] + "…" if len(text) > 40 else text
        report_df["alleles"] = report_df["alleles"].apply(_fmt_alleles)

    cols = [c for c in ["gene", "rel_pos", "consequence", "bio", "frequenz",
                        "vorkommen", "alleles", "status"] if c in report_df]
    st.table(report_df[cols])

    csv = report_df[cols].to_csv().encode("utf-8")
    st.download_button("📥 Bericht als CSV", csv, "ASAV_Report.csv", "text/csv",
                       use_container_width=True)

st.divider()
st.caption("ASAV Pro · SNP-Aware Primer Suite · Daten: Ensembl REST (GRCh38) & ClinVar")
