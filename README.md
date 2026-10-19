# ASAV Pro – SNP-Aware Primer Suite

> PCR-Primer automatisch entwerfen **und** gegen echte Bevölkerungs-Variation
> (gnomAD / ClinVar) absichern – damit kein Primer auf einer polymorphen Stelle sitzt.

**▶️ Live-Demo:** **[asav-primer-suite.streamlit.app](https://asav-primer-suite-84dxtxt3chmka6djabtlnv.streamlit.app/)**
**Stack:** Python · Streamlit · primer3 · Biopython · Ensembl REST · NCBI BLAST

Werkzeug zum **Design und zur Validierung von PCR-Primern** auf Basis echter
Ensembl-Daten (GRCh38). Für jedes Gen werden Primer-Paare generiert,
thermodynamisch bewertet und gegen bekannte Variationen geprüft.

### Was dieses Projekt zeigt
- **Integration mehrerer wissenschaftlicher APIs** (Ensembl REST, Variant Recoder,
  gnomAD-Frequenzen, NCBI BLAST) inkl. Retry/Rate-Limit-Handling und Offline-Cache.
- **Saubere Architektur:** testbare Kernlogik (`asav_core.py`) strikt getrennt von
  der UI (`app.py`).
- **38 Unit-Tests** (pytest) für die wissenschaftliche Logik.
- **Korrekte Bioinformatik:** strangkorrekte Koordinaten, gnomAD-Minor-Allele-
  Frequenzen, Nearest-Neighbor-Tm (primer3), Hardy-Weinberg-Trägerwahrscheinlichkeit.

## Funktionen

- **Gen-Lookup** über die Ensembl-REST-API (GRCh38), inkl. kanonischem
  Transkript, Exon-Struktur und Strang.
- **Ziel flexibel angeben**: per **rsID**, **Exon-Nummer**, **genomischer
  Koordinate** oder **HGVS** (via Ensembl Variant Recoder). Das Amplikon
  umschließt das Ziel mit Mindestabstand zu beiden Primern (saubere Sanger-Reads).
- **Schmelztemperatur + Sekundärstruktur** via **primer3** (SantaLucia-NN;
  Hairpin, Self-/Cross-Dimer), Fallback Biopython `Tm_NN` – nicht die Wallace-
  Faustregel. Plus Heuristiken: GC-Gehalt, 3'-GC-Clamp, PolyX-Runs.
- **SNP-aware Ranking nach Bevölkerungs-Häufigkeit** (**gnomAD**, Population
  wählbar: global / EUR / AFR / EAS / SAS / AMR). Nur *häufige* Varianten
  steuern das Ranking; *seltene* und *pathogene* Stellen werden **angezeigt und
  geflaggt** (nie versteckt), dominieren das Ranking aber nicht.
- **Vorhandensein-Wahrscheinlichkeit** je Variante (Hardy-Weinberg) und
  **Aggregat-Risiko** je Primer: P(Patient trägt ≥1 Variante unter den Primern).
- **Ampel-Empfehlung** (🟢/🟡/🔴) pro Kandidat; strangkorrekte Koordinaten
  auch für Minus-Strang-Gene (z. B. TP53).
- **Spezifität on-demand** via NCBI BLAST (`blastn-short`) gegen das humane Genom.
- **Sequenz-Viewer** mit Primer-/SNP-Highlights; **CSV-Export** der Primer-
  Kandidaten und des SNP-Berichts.
- **Offline-/Demo-Cache** (gespiegelte Ensembl-Antworten) für Präsentationen
  ohne stabiles Netz; **Ensembl-Release** wird angezeigt (Reproduzierbarkeit).

## Installation & Start

```bash
python -m venv venv
venv/Scripts/pip install -r requirements.txt      # Windows
# bzw. source venv/bin/activate && pip install -r requirements.txt  (Linux/macOS)

venv/Scripts/streamlit run app.py
```

## Tests

```bash
venv/Scripts/python -m pytest -q
```

Die wissenschaftliche Kernlogik liegt UI-frei in `asav_core.py` und wird von
`test_core.py` abgedeckt (Koordinaten, Tm, Primer-Qualität, SNP-Bewertung).

## Projektstruktur

| Datei            | Zweck                                                        |
|------------------|--------------------------------------------------------------|
| `app.py`         | Streamlit-Oberfläche, Ensembl-API, Darstellung              |
| `asav_core.py`   | Testbare Kernlogik (keine UI-Abhängigkeiten)                |
| `test_core.py`   | Unit-Tests für `asav_core`                                  |
| `requirements.txt` | Abhängigkeiten (gepinnt)                                   |
| `legacy/`        | Frühe Prototyp-Skripte & Daten-Dumps (nicht mehr benötigt)  |

## Datenquellen

- [Ensembl REST API](https://rest.ensembl.org) – Gene, Sequenzen, Variationen,
  Variant Recoder (HGVS), Populationsfrequenzen (gnomAD/1000G)
- ClinVar-Signifikanz und VEP-Konsequenzen über die Ensembl-Variation-Daten
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov) – Spezifität gegen das humane Genom

> Hinweis: UCSC In-Silico-PCR ist inzwischen durch eine Bot-Challenge
> (Cloudflare) gesperrt; die Spezifitätsprüfung nutzt daher NCBI BLAST.

## Wissenschaftliche Annahmen (Tm)

Zentral in `asav_core.py` dokumentiert: Na⁺ 50 mM, Mg²⁺ 1,5 mM, dNTPs 0,6 mM,
Primer-Konzentration 250 nM. Bei nicht-ACGT-Basen fällt die Berechnung auf das
GC-Modell (`Tm_GC`) zurück.
