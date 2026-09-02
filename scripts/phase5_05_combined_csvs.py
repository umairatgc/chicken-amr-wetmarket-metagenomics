#!/usr/bin/env python3
# =============================================================================
# phase5_05_combined_csvs.py
# Phase 5 — Step 5: Generate combined AMR + MGE + Taxonomy + Metadata CSV files
#
# For each AMR database, joins with each MGE tool on sample + contig,
# then adds GTDB-Tk taxonomy, Bakta gene annotation context, and sample metadata.
#
# AMR databases  : RGI/CARD · Abricate (all DBs) · AMRFinderPlus
# MGE tools      : MOBsuite · ISEScan · IntegronFinder · geNomad
# Taxonomy       : GTDB-Tk (bacterial species per MAG bin)
# Annotation     : Bakta (gene neighbourhood context per contig)
# Metadata       : Sample Type, Source, ONT Run, PCR results
#
# BUG FIXES applied (v2):
#   1. load_integrons:  "ID_integron" → "ID_replicon" column check
#                       (IntegronFinder 2.x has ID_replicon, not ID_integron)
#   2. load_isescan:    Strip BOM/#-prefix from column names before renaming
#                       (handles #seqID → seqID → contig silently failing)
#   3. load_isescan:    Prefer .tsv, fall back to .csv; deduplicate after concat
#                       (avoids double-counting when both file types exist)
#   4. get_samples:     Warn about samples present in AMR dirs but absent from
#                       contig_mag_map (silently excluded from analysis)
#   5. join_and_save:   Warn when many-to-many IS/integron join expands row count
#   6. load_gtdbtk:     Also load gtdbtk.ar53.summary.tsv (archaeal taxonomy)
#                       Fix hard-coded column selection that crashed on missing cols
#   7. load_mobsuite:   mobtyper_results.txt has sample_id = "consensus" (MOBsuite run ID).
#                       Previous code renamed sample_id → sample then only set sample = real
#                       name if column was absent. Since rename succeeded, mt["sample"] stayed
#                       "consensus" and never matched cr["sample"] = actual sample name →
#                       all mobt_* columns (incl. mobt_predicted_mobility) were NaN for all rows.
#                       Fixed by always setting mt["sample"] = sample (real sample name).
#   8. build_bakta_summary: bakta_cds_count checked (x == "cds") but Bakta writes "CDS"
#                       (uppercase) → count was always 0. Fixed with x.str.upper().eq("CDS").
#   9. load_mobsuite:   mger_mge_count used .values from a separate groupby — fragile row
#                       alignment. Fixed by merging count_map explicitly on contig.
#  10. get_samples:     Sample list was derived from contig_mag_map.tsv — any sample with
#                       no MAG bins (failed/empty binning) was silently excluded from all
#                       16 CSVs even though its AMR hits are valid. Fixed by deriving sample
#                       list from amr_results/rgi/ directory instead.
#  11. load_contig_mag_map: contig_mag_map.tsv already has a 'taxonomy' column with the
#                       full GTDB string. Previously ignored; now used to extract gtdb_genus
#                       and gtdb_species directly — species is always available for binned
#                       contigs without depending on GTDB-Tk files. Also added file-not-found
#                       guard (previously crashed with unhandled FileNotFoundError).
#
# USAGE:
#   conda run -p envs/r_env python3 scripts/phase5_05_combined_csvs.py
# =============================================================================

import os
import glob
import pandas as pd

# =============================================================================
# PATHS
# =============================================================================

BASE    = "/data/biol-ioi-onehealth/YOUR_USERNAME/metagenomics-proj-1"
OUTDIR  = os.path.join(BASE, "phase5", "combined_csvs")
PHASE5  = os.path.join(BASE, "phase5")

os.makedirs(OUTDIR, exist_ok=True)

# =============================================================================
# HELPER — get list of all samples from contig_mag_map
# =============================================================================

def get_samples():
    """
    Return sorted list of ALL sample names that have AMR results.

    FIX (Bug 4 revised): Previously derived sample list from contig_mag_map.tsv,
    which only contains samples that had at least one contig binned into a MAG.
    Any sample with a failed or empty binning run (e.g. PK-UF-1070) was silently
    excluded from all 16 combined CSVs even though its AMR hits are valid.

    Fix: derive sample list from amr_results/rgi/ directory (most complete AMR dir).
    Samples without MAG bins still appear in all CSVs — their taxonomy columns
    will simply be NaN (correct: no bin → no species assignment).
    """
    # Primary: use RGI results dir as the authoritative sample list
    rgi_dir = os.path.join(BASE, "amr_results", "rgi")
    if os.path.exists(rgi_dir):
        samples = sorted(
            d for d in os.listdir(rgi_dir)
            if os.path.isdir(os.path.join(rgi_dir, d))
        )
        print(f"  Samples from amr_results/rgi/: {len(samples)}")
    else:
        # Fallback: use contig_mag_map if RGI dir missing
        print("  WARNING: amr_results/rgi/ not found — falling back to contig_mag_map")
        mag_map_path = os.path.join(PHASE5, "contig_mag_map.tsv")
        mag_map = pd.read_csv(mag_map_path, sep="\t")
        samples = sorted(mag_map["sample"].unique().tolist())

    # Inform user about samples that have no MAG binning (taxonomy will be NaN)
    mag_map_path = os.path.join(PHASE5, "contig_mag_map.tsv")
    if os.path.exists(mag_map_path):
        mag_map = pd.read_csv(mag_map_path, sep="\t")
        binned = set(mag_map["sample"].unique())
        no_bins = [s for s in samples if s not in binned]
        if no_bins:
            print(f"  NOTE: {len(no_bins)} sample(s) have no MAG bins — "
                  f"taxonomy columns will be NaN for these samples:")
            for s in sorted(no_bins):
                print(f"    - {s}")

    return samples

# =============================================================================
# SECTION 0 — LOAD SAMPLE METADATA
# =============================================================================

def load_metadata():
    """
    Load sample metadata from samples.txt in the project base directory.

    Confirmed file: /data/.../metagenomics-proj-1/samples.txt
    Confirmed columns (tab-separated):
      Sr, Sr1, Sample IDs, Sample ID, FASTQ Path, sample, ONT Run,
      Sample Type, Sample Type 2, Source 1, Conc.,
      mcr-1 C2, tetX4 C2, fosA3 C2, tmex multi C2,
      NDM C2, OXA-48 C2, KPC C2, Pos Count

    Key outputs added to every combined CSV:
      sample_id   — PK-UF-1002 etc. (join key, dropped after merge)
      sample_type  — Butcher Flies / Farm etc.
      sample_type2 — more specific type
      source       — Shop 1 / Farm location etc.
      ont_run      — Run 1 / Run 2 etc.
      pcr_*        — PCR result columns

    NOTE: File has BOTH 'Sample IDs' (col 3) and 'Sample ID' (col 4).
    Both contain the same PK-* identifier. We use 'Sample IDs' as sample_id
    and drop 'Sample ID' to avoid duplicate column names after renaming.

    NOTE: File has a column literally named 'sample' (contains the word
    'sample' as a category value). Renamed to 'meta_sample_category' to
    avoid collision with the 'sample' join key in AMR/MGE dataframes.
    """
    print("Loading sample metadata...")
    # Confirmed path: samples.txt at project base dir (not in phase5/)
    path = os.path.join(BASE, "samples.txt")
    if not os.path.exists(path):
        # Fallback: try phase5/sample_metadata.csv (older location)
        path = os.path.join(PHASE5, "sample_metadata.csv")
        if not os.path.exists(path):
            print("  WARNING: samples.txt not found — metadata columns will be absent")
            return pd.DataFrame()
        print("  NOTE: Using fallback path phase5/sample_metadata.csv")

    # Tab-separated (values like 'Butcher Flies', 'Strong Positive' contain spaces)
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Strip BOM and whitespace from column names
    df.columns = [c.strip().replace('﻿', '') for c in df.columns]

    # ── Handle duplicate Sample IDs / Sample ID columns ──────────────────────
    # Both columns exist and contain the same PK-* value.
    # Rename 'Sample IDs' → sample_id, then drop 'Sample ID' to avoid conflict.
    if "Sample IDs" in df.columns:
        df = df.rename(columns={"Sample IDs": "sample_id"})
        if "Sample ID" in df.columns:
            df = df.drop(columns=["Sample ID"])
    elif "Sample ID" in df.columns:
        df = df.rename(columns={"Sample ID": "sample_id"})

    # ── Rename all other useful columns to snake_case ─────────────────────────
    rename_map = {}
    for col in df.columns:
        cs = col.strip()
        if cs == "Sample Type":
            rename_map[col] = "sample_type"
        elif cs == "Sample Type 2":
            rename_map[col] = "sample_type2"
        elif cs == "Source 1":
            rename_map[col] = "source"
        elif cs == "ONT Run":
            rename_map[col] = "ont_run"
        elif cs == "sample":
            # Column literally named 'sample' — rename to avoid join key conflict
            rename_map[col] = "meta_sample_category"
        elif cs == "mcr-1 C2":
            rename_map[col] = "pcr_mcr1"
        elif cs == "tetX4 C2":
            rename_map[col] = "pcr_tetX4"
        elif cs == "fosA3 C2":
            rename_map[col] = "pcr_fosA3"
        elif "tmex" in cs.lower():
            rename_map[col] = "pcr_tmex"
        elif cs == "NDM C2":
            rename_map[col] = "pcr_NDM"
        elif cs == "OXA-48 C2":
            rename_map[col] = "pcr_OXA48"
        elif cs == "KPC C2":
            rename_map[col] = "pcr_KPC"
        elif cs == "Pos Count":
            rename_map[col] = "pcr_pos_count"
    df = df.rename(columns=rename_map)

    # ── Keep only useful columns ──────────────────────────────────────────────
    keep = ["sample_id", "sample_type", "sample_type2", "source", "ont_run",
            "pcr_mcr1", "pcr_tetX4", "pcr_fosA3", "pcr_tmex",
            "pcr_NDM", "pcr_OXA48", "pcr_KPC", "pcr_pos_count"]
    keep_present = [c for c in keep if c in df.columns]
    df = df[keep_present].copy()

    # Drop rows without a valid sample ID
    if "sample_id" in df.columns:
        df = df[df["sample_id"].notna() & (df["sample_id"].astype(str).str.strip() != "")]
    # Drop negative controls
    if "sample_type" in df.columns:
        df = df[~df["sample_type"].str.contains("negative|control|blank",
                                                  case=False, na=False)]

    print(f"  Metadata: {len(df)} samples loaded")
    if "sample_type" in df.columns:
        print(f"  Sample types: {df['sample_type'].value_counts().to_dict()}")
    return df

# =============================================================================
# SECTION 1 — LOAD AMR DATABASES
# =============================================================================

def parse_gene_short(aro_name):
    """
    Extract a short gene symbol from an RGI Best_Hit_ARO model description.

    RGI uses full ARO model names like:
      "Escherichia coli AcrAB-TolC with AcrR mutation conferring resistance to..."
    instead of short symbols. This function parses the actual gene name out.

    Examples:
      "...with AcrR mutation..."         → acrR
      "...with MarR mutations..."        → marR
      "...soxR with mutation..."         → soxR
      "...gyrA conferring resistance..." → gyrA
      "...EF-Tu mutants..."              → tufA
      "...PBP3 conferring..."            → ftsI
      "tetR"                             → tetR  (already short)
    """
    import re
    if not aro_name or (isinstance(aro_name, float)):
        return None
    s = str(aro_name).strip()
    if re.match(r'^[a-zA-Z0-9\-]+$', s):
        return s
    if 'EF-Tu' in s:
        return 'tufA'
    if 'PBP3' in s or 'PBP 3' in s:
        return 'ftsI'
    mut = re.search(r'\bwith\s+([A-Z][a-z]{2,4}[A-Z])\s+mutation', s)
    if mut:
        g = mut.group(1)
        return g[0].lower() + g[1:]
    camel = re.findall(r'\b([a-z]{2,5}[A-Z])\b', s)
    if camel:
        return camel[0]
    lower_gene = re.findall(r'\b([a-z]{3,5}[A-Za-z])\s+(?:with|conferring)', s)
    if lower_gene:
        return lower_gene[0]
    return None


def load_rgi(samples):
    """
    Load RGI/CARD output for all samples.
    Keeps only Strict and Perfect hits.
    """
    print("Loading RGI/CARD results...")
    frames = []
    for sample in samples:
        path = os.path.join(BASE, "amr_results", "rgi", sample, "rgi_output.txt")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
            df = df[df["Cut_Off"].isin(["Strict", "Perfect"])].copy()
            if df.empty:
                continue
            df["sample"] = sample
            df = df.rename(columns={
                "Contig"                                    : "contig",
                "Start"                                     : "start",
                "Stop"                                      : "stop",
                "Best_Hit_ARO"                              : "gene_full",
                "Drug Class"                                : "drug_class",
                "Resistance Mechanism"                      : "resistance_mechanism",
                "AMR Gene Family"                           : "amr_gene_family",
                "Cut_Off"                                   : "cut_off",
                "Best_Identities"                           : "pct_identity",
                "Percentage Length of Reference Sequence"   : "pct_length_ref",
                "Best_Hit_Bitscore"                         : "bitscore",
                "ARO"                                       : "aro_accession",
                "Model_type"                                : "model_type",
                "Orientation"                               : "orientation",
                "ORF_ID"                                    : "orf_id"
            })
            df["amr_source"] = "RGI_CARD"
            df["gene_short"] = df["gene_full"].apply(parse_gene_short)
            frames.append(df)
        except Exception as e:
            print(f"  WARNING: Could not load RGI for {sample}: {e}")
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  RGI: {len(result)} hits across {result['sample'].nunique()} samples")
    return result


def load_abricate(samples):
    """
    Load Abricate results for all samples and all databases.
    Databases: resfinder, card, ncbi, argannot, vfdb, plasmidfinder, ecoh
    Filters: %COVERAGE >= 80%, %IDENTITY >= 90%.

    FIX: Previously matched 'COVERAGE' column (coordinate string like "1-800/1500")
         instead of '%COVERAGE' (numeric percentage). All rows were converted to NaN
         by pd.to_numeric() and filtered out → empty dataframe → file skipped.
         Fixed by requiring '%' in the column name when searching for coverage column.
    """
    print("Loading Abricate results...")
    DBS = ["resfinder", "card", "ncbi", "argannot", "vfdb", "plasmidfinder", "ecoh"]
    frames = []
    for sample in samples:
        for db in DBS:
            path = os.path.join(BASE, "amr_results", "abricate", sample,
                                f"abricate_{db}.txt")
            if not os.path.exists(path):
                continue
            try:
                # NOTE: do NOT use comment="#" here — Abricate's header line is
                # "#FILE\tSEQUENCE\t...\t%COVERAGE\t..." and pandas would skip it,
                # promoting the first data row to column names → %COVERAGE never found
                df = pd.read_csv(path, sep="\t", low_memory=False)
                if df.empty:
                    continue
                # Drop visual coverage map column
                if "COVERAGE_MAP" in df.columns:
                    df = df.drop(columns=["COVERAGE_MAP"])
                # Normalise column names — strip BOM and whitespace
                df.columns = [c.strip().replace('﻿', '') for c in df.columns]

                # FIX: match %COVERAGE (numeric %) NOT COVERAGE (coordinate string)
                # Require '%' in column name to distinguish "1-800/1500" from "98.7"
                cov_col = next((c for c in df.columns
                                if "COVERAGE" in c.upper()
                                and "%" in c
                                and c != "COVERAGE_MAP"), None)
                idt_col = next((c for c in df.columns
                                if "IDENTITY" in c.upper()
                                and "%" in c), None)
                if cov_col is None or idt_col is None:
                    continue

                # Apply quality filters
                df = df[
                    (pd.to_numeric(df[cov_col], errors="coerce") >= 80) &
                    (pd.to_numeric(df[idt_col], errors="coerce") >= 90)
                ].copy()
                if df.empty:
                    continue

                df = df.rename(columns={cov_col: "%COVERAGE", idt_col: "%IDENTITY"})
                df["sample"] = sample
                df = df.rename(columns={
                    "#FILE"      : "file",
                    "SEQUENCE"   : "contig",
                    "START"      : "start",
                    "END"        : "end",
                    "STRAND"     : "strand",
                    "GENE"       : "gene",
                    "COVERAGE"   : "coverage_coords",
                    "GAPS"       : "gaps",
                    "%COVERAGE"  : "pct_coverage",
                    "%IDENTITY"  : "pct_identity",
                    "DATABASE"   : "database",
                    "ACCESSION"  : "accession",
                    "PRODUCT"    : "product",
                    "RESISTANCE" : "resistance"
                })
                df["amr_source"] = f"Abricate_{db.upper()}"
                frames.append(df)
            except Exception as e:
                print(f"  WARNING: Could not load Abricate {db} for {sample}: {e}")
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  Abricate: {len(result)} hits across {result['sample'].nunique()} samples")
    return result


def load_amrfinder(samples):
    """
    Load AMRFinderPlus results for all samples.
    Keeps only AMR type hits.
    """
    print("Loading AMRFinderPlus results...")
    frames = []
    for sample in samples:
        path = os.path.join(BASE, "amr_results", "amrfinder", sample,
                            "amrfinder_results.txt")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
            df = df[df["Type"] == "AMR"].copy()
            if df.empty:
                continue
            df["sample"] = sample
            df = df.rename(columns={
                "Contig id"                    : "contig",
                "Start"                        : "start",
                "Stop"                         : "stop",
                "Strand"                       : "strand",
                "Element symbol"               : "gene",
                "Element name"                 : "gene_name",
                "Scope"                        : "scope",
                "Type"                         : "type",
                "Subtype"                      : "subtype",
                "Class"                        : "drug_class",
                "Subclass"                     : "drug_subclass",
                "Method"                       : "method",
                "Target length"                : "target_length",
                "Reference sequence length"    : "ref_length",
                "% Coverage of reference"      : "pct_coverage",
                "% Identity to reference"      : "pct_identity",
                "Alignment length"             : "alignment_length",
                "Closest reference accession"  : "ref_accession",
                "Closest reference name"       : "ref_name",
                "HMM accession"                : "hmm_accession",
                "HMM description"              : "hmm_description",
                "Protein id"                   : "protein_id"
            })
            df["amr_source"] = "AMRFinderPlus"
            frames.append(df)
        except Exception as e:
            print(f"  WARNING: Could not load AMRFinder for {sample}: {e}")
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  AMRFinder: {len(result)} hits across {result['sample'].nunique()} samples")
    return result


def load_resfinder(samples):
    """
    Load standalone ResFinder results for all samples.
    File: amr_results/resfinder_assembly/SAMPLE/ResFinder_results_tab.txt

    ResFinder_results_tab.txt columns:
      Resistance gene, Identity, Alignment Length/Gene Length, Coverage,
      Position in reference, Contig, Position in contig, Phenotype, Accession no.

    Filters: Identity >= 90%, Coverage >= 80%.
    """
    print("Loading ResFinder results...")
    frames = []
    for sample in samples:
        path = os.path.join(BASE, "amr_results", "resfinder_assembly", sample,
                            "ResFinder_results_tab.txt")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
            if df.empty:
                continue
            # Strip BOM and whitespace from column names
            df.columns = [c.strip().replace('﻿', '') for c in df.columns]

            # Apply quality filters
            # Identity column: "Identity" (numeric %)
            # Coverage column: "Coverage" (numeric %)
            if "Identity" in df.columns:
                df = df[pd.to_numeric(df["Identity"], errors="coerce") >= 90].copy()
            if "Coverage" in df.columns:
                df = df[pd.to_numeric(df["Coverage"], errors="coerce") >= 80].copy()
            if df.empty:
                continue

            df["sample"] = sample
            df = df.rename(columns={
                "Resistance gene"              : "gene",
                "Identity"                     : "pct_identity",
                "Alignment Length/Gene Length" : "alignment_over_gene_length",
                "Coverage"                     : "pct_coverage",
                "Position in reference"        : "position_in_reference",
                "Contig"                       : "contig",
                "Position in contig"           : "position_in_contig",
                "Phenotype"                    : "phenotype",
                "Accession no."                : "accession"
            })
            # Parse contig start/stop from "Position in contig" (format: start..end)
            if "position_in_contig" in df.columns:
                pos = df["position_in_contig"].str.extract(r"(\d+)\.\.(\d+)")
                df["start"] = pd.to_numeric(pos[0], errors="coerce")
                df["stop"]  = pd.to_numeric(pos[1], errors="coerce")
            df["amr_source"] = "ResFinder"
            frames.append(df)
        except Exception as e:
            print(f"  WARNING: Could not load ResFinder for {sample}: {e}")
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  ResFinder: {len(result)} hits across {result['sample'].nunique()} samples")
    return result


# =============================================================================
# SECTION 2 — LOAD MGE TOOLS
# =============================================================================

def load_mobsuite(samples):
    """
    Load MOBsuite results for all samples.
    Merges three files per sample into one enriched per-contig table:

      1. contig_report.txt    — one row per contig; join key: contig_id
                                mob_* prefix applied to all non-key columns
      2. mobtyper_results.txt — one row per plasmid; join key: sample + primary_cluster_id
                                mobt_* prefix
      3. mge.report.txt       — one row per MGE hit on a contig; join key: sample + contig_id
                                aggregated to one row per contig before joining
                                mger_* prefix

    Confirmed column headers (from PK-U3-1101P2):

    contig_report.txt:
      contig_id, molecule_type, primary_cluster_id, secondary_cluster_id,
      predicted_mobility, mash_nearest_neighbor, mash_neighbor_distance,
      rep_type(s), relaxase_type(s), ...

    mobtyper_results.txt:
      sample_id, num_contigs, size, gc, md5,
      rep_type(s), rep_type_accession(s),
      relaxase_type(s), relaxase_type_accession(s),
      mpf_type, mpf_type_accession(s),
      orit_type(s), orit_accession(s),
      predicted_mobility,
      mash_nearest_neighbor, mash_neighbor_distance, mash_neighbor_identification,
      primary_cluster_id, secondary_cluster_id,
      predicted_host_range_overall_rank, predicted_host_range_overall_name,
      observed_host_range_ncbi_rank, observed_host_range_ncbi_name,
      reported_host_range_lit_rank, reported_host_range_lit_name,
      associated_pmid(s)

    mge.report.txt:
      sample_id, molecule_type, primary_cluster_id, secondary_cluster_id,
      contig_id, size, gc, md5,
      mge_id, mge_acs, mge_type, mge_subtype,
      mge_length, mge_start, mge_end,
      contig_start, contig_end, length, sstrand,
      qcovhsp, pident, evalue, bitscore
    """
    print("Loading MOBsuite results...")
    frames = []
    mobt_loaded = 0
    mger_loaded = 0

    # Columns to keep from mobtyper (most analytically useful)
    MOBT_KEEP = [
        "num_contigs", "size", "gc",
        "rep_type(s)", "rep_type_accession(s)",
        "relaxase_type(s)", "relaxase_type_accession(s)",
        "mpf_type", "orit_type(s)",
        "predicted_mobility",
        "mash_nearest_neighbor", "mash_neighbor_distance",
        "mash_neighbor_identification",
        "predicted_host_range_overall_rank", "predicted_host_range_overall_name",
        "associated_pmid(s)"
    ]

    # Columns to aggregate from mge.report (MGE-specific — skip cols already in contig_report)
    MGER_AGG = ["mge_id", "mge_acs", "mge_type", "mge_subtype",
                "mge_length", "mge_start", "mge_end", "pident", "evalue", "bitscore"]

    for sample in samples:
        mob_dir = os.path.join(BASE, "mobile_elements", sample, "mobsuite")

        # ── 1. contig_report.txt ─────────────────────────────────────────────
        cr_path = os.path.join(mob_dir, "contig_report.txt")
        if not os.path.exists(cr_path):
            continue
        try:
            cr = pd.read_csv(cr_path, sep="\t", low_memory=False)
            cr["sample"] = sample
            cr = cr.rename(columns={"contig_id": "contig"})
        except Exception as e:
            print(f"  WARNING: Could not load MOBsuite contig_report for {sample}: {e}")
            continue

        # ── 2. mobtyper_results.txt ──────────────────────────────────────────
        mt_path = os.path.join(mob_dir, "mobtyper_results.txt")
        if os.path.exists(mt_path):
            try:
                mt = pd.read_csv(mt_path, sep="\t", low_memory=False)
                mt.columns = [c.strip().replace('﻿', '') for c in mt.columns]
                # BUG FIX (Bug 7): mobtyper_results.txt has sample_id = "consensus"
                # (the MOBsuite run ID, not the actual sample name).
                # Renaming sample_id → sample and then checking if "sample" exists
                # caused mt["sample"] = "consensus" which never matched
                # cr["sample"] = "PK-U3-XXXX" → all mobt_* columns were NaN.
                # Fix: always overwrite with the real sample name.
                mt["sample"] = sample  # always use the actual sample name
                # Keep only useful columns + join keys
                mt_join = ["sample", "primary_cluster_id"]
                keep = mt_join + [c for c in MOBT_KEEP if c in mt.columns]
                mt = mt[keep].copy()
                # Prefix non-join columns with mobt_
                mt = mt.rename(columns={c: f"mobt_{c}" for c in mt.columns
                                        if c not in mt_join})
                if "primary_cluster_id" in cr.columns:
                    cr = cr.merge(mt, on=["sample", "primary_cluster_id"], how="left")
                    mobt_loaded += 1
            except Exception as e:
                print(f"  WARNING: Could not load MOBsuite mobtyper for {sample}: {e}")

        # ── 3. mge.report.txt ────────────────────────────────────────────────
        mge_path = os.path.join(mob_dir, "mge.report.txt")
        if os.path.exists(mge_path):
            try:
                mge = pd.read_csv(mge_path, sep="\t", low_memory=False)
                mge.columns = [c.strip().replace('﻿', '') for c in mge.columns]
                if not mge.empty and "contig_id" in mge.columns:
                    mge = mge.rename(columns={"contig_id": "contig"})
                    # Keep only MGE-specific columns + join key
                    agg_cols = [c for c in MGER_AGG if c in mge.columns]
                    if agg_cols:
                        mge_sub = mge[["contig"] + agg_cols].copy()
                        # Aggregate multiple MGE rows per contig → one row
                        # Most important: mge_type and mge_subtype as pipe-joined lists
                        # Numeric fields (pident, evalue, bitscore): take best hit
                        agg_spec = {}
                        for c in agg_cols:
                            if c in ("pident", "qcovhsp"):
                                agg_spec[c] = (c, "max")
                            elif c in ("evalue",):
                                agg_spec[c] = (c, "min")
                            elif c in ("bitscore",):
                                agg_spec[c] = (c, "max")
                            else:
                                agg_spec[c] = (c, lambda x: " | ".join(
                                    x.dropna().astype(str).unique().tolist()[:5]))
                        mge_agg = mge_sub.groupby("contig").agg(
                            **{k: pd.NamedAgg(column=v[0], aggfunc=v[1])
                               for k, v in agg_spec.items()}
                        ).reset_index()
                        # Add MGE count per contig — computed inside agg to guarantee
                        # row alignment (using .values on a separate groupby is fragile)
                        count_map = mge_sub.groupby("contig").size().rename("mger_mge_count")
                        mge_agg = mge_agg.merge(count_map, on="contig", how="left")
                        # Prefix mge.report columns
                        mge_agg = mge_agg.rename(columns={c: f"mger_{c}"
                                                          for c in mge_agg.columns
                                                          if c != "contig"})
                        cr = cr.merge(mge_agg, on="contig", how="left")
                        mger_loaded += 1
            except Exception as e:
                print(f"  WARNING: Could not load MOBsuite mge.report for {sample}: {e}")

        # ── Prefix remaining contig_report columns ───────────────────────────
        # Skip join keys and already-prefixed columns
        skip = {"sample", "contig", "primary_cluster_id"}
        cr_cols = [c for c in cr.columns
                   if c not in skip
                   and not c.startswith("mobt_")
                   and not c.startswith("mger_")]
        cr = cr.rename(columns={c: f"mob_{c}" for c in cr_cols})

        frames.append(cr)

    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  MOBsuite: {len(result)} contig records "
          f"(mobtyper enriched: {mobt_loaded} samples, "
          f"mge.report enriched: {mger_loaded} samples)")
    return result


def load_isescan(samples):
    """
    Load ISEScan output for all samples.

    FIX (Bug 2): Strip BOM and leading '#' from column names before renaming.
    Some ISEScan versions write '#seqID' instead of 'seqID'. The rename
    {'seqID': 'contig'} silently does nothing if the column is '#seqID',
    causing the join to produce zero matches with no error or warning.

    FIX (Bug 3): Prefer .tsv files; only fall back to .csv if no .tsv found.
    Previously loading both caused duplicate IS element records when a sample
    had both file types (from pipeline reruns or ISEScan version differences).
    After concat, deduplicate on sample + contig + position to be safe.
    """
    print("Loading ISEScan results...")
    frames = []
    for sample in samples:
        base_dir = os.path.join(BASE, "mobile_elements", sample, "isescan", "isescan")

        # Bug 3 fix: prefer .tsv; only use .csv if no .tsv files found
        tsv_files = glob.glob(os.path.join(base_dir, "*.tsv"))
        if not tsv_files:
            tsv_files = glob.glob(os.path.join(base_dir, "*.csv"))
        if not tsv_files:
            continue

        for path in tsv_files:
            try:
                df = pd.read_csv(path, sep="\t", low_memory=False)
                if df.empty:
                    continue
                # If tab-sep read failed (true CSV), retry with comma
                if len(df.columns) == 1:
                    df = pd.read_csv(path, sep=",", low_memory=False)
                if df.empty:
                    continue

                # Bug 2 fix: strip BOM (﻿) and leading '#' from all column names
                # before renaming. Handles '#seqID' → 'seqID' → renamed to 'contig'
                df.columns = [
                    c.strip().replace('﻿', '').lstrip('#')
                    for c in df.columns
                ]

                df["sample"] = sample
                # Now rename — works whether original was 'seqID' or '#seqID'
                df = df.rename(columns={"seqID": "contig", "seq_name": "contig"})

                # Verify the rename worked — warn if 'contig' column is still absent
                if "contig" not in df.columns:
                    print(f"  WARNING: ISEScan file for {sample} has no seqID/seq_name "
                          f"column — actual columns: {list(df.columns[:5])} — skipping")
                    continue

                mge_cols = [c for c in df.columns if c not in ["sample", "contig"]]
                df = df.rename(columns={c: f"is_{c}" for c in mge_cols})
                frames.append(df)
            except Exception as e:
                print(f"  WARNING: Could not load ISEScan for {sample}: {e}")

    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)

    # Bug 3 fix: deduplicate in case both .tsv and .csv were loaded for a sample
    dedup_cols = ["sample", "contig", "is_isBegin", "is_isEnd"]
    dedup_cols_present = [c for c in dedup_cols if c in result.columns]
    if dedup_cols_present:
        before = len(result)
        result = result.drop_duplicates(subset=dedup_cols_present)
        if len(result) < before:
            print(f"  ISEScan: removed {before - len(result)} duplicate IS records "
                  f"(same sample+contig+position from multiple file formats)")

    print(f"  ISEScan: {len(result)} IS element records")
    return result


def load_integrons(samples):
    """
    Load IntegronFinder output for all samples.
    Reads TWO files per sample from Results_Integron_Finder_consensus/:

      1. consensus.integrons — detailed, one row per element inside each integron
                               (integrase gene, attC sites, gene cassettes).
                               Columns: ID_replicon, pos_beg, pos_end, strand,
                               element, type_elt, annotation, model, type, distance_type
                               → join key: ID_replicon (renamed to contig)

      2. consensus.summary   — one row per integron per contig.
                               Columns: ID_replicon, ID_integron, type_integron, default
                               type_integron values: complete | In0 | CALIN
                                 complete = has integrase + attC sites (most mobile/dangerous)
                                 In0      = integrase only, no attC (dormant)
                                 CALIN    = attC cluster, no integrase
                               → aggregated to one row per contig, then merged into main df

    FIX (Bug 1): The previous check 'ID_integron not in df.columns' was wrong.
    IntegronFinder 2.x uses 'ID_replicon' as the replicon/contig identifier column.
    'ID_integron' does not exist as a column in IntegronFinder 2.x output.
    The old check evaluated True for every sample (ID_integron never found),
    silently skipping all samples → all 3 integron CSVs were never created.
    Fixed by checking for 'ID_replicon' instead.
    """
    print("Loading IntegronFinder results...")
    frames = []
    for sample in samples:
        result_dir = os.path.join(BASE, "mobile_elements", sample, "integronfinder",
                                  "Results_Integron_Finder_consensus")

        # ── 1. consensus.integrons (detailed, element-level) ─────────────────
        path = os.path.join(result_dir, "consensus.integrons")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
            if df.empty:
                # Empty after stripping comments = no integrons found — normal result
                continue
            # Bug 1 fix: check for 'ID_replicon' (IntegronFinder 2.x column name)
            # NOT 'ID_integron' (which doesn't exist in 2.x output)
            if "ID_replicon" not in df.columns:
                print(f"  WARNING: IntegronFinder file for {sample} missing expected "
                      f"columns — actual columns: {list(df.columns[:5])} — skipping")
                continue
            df["sample"] = sample
            df = df.rename(columns={"ID_replicon": "contig"})
            mge_cols = [c for c in df.columns if c not in ["sample", "contig"]]
            df = df.rename(columns={c: f"intg_{c}" for c in mge_cols})
        except pd.errors.EmptyDataError:
            # File contains only '#' comment lines = no integrons found — silent skip
            continue
        except Exception as e:
            print(f"  WARNING: Could not load IntegronFinder .integrons for {sample}: {e}")
            continue

        # ── 2. consensus.summary (per-contig integron type counts) ───────────
        # Confirmed column structure (IntegronFinder 2.0.6):
        #   ID_replicon   CALIN   complete   In0   topology   size
        # CALIN/complete/In0 are integer counts per contig (NOT a type_integron text col).
        #   complete = has integrase + attC sites (most concerning — can capture new genes)
        #   In0      = integrase only, no attC (dormant)
        #   CALIN    = attC cluster only, no integrase
        summary_path = os.path.join(result_dir, "consensus.summary")
        if os.path.exists(summary_path):
            try:
                sm = pd.read_csv(summary_path, sep="\t", comment="#", low_memory=False)
                if not sm.empty and "ID_replicon" in sm.columns:
                    sm = sm.rename(columns={"ID_replicon": "contig"})
                    # Build a human-readable type string from the count columns
                    # e.g. complete=1, CALIN=2, In0=0 → "complete | CALIN"
                    type_cols = [c for c in ["complete", "CALIN", "In0"]
                                 if c in sm.columns]
                    def summarise_types(row):
                        return " | ".join(
                            t for t in type_cols
                            if pd.to_numeric(row.get(t, 0), errors="coerce") > 0
                        )
                    sm["intg_integron_types"]    = sm.apply(summarise_types, axis=1)
                    sm["intg_integron_count"]    = sm[type_cols].apply(
                        pd.to_numeric, errors="coerce").sum(axis=1)
                    sm["intg_complete_count"]    = pd.to_numeric(
                        sm.get("complete", 0), errors="coerce").fillna(0).astype(int)
                    sm["intg_calin_count"]       = pd.to_numeric(
                        sm.get("CALIN",    0), errors="coerce").fillna(0).astype(int)
                    sm["intg_in0_count"]         = pd.to_numeric(
                        sm.get("In0",      0), errors="coerce").fillna(0).astype(int)
                    keep = ["contig", "intg_integron_types", "intg_integron_count",
                            "intg_complete_count", "intg_calin_count", "intg_in0_count"]
                    df = df.merge(sm[keep], on="contig", how="left")
            except pd.errors.EmptyDataError:
                pass  # comment-only summary = no integrons, silent skip
            except Exception as e:
                print(f"  WARNING: Could not load IntegronFinder .summary for {sample}: {e}")

        frames.append(df)
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  IntegronFinder: {len(result)} integron records")
    return result


def load_genomad(samples):
    """
    Load geNomad plasmid summary for all samples.
    Filters: plasmid_score >= 0.7.
    """
    print("Loading geNomad results...")
    frames = []
    for sample in samples:
        path = os.path.join(BASE, "mobile_elements", sample, "genomad",
                            "consensus_summary",
                            "consensus_plasmid_summary.tsv")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
            df = df[df["plasmid_score"] >= 0.7].copy()
            if df.empty:
                continue
            df["sample"] = sample
            df = df.rename(columns={"seq_name": "contig"})
            mge_cols = [c for c in df.columns if c not in ["sample", "contig"]]
            df = df.rename(columns={c: f"gnd_{c}" for c in mge_cols})
            frames.append(df)
        except Exception as e:
            print(f"  WARNING: Could not load geNomad for {sample}: {e}")
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  geNomad: {len(result)} plasmid contig records")
    return result


# =============================================================================
# SECTION 3 — LOAD TAXONOMY AND ANNOTATION
# =============================================================================

def load_gtdbtk(samples):
    """
    Load GTDB-Tk taxonomy for all samples.

    FIX (Bug 6): Now loads BOTH bacterial (bac120) AND archaeal (ar53) taxonomy
    files and concatenates them. Previously only bacterial taxonomy was loaded,
    so any archaeal MAG bins received NaN taxonomy.

    FIX: Column selection now uses conditional logic instead of hard-coded names
    that would raise KeyError when optional columns (closest_genome_reference,
    closest_genome_ani) are absent in older GTDB-Tk output versions.
    """
    print("Loading GTDB-Tk taxonomy...")
    frames = []
    for sample in samples:
        classify_dir = os.path.join(BASE, "taxonomy_mags", sample, "classify")
        # Bug 6 fix: load both bacterial and archaeal taxonomy files
        for fname in ["gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"]:
            path = os.path.join(classify_dir, fname)
            if not os.path.exists(path):
                continue
            try:
                df = pd.read_csv(path, sep="\t", low_memory=False)
                df["sample"] = sample
                df = df.rename(columns={
                    "user_genome"    : "bin",
                    "classification" : "gtdb_classification"
                })
                df["gtdb_genus"]   = df["gtdb_classification"].str.extract(r"g__([^;]+)")
                df["gtdb_species"] = df["gtdb_classification"].str.extract(r"s__([^;]+)")

                # Bug 6 fix: build column list conditionally — some GTDB-Tk versions
                # don't produce closest_genome_reference / closest_genome_ani
                keep_cols = ["sample", "bin", "gtdb_classification",
                             "gtdb_genus", "gtdb_species"]
                for opt_col in ["closest_genome_reference", "closest_genome_ani"]:
                    if opt_col in df.columns:
                        keep_cols.append(opt_col)
                df = df[keep_cols]
                frames.append(df)
            except Exception as e:
                print(f"  WARNING: Could not load GTDB-Tk ({fname}) for {sample}: {e}")

    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  GTDB-Tk: {len(result)} MAG taxonomy records")
    return result


def load_contig_mag_map():
    """
    Load contig → MAG bin mapping table.

    contig_mag_map.tsv confirmed columns:
      sample, contig, bin, completeness, contamination, taxonomy

    The 'taxonomy' column already contains the full GTDB classification string,
    e.g. d__Bacteria;p__Bacillota;...;g__Weissella;s__Weissella paramesenteroides

    Extract gtdb_genus and gtdb_species directly from this column so that
    species names are available for ALL binned contigs without depending on
    separate GTDB-Tk summary files being present.
    (GTDB-Tk files still add closest_genome_reference and closest_genome_ani.)
    """
    print("Loading contig → MAG map...")
    path = os.path.join(PHASE5, "contig_mag_map.tsv")
    if not os.path.exists(path):
        print("  WARNING: contig_mag_map.tsv not found — taxonomy will be absent")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Rename 'taxonomy' → 'gtdb_classification' for consistency with GTDB-Tk columns
    if "taxonomy" in df.columns:
        df = df.rename(columns={"taxonomy": "gtdb_classification"})
        df["gtdb_genus"]   = df["gtdb_classification"].str.extract(r"g__([^;]+)")
        df["gtdb_species"] = df["gtdb_classification"].str.extract(r"s__([^;]+)")
    print(f"  contig_mag_map: {len(df)} contig records, "
          f"{df['sample'].nunique()} samples")
    return df


def load_bakta(samples):
    """
    Load Bakta annotation TSV for all quality bins across all samples.
    Bakta TSV has no header — columns: contig, type, start, end, strand,
    locus_tag, gene, product, [db_xrefs optional]
    """
    print("Loading Bakta annotations...")
    BAKTA_COLS = ["contig", "feature_type", "bak_start", "bak_end",
                  "bak_strand", "locus_tag", "bak_gene", "bak_product",
                  "bak_db_xrefs"]
    frames = []
    for sample in samples:
        ann_dir = os.path.join(BASE, "annotation", sample)
        if not os.path.exists(ann_dir):
            continue
        for bin_name in os.listdir(ann_dir):
            tsv_path = os.path.join(ann_dir, bin_name, f"{bin_name}.tsv")
            if not os.path.exists(tsv_path):
                continue
            try:
                df = pd.read_csv(tsv_path, sep="\t", comment="#",
                                 header=None, low_memory=False)
                if df.empty or df.shape[1] < 8:
                    continue
                actual_cols = BAKTA_COLS[:df.shape[1]]
                df.columns = actual_cols
                df["sample"] = sample
                df["bin"]    = bin_name
                frames.append(df)
            except Exception as e:
                print(f"  WARNING: Could not load Bakta for {sample}/{bin_name}: {e}")
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    print(f"  Bakta: {len(result)} annotated features across all bins")
    return result


# =============================================================================
# SECTION 4 — BUILD LOOKUP TABLES
# =============================================================================

def build_taxonomy_lookup(contig_map, gtdbtk):
    """
    Merge contig → bin map with GTDB-Tk taxonomy.

    contig_map already has gtdb_classification, gtdb_genus, gtdb_species extracted
    from its own 'taxonomy' column — so species info is always present for binned
    contigs regardless of whether GTDB-Tk files were found.

    GTDB-Tk files only add the optional ANI columns (closest_genome_reference,
    closest_genome_ani) which are not in contig_mag_map.tsv.
    """
    print("Building taxonomy lookup (contig → species)...")
    if contig_map.empty:
        return pd.DataFrame()
    if gtdbtk.empty:
        return contig_map
    # Only add ANI/reference columns from GTDB-Tk — genus/species already in contig_map
    ani_cols = ["sample", "bin"]
    for opt in ["closest_genome_reference", "closest_genome_ani"]:
        if opt in gtdbtk.columns:
            ani_cols.append(opt)
    if len(ani_cols) > 2:
        # Drop duplicates in gtdbtk (bac120 + ar53 may have same bin)
        gtdbtk_ani = gtdbtk[ani_cols].drop_duplicates(subset=["sample", "bin"])
        lookup = contig_map.merge(gtdbtk_ani, on=["sample", "bin"], how="left")
    else:
        lookup = contig_map  # GTDB-Tk has no ANI cols to add
    return lookup


def build_bakta_summary(bakta):
    """Summarise Bakta annotation per contig — gene count + product list."""
    if bakta.empty:
        return pd.DataFrame()
    print("Summarising Bakta annotation per contig...")
    summary = bakta.groupby(["sample", "bin", "contig"]).agg(
        bakta_total_features = ("feature_type", "count"),
        bakta_cds_count      = ("feature_type", lambda x: x.str.upper().eq("CDS").sum()),
        bakta_gene_list      = ("bak_gene",    lambda x: ",".join(
                                    x.dropna().astype(str).unique().tolist()[:20])),
        bakta_product_list   = ("bak_product", lambda x: " | ".join(
                                    x.dropna().astype(str).unique().tolist()[:10]))
    ).reset_index()
    return summary


# =============================================================================
# SECTION 5 — JOIN AND SAVE
# =============================================================================

def join_and_save(amr_df, mge_df, tax_lookup, bakta_summary, metadata,
                  amr_name, mge_name):
    """
    Join AMR → MGE on sample + contig, then add taxonomy, Bakta context,
    and sample metadata. Saves result as CSV to OUTDIR.

    FIX (Bug 5): Warn when many-to-many join (ISEScan/integrons) expands row
    count significantly. ISEScan can have multiple IS elements per contig, so
    a contig with 3 ARGs and 4 IS elements produces 12 rows. This is correct
    biologically but can be surprising when counting rows.
    """
    out_name = f"{amr_name}_{mge_name}.csv"
    out_path = os.path.join(OUTDIR, out_name)

    print(f"\n  Building {out_name}...")

    if amr_df.empty or mge_df.empty:
        print(f"    SKIPPED — one or both inputs are empty")
        return

    # Step 1: Left join AMR → MGE (preserves all ARG hits)
    n_amr_before = len(amr_df)
    merged = amr_df.merge(mge_df, on=["sample", "contig"], how="left")

    # Bug 5 fix: warn when join expands row count (many-to-many on same contig)
    if len(merged) > n_amr_before:
        print(f"    NOTE: Row count expanded {n_amr_before:,} → {len(merged):,} "
              f"(multiple {mge_name} records per contig — expected for IS/integron joins)")

    # Step 2: Add GTDB-Tk taxonomy
    if not tax_lookup.empty:
        merged = merged.merge(tax_lookup, on=["sample", "contig"], how="left")

    # Step 3: Add Bakta gene context
    if not bakta_summary.empty:
        merged = merged.merge(
            bakta_summary[["sample", "contig", "bakta_total_features",
                           "bakta_cds_count", "bakta_gene_list",
                           "bakta_product_list"]],
            on=["sample", "contig"],
            how="left"
        )

    # Step 4: Add sample metadata
    # Metadata join key: sample_id (PK-UF-1002 etc.) matches 'sample' column in AMR
    if not metadata.empty and "sample_id" in metadata.columns:
        merged = merged.merge(
            metadata,
            left_on="sample",
            right_on="sample_id",
            how="left"
        )
        # Drop the duplicate sample_id column (same value as 'sample')
        if "sample_id" in merged.columns:
            merged = merged.drop(columns=["sample_id"])

    merged.to_csv(out_path, index=False)
    print(f"    Saved: {out_path}  ({len(merged):,} rows, {len(merged.columns)} columns)")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("  Phase 5 — Step 5: Combined AMR + MGE + Taxonomy CSVs")
    print("=" * 60)

    # ── Get sample list ───────────────────────────────────────────
    samples = get_samples()
    print(f"\nSamples found: {len(samples)}")

    # ── Load sample metadata ──────────────────────────────────────
    print("\n── Loading metadata ───────────────────────────────────")
    metadata = load_metadata()

    # ── Load all AMR databases ────────────────────────────────────
    print("\n── Loading AMR databases ──────────────────────────────")
    rgi       = load_rgi(samples)
    abricate  = load_abricate(samples)
    amrfinder = load_amrfinder(samples)
    resfinder = load_resfinder(samples)

    # ── Load all MGE tools ────────────────────────────────────────
    print("\n── Loading MGE tools ──────────────────────────────────")
    mobsuite  = load_mobsuite(samples)
    isescan   = load_isescan(samples)
    integrons = load_integrons(samples)
    genomad   = load_genomad(samples)

    # ── Load taxonomy and annotation ──────────────────────────────
    print("\n── Loading taxonomy and annotation ────────────────────")
    contig_map    = load_contig_mag_map()
    gtdbtk        = load_gtdbtk(samples)
    bakta_raw     = load_bakta(samples)

    # Build lookup tables
    tax_lookup    = build_taxonomy_lookup(contig_map, gtdbtk)
    bakta_summary = build_bakta_summary(bakta_raw)

    # ── Generate all 16 combined CSVs ────────────────────────────
    print("\n── Generating combined CSV files ──────────────────────")

    # RGI × MGE tools
    join_and_save(rgi, mobsuite,  tax_lookup, bakta_summary, metadata, "rgi", "mobsuite")
    join_and_save(rgi, isescan,   tax_lookup, bakta_summary, metadata, "rgi", "isescan")
    join_and_save(rgi, integrons, tax_lookup, bakta_summary, metadata, "rgi", "integrons")
    join_and_save(rgi, genomad,   tax_lookup, bakta_summary, metadata, "rgi", "genomad")

    # Abricate × MGE tools
    join_and_save(abricate, mobsuite,  tax_lookup, bakta_summary, metadata, "abricate", "mobsuite")
    join_and_save(abricate, isescan,   tax_lookup, bakta_summary, metadata, "abricate", "isescan")
    join_and_save(abricate, integrons, tax_lookup, bakta_summary, metadata, "abricate", "integrons")
    join_and_save(abricate, genomad,   tax_lookup, bakta_summary, metadata, "abricate", "genomad")

    # AMRFinderPlus × MGE tools
    join_and_save(amrfinder, mobsuite,  tax_lookup, bakta_summary, metadata, "amrfinder", "mobsuite")
    join_and_save(amrfinder, isescan,   tax_lookup, bakta_summary, metadata, "amrfinder", "isescan")
    join_and_save(amrfinder, integrons, tax_lookup, bakta_summary, metadata, "amrfinder", "integrons")
    join_and_save(amrfinder, genomad,   tax_lookup, bakta_summary, metadata, "amrfinder", "genomad")

    # ResFinder × MGE tools
    join_and_save(resfinder, mobsuite,  tax_lookup, bakta_summary, metadata, "resfinder", "mobsuite")
    join_and_save(resfinder, isescan,   tax_lookup, bakta_summary, metadata, "resfinder", "isescan")
    join_and_save(resfinder, integrons, tax_lookup, bakta_summary, metadata, "resfinder", "integrons")
    join_and_save(resfinder, genomad,   tax_lookup, bakta_summary, metadata, "resfinder", "genomad")

    # ── Summary ───────────────────────────────────────────────────
    created = [f for f in os.listdir(OUTDIR) if f.endswith(".csv")]
    print("\n" + "=" * 60)
    print(f"  CSVs created: {len(created)} / 16")
    print(f"  Output dir  : {OUTDIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
