# Pipeline Scripts

Analysis scripts for *Convergent Antimicrobial Resistance in Intensive Chicken Farming and Wet-markets: Genomic Evidence of Cross-Environment Dissemination*.

Covers the full workflow from raw Nanopore long reads (61 samples + 2 negative controls) through read QC, host depletion, taxonomic classification, decontamination, de novo assembly and polishing, MAG binning/annotation, AMR and MGE detection, and AMR–MGE co-localisation. Figure-generation scripts are not included here.

All scripts were run on the University of Oxford Advanced Research Computing (ARC) SLURM cluster. Paths are hardcoded to the original project directory (`/data/biol-ioi-onehealth/.../metagenomics-proj-1-umair`) and SLURM `#SBATCH` headers reflect the resources requested on ARC — update the `BASE` variable and SLURM directives at the top of each script before running on a different system.

A separate conda environment was created for each major tool across the pipeline (e.g. `envs/qc`, `envs/assembly`, `envs/medaka`, `envs/metabat2`, `envs/checkm2`, `envs/gtdbtk`, `envs/bakta`, `envs/amrfinder`, `envs/resfinder`, `envs/rgi`, `envs/mob_suite`, `envs/isescan`, `envs/integronfinder`, `envs/genomad`) to avoid dependency conflicts between tools with incompatible package requirements. The relevant `install_*.sh` script for each tool documents its environment setup.

Note on step numbering: the original scripts carry step/phase labels (e.g. "Step 2.1", "Phase 4b") from the authors' internal lab notebook, which numbered steps in chronological run order rather than by pipeline category. Those internal labels are not reproduced below, since they don't align with the section numbers in this README; each script's purpose is instead described directly.

## 1. Read QC and adapter trimming

| Script | Purpose |
|---|---|
| `run_qc.sh` | Initial read QC with NanoPlot |
| `run_qc_cont.sh` | Continuation/re-run of NanoPlot QC |
| `run_porechop_abi.sh` | Adapter trimming with Porechop_ABI v0.5.1 (SLURM array, 63 samples/controls) |
| `run_chopper.sh` | Quality and length filtering of adapter-trimmed reads with Chopper (SLURM array) |
| `run_multiqc.sh` | Aggregates all NanoPlot reports into one MultiQC summary |

## 2. Host depletion and decontamination

| Script | Purpose |
|---|---|
| `download_host_genomes.sh` | Downloads/prepares host reference genomes for depletion |
| `run_host_depletion.sh` | Removes host reads via minimap2 + samtools (SLURM array) |
| `run_get_contaminant_taxids.sh` | Looks up NCBI taxIDs for contaminant species flagged by Decontam |
| `run_decontam.sh` / `decontam_analysis.R` | Identifies and flags reagent/lab contaminants using negative controls (AE-buffer, EB-neg) with the `decontam` R package (prevalence method, threshold 0.5) |
| `run_extract_kraken.sh` | Removes contaminant reads from filtered FASTQs, producing the final clean reads (SLURM array) |

## 3. Taxonomic classification

| Script | Purpose |
|---|---|
| `download_kraken2_db.sh` | Downloads the Kraken2 PlusPF database |
| `run_kraken2.sh` | Taxonomic classification of reads with Kraken2 + Bracken (SLURM array, 63 jobs) |
| `run_kraken2_assembly.sh` | Classifies polished metagenomic assemblies (contigs) with Kraken2 |
| `run_bracken_AEbuffer_fix.sh` | One-off re-run of Bracken for the AE-buffer negative control after an initial job failure |
| `run_krakentools.sh` | Merges all 63 Bracken outputs into one combined abundance matrix |
| `run_krona.sh` | Interactive taxonomic composition visualisation with Krona |

## 4. Diversity analysis

| Script | Purpose |
|---|---|
| `run_phyloseq.sh` / `phyloseq_analysis.R` | Alpha and beta diversity analysis (richness, evenness, Bray-Curtis ordination, relative abundance) across the 61 real samples |
| `run_vegan.sh` / `vegan_analysis.R` | PERMANOVA, betadispersion, and NMDS ordination on Bray-Curtis dissimilarities |

## 5. De novo assembly and polishing

| Script | Purpose |
|---|---|
| `run_assembly.sh` | Per-sample metagenome assembly with Flye (SLURM array, 61 samples) |
| `install_medaka.sh` | Installs the GPU-accelerated Medaka polishing environment |
| `run_medaka.sh` | Assembly polishing with Medaka (GPU, SLURM array) |
| `run_minimap2_coverage.sh` | Maps clean reads back to the polished assembly to generate per-contig coverage depth |
| `check_coverage.sh` | Validates that coverage depth/BAM outputs exist and are non-empty for every sample |

## 6. MAG binning, quality control, taxonomy, and annotation

| Script | Purpose |
|---|---|
| `run_metabat2.sh` | Bins contigs into Metagenome-Assembled Genomes (MAGs) with MetaBAT2 |
| `run_checkm2.sh` | Assesses MAG completeness and contamination with CheckM2 |
| `download_gtdbtk.sh` | Downloads the GTDB-Tk reference database |
| `run_gtdbtk.sh` | Assigns taxonomy to medium/high-quality MAGs with GTDB-Tk |
| `run_bakta.sh` | Annotates quality MAGs with Bakta |
| `run_bakta_plasmid_annot.sh` | Annotates individual MOB-suite plasmid FASTAs with Bakta |

## 7. AMR detection

| Script | Purpose |
|---|---|
| `install_amrfinder.sh` | Installs AMRFinderPlus 4.2.7 |
| `install_resfinder.sh` / `build_kma_index.sh` | Installs ResFinder 4.7.2 + KMA aligner and builds the KMA database index |
| `download_amr_databases.sh` | Downloads current AMRFinderPlus, CARD/RGI, and ResFinder databases |
| `rgi_load_card.sh` | Loads the CARD reference database into RGI |
| `run_resfinder.sh` | Per-sample read-based AMR detection with ResFinder 4.7.2 |
| `run_amrfinder.sh` | Per-sample assembly-based AMR detection with AMRFinderPlus 4.2.7 |
| `run_rgi.sh` | Per-sample assembly-based AMR detection with RGI 6.0.5 / CARD 4.0.1 |
| `run_abricate.sh` | Per-sample assembly-based AMR/virulence screening with Abricate 1.0.4 (7 databases) |
| `run_resfinder_assembly.sh` | Per-sample assembly-based AMR detection with ResFinder 4.7.2 |
| `amr_analysis.R` | Combines and analyses results from all 5 AMR detection approaches above |

## 8. MGE detection

| Script | Purpose |
|---|---|
| `install_mob_suite.sh` | Installs MOB-suite |
| `install_mge_mag_envs.sh` | Installs remaining conda environments for the MGE/MAG pipeline |
| `download_databases.sh` | Downloads required reference databases for the MGE/MAG phases |
| `run_mobsuite.sh` | Classifies contigs as chromosome/plasmid and types plasmids with MOB-suite |
| `run_isescan.sh` | Identifies insertion sequence (IS) elements with ISEScan (SLURM array) |
| `run_isescan_array30_24hr.sh` | One-off re-run of ISEScan for a single sample (array index 30) with extended walltime/CPUs |
| `run_integronfinder.sh` | Detects integrons with IntegronFinder 2 |
| `run_genomad.sh` | Classifies assemblies into plasmid/virus/chromosome with geNomad |

## 9. AMR–MGE co-localisation and combined files

| Script | Purpose |
|---|---|
| `run_build_contig_mag_map.sh` | Builds a contig-to-MAG mapping table from quality-filtered CheckM2 bins, joined with GTDB-Tk taxonomy |
| `phase5_01_combine_args.R` | Combines ARG calls from RGI, Abricate, and AMRFinderPlus into one table |
| `phase5_02_combine_mges.R` | Combines MGE calls from ISEScan, IntegronFinder, geNomad, and MOB-suite into one table |
| `phase5_03_colocalisation.R` | For each ARG hit, finds MGEs on the same contig within 10 kb; joins with MAG taxonomy and sample metadata |
| `phase5_05_combined_csvs.py` / `run_phase5_05.sh` | Generates the combined AMR + MGE + taxonomy + metadata CSVs (one per AMR database × MGE tool pairing), used to produce the supplementary data files |

## Suggested execution order

1. QC and trimming (§1)
2. Host depletion (§2, `download_host_genomes.sh` → `run_host_depletion.sh`)
3. Taxonomic classification (§3) → decontamination (§2, `run_decontam.sh`) → diversity analysis (§4) → clean-read extraction (§2, `run_extract_kraken.sh`)
4. Assembly and polishing (§5)
5. MAG binning, QC, taxonomy, annotation (§6)
6. AMR detection (§7) and MGE detection (§8), run in parallel on polished assemblies
7. Co-localisation and combined files (§9)
