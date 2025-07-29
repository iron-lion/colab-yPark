import sys
from src import *


SPECIES = "Killifish"
#SPECIES = "C.elegans"
GENOME_FASTA = "Nothobranchius_furzeri.Nfu_20140520.dna_sm.toplevel.fa"
#GENOME_FASTA = "Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa"
REFSEQ_GTF = "Nothobranchius_furzeri.Nfu_20140520.114.gtf"
#REFSEQ_GTF = "Caenorhabditis_elegans.WBcel235.114.gtf"
JASPAR_PFM_FILE = "JASPAR2022_combined_matrices_Vertebrata_pfm.txt"
#JASPAR_PFM_FILE = "20250728121304_JASPAR2024_combined_matrices_220990_pfm.txt"
GENES_OF_INTEREST = ["ENSNFUG00015011653" ]
#GENES_OF_INTEREST = ["WBGene00022842", "WBGene00014148", ]
PROMOTER_LENGTH = 2000
PWM_SCORE_THRESHOLD = 0.8

# --- Main Execution ---
if __name__ == "__main__":
    # --- 1. Load Genome ---
    genome_data = load_genome(GENOME_FASTA)
    if not genome_data:
        print("Error: Genome data could not be loaded. Exiting.")
        sys.exit(1)

    # --- 2. Load JASPAR Motifs ---
    jaspar_motifs_obj = load_jaspar_pfm(JASPAR_PFM_FILE)
    if not jaspar_motifs_obj:
        print("Error: No JASPAR motifs loaded from the PFM file. Check file path and format. Exiting.")
        sys.exit(1)

    # --- 3. Get Promoter Coordinates for TARGET Genes ---
    print(f"\n--- Processing Target Genes: {len(GENES_OF_INTEREST)} genes specified ---")
    target_promoter_coords = get_promoter_coordinates(REFSEQ_GTF, GENES_OF_INTEREST, PROMOTER_LENGTH)

    if not target_promoter_coords:
        print("Warning: No promoter coordinates found for target genes. Cannot scan for motifs. Exiting.")
        sys.exit(1)

    genes_to_process_target = {gene_id: coords for gene_id, coords in target_promoter_coords.items()
                               if coords['chr'] in genome_data}

    if not genes_to_process_target:
        print("Warning: No target genes to process after filtering by genome presence. Exiting.")
        sys.exit(1)

    # --- 4. Extract Promoter Sequences for TARGET Genes ---
    target_promoter_seqs = extract_promoter_sequences(genome_data, genes_to_process_target)
    if not target_promoter_seqs:
        print("Warning: No promoter sequences could be extracted for target genes. Exiting.")
        sys.exit(1)

    # --- 5. Scan Target Promoter Regions for Known Motifs ---
    print("\n--- Scanning Target Promoter Regions for Known Motifs ---")
    target_motif_hits = scan_for_hits(target_promoter_seqs, jaspar_motifs_obj, PWM_SCORE_THRESHOLD)
    if target_motif_hits:
        sorted_target_hits = sorted(target_motif_hits.items(), key=lambda item: item[1]['count'], reverse=True)
        for motif_name, data in sorted_target_hits:
            print(f"Motif: {motif_name}, Hits: {data['count']}, Genes with hits: {len(data['genes'])}")
    else:
        print("No motifs found in target promoter regions with the given threshold.")

    # --- Step 6: Perform Statistical Enrichment Analysis ---
    print("\n--- Performing Motif Enrichment Analysis (Fisher's Exact Test) ---")
    enrichment_results = []

    # --- Run again for background --- 
    all_gene_promoter_coords = get_promoter_coordinates(REFSEQ_GTF, [], PROMOTER_LENGTH)
    valid_promoter_coords = {gene_id: coords for gene_id, coords in all_gene_promoter_coords.items()
                                 if coords['chr'] in genome_data}
    background_promoter_seqs = extract_promoter_sequences(genome_data, valid_promoter_coords)
    background_motif_hits = scan_for_hits(background_promoter_seqs, jaspar_motifs_obj, PWM_SCORE_THRESHOLD)

    total_target_genes = len(GENES_OF_INTEREST)
    total_background_genes = len(valid_promoter_coords)

    import dill, pickle
    with open(f"background_{SPECIES}_PWM_{PWM_SCORE_THRESHOLD}.pkl", 'wb') as outp:
        dill.dump(background_motif_hits, outp, pickle.HIGHEST_PROTOCOL)
