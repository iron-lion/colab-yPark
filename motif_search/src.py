import sys, math
from collections import defaultdict, Counter
from Bio import SeqIO 
from Bio import motifs
from Bio.Seq import Seq
import tqdm


def calculate_pseudocounts(motif):
    """Calculate pseudocounts.

    Computes the root square of the total number of sequences multiplied by
    the background nucleotide.
    """
    alphabet = motif.alphabet                                                                                                                                                             
    background = motif.background

    # It is possible to have unequal column sums so use the average
    # number of instances.
    total = 0
    for i in range(motif.length):
        total += sum(motif.counts[letter][i] for letter in alphabet)

    avg_nb_instances = total / motif.length
    sq_nb_instances = math.sqrt(avg_nb_instances)

    if background:
        background = dict(background)
    else:
        background = dict.fromkeys(sorted(alphabet), 1.0)

    total = sum(background.values())
    pseudocounts = {}

    for letter in alphabet:
        background[letter] /= total
        pseudocounts[letter] = sq_nb_instances * background[letter]

    return pseudocounts



# --- 1. Load Genome ---
def load_genome(fasta_file):
    """
    Loads the genome sequence from a FASTA file into a dictionary.
    Keys are chromosome names, values are SeqIO objects.
    """
    print(f"Loading genome from {fasta_file}...")
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    print("Genome loaded.")
    return genome


# --- 2. Parse Gene Annotations and Get Promoter Coordinates ---
def get_promoter_coordinates(gtf_file, genes_of_interest, promoter_length):
    """
    Parses a GTF/GFF file to get gene coordinates and calculate promoter regions.
    Returns a dictionary: {gene_id: {'chr': str, 'start': int, 'end': int, 'strand': str}}
    Note: This is a simplified parser. Real GTF/GFF parsing can be complex.
    """
    print(f"Parsing annotations from {gtf_file} and identifying promoter regions...")
    gene_promoters = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            # Assuming GTF/GFF format: seqname source feature start end score strand frame attributes
            seqname = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            if feature_type == "gene": # Or "CDS" depending on what you define as gene start
                # Extract gene ID from attributes (this part is highly dependent on GTF/GFF format)
                # Example for Ensembl GTF: gene_id "WBGene0000001"; gene_version "1";
                #gene_id_match = [attr.strip().split(' ')[1].strip('"') for attr in attributes.split(';') if "gene_id" in attr]
                gene_id_match = attributes.split(';')[0].replace('gene_id ','').replace('"','')
                if not gene_id_match:
                    continue
                gene_id = gene_id_match

                if genes_of_interest and gene_id not in genes_of_interest:
                    continue

                promoter_start = -1
                promoter_end = -1

                if strand == '+':
                    # Promoter is upstream of the start site
                    promoter_start = start - promoter_length
                    promoter_end = start - 1 # Or just 'start' if you want to include the TSS base
                elif strand == '-':
                    # Promoter is upstream of the end site (for reverse strand)
                    promoter_start = end + 1 # Or just 'end'
                    promoter_end = end + promoter_length

                # Ensure promoter coordinates are not negative
                promoter_start = max(0, promoter_start)

                if promoter_start != -1 and promoter_end != -1:
                    gene_promoters[gene_id] = {
                        'chr': seqname,
                        'start': promoter_start,
                        'end': promoter_end,
                        'strand': strand
                    }
    print(f"Found promoter regions for {len(gene_promoters)} genes.")
    return gene_promoters


# --- 3. Extract Promoter Sequences ---
def extract_promoter_sequences(genome, gene_promoters):
    """
    Extracts the DNA sequence for each promoter region.
    Returns a dictionary: {gene_id: 'sequence'}
    """
    print("Extracting promoter sequences...")
    promoter_sequences = {}
    for gene_id, coords in gene_promoters.items():
        chrom = coords['chr']
        start = coords['start']
        end = coords['end']
        strand = coords['strand']

        if chrom not in genome:
            print(f"Warning: Chromosome '{chrom}' not found in genome for gene {gene_id}. Skipping.")
            continue

        # Biopython Seq objects are 0-indexed, slice is [start:end]
        # Adjusting for 1-based GTF to 0-based Python slicing
        seq_obj = genome[chrom].seq[start:end]

        if strand == '-':
            seq_obj = seq_obj.reverse_complement()
        promoter_sequences[gene_id] = str(seq_obj).upper() # Convert to string and uppercase
    print(f"Extracted sequences for {len(promoter_sequences)} promoters.")
    return promoter_sequences


# --- Revised Common K-mer to Known Motif Comparison (using PFM/PWM scanning) ---
def compare_common_kmers_with_pfms(common_kmers, jaspar_motifs, promoter_sequences, score_threshold):
    """
    Compares found common k-mers by scanning promoter regions with JASPAR PWMs.
    This is a more robust way to link found k-mers to known motifs.
    It identifies which promoter regions contain matches to known motifs.
    """
    print(f"\n--- Scanning promoter regions for known JASPAR motifs (score threshold: {score_threshold*100:.1f}%) ---")
    motif_hits_in_promoters = defaultdict(lambda: defaultdict(list)) # {motif_name: {gene_id: [(pos, score, seq, strand), ...]}}
    overall_motif_counts = defaultdict(int) # {motif_name: total occurrences}

    for motif_name, motif_obj in jaspar_motifs.items():
        print(f"Scanning for motif: {motif_name} (Length: {len(motif_obj)})")
        for gene_id, promoter_seq in promoter_sequences.items():
            # Ensure promoter sequence is long enough to contain the motif
            if len(promoter_seq) >= len(motif_obj):
                matches = scan_sequence_with_pwm(promoter_seq, motif_obj, score_threshold)
                if matches:
                    motif_hits_in_promoters[motif_name][gene_id].extend(matches)
                    overall_motif_counts[motif_name] += len(matches)

    print("\n--- Summary of Known Motif Hits in Promoters ---")
    if overall_motif_counts:
        for motif_name, count in sorted(overall_motif_counts.items(), key=lambda item: item[1], reverse=True):
            continue#print(f"Motif: {motif_name} (Total hits: {count})")
            # You can uncomment to print individual hits:
            # for gene_id, hits in motif_hits_in_promoters[motif_name].items():
            #     for pos, score, seq, strand in hits:
            #         print(f"  Gene {gene_id}: Match '{seq}' at position {pos} (strand {strand}, score {score:.2f})")
    else:
        print("No known motifs found in promoter regions with the given threshold.")

    # You could also link your *found k-mers* directly to the motifs
    # by checking if your exact k-mers (from common_kmers) achieve a high score
    # against the PWMs, but the primary use of PWMs is to scan a sequence.
    # For now, we'll just report the found motif occurrences.
    return motif_hits_in_promoters


# --- scan_sequence_with_pwm (The latest corrected version with .min_score()/.max_score())
def scan_sequence_with_pwm(sequence, motif_obj, threshold_percent):
    """
    Scans a DNA sequence for matches to a motif using its PWM/PSSM.
    Returns a list of (start_position, score, matching_sequence) tuples
    for matches exceeding the threshold, considering both strands.
    """
    matches = []

    pssm = motif_obj.pssm

    min_score = pssm.min
    max_score = pssm.max

    if (max_score - min_score) > 0:
        score_cutoff = (max_score - min_score) * threshold_percent + min_score
    else:
        # print(f"Warning: Motif '{motif_obj.name}' has min_score == max_score. Skipping scanning.", file=sys.stderr)
        return [] # Return empty if no valid score range

    # Scan forward strand
    for position, score in pssm.search(str(sequence), threshold=score_cutoff):
        if position + len(motif_obj) <= len(sequence):
            matched_seq = sequence[position : position + len(motif_obj)]
            matches.append((position, score, matched_seq, '+'))

    # Scan reverse complement strand
    rev_comp_sequence_obj = Seq(sequence).reverse_complement()
    rev_comp_sequence_str = str(rev_comp_sequence_obj)

    for position, score in pssm.search(rev_comp_sequence_str, threshold=score_cutoff):
        if position + len(motif_obj) <= len(rev_comp_sequence_str):
            matched_seq = rev_comp_sequence_str[position : position + len(motif_obj)]
            matches.append((position, score, matched_seq, '-'))

    return matches


def scan_for_hits(promoter_seqs, jaspar_motifs, threshold):
    """Helper function to scan a set of promoters and return motif hits."""
    motif_hits = defaultdict(lambda: {'count': 0, 'genes': set()})
    for motif_name, motif_obj in tqdm.tqdm(jaspar_motifs.items()):
        if not motif_obj or len(motif_obj) == 0:
            continue
        try:
            _ = motif_obj.pssm
        except Exception as e:
            continue # Already warned in parsing step

        for gene_id, promoter_seq in promoter_seqs.items():
            if len(promoter_seq) >= len(motif_obj):
                matches = scan_sequence_with_pwm(promoter_seq, motif_obj, threshold)
                if matches:
                    motif_hits[motif_name]['count'] += len(matches)
                    motif_hits[motif_name]['genes'].add(gene_id)
    return motif_hits
  

def load_jaspar_pfm(file_path):
    """
    Parses a JASPAR-like PFM file with a custom floating-point format.
    Assumes 4 rows (A, C, G, T) for each motif.
    Converts float counts to integers by rounding for Bio.motifs.Motif object.
    Crucially, passes PSSM_PSEUDOCOUNTS to the Motif constructor.
    """
    print(f"Starting custom parsing of PFM file: {file_path}")
    jaspar_motifs = {}
    
    current_motif_name = None
    current_pfm_rows = [] # Will temporarily store [A_counts, C_counts, G_counts, T_counts]

    try:
        with open(file_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                original_line = line # Keep original for error messages
                line = line.strip()

                if not line: # Skip empty lines
                    continue

                if line.startswith('>'):
                    # A new motif is starting. First, try to process the previous one.
                    if current_motif_name: # Check if there was a previous motif being built
                        if len(current_pfm_rows) == 4:
                            try:
                                integer_pfm_dict = {
                                    'A': [int(round(x)) for x in current_pfm_rows[0]],
                                    'C': [int(round(x)) for x in current_pfm_rows[1]],
                                    'G': [int(round(x)) for x in current_pfm_rows[2]],
                                    'T': [int(round(x)) for x in current_pfm_rows[3]],
                                }
                                
                                # Check if any row is empty (e.g., due to parsing errors resulting in empty lists)
                                if any(not row for row in integer_pfm_dict.values()):
                                    raise ValueError("One or more PFM rows are empty after conversion.")

                                # THE FIX HERE: Pass pseudocounts to the Motif constructor
                                motif_obj = motifs.Motif(alphabet="ACGT", counts=integer_pfm_dict)
                                motif_obj.pseudocounts = calculate_pseudocounts(motif_obj)

                                motif_obj.name = current_motif_name
                                jaspar_motifs[current_motif_name] = motif_obj
                                # print(f"Successfully parsed motif: {current_motif_name}") # Debugging line
                            except Exception as e:
                                print(f"Error processing motif '{current_motif_name}' (started around line {line_num - len(current_pfm_rows) -1}): Incomplete or malformed data preventing motif creation: {e}", file=sys.stderr)
                        else:
                            print(f"Warning: Motif '{current_motif_name}' (started around line {line_num - len(current_pfm_rows) -1}) has {len(current_pfm_rows)} data rows instead of 4. Skipping.", file=sys.stderr)

                    # Reset for the new motif
                    parts = line[1:].split(' ')
                    if len(parts) > 1:
                        current_motif_name = parts[1]
                    else:
                        current_motif_name = parts[0]
                    current_pfm_rows = [] # Clear rows for the new motif
                    # print(f"Detected new motif header: {current_motif_name}") # Debugging line

                else:
                    # This is a data line (A, C, G, T counts)
                    if current_motif_name is None:
                        print(f"Warning: Data line found before any motif header at line {line_num}: '{original_line.strip()}'. Skipping.", file=sys.stderr)
                        continue

                    if len(current_pfm_rows) >= 4:
                        print(f"Warning: Motif '{current_motif_name}' at line {line_num}: More than 4 data rows detected. Skipping extra row: '{original_line.strip()}'.", file=sys.stderr)
                        continue
                    
                    try:
                        counts = [float(x) for x in line.split()]
                        # Essential check: ensure 'counts' list is not empty after splitting
                        if not counts:
                            print(f"Warning: Motif '{current_motif_name}' at line {line_num}: Data line is empty after splitting. Skipping: '{original_line.strip()}'.", file=sys.stderr)
                            continue
                        
                        current_pfm_rows.append(counts)
                    except ValueError as e:
                        print(f"Error parsing data line for motif '{current_motif_name}' at line {line_num}: '{original_line.strip()}'. Error: {e}. Skipping this line.", file=sys.stderr)
                    
        # After the loop, process the very last motif in the file
        if current_motif_name: # Check if there was any motif at all
            if len(current_pfm_rows) == 4:
                try:
                    integer_pfm_dict = {
                        'A': [int(round(x)) for x in current_pfm_rows[0]],
                        'C': [int(round(x)) for x in current_pfm_rows[1]],
                        'G': [int(round(x)) for x in current_pfm_rows[2]],
                        'T': [int(round(x)) for x in current_pfm_rows[3]],
                    }
                    if any(not row for row in integer_pfm_dict.values()):
                        raise ValueError("One or more PFM rows are empty after conversion in final motif.")

                    # THE FIX HERE: Pass pseudocounts to the Motif constructor for the last motif
                    motif_obj = motifs.Motif(alphabet="ACGT", counts=integer_pfm_dict)
                    motif_obj.pseudocounts = calculate_pseudocounts(motif_obj)

                    motif_obj.name = current_motif_name
                    jaspar_motifs[current_motif_name] = motif_obj
                    # print(f"Successfully parsed final motif: {current_motif_name}") # Debugging line
                except Exception as e:
                    print(f"Error processing final motif '{current_motif_name}': Incomplete or malformed data preventing motif creation: {e}", file=sys.stderr)
            else:
                print(f"Warning: Final motif '{current_motif_name}' has {len(current_pfm_rows)} data rows instead of 4. Skipping.", file=sys.stderr)

    except FileNotFoundError:
        print(f"Critical Error: PFM file not found at {file_path}", file=sys.stderr)
        return {} # Return empty dict on critical file error
    except Exception as e:
        print(f"An unexpected critical error occurred during file reading: {e}", file=sys.stderr)
        return {} # Return empty dict on critical file error

    print(f"Finished custom PFM parsing. Loaded {len(jaspar_motifs)} motifs.")
    return jaspar_motifs
