from collections import defaultdict
import pandas as pd
from pathlib import Path
import json
from tqdm import tqdm
import math
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

# Threshold for Chi-square filtering
CHI_SQUARE_THRESHOLD = 3.841

# Global background frequency distribution for amino acids
AA_BACKGROUND = {
    "X": 0.0,
    "S": 0.08621733200464735,
    "L": 0.07651555159603048,
    "P": 0.07203973199729977,
    "T": 0.06922675985837555,
    "A": 0.06405039516162363,
    "G": 0.0596202682782342,
    "V": 0.06135659146120252,
    "R": 0.05800586887338755,
    "D": 0.05276268719088534,
    "N": 0.051973971225069915,
    "Y": 0.0454190182724939,
    "H": 0.05010516212877539,
    "F": 0.04232959372517324,
    "I": 0.0396433672087032,
    "M": 0.033543274904825976,
    "K": 0.03183794929256655,
    "Q": 0.03255870021445727,
    "E": 0.03274491525034556,
    "W": 0.030877254212225442,
    "C": 0.009171607143677185
}


def create_kmer_dict(custom_db_path: str, kmer_length: int, wild_cards=False, wild_cards_pattern=[]):
    '''
    Creates a temporary dictionary to act as a database using a file as well as a protein table to reference
    the dictionary. Computes tiling for any kmer size between 4 and 8.
    INPUT
    custom_db_path: A string of a filepath where the custom proteome file is saved (Fasta file format)
    kmer_length: an integer of the length of the kmer to be created
    wild_cards: a boolean that describes whether to generate wildcards in the output
    wild_cards_pattern: a list of indexes (0-indexed) for the wildcards in the kmer
    OUTPUT
    protein_table: A 2d array representing each protein in the custom proteom table. metadata [proteinID,header,len(sequence)]
    seq_dict: A dictionary of key = kmerString and value is a list of tuples in format (proteinId, positionOfKmer)
    '''
    if kmer_length < 4 or kmer_length > 8:
        raise ValueError("Size must be between 4 and 8.")

    seq_dict = collections.defaultdict(list)
    protein_table = []
    protein_id = 0

    def process_sequence(sequence, protein_id, seq_dict, wild_cards):
        # Process the sequence to find xmers of the given size
        for i in range(len(sequence) - kmer_length + 1):
            cur_seq = sequence[i:i+kmer_length]
            seq_dict[cur_seq].append([protein_id, i])
            if wild_cards:
                add_wildcards(cur_seq, protein_id, i, seq_dict, kmer_length, wild_cards_pattern)

    def add_wildcards(cur_seq, protein_id, i, seq_dict, size, pattern=None):
        """
        Adds wildcard sequences to the dictionary based on the given pattern or all possible combinations.

        Parameters:
        - ccur_seq: The current sequence being processed.
        - protein_id: The ID of the protein being processed.
        - i: The position in the sequence.
        - seq_dict: The dictionary to store sequences.
        - size: The size of the sequence.
        - pattern: Optional list of indexes for wildcards (e.g., [3] for "ABCXE").
        """
        num_wildcards = size - 4  # Minimum of 4 fixed positions
        positions = range(1, size - 1)  # Wildcards cannot be at the first or last position

        if pattern:
            # Replace the specified indexes in the pattern with 'X'
            wild_card_seq = list(cur_seq)  # Convert to list for mutability
            for pos in pattern:
                if 0 < pos < size - 1:  # Ensure the index is within bounds
                    wild_card_seq[pos] = 'X'
            wild_card_seq = ''.join(wild_card_seq)  # Convert back to string
            seq_dict[wild_card_seq].append([protein_id, i])
        else:
            # Generate all combinations of positions to place wildcards
            for wildcard_positions in combinations(positions, num_wildcards):
                wild_card_seq = list(cur_seq)
                for pos in wildcard_positions:
                    wild_card_seq[pos] = 'X'
                wild_card_seq = ''.join(wild_card_seq)
                seq_dict[wild_card_seq].append([protein_id, i])

    for file in os.listdir(custom_db_path):
        with open(os.path.join(custom_db_path, file), 'r') as f:
            sequence = ''
            header = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):  # This is a header line
                    if sequence:  # Save the previous protein's data
                        protein_table.append([protein_id, header, len(sequence)])
                        process_sequence(sequence, protein_id, seq_dict, wild_cards)
                    protein_id += 1
                    header = line
                    sequence = ''
                else:
                    sequence += line.replace('X', 'Q')  # Build the sequence string

            if sequence:  # Don't forget to process the last sequence in the file
                protein_table.append([protein_id, header, len(sequence)])
                process_sequence(sequence, protein_id, seq_dict, wild_cards)

    return protein_table, seq_dict


def split_input_by_group(input_path: Path):
    """
    Splits a single Excel/CSV file containing mixed AD_ and NC_ columns into two separate files:
    one for AD (Diseased) and one for NC (Normal Control).
    
    Parameters:
    - input_path: Path to the original .csv or .xlsx file
    
    Output:
    - Two files are saved with suffixes `_AD` and `_NC` respectively, preserving the original extension.
    """
    if input_path.suffix == '.xlsx':
        df = pd.read_excel(input_path)
    elif input_path.suffix == '.csv':
        df = pd.read_csv(input_path)
    else:
        raise ValueError("Unsupported file format.")

    # Collect column pairs
    ad_cols = []
    nc_cols = []

    i = 0
    while i < df.shape[1] - 1:  # ensure i+1 exists
        header = str(df.columns[i]).strip()

        # Skip empty or auto-generated column
        if header == "" or header.startswith("Unnamed:"):
            print(f"Skipping empty or invalid column at index {i}: {header}")
            i += 1  # Only move by 1, so we don't skip the next one
            continue

        if header.startswith("AD_"):
            ad_cols.extend(df.columns[i:i+2])
            i += 2
        elif header.startswith("NC_"):
            nc_cols.extend(df.columns[i:i+2])
            i += 2
        else:
            print(f"Skipping unrecognized prefix: {header}")
            i += 1  # Unknown prefix — move forward by 1 and keep looking

    # Build new DataFrames
    df_ad = df[ad_cols]
    df_nc = df[nc_cols]
    for df_split in [df_ad, df_nc]:
        for j in range(1, df_split.shape[1], 2):
            df_split.columns.values[j] = "count"
    # Construct output paths
    base_name = input_path.stem
    suffix = input_path.suffix
    output_ad = input_path.parent / f"{base_name}_AD{suffix}"
    output_nc = input_path.parent / f"{base_name}_NC{suffix}"

    # Save
    if suffix == '.xlsx':
        df_ad.to_excel(output_ad, index=False)
        df_nc.to_excel(output_nc, index=False)
    else:
        df_ad.to_csv(output_ad, index=False)
        df_nc.to_csv(output_nc, index=False)

    print(f"AD file saved to {output_ad}")
    print(f"NC file saved to {output_nc}")
    return output_ad, output_nc

def apply_chi_square_filter(kmer_dict, aa_background, threshold, output_csv=False, filename=None, normalize_expected=True):
    """
    Applies chi-square filtering to a k-mer frequency dictionary.

    Parameters:
    - kmer_dict: Dictionary with k-mer strings as keys and their total frequency as values.
    - aa_background: Dictionary of expected amino acid frequencies (e.g., AA_BG global metadata).
    - threshold: Chi-square threshold value; default is 3.841 for p=0.05, df=1.
    - output_csv: If True, save filtered results to a CSV file.
    - filename: Output filename for saving the filtered CSV (required if output_csv is True).
    - normalize_expected: If True, divide expected frequency by total number of unique k-mers.

    Returns:
    - A filtered dictionary containing only k-mers that pass the chi-square test.
    """

    filtered = {}
    total_count = sum(kmer_dict.values())  # Total observed frequency

    for kmer, observed in kmer_dict.items():
        expected = 0
        # Sum up expected frequency based on amino acid background
        for aa in kmer:
            expected += aa_background.get(aa, 0) * total_count
        # Optionally normalize expected value by number of unique k-mers
        if normalize_expected:
            expected /= len(kmer_dict)
        if expected == 0:
            continue  # Avoid division by zero
        chi2 = (observed - expected) ** 2 / expected
        if chi2 >= threshold:
            filtered[kmer] = observed
    if output_csv and filename:
        pd.DataFrame(sorted(filtered.items()), columns=["Sequence", "Frequency"]).to_csv(filename, index=False)

    return filtered

# Step 1: Tile and count k-mers from input Excel or CSV
def tile_patient_file(path, 
                      kmer_length=4, 
                      wildcard_positions=[], 
                      return_individual_dicts=False, 
                      apply_chi_square=True, 
                      aa_background=AA_BACKGROUND, 
                      chi_square_threshold=3.841, 
                      normalize_expected=True):
    if path.suffix == '.xlsx':
        df = pd.read_excel(path)
    elif path.suffix == '.csv':
        df = pd.read_csv(path)
    else:
        raise ValueError("Unsupported file format.")

    patient_dicts = {}
    total_kmer_dict = defaultdict(int)

    for col in tqdm(range(0, df.shape[1], 2), total=df.shape[1] // 2, desc=f"Processing {Path(path).stem} patients"):
        seq_col = df.columns[col]
        count_col = df.columns[col + 1]
        patient_name = seq_col
        single_patient_dict = defaultdict(int)

        for seq, count in zip(df[seq_col], df[count_col]):
            if pd.isna(seq) or pd.isna(count):
                continue
            try:
                count = int(count)
                for i in range(len(seq) - kmer_length + 1):
                    kmer = list(seq[i:i + kmer_length])
                    for pos in wildcard_positions:
                        if 0 <= pos < len(kmer):
                            kmer[pos] = 'X'
                    kmer_str = ''.join(kmer)
                    single_patient_dict[kmer_str] += count
                    total_kmer_dict[kmer_str] += count
            except Exception as e:
                print(f"Skipping invalid row: {seq}, {count} ({e})")
                continue

        patient_dicts[patient_name] = single_patient_dict

    if apply_chi_square:
        if aa_background is None:
            raise ValueError("aa_background must be provided.")
        total_kmer_dict = apply_chi_square_filter(
            total_kmer_dict,
            aa_background=aa_background,
            threshold=chi_square_threshold,
            output_csv=False,
            normalize_expected=normalize_expected
        )
    
    # if output_csv:
    #     df_out = pd.DataFrame(sorted(total_kmer_dict.items()), columns=["Sequence", "Frequency"])
    #     wildcard_str = ''.join(str(pos) for pos in wildcard_positions)
    #     suffix = "_filtered" if apply_chi_square else "_raw"
    #     output_filename = str(path.parent / f"{Path(path).stem}_{kmer_length}mers_[{wildcard_str}]{suffix}.csv")
    #     df_out.to_csv(output_filename, index=False)

    # Return both patient_dicts and the filtered kmer set
    if return_individual_dicts:
        return patient_dicts, set(total_kmer_dict.keys())
    else:
        return total_kmer_dict


def build_kmer_matrix(patient_dicts, normalize=True, output_path=None):
    """
    Builds a patient-by-kmer matrix with optional normalization.
    If output_path is provided, writes the matrix to that file.
    
    Parameters:
        patient_dicts (dict): Dictionary of patient_name → {kmer: count}
        normalize (bool): If True, normalize each patient's counts to proportions
        output_
        path (str or Path): If set, writes the matrix to CSV at this path

    Returns:
        matrix (pd.DataFrame): kmer x patient matrix
    """
    all_kmers = set()
    for pdict in patient_dicts.values():
        all_kmers.update(pdict.keys())

    all_kmers = sorted(all_kmers)
    columns_dict = {}

    for patient, kmers in tqdm(patient_dicts.items(), desc="Building k-mer matrix"):
        total = sum(kmers.values()) if normalize else 1
        columns_dict[patient] = pd.Series({kmer: freq / total for kmer, freq in kmers.items()})

    matrix = pd.DataFrame(columns_dict).fillna(0)
    matrix = matrix.reindex(sorted(matrix.index))  # Ensure consistent row order

    if output_path:
        matrix.to_csv(output_path)
        print(f"K-mer matrix saved to {output_path}")

    return matrix


def run_mannwhitney(matrix, pos_cols, neg_cols, output_filename="mannwhitney_results.csv", fdr_method='fdr_tsbky'):
    """
    Performs Mann–Whitney U test for each k-mer across positive and negative samples.
    Computes Mean Rank Difference (MRD) and splits results into:
      - Full result CSV
      - AD-elevated (MRD > 0)
      - NC-elevated (MRD < 0)
    """

    pvals = []
    mean_rank_diffs = []
    kmers = []

    for kmer in tqdm(matrix.index, desc="Running Mann–Whitney U test"):
        pos_vals = matrix.loc[kmer, pos_cols].values
        neg_vals = matrix.loc[kmer, neg_cols].values
        try:
            # Combine and rank values
            combined = list(pos_vals) + list(neg_vals)
            ranks = rankdata(combined)

            rank_pos = ranks[:len(pos_vals)]
            rank_neg = ranks[len(pos_vals):]

            mean_rank_diff = rank_pos.mean() - rank_neg.mean()

            # MWU test
            stat, p = mannwhitneyu(pos_vals, neg_vals, alternative='two-sided')

        except Exception as e:
            print(f"Skipping {kmer}: {e}")
            p = 1.0
            mean_rank_diff = 0.0

        kmers.append(kmer)
        pvals.append(p)
        mean_rank_diffs.append(mean_rank_diff)

    # FDR correction
    fdr_adjusted_pvals = multipletests(pvals, method=fdr_method)[1]

    # Full result dataframe
    result_df = pd.DataFrame({
        "kmer": kmers,
        "p_value": pvals,
        "mean_rank_diff": mean_rank_diffs,
        "fdr_corrected_p": fdr_adjusted_pvals
    })

    result_df.sort_values("p_value", inplace=True)
    result_df.to_csv(output_filename, index=False)
    print(f"Full Mann–Whitney U test results saved to {output_filename}")

    # Split into AD- and NC-elevated
    ad_elevated = result_df[result_df["mean_rank_diff"] > 0]
    nc_elevated = result_df[result_df["mean_rank_diff"] < 0]

    ad_filename = output_filename.replace(".csv", "_AD.csv")
    nc_filename = output_filename.replace(".csv", "_NC.csv")

    ad_elevated.to_csv(ad_filename, index=False)
    nc_elevated.to_csv(nc_filename, index=False)

    print(f"AD-elevated k-mers saved to {ad_filename}")
    print(f"NC-elevated k-mers saved to {nc_filename}")

    return result_df

# Step 4: Wrapper function for end-to-end analysis from file intake to U-test
def analyze_groups(input_path, pos_file, neg_file, k=7, wildcard_positions=[], normalize=True):
    """
    Top-level function:
    1. Tiles and processes patient data from positive and negative Excel/CSV files
    2. Builds the patient-by-kmer frequency matrix
    3. Runs Mann–Whitney U test to compare frequencies between groups
    4. Outputs a ranked CSV of p-values and effect sizes
    """
    base_stem = Path(input_path).stem
    wildcard_str = ''.join(str(i) for i in wildcard_positions)
    
    pos_dicts, filtered_kmers_pos = tile_patient_file(Path(pos_file), kmer_length=k, wildcard_positions=wildcard_positions, return_individual_dicts=True)
    neg_dicts, filtered_kmers_neg = tile_patient_file(Path(neg_file), kmer_length=k, wildcard_positions=wildcard_positions, return_individual_dicts=True)

    # Only keep kmers that passed chi-square in BOTH groups (optional, or take union instead)
    filter_kmers = filtered_kmers_pos & filtered_kmers_neg

    all_dicts = {**pos_dicts, **neg_dicts}
    output_path = Path(pos_file).parent / f"{base_stem}_matrix_{k}mers_[{''.join(str(i) for i in wildcard_positions)}].csv"
    matrix = build_kmer_matrix(all_dicts, normalize=normalize, output_path = output_path)

    # Subset matrix to filtered kmers only
    matrix = matrix.loc[matrix.index.intersection(filter_kmers)]

    pos_cols = list(pos_dicts.keys())
    neg_cols = list(neg_dicts.keys())

    
    output_path = Path(pos_file).parent / f"{base_stem}_U_test_{k}mers_[{wildcard_str}].csv"
    run_mannwhitney(matrix, pos_cols, neg_cols, output_filename=str(output_path))

if __name__ ==  "__main__":
    #tile_patient_file(Path("/Users/ciao/Downloads/SI_file1NC.xlsx"), 5, [1, 3, 5])
    input_path = Path("/Users/ciao/Desktop/module2_test_FDR/SI_file1.csv")
    pos_file, neg_file = split_input_by_group(input_path)
    for i in range(4, 8):
        analyze_groups(input_path, pos_file, neg_file, k=i, wildcard_positions=[], normalize=True)