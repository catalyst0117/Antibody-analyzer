import os
import re
import gzip
import shutil
from pathlib import Path
from collections import Counter, defaultdict
from Bio.Seq import Seq
import pandas as pd
#just a quick example of not including bg file
# def process_multiple_fastq(input_dir, background_file=None, output_name='sequence_matrix.xlsx'):
#     input_path = Path(input_dir)
#     pattern = re.compile(r'TTCTCACTCT.{0,36}')
#     sample_counters = {}
#     all_sequences = set()

#     # Load background sequences if provided
#     background_set = set()
#     if background_file and Path(background_file).exists():
#         print(f"⏳ Processing background FASTQ: {background_file}")
#         background_set = extract_peptides_from_fastq(background_file)
#         print(f"✅ Loaded {len(background_set)} background peptide sequences.")

#         # ✅ Save background peptides to a file for reference
#         bg_out_path = Path(background_file).with_name("processed_background_peptides.txt")
#         with open(bg_out_path, 'w') as f:
#             for seq in sorted(background_set):
#                 f.write(seq + '\n')
#         print(f"📄 Background peptides saved to: {bg_out_path}")
    
#     # def extract_peptides_from_fastq(fastq_path):
#     #     peptides = set()
#     #     with open(fastq_path, 'r') as f:
#     #         for line in f:
#     #             if line.startswith('@') or line.startswith('+'):
#     #                 continue  # Skip headers and quality markers
#     #             match = pattern.search(line)
#     #             if match:
#     #                 nt_seq = match.group(0)[10:]
#     #                 if len(nt_seq) == 36 and set(nt_seq).issubset({'A', 'T', 'C', 'G'}):
#     #                     aa = str(Seq(nt_seq).translate()).replace('*', 'X')
#     #                     if len(aa) >= 12:
#     #                         peptides.add(aa)
#     #     return peptides
#     # background_set = extract_peptides_from_fastq(background_file)
#     # print(f"✅ Loaded {len(background_set)} background peptide sequences.")
    

#     # Process each FASTQ or FASTQ.GZ file
#     for file in sorted(input_path.glob("*.fastq*")):
#         sample_name = file.stem.split(".")[0]
#         print(f"\nProcessing sample: {sample_name}")

#         # Unzip if necessary
#         if file.suffix == '.gz':
#             unzipped_path = file.with_suffix('')
#             with gzip.open(file, 'rb') as f_in, open(unzipped_path, 'wb') as f_out:
#                 shutil.copyfileobj(f_in, f_out)
#             file = unzipped_path

#         # Step 1: Extract variable regions
#         with open(file, 'r') as f:
#             sequences = [
#                 match.group(0)[10:]  # trim TTCTCACTCT
#                 for line in f
#                 if (match := pattern.search(line))
#             ]

#         # Step 2: Translate valid DNA to protein
#         translated = []
#         for seq in sequences:
#             if len(seq) == 36 and set(seq).issubset({'A', 'T', 'C', 'G'}):
#                 aa = str(Seq(seq).translate()).replace('*', 'X')
#                 if len(aa) >= 12 and aa not in background_set:
#                     translated.append(aa)

#         # Step 3: Count and store
#         counter = Counter(translated)
#         sample_counters[sample_name] = counter
#         all_sequences.update(counter.keys())
#         print(f"  {len(counter)} unique sequences kept")

#     # Step 4: Build multi-column sample layout
#     sorted_sequences = sorted(all_sequences)
#     # Step 4: Build multi-column sample layout with sorting by frequency
#     df_parts = []

#     for sample, counter in sample_counters.items():
#         # Sort the sequences by count descending for this sample
#         sorted_pairs = sorted(counter.items(), key=lambda x: -x[1])

#         # Split into two columns
#         sequences = [seq for seq, _ in sorted_pairs]
#         counts = [count for _, count in sorted_pairs]

#         sample_df = pd.DataFrame({
#             f"{sample}_Sequence": sequences,
#             f"{sample}_Count": counts
#         })

#         df_parts.append(sample_df)

#     # Pad shorter dataframes with NaNs so concat works
#     max_len = max(len(df) for df in df_parts)
#     df_parts = [df.reindex(range(max_len)) for df in df_parts]

#     # Combine all sample columns side-by-side
#     df = pd.concat(df_parts, axis=1)

#     # Always save both XLSX and CSV
#     xlsx_path = input_path / (output_name if output_name.endswith('.xlsx') else output_name + ".xlsx")

#     df.to_excel(xlsx_path, index=False)

#     print(f"\n✅ Output saved as:\n  Excel: {xlsx_path}")

# if __name__ == "__main__":
#     path = "/Users/ciao/Desktop/fastq_data"
#     process_multiple_fastq(input_dir=path, background_file=None, output_name="result_without_bead")

import os
import re
import gzip
import shutil
from pathlib import Path
from collections import Counter
from Bio.Seq import Seq
import pandas as pd

def process_multiple_fastq(input_dir, background_file=None, output_name='sequence_matrix.xlsx'):
    input_path = Path(input_dir)
    pattern = re.compile(r'TTCTCACTCT.{0,36}')
    sample_counters = {}
    all_sequences = set()

    # ✅ Helper: Convert background FASTQ to peptide set
    def extract_peptides_from_fastq(fastq_path):
        peptides = set()
        with open(fastq_path, 'r') as f:
            for line in f:
                if line.startswith('@') or line.startswith('+'):
                    continue  # Skip headers and quality markers
                match = pattern.search(line)
                if match:
                    nt_seq = match.group(0)[10:]
                    if len(nt_seq) == 36 and set(nt_seq).issubset({'A', 'T', 'C', 'G'}):
                        aa = str(Seq(nt_seq).translate()).replace('*', 'X')
                        if len(aa) >= 12:
                            peptides.add(aa)
        return peptides

    # ✅ Load and process background sequences if provided
    background_set = set()
    if background_file and Path(background_file).exists():
        print(f"⏳ Processing background FASTQ: {background_file}")
        background_set = extract_peptides_from_fastq(background_file)
        print(f"✅ Loaded {len(background_set)} background peptide sequences.")
        bg_out_path = Path(background_file).with_name("processed_background_peptides.txt")
        with open(bg_out_path, 'w') as f:
            for seq in sorted(background_set):
                f.write(seq + '\n')
        print(f"📄 Background peptides saved to: {bg_out_path}")

    # 🔁 Process each FASTQ or FASTQ.GZ file
    for file in sorted(input_path.glob("*.fastq*")):
        sample_name = file.stem.split(".")[0]
        print(f"\n📦 Processing sample: {sample_name}")

        # Unzip if necessary
        if file.suffix == '.gz':
            unzipped_path = file.with_suffix('')
            with gzip.open(file, 'rb') as f_in, open(unzipped_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            file = unzipped_path

        # Step 1: Extract variable regions
        with open(file, 'r') as f:
            sequences = [
                match.group(0)[10:]
                for line in f
                if (match := pattern.search(line))
            ]

        # Step 2: Translate valid DNA to aminoacids
        translated = []
        filtered_out = []

        for seq in sequences:
            #complete codons 36 bases in length
            if len(seq) == 36 and set(seq).issubset({'A', 'T', 'C', 'G'}):
                aa = str(Seq(seq).translate()).replace('*', 'X')
                if len(aa) >= 12:
                    if aa in background_set:
                        filtered_out.append(aa)
                    else:
                        translated.append(aa)

        # Optional: log some filtered sequences
        if filtered_out:
            print(f"🛑 {len(filtered_out)} sequences filtered from {sample_name} due to background overlap.")

            # Save to file
            filtered_path = input_path / f"{sample_name}_filtered_out.txt"
            with open(filtered_path, 'w') as f:
                for seq in sorted(set(filtered_out)):
                    f.write(seq + '\n')
            print(f"   📄 Saved filtered sequences to: {filtered_path}")

            # Step 3: Count and store
            counter = Counter(translated)
            sample_counters[sample_name] = counter
            all_sequences.update(counter.keys())
            print(f"✅ {len(counter)} unique sequences kept after background filtering")

        # Step 4: Build per-sample columns sorted by frequency
        df_parts = []
        for sample, counter in sample_counters.items():
            sorted_pairs = sorted(counter.items(), key=lambda x: -x[1])
            sequences = [seq for seq, _ in sorted_pairs]
            counts = [count for _, count in sorted_pairs]
            sample_df = pd.DataFrame({
                f"{sample}_Sequence": sequences,
                f"{sample}_Count": counts
            })
            df_parts.append(sample_df)

        # Step 5: Pad shorter DataFrames so they align
        max_len = max(len(df) for df in df_parts)
        df_parts = [df.reindex(range(max_len)) for df in df_parts]

        # Step 6: Combine side by side and save
        df = pd.concat(df_parts, axis=1)
        xlsx_path = input_path / (output_name if output_name.endswith('.xlsx') else output_name + ".xlsx")
        
        df.to_excel(xlsx_path, index=False)
        

        print(f"\n📁 Output saved as:\n  Excel: {xlsx_path}")

# ✅ Example usage
if __name__ == "__main__":
    path = "/Users/ciao/Desktop/fastq_data"
    background = "/Users/ciao/Desktop/beads_no_phage_3rd_UA.fastq"
    process_multiple_fastq(input_dir=path, background_file=background, output_name="result")