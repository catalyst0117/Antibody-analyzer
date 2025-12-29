export type FastqSampleSummary = {
  sample_name: string;
  total_sequences: number;
  unique_sequences: number;
  filtered_sequences: number;
  output_columns: string[];
};

export type FastqResponse = {
  result_id: string;
  summary: FastqSampleSummary[];
  excel_filename: string;
  background_filename: string | null;
  filtered_files: string[];
};

export type KmerResultSummary = {
  k: number;
  total_kmers: number;
  ad_elevated: number;
  nc_elevated: number;
  result_filename: string;
  ad_filename: string;
  nc_filename: string;
  matrix_filename: string;
};

export type KmerResponse = {
  result_id: string;
  runs: KmerResultSummary[];
};
