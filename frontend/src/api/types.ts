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
  volcano_filename?: string;
};

export type KmerResponse = {
  result_id: string;
  runs: KmerResultSummary[];
};

export type KmerTaskCreatedResponse = {
  task_id: string;
};

export type KmerTaskStatusResponse = {
  task_id: string;
  status: "queued" | "running" | "succeeded" | "failed";
  progress: number;
  message: string;
  result: KmerResponse | null;
  error: string | null;
};

export type Module3Response = {
  result_id: string;
  output_folder_name: string;
  positive_mapping_filename: string;
  negative_mapping_filename: string;
  positive_clean_filename: string;
  negative_clean_filename: string;
  positive_manhattan_filename: string;
  negative_manhattan_filename: string;
  top_n: number | null;
  wildcards: boolean;
  q_cutoff: number;
};
