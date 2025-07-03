#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

// Structure to hold nucleotide counts
typedef struct {
    int A, T, G, C;
    int INS, DEL;
    int total_reads;
    int coverage;
    char ref_base;
} nuc_counts_t;

// Structure for position information
typedef struct {
    int original_index;
    const char *chr;
    int pos;
} pos_info_t;

// Structure for cancer hotspot information
typedef struct {
    char chr[32];
    int pos;
    char ref[16];
    char alt[16];
    char gene[64];
    char variant_class[32];
    char protein_change[128];
    char meta[256];
} hotspot_info_t;

// Structure for cancer hotspot processing parameters
typedef struct {
    double vaf_threshold;
    int min_depth;
    int min_alt_reads;
    const char *output_prefix;
} hotspot_params_t;

// Function to initialize nucleotide counts
void init_counts(nuc_counts_t *counts) {
    counts->A = 0;
    counts->T = 0;
    counts->G = 0;
    counts->C = 0;
    counts->INS = 0;
    counts->DEL = 0;
    counts->total_reads = 0;
    counts->coverage = 0;
    counts->ref_base = 'N';
}

// Function to parse CIGAR and count nucleotides at a specific position
void parse_cigar_and_count(bam1_t *read, int target_pos, nuc_counts_t *counts) {
    uint32_t *cigar = bam_get_cigar(read);
    uint32_t n_cigar = read->core.n_cigar;
    uint8_t *seq = bam_get_seq(read);
    
    int ref_pos = read->core.pos;
    int read_pos = 0;
    
    for (uint32_t i = 0; i < n_cigar; i++) {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        
        switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                if (ref_pos <= target_pos && target_pos < ref_pos + len) {
                    int offset = target_pos - ref_pos;
                    int seq_pos = read_pos + offset;
                    uint8_t base = bam_seqi(seq, seq_pos);
                    
                    switch (base) {
                        case 1: counts->A++; break;  // A
                        case 2: counts->C++; break;  // C
                        case 4: counts->G++; break;  // G
                        case 8: counts->T++; break;  // T
                    }
                    counts->coverage++;
                }
                ref_pos += len;
                read_pos += len;
                break;
                
            case BAM_CINS:
                if (ref_pos == target_pos) {
                    counts->INS++;
                }
                read_pos += len;
                break;
                
            case BAM_CDEL:
                if (ref_pos <= target_pos && target_pos < ref_pos + len) {
                    counts->DEL++;
                }
                ref_pos += len;
                break;
                
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
                read_pos += len;
                break;
                
            case BAM_CREF_SKIP:
                ref_pos += len;
                break;
                
            case BAM_CPAD:
                break;
        }
    }
}

// Function to count nucleotides at a specific position
void count_at_position(samFile *fp, hts_idx_t *idx, bam_hdr_t *hdr, 
                      const char *chr, int pos, int mapq, int flag, 
                      faidx_t *fai, nuc_counts_t *counts) {
    
    init_counts(counts);
    
    // Get reference base if fasta index is provided
    if (fai) {
        int len;
        char *ref_seq = faidx_fetch_seq(fai, chr, pos, pos, &len);
        if (ref_seq && len > 0) {
            counts->ref_base = ref_seq[0];
            free(ref_seq);
        }
    }
    
    char bam_region[256];
    snprintf(bam_region, sizeof(bam_region), "%s:%d-%d", chr, pos + 1, pos + 1);
    
    hts_itr_t *iter = sam_itr_querys(idx, hdr, bam_region);
    if (!iter) return;
    
    bam1_t *read = bam_init1();
    if (!read) {
        hts_itr_destroy(iter);
        return;
    }
    
    while (sam_itr_next(fp, iter, read) >= 0) {
        counts->total_reads++;
        
        // Filter by mapping quality
        if (read->core.qual < mapq) continue;
        
        // Filter by flags
        if (read->core.flag & flag) continue;
        
        // Parse CIGAR and count nucleotides
        parse_cigar_and_count(read, pos, counts);
    }
    
    bam_destroy1(read);
    hts_itr_destroy(iter);
}

// Chromosome-based batch processing for improved performance
void process_chromosome_batch(samFile *fp, hts_idx_t *idx, bam_hdr_t *hdr,
                             const char *chr, pos_info_t *chr_positions, int n_chr_pos,
                             int min_mapq, int filter_flag, faidx_t *fai,
                             SEXP result_chr, SEXP result_pos, SEXP result_A, SEXP result_T, 
                             SEXP result_G, SEXP result_C, SEXP result_INS, SEXP result_DEL,
                             SEXP result_total, SEXP result_cov, SEXP result_ref) {
    
    // Sort positions within chromosome for sequential access
    for (int i = 0; i < n_chr_pos - 1; i++) {
        for (int j = i + 1; j < n_chr_pos; j++) {
            if (chr_positions[i].pos > chr_positions[j].pos) {
                pos_info_t temp = chr_positions[i];
                chr_positions[i] = chr_positions[j];
                chr_positions[j] = temp;
            }
        }
    }
    
    // Process positions sequentially for this chromosome
    for (int i = 0; i < n_chr_pos; i++) {
        int pos = chr_positions[i].pos;
        int orig_idx = chr_positions[i].original_index;
        
        nuc_counts_t counts;
        count_at_position(fp, idx, hdr, chr, pos, min_mapq, filter_flag, fai, &counts);
        
        // Store results in original order
        SET_STRING_ELT(result_chr, orig_idx, mkChar(chr));
        INTEGER(result_pos)[orig_idx] = pos + 1; // Convert back to 1-based
        INTEGER(result_A)[orig_idx] = counts.A;
        INTEGER(result_T)[orig_idx] = counts.T;
        INTEGER(result_G)[orig_idx] = counts.G;
        INTEGER(result_C)[orig_idx] = counts.C;
        INTEGER(result_INS)[orig_idx] = counts.INS;
        INTEGER(result_DEL)[orig_idx] = counts.DEL;
        INTEGER(result_total)[orig_idx] = counts.total_reads;
        INTEGER(result_cov)[orig_idx] = counts.coverage;
        
        // Store reference base
        char ref_str[2] = {counts.ref_base, '\0'};
        SET_STRING_ELT(result_ref, orig_idx, mkChar(ref_str));
    }
}

// Calculate VAF for cancer hotspot analysis
double calculate_hotspot_vaf(const hotspot_info_t *hotspot, const nuc_counts_t *counts, int *passes_filters, const hotspot_params_t *params) {
    *passes_filters = 0;
    
    int total_reads = counts->A + counts->T + counts->G + counts->C + counts->INS + counts->DEL;
    if (total_reads < params->min_depth) {
        return -1.0;
    }
    
    int alt_count = 0;
    double vaf = 0.0;
    
    // Determine variant type and calculate VAF
    if (strcmp(hotspot->variant_class, "Nonsense_Mutation") == 0 || 
        strcmp(hotspot->variant_class, "Missense_Mutation") == 0 ||
        strcmp(hotspot->variant_class, "Silent") == 0) {
        // SNV - get alt allele count
        if (hotspot->alt[0] == 'A') alt_count = counts->A;
        else if (hotspot->alt[0] == 'T') alt_count = counts->T;
        else if (hotspot->alt[0] == 'G') alt_count = counts->G;
        else if (hotspot->alt[0] == 'C') alt_count = counts->C;
        
        if (alt_count >= params->min_alt_reads) {
            vaf = (double)alt_count / total_reads;
            if (vaf >= params->vaf_threshold) {
                *passes_filters = 1;
            }
        }
    } else if (strstr(hotspot->variant_class, "Ins") != NULL) {
        // Insertion
        alt_count = counts->INS;
        if (alt_count >= params->min_alt_reads) {
            vaf = (double)alt_count / total_reads;
            if (vaf >= params->vaf_threshold) {
                *passes_filters = 1;
            }
        }
    } else if (strstr(hotspot->variant_class, "Del") != NULL) {
        // Deletion
        alt_count = counts->DEL;
        if (alt_count >= params->min_alt_reads) {
            vaf = (double)alt_count / total_reads;
            if (vaf >= params->vaf_threshold) {
                *passes_filters = 1;
            }
        }
    }
    
    return vaf;
}

// HTML generation functions for cancer hotspots
void write_html_header(FILE *fp, const char *bam_name) {
    fprintf(fp, "<!DOCTYPE html>\n");
    fprintf(fp, "<html lang=\"en\">\n");
    fprintf(fp, "<head>\n");
    fprintf(fp, "<meta charset=\"UTF-8\">\n");
    fprintf(fp, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
    fprintf(fp, "<title>%s | Cancer Hotspots Analysis</title>\n", bam_name);
    
    // CSS styling
    fprintf(fp, "<style>\n");
    fprintf(fp, "body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background-color: #f8f9fa; }\n");
    fprintf(fp, ".container { max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }\n");
    fprintf(fp, "h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }\n");
    fprintf(fp, "h2 { color: #34495e; margin-top: 30px; }\n");
    fprintf(fp, "table { width: 100%%; border-collapse: collapse; margin: 20px 0; font-size: 14px; }\n");
    fprintf(fp, "th { background-color: #3498db; color: white; padding: 12px 8px; text-align: left; position: sticky; top: 0; }\n");
    fprintf(fp, "td { padding: 8px; border-bottom: 1px solid #ecf0f1; }\n");
    fprintf(fp, "tr:nth-child(even) { background-color: #f8f9fa; }\n");
    fprintf(fp, "tr:hover { background-color: #e8f4fd; }\n");
    fprintf(fp, ".high-vaf { background-color: #e74c3c; color: white; font-weight: bold; }\n");
    fprintf(fp, ".medium-vaf { background-color: #f39c12; color: white; font-weight: bold; }\n");
    fprintf(fp, ".low-vaf { background-color: #27ae60; color: white; }\n");
    fprintf(fp, ".gene-link { color: #3498db; text-decoration: none; font-weight: bold; }\n");
    fprintf(fp, ".gene-link:hover { text-decoration: underline; }\n");
    fprintf(fp, ".summary-table { width: auto; margin: 20px 0; }\n");
    fprintf(fp, ".summary-table th, .summary-table td { padding: 8px 15px; }\n");
    fprintf(fp, "</style>\n");
    fprintf(fp, "</head>\n");
    fprintf(fp, "<body>\n");
    fprintf(fp, "<div class=\"container\">\n");
    fprintf(fp, "<h1>Cancer Hotspots Analysis Report</h1>\n");
    fprintf(fp, "<p><strong>Sample:</strong> %s</p>\n", bam_name);
    
    // Add current timestamp
    time_t now = time(NULL);
    struct tm *tm_info = localtime(&now);
    char timestamp[64];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", tm_info);
    fprintf(fp, "<p><strong>Generated:</strong> %s</p>\n", timestamp);
}

void write_html_summary(FILE *fp, int total_variants, int passing_variants, double avg_depth, const hotspot_params_t *params) {
    fprintf(fp, "<h2>Analysis Summary</h2>\n");
    fprintf(fp, "<table class=\"summary-table\">\n");
    fprintf(fp, "<tr><th>Parameter</th><th>Value</th></tr>\n");
    fprintf(fp, "<tr><td>Total variants analyzed</td><td>%d</td></tr>\n", total_variants);
    fprintf(fp, "<tr><td>Variants passing filters</td><td>%d (%.1f%%)</td></tr>\n", 
           passing_variants, total_variants > 0 ? (double)passing_variants / total_variants * 100 : 0.0);
    fprintf(fp, "<tr><td>Average depth of coverage</td><td>%.1f</td></tr>\n", avg_depth);
    fprintf(fp, "<tr><td>VAF threshold</td><td>%.3f</td></tr>\n", params->vaf_threshold);
    fprintf(fp, "<tr><td>Minimum depth</td><td>%d</td></tr>\n", params->min_depth);
    fprintf(fp, "<tr><td>Minimum alt reads</td><td>%d</td></tr>\n", params->min_alt_reads);
    fprintf(fp, "</table>\n");
}

void write_html_table_header(FILE *fp) {
    fprintf(fp, "<h2>Detected Variants</h2>\n");
    fprintf(fp, "<table>\n");
    fprintf(fp, "<thead>\n");
    fprintf(fp, "<tr>\n");
    fprintf(fp, "<th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th>\n");
    fprintf(fp, "<th>Gene</th><th>Variant Type</th><th>Protein Change</th>\n");
    fprintf(fp, "<th>VAF</th><th>A</th><th>T</th><th>G</th><th>C</th><th>Ins</th><th>Del</th>\n");
    fprintf(fp, "</tr>\n");
    fprintf(fp, "</thead>\n");
    fprintf(fp, "<tbody>\n");
}

void write_html_variant_row(FILE *fp, const hotspot_info_t *hotspot, double vaf, const nuc_counts_t *counts) {
    // Determine VAF class for coloring
    const char *vaf_class = "";
    if (vaf >= 0.3) vaf_class = "high-vaf";
    else if (vaf >= 0.1) vaf_class = "medium-vaf";
    else if (vaf >= 0.05) vaf_class = "low-vaf";
    
    fprintf(fp, "<tr>\n");
    fprintf(fp, "<td>%s</td>\n", hotspot->chr);
    fprintf(fp, "<td>%d</td>\n", hotspot->pos);
    fprintf(fp, "<td>%s</td>\n", hotspot->ref);
    fprintf(fp, "<td>%s</td>\n", hotspot->alt);
    fprintf(fp, "<td><a href=\"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=%s#variants\" class=\"gene-link\" target=\"_blank\">%s</a></td>\n", 
           hotspot->gene, hotspot->gene);
    fprintf(fp, "<td>%s</td>\n", hotspot->variant_class);
    fprintf(fp, "<td>%s</td>\n", hotspot->protein_change);
    fprintf(fp, "<td class=\"%s\">%.3f</td>\n", vaf_class, vaf);
    fprintf(fp, "<td>%d</td>\n", counts->A);
    fprintf(fp, "<td>%d</td>\n", counts->T);
    fprintf(fp, "<td>%d</td>\n", counts->G);
    fprintf(fp, "<td>%d</td>\n", counts->C);
    fprintf(fp, "<td>%d</td>\n", counts->INS);
    fprintf(fp, "<td>%d</td>\n", counts->DEL);
    fprintf(fp, "</tr>\n");
}

void write_html_footer(FILE *fp) {
    fprintf(fp, "</tbody>\n");
    fprintf(fp, "</table>\n");
    fprintf(fp, "<p style=\"margin-top: 30px; font-size: 12px; color: #7f8c8d;\">\n");
    fprintf(fp, "Generated by <a href=\"https://github.com/PoisonAlien/maftools\">maftools::cancerhotspots()</a>\n");
    fprintf(fp, "</p>\n");
    fprintf(fp, "</div>\n");
    fprintf(fp, "</body>\n");
    fprintf(fp, "</html>\n");
}

// Main C function called from R
SEXP bamrc_c(SEXP bam_path, SEXP chr_vec, SEXP pos_vec, SEXP mapq, SEXP flag, SEXP fasta_path, SEXP include_idxstats, SEXP verbose, SEXP hotspot_data) {
    
    // Get input parameters
    const char *bam_file = CHAR(STRING_ELT(bam_path, 0));
    int n_pos = LENGTH(chr_vec);
    int min_mapq = INTEGER(mapq)[0];
    int filter_flag = INTEGER(flag)[0];
    int calc_idxstats = 0;
    int is_verbose = 0;
    int is_hotspot_mode = 0;
    
    // Check if idxstats should be calculated
    if (include_idxstats != R_NilValue && LENGTH(include_idxstats) > 0) {
        calc_idxstats = LOGICAL(include_idxstats)[0];
    }
    
    // Check verbose parameter
    if (verbose != R_NilValue && LENGTH(verbose) > 0) {
        is_verbose = LOGICAL(verbose)[0];
    }
    
    // Check if cancer hotspot mode is enabled
    hotspot_info_t *hotspots = NULL;
    hotspot_params_t hotspot_params = {0};
    if (hotspot_data != R_NilValue && LENGTH(hotspot_data) > 0) {
        is_hotspot_mode = 1;
        
        // Parse hotspot data (list with: data.frame, vaf_threshold, min_depth, min_alt_reads, output_prefix)
        SEXP hotspot_df = VECTOR_ELT(hotspot_data, 0);
        hotspot_params.vaf_threshold = REAL(VECTOR_ELT(hotspot_data, 1))[0];
        hotspot_params.min_depth = INTEGER(VECTOR_ELT(hotspot_data, 2))[0];
        hotspot_params.min_alt_reads = INTEGER(VECTOR_ELT(hotspot_data, 3))[0];
        hotspot_params.output_prefix = CHAR(STRING_ELT(VECTOR_ELT(hotspot_data, 4), 0));
        
        // Allocate hotspot array
        hotspots = (hotspot_info_t *)malloc(n_pos * sizeof(hotspot_info_t));
        if (!hotspots) {
            error("Memory allocation failed for hotspot data");
        }
        
        // Parse hotspot dataframe (columns: chr, pos, ref, alt, gene, variant_class, protein_change, meta)
        SEXP ref_col = VECTOR_ELT(hotspot_df, 2);
        SEXP alt_col = VECTOR_ELT(hotspot_df, 3);
        SEXP gene_col = VECTOR_ELT(hotspot_df, 4);
        SEXP vc_col = VECTOR_ELT(hotspot_df, 5);
        SEXP pc_col = VECTOR_ELT(hotspot_df, 6);
        SEXP meta_col = VECTOR_ELT(hotspot_df, 7);
        
        for (int i = 0; i < n_pos; i++) {
            strncpy(hotspots[i].chr, CHAR(STRING_ELT(chr_vec, i)), sizeof(hotspots[i].chr) - 1);
            hotspots[i].chr[sizeof(hotspots[i].chr) - 1] = '\0';
            hotspots[i].pos = INTEGER(pos_vec)[i];
            strncpy(hotspots[i].ref, CHAR(STRING_ELT(ref_col, i)), sizeof(hotspots[i].ref) - 1);
            hotspots[i].ref[sizeof(hotspots[i].ref) - 1] = '\0';
            strncpy(hotspots[i].alt, CHAR(STRING_ELT(alt_col, i)), sizeof(hotspots[i].alt) - 1);
            hotspots[i].alt[sizeof(hotspots[i].alt) - 1] = '\0';
            strncpy(hotspots[i].gene, CHAR(STRING_ELT(gene_col, i)), sizeof(hotspots[i].gene) - 1);
            hotspots[i].gene[sizeof(hotspots[i].gene) - 1] = '\0';
            strncpy(hotspots[i].variant_class, CHAR(STRING_ELT(vc_col, i)), sizeof(hotspots[i].variant_class) - 1);
            hotspots[i].variant_class[sizeof(hotspots[i].variant_class) - 1] = '\0';
            strncpy(hotspots[i].protein_change, CHAR(STRING_ELT(pc_col, i)), sizeof(hotspots[i].protein_change) - 1);
            hotspots[i].protein_change[sizeof(hotspots[i].protein_change) - 1] = '\0';
            strncpy(hotspots[i].meta, CHAR(STRING_ELT(meta_col, i)), sizeof(hotspots[i].meta) - 1);
            hotspots[i].meta[sizeof(hotspots[i].meta) - 1] = '\0';
        }
    }
    
    // Open fasta index if provided
    faidx_t *fai = NULL;
    if (fasta_path != R_NilValue) {
        const char *fasta_file = CHAR(STRING_ELT(fasta_path, 0));
        fai = fai_load(fasta_file);
        if (!fai) {
            error("Cannot load fasta index for file: %s", fasta_file);
        }
    }
    
    // Open BAM file
    samFile *fp = sam_open(bam_file, "r");
    if (!fp) {
        error("Cannot open BAM file: %s", bam_file);
    }
    
    // Load header
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        error("Cannot read BAM header");
    }
    
    // Load index
    hts_idx_t *idx = sam_index_load(fp, bam_file);
    if (!idx) {
        bam_hdr_destroy(hdr);
        sam_close(fp);
        if (fai) fai_destroy(fai);
        error("Cannot load BAM index");
    }
    
    // Calculate total mapped reads if requested
    uint64_t total_mapped = 0;
    if (calc_idxstats) {
        int32_t n_contigs = hdr->n_targets;
        for (int i = 0; i < n_contigs; i++) {
            uint64_t n_mapped = 0, n_unmapped = 0;
            if (hts_idx_get_stat(idx, i, &n_mapped, &n_unmapped) == 0) {
                total_mapped += n_mapped;
            }
        }
    }
    
    // Allocate result vectors
    SEXP result_chr = PROTECT(allocVector(STRSXP, n_pos));
    SEXP result_pos = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_A = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_T = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_G = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_C = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_INS = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_DEL = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_total = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_cov = PROTECT(allocVector(INTSXP, n_pos));
    SEXP result_ref = PROTECT(allocVector(STRSXP, n_pos));
    
    // Create position info structure for efficient processing
    pos_info_t *pos_info = (pos_info_t *)malloc(n_pos * sizeof(pos_info_t));
    if (!pos_info) {
        error("Memory allocation failed for position info");
    }
    
    // Populate position info
    for (int i = 0; i < n_pos; i++) {
        pos_info[i].original_index = i;
        pos_info[i].chr = CHAR(STRING_ELT(chr_vec, i));
        pos_info[i].pos = INTEGER(pos_vec)[i] - 1; // Convert to 0-based
    }
    
    if (is_verbose) {
        Rprintf("Processing %d positions:\n", n_pos);
    }
    
    // Group positions by chromosome for efficient processing
    // Find unique chromosomes
    char **unique_chrs = (char **)malloc(n_pos * sizeof(char *));
    int n_unique_chrs = 0;
    
    for (int i = 0; i < n_pos; i++) {
        const char *chr = pos_info[i].chr;
        int found = 0;
        for (int j = 0; j < n_unique_chrs; j++) {
            if (strcmp(unique_chrs[j], chr) == 0) {
                found = 1;
                break;
            }
        }
        if (!found) {
            unique_chrs[n_unique_chrs] = (char *)malloc((strlen(chr) + 1) * sizeof(char));
            strcpy(unique_chrs[n_unique_chrs], chr);
            n_unique_chrs++;
        }
    }
    
    int processed_positions = 0;
    
    // Process each chromosome
    for (int chr_idx = 0; chr_idx < n_unique_chrs; chr_idx++) {
        const char *current_chr = unique_chrs[chr_idx];
        
        // Count positions for this chromosome
        int chr_pos_count = 0;
        for (int i = 0; i < n_pos; i++) {
            if (strcmp(pos_info[i].chr, current_chr) == 0) {
                chr_pos_count++;
            }
        }
        
        if (chr_pos_count == 0) continue;
        
        // Collect positions for this chromosome
        pos_info_t *chr_positions = (pos_info_t *)malloc(chr_pos_count * sizeof(pos_info_t));
        int chr_pos_idx = 0;
        for (int i = 0; i < n_pos; i++) {
            if (strcmp(pos_info[i].chr, current_chr) == 0) {
                chr_positions[chr_pos_idx++] = pos_info[i];
            }
        }
        
        // Process this chromosome batch
        process_chromosome_batch(fp, idx, hdr, current_chr, chr_positions, chr_pos_count,
                               min_mapq, filter_flag, fai,
                               result_chr, result_pos, result_A, result_T, result_G, result_C,
                               result_INS, result_DEL, result_total, result_cov, result_ref);
        
        processed_positions += chr_pos_count;
        
        if (is_verbose) {
            int bar_width = 50;
            double progress = (double)processed_positions / n_pos;
            int pos = bar_width * progress;
            
            Rprintf("\r[");
            for (int i = 0; i < bar_width; ++i) {
                if (i < pos) Rprintf("=");
                else if (i == pos) Rprintf(">");
                else Rprintf(" ");
            }
            Rprintf("] %.1f%%", progress * 100);
        }
        
        free(chr_positions);
        R_CheckUserInterrupt(); // Allow user to interrupt
    }
    
    if(is_verbose){
        Rprintf("\n");
    }
    
    // Cleanup
    for (int i = 0; i < n_unique_chrs; i++) {
        free(unique_chrs[i]);
    }
    free(unique_chrs);
    free(pos_info);
    
    // Create data.frame
    int n_cols = fai ? 11 : 10;
    SEXP result_df = PROTECT(allocVector(VECSXP, n_cols));
    SET_VECTOR_ELT(result_df, 0, result_chr);
    SET_VECTOR_ELT(result_df, 1, result_pos);
    SET_VECTOR_ELT(result_df, 2, result_A);
    SET_VECTOR_ELT(result_df, 3, result_T);
    SET_VECTOR_ELT(result_df, 4, result_G);
    SET_VECTOR_ELT(result_df, 5, result_C);
    SET_VECTOR_ELT(result_df, 6, result_INS);
    SET_VECTOR_ELT(result_df, 7, result_DEL);
    SET_VECTOR_ELT(result_df, 8, result_total);
    SET_VECTOR_ELT(result_df, 9, result_cov);
    if (fai) {
        SET_VECTOR_ELT(result_df, 10, result_ref);
    }
    
    // Set column names
    SEXP col_names = PROTECT(allocVector(STRSXP, n_cols));
    SET_STRING_ELT(col_names, 0, mkChar("chr"));
    SET_STRING_ELT(col_names, 1, mkChar("pos"));
    SET_STRING_ELT(col_names, 2, mkChar("A"));
    SET_STRING_ELT(col_names, 3, mkChar("T"));
    SET_STRING_ELT(col_names, 4, mkChar("G"));
    SET_STRING_ELT(col_names, 5, mkChar("C"));
    SET_STRING_ELT(col_names, 6, mkChar("INS"));
    SET_STRING_ELT(col_names, 7, mkChar("DEL"));
    SET_STRING_ELT(col_names, 8, mkChar("total_reads"));
    SET_STRING_ELT(col_names, 9, mkChar("coverage"));
    if (fai) {
        SET_STRING_ELT(col_names, 10, mkChar("ref_base"));
    }
    setAttrib(result_df, R_NamesSymbol, col_names);
    
    // Set class to data.frame
    SEXP row_names = PROTECT(allocVector(INTSXP, 2));
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -n_pos;
    setAttrib(result_df, R_RowNamesSymbol, row_names);
    
    SEXP class_name = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(class_name, 0, mkChar("data.frame"));
    setAttrib(result_df, R_ClassSymbol, class_name);
    
    // Add idxstats as attribute if calculated
    if (calc_idxstats) {
        SEXP idxstats_attr = PROTECT(allocVector(REALSXP, 1));
        REAL(idxstats_attr)[0] = (double)total_mapped;
        setAttrib(result_df, install("total_mapped_reads"), idxstats_attr);
        UNPROTECT(1);
    }
    
    // Process cancer hotspots if enabled
    if (is_hotspot_mode && hotspots) {
        if (is_verbose) {
            Rprintf("Generating cancer hotspot report...\n");
        }
        
        // Open output files
        char html_file[512], tsv_file[512];
        snprintf(html_file, sizeof(html_file), "%s.html", hotspot_params.output_prefix);
        snprintf(tsv_file, sizeof(tsv_file), "%s.tsv", hotspot_params.output_prefix);
        
        FILE *html_fp = fopen(html_file, "w");
        FILE *tsv_fp = fopen(tsv_file, "w");
        
        if (html_fp && tsv_fp) {
            // Get BAM basename for report
            const char *bam_basename = strrchr(bam_file, '/');
            bam_basename = bam_basename ? bam_basename + 1 : bam_file;
            
            // Write HTML header
            write_html_header(html_fp, bam_basename);
            
            // Write TSV header
            fprintf(tsv_fp, "loci\tfa_ref\tNT_change\tHugo_Symbol\tVariant_Classification\tAA_change\tMeta\tVAF\tA\tT\tG\tC\tIns\tDel\n");
            
            // Process each hotspot
            int passing_variants = 0;
            double total_depth = 0.0;
            
            for (int i = 0; i < n_pos; i++) {
                // Get nucleotide counts for this position
                nuc_counts_t counts;
                int pos_idx = -1;
                
                // Find this position in our results
                for (int j = 0; j < n_pos; j++) {
                    if (strcmp(CHAR(STRING_ELT(result_chr, j)), hotspots[i].chr) == 0 &&
                        INTEGER(result_pos)[j] == hotspots[i].pos) {
                        pos_idx = j;
                        break;
                    }
                }
                
                if (pos_idx >= 0) {
                    counts.A = INTEGER(result_A)[pos_idx];
                    counts.T = INTEGER(result_T)[pos_idx];
                    counts.G = INTEGER(result_G)[pos_idx];
                    counts.C = INTEGER(result_C)[pos_idx];
                    counts.INS = INTEGER(result_INS)[pos_idx];
                    counts.DEL = INTEGER(result_DEL)[pos_idx];
                    counts.ref_base = CHAR(STRING_ELT(result_ref, pos_idx))[0];
                    
                    // Calculate VAF
                    int passes_filters = 0;
                    double vaf = calculate_hotspot_vaf(&hotspots[i], &counts, &passes_filters, &hotspot_params);
                    
                    total_depth += (counts.A + counts.T + counts.G + counts.C + counts.INS + counts.DEL);
                    
                    // Write to TSV (all variants)
                    fprintf(tsv_fp, "%s:%d\t%c\t%s>%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\n",
                           hotspots[i].chr, hotspots[i].pos, counts.ref_base,
                           hotspots[i].ref, hotspots[i].alt, hotspots[i].gene,
                           hotspots[i].variant_class, hotspots[i].protein_change, hotspots[i].meta,
                           vaf, counts.A, counts.T, counts.G, counts.C, counts.INS, counts.DEL);
                    
                    // Write to HTML (only passing variants)
                    if (passes_filters) {
                        if (passing_variants == 0) {
                            write_html_table_header(html_fp);
                        }
                        write_html_variant_row(html_fp, &hotspots[i], vaf, &counts);
                        passing_variants++;
                    }
                }
            }
            
            // Write HTML summary and footer
            double avg_depth = n_pos > 0 ? total_depth / n_pos : 0.0;
            write_html_summary(html_fp, n_pos, passing_variants, avg_depth, &hotspot_params);
            write_html_footer(html_fp);
            
            fclose(html_fp);
            fclose(tsv_fp);
            
            if (is_verbose) {
                Rprintf("Cancer hotspot analysis complete:\n");
                Rprintf("  Total variants: %d\n", n_pos);
                Rprintf("  Variants passing filters: %d\n", passing_variants);
                Rprintf("  HTML report: %s\n", html_file);
                Rprintf("  TSV file: %s\n", tsv_file);
            }
        } else {
            if (html_fp) fclose(html_fp);
            if (tsv_fp) fclose(tsv_fp);
            warning("Could not create output files for cancer hotspot analysis");
        }
        
        free(hotspots);
    }
    
    // Cleanup
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    if (fai) fai_destroy(fai);
    
    UNPROTECT(15);
    return result_df;
}