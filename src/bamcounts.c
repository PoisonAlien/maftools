#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <math.h>
#include <time.h>

#include <Rdefines.h>
#include <Rinternals.h>

#include <R.h>

// Constants for optimization
#define MAX_LOCI_STR 256
#define BUFFER_SIZE 1024
#define MAX_READ_LEN 1000
#define MAX_PATH 1000

// Processing modes
typedef enum {
    MODE_SNPCOUNTS,     // Simple nucleotide counting with BAM stats
    MODE_NTCOUNTS,      // Simple nucleotide counting
    MODE_SOMATICFREQ    // Cancer hotspot analysis with VAF calculation and HTML
} count_mode_t;

// Input format types
typedef enum {
    INPUT_SIMPLE,       // 2-column: chr, pos
    INPUT_VARIANT       // 8-column: chr, pos, ref, alt, gene, vc, pc, cosid
} input_format_t;

// Nucleotide to index mapping for faster lookup
static const int NT_TO_INDEX[256] = {
    ['A'] = 0, ['a'] = 0,
    ['T'] = 1, ['t'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['C'] = 3, ['c'] = 3,
    ['N'] = -1, ['n'] = -1
};

char *VERSION = "0.2.0";

// Function prototypes for HTML output (only used in somaticfreq mode)
void printhead(FILE *fn, char *bam_fn);
void printrow(FILE *fn, char *chrom, char *start, char *ref, char *alt, char *gene, char *vc, char *pc, char *cosid, float vaf, float counts[]);
void printargtbl(FILE *fn, int mapq, char *bam, int tot_var, int som_vars, float vaf, int cov, float doc, int treads);
void printfooter(FILE *fn);
char *basename(const char *path);
char *removeExt(char* myStr);

// Variant information structure
typedef struct {
    char *chrom;
    char *start;
    char *ref;
    char *alt;
    char *gene;
    char *vc;
    char *pc;
    char *cosid;
} variant_info_t;

// Processing parameters structure
typedef struct {
    uint32_t mapq;
    uint32_t flag;
    float vaf_threshold;
    int min_depth;
    int min_alt_reads;
    const char *fafile;
    const char *output_prefix;
    count_mode_t mode;
} processing_params_t;

// Optimized nucleotide counting using improved CIGAR parsing
void count_nucleotides_at_position(samFile *fp_in, hts_idx_t *fp_idx, bam_hdr_t *bamHdr, 
                                   const char *chrom, int32_t pos, uint32_t mapq, uint32_t flag,
                                   float *nt_counts, int *total_reads) {
    
    // Initialize counts
    memset(nt_counts, 0, 6 * sizeof(float));
    *total_reads = 0;
    
    // Create region string efficiently
    char region[MAX_LOCI_STR];
    snprintf(region, sizeof(region), "%s:%d-%d", chrom, pos + 1, pos + 1);
    
    // Query reads in region
    hts_itr_t *iter = sam_itr_querys(fp_idx, bamHdr, region);
    if (iter == NULL) {
        return;
    }
    
    bam1_t *aln = bam_init1();
    
    // Process each alignment
    while (sam_itr_next(fp_in, iter, aln) >= 0) {
        // Apply quality filters early
        if (aln->core.qual < mapq || (aln->core.flag & flag)) {
            continue;
        }
        
        (*total_reads)++;
        
        // Get alignment start and sequence
        int32_t aln_start = aln->core.pos;
        int32_t aln_end = bam_endpos(aln);
        
        // Check if alignment covers target position
        if (pos < aln_start || pos >= aln_end) {
            continue;
        }
        
        // Get reference position in read coordinates
        uint32_t *cigar = bam_get_cigar(aln);
        uint8_t *seq = bam_get_seq(aln);
        
        int32_t ref_pos = aln_start;
        int32_t read_pos = 0;
        
        // Walk through CIGAR to find target position
        for (int i = 0; i < aln->core.n_cigar; i++) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            
            switch (op) {
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF:
                    if (ref_pos <= pos && pos < ref_pos + len) {
                        // Target position found in match region
                        int offset = pos - ref_pos;
                        char base = seq_nt16_str[bam_seqi(seq, read_pos + offset)];
                        
                        int idx = NT_TO_INDEX[(unsigned char)base];
                        if (idx >= 0 && idx < 4) {
                            nt_counts[idx]++;
                        }
                        goto next_read;
                    }
                    ref_pos += len;
                    read_pos += len;
                    break;
                    
                case BAM_CINS:
                    if (ref_pos == pos) {
                        nt_counts[4]++; // Insertion count
                        goto next_read;
                    }
                    read_pos += len;
                    break;
                    
                case BAM_CDEL:
                    if (ref_pos <= pos && pos < ref_pos + len) {
                        nt_counts[5]++; // Deletion count
                        goto next_read;
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
            }
        }
        
        next_read:
        continue;
    }
    
    bam_destroy1(aln);
    hts_itr_destroy(iter);
}

// Calculate VAF based on variant type and alleles (for somaticfreq mode)
float calculate_vaf(const char *ref, const char *alt, const char *vc, float *nt_counts, 
                   int total_reads, int min_alt_reads, int *pass_reads) {
    
    *pass_reads = 1;
    float vaf = 0.0;
    
    if (total_reads == 0) {
        *pass_reads = 0;
        return -1.0;
    }
    
    if (strcmp(vc, "INS") == 0) {
        if (nt_counts[4] < min_alt_reads) {
            *pass_reads = 0;
        }
        vaf = nt_counts[4] / total_reads;
    } else if (strcmp(vc, "DEL") == 0) {
        if (nt_counts[5] < min_alt_reads) {
            *pass_reads = 0;
        }
        vaf = nt_counts[5] / total_reads;
    } else {
        // SNV
        int alt_idx = NT_TO_INDEX[(unsigned char)alt[0]];
        if (alt_idx >= 0 && alt_idx < 4 && strcmp(ref, "-") != 0) {
            if (nt_counts[alt_idx] < min_alt_reads) {
                *pass_reads = 0;
            }
            vaf = nt_counts[alt_idx] / total_reads;
        } else {
            vaf = -1.0;
        }
    }
    
    return vaf;
}

// Parse input line based on format
int parse_input_line(char *buffer, input_format_t format, variant_info_t *variant) {
    // Remove newline
    size_t len = strlen(buffer);
    if (len > 0 && buffer[len-1] == '\n') {
        buffer[len-1] = '\0';
        len--;
    }
    
    if (len == 0) return 0;
    
    if (format == INPUT_SIMPLE) {
        // Simple format: chr, pos
        variant->chrom = strtok(buffer, "\t");
        variant->start = strtok(NULL, "\t");
        
        if (!variant->chrom || !variant->start) {
            return 0;
        }
        
        // Set defaults for unused fields
        variant->ref = "";
        variant->alt = "";
        variant->gene = "";
        variant->vc = "";
        variant->pc = "";
        variant->cosid = "";
        
    } else if (format == INPUT_VARIANT) {
        // Variant format: chr, pos, ref, alt, gene, vc, pc, cosid
        variant->chrom = strtok(buffer, "\t");
        variant->start = strtok(NULL, "\t");
        variant->ref = strtok(NULL, "\t");
        variant->alt = strtok(NULL, "\t");
        variant->gene = strtok(NULL, "\t");
        variant->vc = strtok(NULL, "\t");
        variant->pc = strtok(NULL, "\t");
        variant->cosid = strtok(NULL, "\t");
        
        if (!variant->chrom || !variant->start || !variant->ref || !variant->alt || 
            !variant->gene || !variant->vc || !variant->pc || !variant->cosid) {
            return 0;
        }
    }
    
    return 1; // Success
}

// Get reference base from FASTA
char get_reference_base(faidx_t *fa, const char *chrom, const char *start_str) {
    if (!fa) {
        return 'N';
    }
    
    char loci_str[MAX_LOCI_STR];
    snprintf(loci_str, sizeof(loci_str), "%s:%s-%s", chrom, start_str, start_str);
    
    int seq_len;
    char *seq = fai_fetch(fa, loci_str, &seq_len);
    char ref_base = 'N';
    
    if (seq && seq_len > 0) {
        ref_base = seq[0];
        free(seq);
    }
    
    return ref_base;
}

// Unified BAM processing engine
void bamcounts_engine(const char *bam, const char *input_file, processing_params_t *params) {
    
    int vars_processed = 0;
    int som_vars = 0;        // For somaticfreq mode
    float mean_doc = 0.0;    // For somaticfreq mode
    int novaf = 0;           // For somaticfreq mode
    uint64_t total_mapped = 0; // For snpcounts mode
    
    // Determine input format and output files
    input_format_t input_format = (params->mode == MODE_SOMATICFREQ) ? INPUT_VARIANT : INPUT_SIMPLE;
    
    char html_file[MAX_PATH] = "";
    char tsv_file[MAX_PATH];
    
    if (params->mode == MODE_SOMATICFREQ) {
        // Create both HTML and TSV files
        char *bn = basename(bam);
        char *bnn = NULL;
        if (bn) {
            bnn = removeExt(bn);
        }
        
        if (params->output_prefix == NULL && bnn) {
            snprintf(html_file, sizeof(html_file), "%s.html", bnn);
            snprintf(tsv_file, sizeof(tsv_file), "%s.tsv", bnn);
        } else {
            snprintf(html_file, sizeof(html_file), "%s.html", params->output_prefix);
            snprintf(tsv_file, sizeof(tsv_file), "%s.tsv", params->output_prefix);
        }
        
        if (bn) free(bn);
        if (bnn) free(bnn);
    } else {
        // Only TSV file
        snprintf(tsv_file, sizeof(tsv_file), "%s.tsv", params->output_prefix);
    }
    
    // Suppress HTSlib warnings
    hts_verbose = 0;
    
    // Open input file
    FILE *input_fp = fopen(input_file, "r");
    if (!input_fp) {
        Rf_error("Cannot open input file: %s", input_file);
        return;
    }
    
    // Open output files
    FILE *html_fp = NULL;
    if (params->mode == MODE_SOMATICFREQ) {
        html_fp = fopen(html_file, "w");
        if (!html_fp) {
            Rf_error("Cannot create HTML file: %s", html_file);
            fclose(input_fp);
            return;
        }
    }
    
    FILE *tsv_fp = fopen(tsv_file, "w");
    if (!tsv_fp) {
        Rf_error("Cannot create TSV file: %s", tsv_file);
        fclose(input_fp);
        if (html_fp) fclose(html_fp);
        return;
    }
    
    // Open BAM file and index
    samFile *fp_in = hts_open(bam, "r");
    if (!fp_in) {
        Rf_error("Cannot open BAM file: %s", bam);
        fclose(input_fp);
        fclose(tsv_fp);
        if (html_fp) fclose(html_fp);
        return;
    }
    
    hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
    if (!fp_idx) {
        Rf_error("Cannot load BAM index for: %s", bam);
        sam_close(fp_in);
        fclose(input_fp);
        fclose(tsv_fp);
        if (html_fp) fclose(html_fp);
        return;
    }
    
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    if (!bamHdr) {
        Rf_error("Cannot read BAM header from: %s", bam);
        hts_idx_destroy(fp_idx);
        sam_close(fp_in);
        fclose(input_fp);
        fclose(tsv_fp);
        if (html_fp) fclose(html_fp);
        return;
    }
    
    // Load FASTA index if provided
    faidx_t *fa = NULL;
    if (strcmp(params->fafile, "NULL") != 0) {
        fa = fai_load(params->fafile);
    }
    
    // Calculate total mapped reads for snpcounts mode
    if (params->mode == MODE_SNPCOUNTS) {
        int32_t n_contigs = bamHdr->n_targets;
        int i;
        for (i = 0; i < n_contigs; i++) {
            uint64_t n_mapped = 0, n_unmapped = 0;
            if (hts_idx_get_stat(fp_idx, i, &n_mapped, &n_unmapped) == 0) {
                total_mapped += n_mapped;
            }
        }
    }
    
    // Print headers
    if (html_fp) {
        char *bn = basename(bam);
        char *bnn = NULL;
        if (bn) {
            bnn = removeExt(bn);
        }
        if (bnn) {
            printhead(html_fp, bnn);
        }
        if (bn) free(bn);
        if (bnn) free(bnn);
    }
    
    // Write TSV header based on mode
    if (params->mode == MODE_SNPCOUNTS) {
        fprintf(tsv_fp, "#idxstats_mapped_reads\t%llu\n", total_mapped);
        fprintf(tsv_fp, "loci\tfa_ref\tA\tT\tG\tC\tIns\tDel\n");
    } else if (params->mode == MODE_NTCOUNTS) {
        fprintf(tsv_fp, "loci\tfa_ref\tA\tT\tG\tC\tIns\tDel\n");
    } else if (params->mode == MODE_SOMATICFREQ) {
        fprintf(tsv_fp, "loci\tfa_ref\tNT_change\tHugo_Symbol\tVariant_Classification\tAA_change\tMeta\tVAF\tA\tT\tG\tC\tIns\tDel\n");
        
        // Print processing information
        Rprintf("Input BAM file           :  %s\n", bam);
        Rprintf("Variants                 :  %s\n", input_file);
        Rprintf("VAF filter               :  %.3f\n", params->vaf_threshold);
        Rprintf("min reads for t_allele   :  %d\n", params->min_alt_reads);
        Rprintf("MAPQ filter              :  %d\n", params->mapq);
        Rprintf("FLAG filter              :  %d\n", params->flag);
        Rprintf("Coverage filter          :  %d\n", params->min_depth);
        Rprintf("HTSlib version           :  %s\n\n", hts_version());
    }
    
    // Process each line in the input file
    char buffer[BUFFER_SIZE];
    variant_info_t variant;
    
    while (fgets(buffer, sizeof(buffer), input_fp)) {
        if (!parse_input_line(buffer, input_format, &variant)) {
            continue;
        }
        
        int32_t target_pos = atoi(variant.start) - 1; // Convert to 0-based
        
        // Get reference base
        char ref_base = get_reference_base(fa, variant.chrom, variant.start);
        
        // Count nucleotides at this position
        float nt_counts[6];
        int total_reads;
        count_nucleotides_at_position(fp_in, fp_idx, bamHdr, variant.chrom, target_pos, 
                                     params->mapq, params->flag, nt_counts, &total_reads);
        
        vars_processed++;
        
        if (params->mode == MODE_SOMATICFREQ) {
            mean_doc += total_reads;
            
            if (vars_processed % 1000 == 0) {
                Rprintf("Processed %d entries..\n", vars_processed);
            }
            
            // Calculate VAF and determine if variant passes filters
            int pass_reads;
            float vaf = calculate_vaf(variant.ref, variant.alt, variant.vc, nt_counts, 
                                    total_reads, params->min_alt_reads, &pass_reads);
            
            // Write TSV output
            fprintf(tsv_fp, "%s:%s\t%c\t%s>%s\t%s\t%s\t%s\t%s\t%.3f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",
                    variant.chrom, variant.start, ref_base, variant.ref, variant.alt, 
                    variant.gene, variant.vc, variant.pc, variant.cosid, vaf,
                    nt_counts[0], nt_counts[1], nt_counts[2], nt_counts[3], nt_counts[4], nt_counts[5]);
            
            // Write HTML output for variants that pass filters
            if (!isnan(vaf) && vaf >= 0) {
                if (vaf >= params->vaf_threshold && total_reads > params->min_depth && pass_reads == 1) {
                    som_vars++;
                    printrow(html_fp, variant.chrom, variant.start, variant.ref, variant.alt, 
                           variant.gene, variant.vc, variant.pc, variant.cosid, vaf, nt_counts);
                }
            } else {
                novaf++;
            }
            
        } else {
            // Simple modes (SNPCOUNTS, NTCOUNTS)
            if (vars_processed % 1000 == 0) {
                Rprintf("Processed %d entries..\n", vars_processed);
            }
            
            fprintf(tsv_fp, "%s:%s\t%c\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",
                    variant.chrom, variant.start, ref_base,
                    nt_counts[0], nt_counts[1], nt_counts[2], nt_counts[3], 
                    nt_counts[4], nt_counts[5]);
        }
    }
    
    // Finalize outputs
    if (params->mode == MODE_SOMATICFREQ) {
        mean_doc = mean_doc / vars_processed;
        
        char *bn = basename(bam);
        if (bn) {
            printargtbl(html_fp, params->mapq, bn, vars_processed, som_vars, 
                       params->vaf_threshold, params->min_depth, mean_doc, params->min_alt_reads);
            free(bn);
        }
        printfooter(html_fp);
        
        Rprintf("Done!\n\n");
        Rprintf("Summary:\n");
        Rprintf("Total variants processed :  %d\n", vars_processed);
        Rprintf("Variants > %.2f threshold:  %d\n", params->vaf_threshold, som_vars);
        Rprintf("Avg. depth of coverage   :  %.2f\n", mean_doc);
        Rprintf("Output html report       :  %s\n", html_file);
        Rprintf("Output TSV file          :  %s\n", tsv_file);
        
        if (novaf > 0) {
            Rf_warning("Could not estimate VAF for %d variants. VAF has been set to -1 in %s.\nManually inspect them in IGV.\n", novaf, tsv_file);
        }
    }
    
    // Cleanup
    if (fa) {
        fai_destroy(fa);
    }
    bam_hdr_destroy(bamHdr);
    hts_idx_destroy(fp_idx);
    sam_close(fp_in);
    fclose(input_fp);
    fclose(tsv_fp);
    if (html_fp) fclose(html_fp);
}

// Utility functions (from original somaticfreq.c)
char *basename(const char *path) {
    const char *s = strrchr(path, '/');
    if (!s) {
        return strdup(path);
    } else {
        return strdup(s + 1);
    }
}

char *removeExt(char* myStr) {
    char *retStr;
    char *lastExt;
    if (myStr == NULL) return NULL;
    if ((retStr = malloc(strlen(myStr) + 1)) == NULL) return NULL;
    strcpy(retStr, myStr);
    lastExt = strrchr(retStr, '.');
    if (lastExt != NULL)
        *lastExt = '\0';
    return retStr;
}

void printhead(FILE *fn, char *bam_fn) {
    fprintf(fn, "<!DOCTYPE html>\n");
    fprintf(fn, "<html lang=\"en\">");

    fprintf(fn, "<head>\n");
    fprintf(fn, "<meta charset=\"UTF-8\">\n");
    fprintf(fn, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
    fprintf(fn, "<title>%s | somatic variants</title>\n", bam_fn);
    fprintf(fn, "<link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css\">\n");
    fprintf(fn, "<script src=\"https://code.jquery.com/jquery-3.5.1.js\"></script>\n");
    fprintf(fn, "<script src=\"https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js\"></script>\n");
    fprintf(fn, "<style>\n");
    fprintf(fn, "body {padding: 0px; margin: 0; font-family: Verdana, Geneva, Tahoma, sans-serif;}\n");
    fprintf(fn, "tr {transition: all .2s ease-in; cursor: pointer; color: #2c3e50;}\n");
    fprintf(fn, "th, td {padding: 12px; text-align: left; border-bottom: 1px solid #ddd;}\n");
    fprintf(fn, "#header {background-color: #2c3e50; color: #fff;}\n");
    fprintf(fn, "h1 {font-weight: 600; text-align: left; color: #2c3e50; padding: 10px 0px;}\n");
    fprintf(fn, "h3 {text-align: left; color: #c0392b; padding: 5px}\n");
    fprintf(fn, "tr:hover {background-color: #f5f5f5; transform: scale(1.02); box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n");
    fprintf(fn, "@media only screen and (max-width: 768px) {table {width: 90%%;}}\n");
    fprintf(fn, ".details tr { line-height: 15px; }\n");
    fprintf(fn, ".details th, td {padding: 5px; text-align: left; border-bottom: 1px solid #ddd; font-family:'Courier New', Courier, monospace;}\n");
    fprintf(fn, "</style>\n");
    fprintf(fn, "<script>\n");
    fprintf(fn, "$(document).ready(function() {\n");
    fprintf(fn, "$('#cosmic').DataTable();});\n");
    fprintf(fn, "</script>\n");
    fprintf(fn, "</head>\n");

    fprintf(fn, "<h3> <a href=\"https://www.cancerhotspots.org/\">Cancer hotspots</a> </h3>");

    time_t t = time(NULL);
    struct tm cutime = *localtime(&t);
    fprintf(fn, "<table class=\"details\">\n");
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Date</td><td>%d-%02d-%02d %02d:%02d:%02d</td></tn>", cutime.tm_year + 1900, cutime.tm_mon + 1, cutime.tm_mday, cutime.tm_hour, cutime.tm_min, cutime.tm_sec);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Sample</td><td>%s</td></tn>", bam_fn);
    fprintf(fn, "</table>\n");

    fprintf(fn, "<p style=\"margin-top:2.5em\"> \n</p>\n");
    fprintf(fn, "<p style=\"margin-top:2.5em\"> \n</p>\n");
    fprintf(fn, "<p style=\"margin-top:2.5em\"> \n</p>\n");

    fprintf(fn, "<table id=\"cosmic\" class=\"display\">\n");
    fprintf(fn, "<thead><tr id=\"header\"><th>Chr</th><th>Pos</th><th>Ref</th><th>Alt</th><th>Gene</th><th>Type</th><th>AA change</th><th>Meta</th><th>VAF</th><th>A</th><th>T</th><th>G</th><th>C</th><th>Ins</th><th>Del</th></tr></thead>\n");
    fprintf(fn, "<tbody>\n");
}

void printargtbl(FILE *fn, int mapq, char *bam, int tot_var, int som_vars, float vaf, int cov, float doc, int treads) {
    fprintf(fn, "</table>\n"); // end previous table before printing summary
    fprintf(fn, "<h3 >Summary </h3>");
    fprintf(fn, "<table class=\"details\">\n");
    fprintf(fn, "<tr><td style=\"font-weight:bold\">BAM</td><td>%s</td></tn>", bam);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Variants queried</td><td>%d</td></tn>", tot_var);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Variants with somatic evidance</td><td>%d [VAF > %.2f]</td></tn>", som_vars, vaf);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Avg. depth of coverage</td><td>%.2f</td></tn>", doc);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">MAPQ filter</td><td>%d</td></tn>", mapq);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Coverage filter</td><td>%d</td></tn>", cov);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">VAF filter</td><td>%.2f</td></tn>", vaf);
    fprintf(fn, "<tr><td style=\"font-weight:bold\">Min. number of reads supporting tumor allele </td><td>%d</td></tn>", treads);
    fprintf(fn, "</table>\n");
    fprintf(fn, "</tbody>\n");
    fprintf(fn, "</body>\n");
    fprintf(fn, "</html>\n");
}

void printrow(FILE *fn, char *chrom, char *start, char *ref, char *alt, char *gene, char *vc, char *pc, char *cosid, float vaf, float counts[]) {
    fprintf(fn, "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><a target=\"_blank\" href=\"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=%s#variants\"> %s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td><td>%.0f</td><td>%.0f</td><td>%.0f</td><td>%.0f</td><td>%.0f</td><td>%.0f</td></tr>\n",
        chrom, start, ref, alt, gene, gene, vc, pc, cosid, vaf, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]);
}

void printfooter(FILE *fn) {
    fprintf(fn, "<p style=\"margin-top:2.5em\"> \n</p>\n");
    fprintf(fn, "<div id=\"footer\">\n");
    fprintf(fn, "<span style=\"float:left;font-family:'Courier New', Courier, monospace; padding: 5px;\" >Generated by <a href=\"https://github.com/PoisonAlien/maftools\">maftools::cancerhotspots()</a> </span>");
    fprintf(fn, "</div>");
}

// R interface compatibility wrapper functions
SEXP snpc(SEXP filename, SEXP bedname, SEXP qual, SEXP flag, SEXP fa, SEXP op_file) {
    processing_params_t params = {
        .mapq = Rf_asInteger(qual),
        .flag = Rf_asInteger(flag),
        .vaf_threshold = 0.0,
        .min_depth = 0,
        .min_alt_reads = 0,
        .fafile = Rf_translateChar(Rf_asChar(fa)),
        .output_prefix = Rf_translateChar(Rf_asChar(op_file)),
        .mode = MODE_SNPCOUNTS
    };
    
    bamcounts_engine(
        Rf_translateChar(Rf_asChar(filename)), 
        Rf_translateChar(Rf_asChar(bedname)),
        &params
    );
    
    return R_NilValue;
}

SEXP ntc(SEXP filename, SEXP bedname, SEXP qual, SEXP flag, SEXP fa, SEXP op_file) {
    processing_params_t params = {
        .mapq = Rf_asInteger(qual),
        .flag = Rf_asInteger(flag),
        .vaf_threshold = 0.0,
        .min_depth = 0,
        .min_alt_reads = 0,
        .fafile = Rf_translateChar(Rf_asChar(fa)),
        .output_prefix = Rf_translateChar(Rf_asChar(op_file)),
        .mode = MODE_NTCOUNTS
    };
    
    bamcounts_engine(
        Rf_translateChar(Rf_asChar(filename)), 
        Rf_translateChar(Rf_asChar(bedname)),
        &params
    );
    
    return R_NilValue;
}

SEXP readb(SEXP filename, SEXP bedname, SEXP depth, SEXP t_depth, SEXP vaf, SEXP fa, SEXP op_file, SEXP qual, SEXP flag) {
    processing_params_t params = {
        .mapq = Rf_asInteger(qual),
        .flag = Rf_asInteger(flag),
        .vaf_threshold = Rf_asReal(vaf),
        .min_depth = Rf_asInteger(depth),
        .min_alt_reads = Rf_asInteger(t_depth),
        .fafile = Rf_translateChar(Rf_asChar(fa)),
        .output_prefix = Rf_translateChar(Rf_asChar(op_file)),
        .mode = MODE_SOMATICFREQ
    };
    
    bamcounts_engine(
        Rf_translateChar(Rf_asChar(filename)), 
        Rf_translateChar(Rf_asChar(bedname)),
        &params
    );
    
    return R_NilValue;
}
