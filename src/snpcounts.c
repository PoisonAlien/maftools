#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <math.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <R.h>

int countlines(char *filename){

  FILE *fp = fopen(filename,"r");
  int ch=0;
  int lines=0;

  if (fp == NULL){
    return 0;
  }

  while(!feof(fp)){
    ch = fgetc(fp);
    if(ch == '\n'){
      lines++;
    }
  }
  fclose(fp);
  return lines;
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// Progress bar modified from: https://stackoverflow.com/a/36315819
void printProgress(double percentage, int done, int tot) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  fprintf(stderr, "\r%3d%% [%.*s%*s] %d/%d", val, lpad, PBSTR, rpad, "", done, tot);
  fflush(stdout);
}

void snpcounts(const char *bam, const char *bedfile, uint32_t q, uint32_t F, const char *fafile, const char *op){

  int vars_gt = 1; //No. of BED entries

  hts_verbose = 0; //suppresses htslib warnings

  char tsv_file[1000];
  strcpy(tsv_file, op);
  strcat(tsv_file, ".tsv");

  //Open bed file
  //int nloci = countlines(bedfile);

  FILE *bed_fp;
  bed_fp = fopen(bedfile, "r");
  char buff[1000];

  //Open TSV report file
  FILE *tsv_fp;
  tsv_fp = fopen(tsv_file, "w" );

  //fasta file
  char *seq;
  faidx_t *fa = fai_load(fafile);

  //BAM file
  samFile *fp_in = hts_open(bam,"r"); //open bam file
  hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
  bam1_t *aln = bam_init1(); //initialize an alignment

  uint64_t n_mapped = 0;
  uint64_t n_unmapped = 0;
  uint64_t tot_mapped = 0;
  int res = 0;

  int32_t n_contigs =bamHdr->n_targets;
  for(int i = 0; i < n_contigs; i++){
    res = hts_idx_get_stat(fp_idx, i, &n_mapped, &n_unmapped);
    if(res == 0){
      tot_mapped = tot_mapped + n_mapped;
      }
    }

  fprintf(tsv_fp, "#idxstats_mapped_reads\t%llu\n", tot_mapped);
  fprintf(tsv_fp, "loci\tfa_ref\tA\tT\tG\tC\tIns\tDel\n");

  //For every loci in the BED file
  while(fgets(buff,1000,bed_fp) != NULL){

    //Remove trailing new line chars
    int len = strlen(buff);
    if(buff[len-1] == '\n' ){
      buff[len-1] = 0;
    }

    char *chrom = strtok(buff,"\t");
    char *start = strtok(NULL,"\t");

    char loci[250] = "";
    strcat(loci, chrom); strcat(loci, ":"); strcat(loci, start); strcat(loci, "-"); strcat(loci, start);

    //Fetch base at target loci from fasta file
    if(fa != NULL){
      int templen = 100;
      seq = fai_fetch(fa, loci, &templen);
      fprintf(tsv_fp, "%s:%s\t%s", chrom, start, seq);
      free(seq);
    }else{
      fprintf(tsv_fp, "%s:%s\tNA", chrom, start);
    }

    int32_t target_pos = atoi(start) -1; //input position are 1 based

    //load reads in target loci
    hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, loci);

    //Keep track of total reads and nt counts per loci
    int32_t tot_reads = 0;
    float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    //printProgress(vars_gt/(double)nloci, vars_gt, nloci);
    vars_gt = vars_gt + 1;


//if(vars_gt % 1000 == 0){
  //    Rprintf("%d | ", vars_gt); //verbose processed loci
    //}

    //For every read in the BAM file of target region
    while(sam_itr_next(fp_in, samitr, aln) > 0){

      int32_t pos = aln->core.pos ; //left most position of alignment in zero based coordinate (0-based)
      uint32_t len = aln->core.l_qseq; //length of the read.
      uint32_t* cig = bam_get_cigar(aln);
      uint8_t *qs = bam_get_seq(aln); //quality string
      //char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
      //uint16_t samflag = aln->core.flag; // flag
      //uint32_t q2 = aln->core.qual ; //mapping quality


      //MAPQ and FLAG filter
      if(aln->core.qual <= q){
        continue;
      }

      if(aln->core.flag >= F){
        continue;
      }

      tot_reads = tot_reads +1;
      //char ins_seq[100]; char del_seq[100];

      //get nucleotide id and converts them into IUPAC id.
      char *qseq = (char *)malloc(len);
      int i = 0;
      for(i=0; i< len ; i++){
        qseq[i] = seq_nt16_str[bam_seqi(qs,i)];
      }

      //target position on the read
      int32_t pos_onread = 0;

      //For every CIGAR string
      int k = 0;
      for(k=0;k< aln->core.n_cigar ;++k){
        int cop =cig[k] & BAM_CIGAR_MASK; // CIGAR string
        int cl = cig[k] >> BAM_CIGAR_SHIFT; // CIGAR length

        if(BAM_CIGAR_STR[cop] == 'M'){
          pos_onread = pos_onread + cl;
          pos = pos + cl;
        }else if(BAM_CIGAR_STR[cop] == 'S'){
          pos_onread = pos_onread + cl;
        }else if(BAM_CIGAR_STR[cop] == 'I'){
          pos_onread = pos_onread + cl;
        }else if(BAM_CIGAR_STR[cop] == 'D'){
          pos = pos + cl;
        }

        if(pos > target_pos){
          if(BAM_CIGAR_STR[cop] == 'M'){
            pos_onread = pos_onread - (pos - target_pos);
            if(qseq[pos_onread] == 'A'){
              nt[0] = nt[0] + 1;
            }else if(qseq[pos_onread] == 'T'){
              nt[1] = nt[1] + 1;
            }else if(qseq[pos_onread] == 'G'){
              nt[2] = nt[2] + 1;
            }else if(qseq[pos_onread] == 'C'){
              nt[3] = nt[3] + 1;
            }
            break;
          }
        }else if(pos == target_pos){
          if(BAM_CIGAR_STR[cop] == 'I'){
            nt[4] = nt[4] + 1;
            // insertion sequence
            // for(int i = 0; i < cl; i++){
              //     //strcat(ins_seq, &qseq[pos_onread+i]);
              // }
            break;
          }else if(BAM_CIGAR_STR[cop] == 'D'){
            nt[5] = nt[5] + 1;
            // deletion sequence
            // for(int i = 0; i < cl; i++){
              //     strcat(del_seq, qseq[pos_onread+i]);
              // }
            break;
          }
        }
      }

      free(qseq);
    }

    hts_itr_destroy(samitr);
    fprintf(tsv_fp, "\t%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n",  nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);
  }

  //fprintf(tsv_fp, "\t%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n",  nt[0], nt[1], nt[2], nt[3], nt[4], nt[5]);

  bam_destroy1(aln);
  bam_hdr_destroy(bamHdr);
  fai_destroy(fa);
  sam_close(fp_in);
  fclose(bed_fp);
  fclose(tsv_fp);
  //Rprintf("\n Done!");
}

SEXP snpc(SEXP filename, SEXP bedname, SEXP qual, SEXP flag, SEXP fa, SEXP op_file){
  //SEXP tbl = PROTECT(Rf_allocVector(STRSXP, 1));
  //char *parse_bam (const char *bam, const char *bedfile, int d, int t, float v, const char *fafile, const char *op, uint32_t q, uint32_t F){
  snpcounts(Rf_translateChar(Rf_asChar(filename)), Rf_translateChar(Rf_asChar(bedname)),Rf_asInteger(qual), Rf_asInteger(flag),
             Rf_translateChar(Rf_asChar(fa)), Rf_translateChar(Rf_asChar(op_file)));
    //nt[0] = 1;
    return 0;
  }
