#include "array-list.h"
#include "momo.h"

// Structure for tracking momo command line algorithm.
typedef enum motifx_status {
  ACTIVE,
  INACTIVE,
  DELETED
} MOTIFX_STATUS_T;

typedef struct regexmotif {
  bool* conserved;
  MATRIX_T* residues;
//  bool* representation;
//  ARRAYLST_T* residues;
}  REGEX_MOTIF_T;

MATRIX_T* get_count_matrix(MATRIX_T* count,
                           ARRAYLST_T* sequences,
                           MOTIFX_STATUS_T** status,
                           MOMO_OPTIONS_T* options,
                           SUMMARY_T* summary);

void print_matrix_to_terminal(MATRIX_T* matrix, MOMO_OPTIONS_T* options, SUMMARY_T* summary);

void create_motifs(MOMO_OPTIONS_T* options,
                   SUMMARY_T* summary,
                   SEQ_T** all_sequences,
                   int num_sequences);
