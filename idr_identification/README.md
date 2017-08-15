Holds data and scripts relevant to the part of the project where Tuberculosis
virulence factors were searched for disordered regions.

## Files

- disordered_pdbids_no-partials_no-specificforms.tsv contains selected RCSB PDB
  IDs of virulence factor structures.

- disordered_pdbids.tsv is similar to
  disordered_pdbids_no-partials_no-specificforms.tsv but it is less selective.

- disorder_info.tsv is a table containing information on the longest IDRs in
  each virulence factor, as well as the overall % disorder for each VF.

- query_pondr.py is a python3 script that automatically runs the PONDR VLXT
  disorder predictor on VF amino acid sequences.

- sorted_disordered.tsv is a table of all VFs sorted by "how disordered" they
  are.

- tb_vf_aaseqs.fa is a FASTA file containing amino acid sequences of all
  virulence factors.

- top150_disordered_with-pdbids.tsv is a table containing the "150 most
  disordered" virulence factors, along with what PDB entries were found for them
  (if any).

- TubercuList_R27_FASTA.txt is a FASTA file containing the amino acid sequences
  of all TB proteins in the R27 TubercuList database.

## Subdirectories

- nonvf/ contains disorder information on non-VF TB proteins.
