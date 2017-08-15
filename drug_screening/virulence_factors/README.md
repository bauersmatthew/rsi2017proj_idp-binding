Contains virulence factor structure and search space data.

## Files

- disordered_pdbids_no-partials_no-specificforms.tsv contains selected RCSB PDB
  IDs of virulence factor structures.

- get_pdbs.py is a python3 script that downloads PDB files from the RSCB PDB
  when given PDB IDs.

- idr_locs.tsv is a table containing the locations of the largest IDR in each
  virulence factor.

- model_quality_order.txt is a list of selected structures ranked in decreasing
  order by quality.

- pdb_parser.py is a python3 module containing routines for parsing and editing
  PDB files.

- relevant_chains_1only.tsv is a table containing which "chain" in each VF
  structure is actually the VF molecule that we are interested in.

- search_spaces.tsv is a table containing the C-side and N-side search spaces on
  each VF structure for which the IDR was found in the structure.

- select_chains.py is a python3 script that edits PDB files to contain only the
  selected chains.

- select_search_spaces.py is a python3 script that finds the coordinates of 5 aa
  long search spaces on the N-side and C-side of the long IDR in VF structures.

- tb_vf_aaseqs.fa is a FASTA file containing the amino acid sequences of all TB
  virulence factors.

## Subdirectories

- pdbqt/ contains PDBQT-formatted virulence factor structures.
