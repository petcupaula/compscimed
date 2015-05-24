Implementation of some of the algorithms used in bioinformatics.


---


**Get a local copy** of the _compscimed_ repository with the command:

`hg clone https://code.google.com/p/compscimed/ path/to/your/local/copy`


---


**Example**: Performing a global sequence alignment of two amino-acid sequences (HEAGAWGHEE and PAWHEAE), using the Needleman-Wunsch algorithm, with a linear gap model with gap penalty of 8, using the BLOSUM50 scoring matrix:

`python align.py global blosum50 8 HEAGAWGHEE PAWHEAE >> alignment.txt`

**Example**: Performing a local sequence alignment of two amino-acid sequences (HEAGAWGHEE and PAWHEAE), using the Smith-Waterman algorithm, with an affine gap model with gap open penalty of 10 and gap extend penalty of 1, using the BLOSUM50 scoring matrix:

`python align_affine.py local blosum50 10 1 HEAGAWGHEE PAWHEAE >> alignment.txt`


---


**Example**: Printing some statistics about a 23andMe data file:

```
python txt2p.py genome_you.txt genome_you.p #read data file
python snpstats.py genome_you.p >> snpstats.txt
```

**Example**: Extract the SNP data corresponding to chromosome 15 from your 23andMe data file:

```
python txt2p.py genome_you.txt genome_you.p #read data file
python extract15.py genome_you.p chr15.p
```