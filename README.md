# Telomere_2023

Code for TeloReader presented in the paper 'Telomerase-independent survival leads to a mosaic of complex subtelomere rearrangements in Chlamydomonas reinhardtii' Frederic Chaux, Nicolas Agier, Clotilde Garrido, Gilles Fischer, Stephan Eberhard, Zhou Xu

@author: Clotilde Garrido, Sorbonne Université - CNRS

# Dependencies
Python 3.9.12 
pandas 1.4.3
numpy 1.19.5
matplotlib 3.5.1

# Installation
Download directory TELOREADER and install the dependencies.

# Usage
python Teloreader.py <strain> <path> <fasta>

# Help
python -h Teloreader.py

# Output
  Two text files corresponding to the G-rich- and C-rich-specific 8-mer scores calculated from the fasta file.
  
  Output directory can by specify with option -o.
  
  One fasta file coresponding to all found telomeric sequences.
  
  One csv file summarizing the results.
