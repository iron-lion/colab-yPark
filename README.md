# Genome exploration

- [Enriched Motif Search (C. elegans, Mouse, Killifish)](https://colab.research.google.com/github/iron-lion/colab-yPark/blob/main/notebooks/promoter_k_mer.ipynb)
This program identifies the most common DNA sequence motifs within the promoter regions of a given list of genes. It extracts upstream promoter sequences (typically a fixed length upstream of the transcription start site) from a reference genome and search motifs (with [JASPAR](https://jaspar2022.genereg.net/) database)
Each motif is then counted based on the number of occurrences across all input gene regions. Finally, the motifs are ranked and listed according to their frequency, reflecting how often each motif appears within the gene set.
(Similar to [FIMO](https://meme-suite.org/meme/tools/fimo), but easy to use)

# Fluorescent microscopy image process

- [C. elegans Head signal counter](https://colab.research.google.com/github/iron-lion/colab-yPark/blob/main/notebooks/worm_head_signal_detect.ipynb)
Segment head part of Worm and calculate peaks of fluorescent signals.

- [C. elegans segmentation and signal counter](https://colab.research.google.com/github/iron-lion/colab-yPark/blob/main/notebooks/micro_sam_worm_body_signal_detect.ipynb)
Segment the entire body of the worm (or cells) and calculate peaks of fluorescent signals.
