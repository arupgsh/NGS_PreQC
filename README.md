# PreQC

PreQC is a collection of Python3 scripts to extract informations like instrumentused for sequencing, flowcell name, read length distribution and library strandedness from raw sequencing data.

## Scripts currently available

` fast_illumina_seq_detector.py ` : This script detects the name of sequencer and flowcell from Illumina sequencing data.

**Usage** 

```

$ python3 fast_illumina_seq_detector.py file.fastq.gz

{'seq_name': 'NovaSeq', 'qual': 'High', 'flowcell_name': 'S4 flow cell'}

```


