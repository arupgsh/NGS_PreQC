#!/usr/bin/python3

import gzip
import sys
import re

# Code source: https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py
# Patterns from : https://github.com/stjudecloud/ngsderive/blob/master/ngsderive/commands/instrument.py


# dictionary of instrument id regex: [platform(s)]
sequencer_ids = {
        "HWI-M[0-9]{4}$" : ["MiSeq"],
        "HWUSI" : ["Genome Analyzer IIx"],
        "M[0-9]{5}$" : ["MiSeq"],
        "HWI-C[0-9]{5}$" : ["HiSeq 1500"],
        "C[0-9]{5}$" : ["HiSeq 1500"],
        "HWI-D[0-9]{5}$" : ["HiSeq 2500"],
        "D[0-9]{5}$" : ["HiSeq 2500"],
        "J[0-9]{5}$" : ["HiSeq 3000"],
        "K[0-9]{5}$" : ["HiSeq 3000","HiSeq 4000"],
        "E[0-9]{5}$" : ["HiSeq X"],
        "NB[0-9]{6}$": ["NextSeq"],
        "NS[0-9]{6}$" : ["NextSeq"],
        "MN[0-9]{5}$" : ["MiniSeq"],
        "A[0-9]{5}$" :  ["NovaSeq"]
        }

# dictionary of flow cell id regex: ([platform(s)], flow cell version and yeild)
flowcell_ids = {
    "^C[A-Z0-9]{4}ANXX$": (["HiSeq 1500","HiSeq 2000",  "HiSeq 2500",], "High Output (8-lane) v4 flow cell"),
    "^C[A-Z0-9]{4}ACXX$": (["HiSeq 1000","HiSeq 1500",  "HiSeq 2000","HiSeq 2500"], "High Output (8-lane) v3 flow cell"),
    "^D[A-Z0-9]{4}ACXX$": (["HiSeq 1000","HiSeq 1500",  "HiSeq 2000","HiSeq 2500"], "High Output (8-lane) v3 flow cell"),
    "^H[A-Z0-9]{4}ADXX$": (["HiSeq 1500","HiSeq 2000",  "HiSeq 2500"], "Rapid Run (2-lane) v1 flow cell"),
    "^H[A-Z0-9]{4}BCXX$": (["HiSeq 1500","HiSeq 2500"], "Rapid Run (2-lane) v2 flow cell"),
    "^H[A-Z0-9]{4}BCXY$": (["HiSeq 1500","HiSeq 2500"], "Rapid Run (2-lane) v2 flow cell"),
    "^H[A-Z0-9]{4}BBXX$": (["HiSeq 4000"],  "(8-lane) v1 flow cell"),
    "^H[A-Z0-9]{4}BBXY$": (["HiSeq 4000"], "(8-lane) v1 flow cell"),
    "^H[A-Z0-9]{4}CCXX$": (["HiSeq X"], "(8-lane) flow cell"),
    "^H[A-Z0-9]{4}CCXY$": (["HiSeq X"], "(8-lane) flow cell"),
    "^H[A-Z0-9]{4}ALXX$": (["HiSeq X"], "(8-lane) flow cell"),
    "^H[A-Z0-9]{4}BGX[A-Z,0-9]$": (["NextSeq"], "High output flow cell"),
    "^H[A-Z0-9]{4}AFXX$": (["NextSeq"], "Mid output flow cell"),
    "^H[A-Z0-9]{5}RXX$": (["NovaSeq"], "S1 flow cell"),
    "^H[A-Z0-9]{5}RXX$": (["NovaSeq"], "SP flow cell"),
    "^H[A-Z0-9]{5}MXX$": (["NovaSeq"], "S2 flow cell"),
    "^H[A-Z0-9]{5}SXX$": (["NovaSeq"], "S4 flow cell"),
    "^A[A-Z0-9]{4}$": (["MiSeq"], "MiSeq flow cell"),
    "^B[A-Z0-9]{4}$": (["MiSeq"], "MiSeq flow cell"),
    "^D[A-Z0-9]{4}$": (["MiSeq"], "MiSeq nano flow cell"),
    "^G[A-Z0-9]{4}$": (["MiSeq"], "MiSeq micro flow cell")
     # "^D[A-Z0-9]{4}$" : ["HiSeq 2000", "HiSeq 2500"],  # Unknown HiSeq flow cell examined in SJ data
}

# do intersection of lists
def intersect(a, b):
    return list(set(a) & set(b))

def union(a, b):
    return list(set(a) | set(b))

def instrument_data(instrument_id):
    #instrument_name = []
    for instrument,instrument_name in sequencer_ids.items():
        if re.search(instrument,instrument_id):
            return instrument_name

def flowcell_data(flowcell_id):
    for flowcell,details in flowcell_ids.items():
        if re.search(flowcell,flowcell_id):
            fc_instrument_name = details[0]
            flowcell_name = details[1]
            return fc_instrument_name,flowcell_name

def info_validator(instrument_id, flowcell_id):
    instrument_name = instrument_data(instrument_id)
    if flowcell_data(flowcell_id):
        fc_instrument_name,flowcell_name = flowcell_data(flowcell_id)
    else:
        fc_instrument_name,flowcell_name = ["",""]

    sequencer_details = {}

    if not flowcell_name:
        flowcell_name = "Failed"

    if not instrument_name and not fc_instrument_name:
        sequencer_details["seq_name"] = "Failed"
        sequencer_details["qual"] = "NA"
    if not instrument_name:
        sequencer_details["seq_name"] = fc_instrument_name
        sequencer_details["qual"] = "Medium"
    if fc_instrument_name:
        sequencer_details["seq_name"] = instrument_name
        sequencer_details["qual"] = "Medium"
    
    #when instrument name supported by both datasets
    instrument = intersect(instrument_name, fc_instrument_name)
    if instrument:
        sequencer_details["seq_name"] = instrument
        sequencer_details["qual"] = "High"
    else:
        instrument = union(instrument_name, fc_instrument_name)
        sequencer_details["seq_name"] = instrument
        sequencer_details["qual"] = "Low"
    sequencer_details["seq_name"] = '/'.join(sequencer_details["seq_name"])
    sequencer_details["flowcell_name"] = flowcell_name
    print(sequencer_details)


def detect_instument(fastq):
    with gzip.open(fastq, 'rt') as f:
        for i in range(4):
            line = f.readline()
            if line.startswith('@'):
                seq_header = line.split(':')
                instrument_id = seq_header[0].strip('@')
                flowcell_id = seq_header[2]
                info_validator(instrument_id,flowcell_id)

detect_instument(sys.argv[1]) # this will take fastq.gz as input and return a dictionary
