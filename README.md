# PCR-Sim

This tool aims to predict the results of a given set of diagnostic TaqMan RT-qPCRs on a set of influenza genome sequences through assessing the level of complementarity between each target sequence and the sequence of the primers and probes of each RT-qPCR. The complementarity between each primer/probe sequence and the target is assessed via the predicted free energy of hybridisation associated with duplex formation between each pair of sequences, calculated via system calls to UNAFOLD. These values are normalised to the binding energy expected of a perfect match for each primer/probe sequence to allow comparison between different RT-qPCRs. Mismatches in the 3 bases at the 3' (primers) and 5' (probe) are weighted more heavily. Using primer/probe binding values of RT-qPCRs on viruses with known results, thresholds for the normalised primer/probe binding energy values were set to allow prediction of RT-qPCR results for unknown RT-qPCRs and unknown viruses.

USAGE:
python pcr_sim.py [arguments_file]

An example arguments file is included in this repository. The arguments specified are:
 - genomes_file_path - the path to a FASTA file containined the viral genome sequences to be assessed. The format of the FASTA headers should follow the standard GISAID format ('>[strain name] | [segment number]'). A helper script, reformat.py, is included in this repository to convert FASTA files from the FluDB format for use with this tool.
 
 - pcr_conditions_file_path - the path to a .csv file containing the sequences of each RT-qPCR primer/probe (in 5'-3' orientation) as well as some other information regarding the reaction conditions of each RT-qPCR. An example can be found in this repository as pcr_conditions.csv.
 
 - output_file - the path specifying where the main output file should be created.
 
 - max_amp_size - the upper limit for potential amplicon lengths (in nucleotides) to be considered. Amplicons above this limit will automatically be predicted as negative.
 
 - cores - the number of CPU cores to be used.
 
DEPENDENCIES:
 - python v2.7.15
 - pandas v0.24.1
 - numpy v1.15.4
 - BioPython v1.72
 - OligoArrayAux (UNAFOLD)
