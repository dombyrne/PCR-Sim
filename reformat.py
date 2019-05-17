'''
Reformats FluDB-format FASTA files to the standard GISAID format for 
use with pcr_sim.py.
'''

import sys

input_fname = sys.argv[1]
output_fname = sys.argv[2]

with open(output_fname, "w") as fh:
	for line in open(input_fname).readlines():
		if line.startswith(">"):
			f = line.split("|")
			for i in f:
				if "Strain Name:" in i:
					strain_name = i.split(":")[-1].strip()
				if "Segment:" in i:
					segment_number = i.split(":")[-1].strip()
			new_header = ">" + strain_name + " | " + segment_number + "\n"
			fh.write(new_header)
		else:
			fh.write(line)
