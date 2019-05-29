'''
Dominic Byrne - May 2019
Animal and Plant Health Agency

This tool was created in an attempt to predict the results of
diagnostic RT-qPCRs based on the degree of primer/probe - target
complementarity between given influenza genome sequences and RT-qPCRs.

System calls to UNAFOLD are used to predict the free energy of hybridisation
between a given primer/probe sequence and the target. These values
are normalised to the binding energy associated with a perfect match for 
each primer/probe sequence. These binding index values are compared to thresholds
(set by analysing the binding of RT-qPCR reactions with 
experimentally-determined results) to predict the outcome of each RT-qPCR.

The main output is a .csv file with the predicted results (P or N) for each
RT-qPCR - virus combination. The actual binding index values are saved
to another .csv file named 'dg_vals.csv' within the temporary directory
created during the analysis, however this is deleted after each run by default.
'''

import sys, os, subprocess, time, re
import pandas as pd
import numpy as np
from Bio import Seq
from itertools import product
from multiprocessing.dummy import Pool as ThreadPool

args = sys.argv[1]
for line in open(args).readlines():
	exec(line.strip())

def rev_comp(sequence):
	
	complement = {"A":"T","C":"G","G":"C","T":"A",
				  "R":"Y","Y":"R","M":"K","K":"M",
				  "S":"S","W":"W","H":"D","B":"V",
				  "V":"B","D":"H","N":"N","-":"-"}
	rev_seq = sequence[::-1]
	rev_comp_seq = "".join([complement[base] for base in rev_seq])
	
	return rev_comp_seq


def disambiguate(sequence):
	
	d = Seq.IUPAC.IUPACData.ambiguous_dna_values
	d.update({"-":"-"})
	seq_list = list(map("".join, product(*map(d.get, sequence))))
	
	return seq_list


def get_conditions(conditions_file):
	
	cond = pd.read_csv(pcr_conditions_file_path)
	cond.set_index("PCR Name", inplace=True)
	
	return cond
	

def correct_primpro_orientation(cond_frame):
	
	for pcr in cond_frame.index:
		rev_seq = cond_frame.loc[pcr,"Rev primer seq"]
		cond_frame.loc[pcr,"Rev primer seq"] = rev_comp(rev_seq)
		
		if cond_frame.loc[pcr,"Probe target strand"] == "+":
			pro_seq = cond_frame.loc[pcr,"Probe seq"]
			cond_frame.loc[pcr,"Probe seq"] = rev_comp(pro_seq)


def parse_genomes_file(genomes_file, cond_frame):
	
	segments_wanted=sorted([str(i) for i in list(cond_frame["Target segment number"].unique())])
	data = open(genomes_file).read()
	seqs = [seq for seq in data.split(">") if len(seq.strip()) > 0]
	seq_records = {seq.split("\n")[0] : "".join(seq.split("\n")[1:]).upper() 
				   for seq in seqs}
	seq_df = pd.DataFrame()
	
	for rec in seq_records.keys():
		seg_name = rec.split("|")[1].strip()
		seg_id = "Seg" + seg_name
		strain = rec.split("|")[0].strip().replace("/", "_").replace(" ", "_")
		seq = seq_records[rec].strip().replace("\r", "")
		
		if seg_name in segments_wanted:
			seq_df.loc[strain, seg_id] = seq

	return seq_df
	


def make_temp_seq_files(gen_frame, cond_frame, temp_dir_path):
	
	strains = gen_frame.index
	segments = gen_frame.columns
	pcrs = cond_frame.index
	
	file_dict = {}
	indi_segs = [[strain,seg] for strain in strains for seg in segments]
	for [strain, seg] in indi_segs:
		header = ">" + strain + "_" + seg
		fname = temp_dir_path + "/" + strain + "_" + seg + ".fa"
		file_dict.update({(strain, seg): fname})
		seq = gen_frame.loc[strain, seg]
		if type(seq) == str:
			with open(fname, "w") as fh:
				fh.write("\n".join([header, seq]))
	
	for pcr in pcrs:
		fwd_seq = cond_frame.loc[pcr, "Fwd primer seq"]
		rev_seq = cond_frame.loc[pcr, "Rev primer seq"]
		pro_seq = cond_frame.loc[pcr, "Probe seq"]
		seqs = {"fwd": fwd_seq, "rev": rev_seq, "pro": pro_seq}
		for i in seqs.keys():
			fname = temp_dir_path + "/" + pcr + "_" + i + ".fa"
			file_dict.update({(pcr, i): fname})
			header = ">" + pcr + "_" + i
			with open(fname, "w") as fh:
				fh.write("\n".join([header, seqs[i]]))
	
	return file_dict


def align_seqs(gene_file, pcr_file):
	
	gene_id = ".".join(gene_file.split(".")[:-1])
	pcr_id = ".".join(pcr_file.split("/")[-1].split(".")[:-1])
	outfile = gene_id + "_" + pcr_id + ".aln"
	subprocess.check_output("water -asequence %s -bsequence %s -gapopen 1 -gapextend 10 -outfile %s -aformat markx10" 
			  % (gene_file, pcr_file, outfile), stderr = subprocess.STDOUT, shell=True)
	
	return outfile


def get_binding_site(alignment_file):
	
	aln = open(alignment_file).read()
	start = int(re.findall(r'al_start: (\d+)', aln)[0]) - 1
	stop = int(re.findall(r'al_stop: (\d+)', aln)[0]) - 1
	
	return [start, stop]


def check_binding_positions(binding_sites):
	
	all_pos = sum(binding_sites, [])
	uniq_sorted_all_pos = sorted(set(all_pos))
	
	return all_pos == uniq_sorted_all_pos


def disambiguate_amplicons(gen_frame, cond_frame, file_dict, temp_dir_path):
	
	strains = gen_frame.index
	segments = gen_frame.columns
	pcrs = cond_frame.index
	target_segs = {pcr: "Seg" + str(cond_frame.loc[pcr, "Target segment number"]) for pcr in pcrs}
	combos = [[strain, seg, pcr]
			   for strain in strains
			   for seg in segments
			   for pcr in pcrs if seg == target_segs[pcr]]
	pcr_components = ["fwd", "pro", "rev"]

	amp_dict = {pcr: {} for pcr in pcrs}
	for strain, seg, pcr in combos:
		gene_file = file_dict[strain, seg]
		if not os.path.isfile(gene_file):
			continue
		pcr_files = [file_dict[pcr, comp] for comp in pcr_components]
		alignments = [align_seqs(gene_file, pcr_file)
					  for pcr_file in pcr_files]
		binding_sites = [get_binding_site(aln) for aln in alignments]
		start = min(binding_sites)[0]
		end = max(binding_sites)[1]
		gene = gen_frame.loc[strain, seg]
		amplicon = gene[start : end + 1]	
		amp_id = "_".join([strain, seg])
		
		if len(amplicon) <= max_amp_size:
			dis_amps = disambiguate(amplicon.replace("N", "-"))
			amp_dict[pcr].update({amp_id : dis_amps})
	
	return amp_dict


def assess_pcr_binding(amp_dict, cond_frame, temp_dir_path, n_cores, empty_out_frame):
	
	pcrs = amp_dict.keys()
	pcr_comps = ["Fwd primer seq", "Probe seq", "Rev primer seq"]
	for pcr in pcrs:
		temp = str(cond_frame.loc[pcr, "Annealing temp"])
		Mg = cond_frame.loc[pcr, "Mg conc"]
		if not np.isnan(Mg):
			melt_args = "-M " + str(Mg) + " "
		else:
			melt_args = ""
		melt_args += "-t %s -C 1 -n DNA --allpairs" % temp 
		
		for comp in pcr_comps:
			comp_id = comp[:3].lower()
			comp_seq = rev_comp(cond_frame.loc[pcr, comp])
			len_comp_seq = len(comp_seq)
			
			dis_comp_seqs = disambiguate(comp_seq)
			r_dis_comp_seqs = [rev_comp(seq) for seq in dis_comp_seqs]
			
			pcr_seqs = zip(dis_comp_seqs, r_dis_comp_seqs)
			arg_tups = [(temp_dir_path, pcr, melt_args, comp_id, pcr_seqs, amp_dict, amp_id, empty_out_frame)
					for amp_id in amp_dict[pcr].keys()]
			
			
			os.chdir(temp_dir_path)
			pool = ThreadPool(cores)
			pool.map(get_dg_vals, arg_tups)
			os.chdir("..")
			
			
	return empty_out_frame


def get_dg_vals(args):
	
	temp_dir_path = args[0]
	pcr = args[1]
	melt_args = args[2]
	comp_id = args[3]
	pcr_seqs = args[4]
	amp_dict = args[5]
	amp_id = args[6]
	empty_out_frame = args[7]
	
	fname_stem = temp_dir_path + "/" + "_".join([pcr, comp_id, amp_id])
	comp_fname = fname_stem + "_1.fa"
	amp_fname = fname_stem + "_2.fa"
	r_comp_fname = fname_stem + "_3.fa"

	
	with open(comp_fname, "w") as comp_fh, \
		 open(amp_fname, "w") as amp_fh, \
		 open(r_comp_fname, "w") as r_comp_fh:
		
		for amp_i, amp_seq in enumerate(amp_dict[pcr][amp_id]):
			amp_header = ">" + "_".join([amp_id, str(amp_i + 1)])
			for comp_i, seqs in enumerate(pcr_seqs):
				comp_seq = seqs[0]
				r_comp_seq = seqs[1]
				
				comp_header = ">" + "_".join([pcr, comp_id, str(comp_i + 1)])
				r_comp_header = comp_header + "_rev_comp"
				
				amp_fh.write(amp_header + "\n")
				amp_fh.write(amp_seq + "\n")
				
				comp_fh.write(comp_header + "\n")
				comp_fh.write(comp_seq + "\n")
				
				r_comp_fh.write(r_comp_header + "\n")
				r_comp_fh.write(r_comp_seq + "\n")
				
	
	file1 = comp_fname.split("/")[-1]
	file2 = amp_fname.split("/")[-1]
	file3 = r_comp_fname.split("/")[-1]
	
	target_out = subprocess.check_output(" ".join(["melt.pl", melt_args, file1, file2]), shell=True).split("\n")
	max_out = subprocess.check_output(" ".join(["melt.pl", melt_args, file1, file3]), shell=True).split("\n")
	
	dG_vals = [float(val) for val in [line.split("\t")[0].strip() for line in target_out if "\t" in line][1:]]
	max_dG_vals = [float(val) for val in [line.split("\t")[0].strip() for line in max_out if "\t" in line][1:]]
	
	dG_final = [round(vals[0]/vals[1], 3) for vals in zip(dG_vals, max_dG_vals)]
	
	strain = "_".join(amp_id.split("_")[:-1])
	pcr_comp = pcr + "_" + comp_id
	empty_out_frame.loc[strain, pcr_comp] = max(dG_final)			



conditions = get_conditions(pcr_conditions_file_path)
genomes = parse_genomes_file(genomes_file_path, conditions)
correct_primpro_orientation(conditions)

cwd = os.getcwd()
temp_dir = cwd + "/temp_pcr_dir_" + str(time.time())
os.mkdir(temp_dir)

seq_files = make_temp_seq_files(genomes, conditions, temp_dir)
amplicons =  disambiguate_amplicons(genomes, conditions, seq_files, temp_dir)
comps = ["fwd", "rev", "pro"]
cols = [pcr + "_" + comp for pcr in conditions.index for comp in comps]
empty_frame = pd.DataFrame(index = genomes.index, columns = cols)

dg_vals = assess_pcr_binding(amplicons, conditions, temp_dir, cores, empty_frame).fillna(0)
dg_vals.to_csv(temp_dir + "/dg_vals.csv")

results = pd.DataFrame(index = genomes.index, columns = conditions.index)
positives = {pcr : [] for pcr in conditions.index}
negatives = {pcr : [] for pcr in conditions.index}

'''Outputs'''
for strain in dg_vals.index:
	for pcr in conditions.index:
		fwd_val = dg_vals.loc[strain, pcr + "_fwd"]
		rev_val = dg_vals.loc[strain, pcr + "_rev"]
		pro_val = dg_vals.loc[strain, pcr + "_pro"]
		
		prim = np.mean([fwd_val, rev_val])
		
		if prim >= 0.7 and pro_val >= 0.7:
			results.loc[strain, pcr] = "P"
			positives[pcr].append(strain)
		else:
			results.loc[strain, pcr] = "N"
			negatives[pcr].append(strain)

for pcr in conditions.index:
	print "%s predicted positives: %d" % (pcr, len(positives[pcr]))
	print "%s predicted negatives: %d" % (pcr, len(negatives[pcr]))

results.to_csv(output_file)
#os.system("rm -r %s" % temp_dir)
