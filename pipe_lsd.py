#!/usr/bin/env python3
# Create input for lsd tool.
# 1-based coordinates
# Author: Gennady Khvorykh, info@inzilico.com 

import sys
import re
import os
import argparse
import func
from multiprocessing import Pool
import subprocess

# Initiate
parser = argparse.ArgumentParser(description = 'Create input files for lsd tool')
parser.add_argument('--cpu', type = int, help = 'The number of CPU to be used (default: 12)', default = 12)
parser.add_argument('--out', type = str, help = 'path/to/folder to save output files', default = 'out')
parser.add_argument('vcf1', type = str, help = 'path/to/file.vcf.gz with population variants')
parser.add_argument('vcf2', type = str, help = 'path/to/file.vcf.gz with haploid ancestor variants')
parser.add_argument('gff', type = str, help = 'path/to/file.gff with the list of windows')
parser.add_argument('pops', type = str, help = 'path/to/file.txt with the population id')
args = parser.parse_args()
resources = "res.cfg"
res = func.read_resources(resources)

vcf1 = args.vcf1
vcf2 = args.vcf2
gff = args.gff
pops = args.pops
out = args.out
cpu = args.cpu

def check_file(x):
	if not os.path.isfile(x):
		print(f'{x} doesn\'t exist!')
		sys.exit(1)

def read_gff(gff):
	out = []
	print(f'Loading {gff}...')
	f = open(gff, 'r')
	for l in f:
		l = l.rstrip('\r\n')
		x = l.split()
		try:
			id = re.search('ID=(.*)', x[8]).group(1)
		except AttributeError:
			print('ID attribute is not provided in gff file')
			sys.exit(1)
		out.append([x[0], x[3], x[4], id])
	f.close()
	size = len(out)
	print(f'Regions read: {size}')
	return(out)

def make_input_for_ancestor_sequences(r):
	reg = f'{r[0]}:{r[1]}-{r[2]}'
	# Split reference fasta into regions
	output1 = f'{out}/{r[3]}.ref.fa'
	subprocess.run([res['samtools'], 'faidx', res['ref'], reg, '-o', output1], check = True)
	# Split haploid ancestor variants into regions
	output2 = f'{out}/{r[3]}.anc.vcf'
	subprocess.run([res['bcftools'], 'view', vcf2, '--regions', reg, '-Ov', '-o', output2], check = True)

def conc_fasta(r):
	win = r[3]
	f1 = f'{out}_win/{win}.fas'
	f2 = f'{out}/{win}.anc.fa'
	# Check input
	for fn in [f1, f2]: check_file(fn)
	output = f'{out}/{win}.fa'  
	f = open(output, 'w')
	subprocess.run(['cat', f1, f2], check = True, stdout = f)
	f.close()

def build_trees(regions):
	for r in regions:
		win = r[3]
		fa = f'{out}/{win}.fa'
		check_file(fa)
		output = f'{out}/{win}.nwk'
		subprocess.run([res['ft'], '-gtr', '-nt', '-nosupport', '-quiet', '-out', output, fa], check = True)

def get_size(x):
	p = subprocess.Popen([res['bcftools'], 'query', '-l', x], stdout = subprocess.PIPE)
	s = subprocess.run(['wc', '-l'], stdin = p.stdout, capture_output = True, text = True).stdout.strip("\n")
	p.wait()
	return(int(s))

def lsd(out, n):
	# Step in the folder with win*.nwk files
	cwd = os.getcwd()
	os.chdir(out)

	# Concatenate all trees into one file
	fout = "trees.nwk"
	if os.path.exists(fout):
		os.remove(fout)
	cmd = f'ls win*.nwk | sort -V | xargs cat > {fout}'
	os.system(cmd)

	# Correct ids
	fin = fout
	fout = "trees.cor.nwk"
	if os.path.exists(fout):
		os.remove(fout)
	cmd = f'sed \'s/_0/A/g\' {fin} | sed \'s/_1/B/g\' > {fout}'
	os.system(cmd)
	
	# run LSD
	fin = fout
	fout = "data"
	size = str(2*n + 1) 
	subprocess.run([res['lsd'], '-t', fin, '-n', size, '-r', 'ANC', '-g', pops, '-o', fout], check = True)
	os.chdir(cwd)

def main():
	# Check input
	for fn in [vcf1, vcf2, gff, pops]: check_file(fn)

	# Get the number of samples in the population data
	n = get_size(vcf1)
	
	# Show input
	print(f'Population data: {vcf1}')
	print(f'Ancestral data: {vcf2}')
	print(f'gff: {gff}')
	print(f'Populations: {pops}')
	print(f'Sample size: {n}')
	print(f'Output: {out}')
	print(f'CPU: {cpu}')

	# Create output folder
	if not os.path.exists(out):
		os.makedirs(out)

	# Read gff
	reg = read_gff(gff)

	# Make fasta and vcf to create ancestor sequences for each region
	with Pool(cpu) as p:
		p.map(make_input_for_ancestor_sequences, reg)
	# Create ancestor sequences
	subprocess.run(['Rscript', res['vcf2fasta'], out, gff, str(cpu)], check = True)
	# Create population sequences
	subprocess.run(
		[res['vcf2fa'], '--fasta', res['ref'], '--vcf', vcf1, '--gff', gff, '--out', out, '--feat', "win"], 
		check = True)
	# Concatenate population and ancestral sequences
	with Pool(cpu) as p:
		p.map(conc_fasta, reg)
	# Build trees
	build_trees(reg)

	# Estimate LSD
	lsd(out, n)

if __name__ == "__main__":
	main()
