#!/usr/bin/python

import os,sys
import random
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna

class bcolors:
	FAIL = '\033[91m'
	WARNING = '\033[93m'
	OKGREEN = '\033[92m'
	ENDC = '\033[0m'

def remove_stop_codons(self):
	tmp_seq = self
	counter = 0
	codons_list = []
	while counter < len(tmp_seq):
		codons_list.append(tmp_seq[counter:counter+3])
		counter += 3
	codons_list_no_start = []
	for i, codon in enumerate(codons_list):
		if (codon == 'TAG' or codon == 'TAA' or codon == 'TGA'):
			codons_list_no_start.append('TGG')
		else:
			codons_list_no_start.append(codon)
	new_seq = ''.join(codons_list_no_start)
	return new_seq

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage)
	parser.add_option("-n","--noncoding_transcripts",action="store",dest="noncoding",help="Path to FASTA file that contains noncoding transcripts.")
	parser.add_option("-c","--coding_orfs",action="store",dest="orfs",help="Path to CPAT output that contains the longest ORFs of coding transcripts")
	parser.add_option("-o","--output",action="store",dest="output",help="Path to output FASTA file (must have extension '.fasta')")

	(options,args)=parser.parse_args()

	if not options.noncoding:
		print >>sys.stderr, bcolors.FAIL + "\nError: no FASTA file containing non-coding transcripts defined\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	if not options.orfs:
		print >>sys.stderr, bcolors.FAIL + "\nError: no CPAT output defined\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	if not options.output:
		print >>sys.stderr, bcolors.FAIL + "\nError: no output file defined\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	fasta_sequences = SeqIO.parse(open(options.noncoding),'fasta')
	sequences = []
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		sequences.append(sequence)

	orfs_file = open(options.orfs)
	fasta_records = []
	count = 0
	for i, line in enumerate(orfs_file):
		if i == 0:
			continue
		else:
			count += 1
			line = line.rstrip()
			fields = line.split('\t')
			orf_length = int(fields[3]) - 6
			mRNA_size = int(fields[1])
			diff = mRNA_size - int(fields[3])
			rand_seq = random.choice(sequences)
			while (len(rand_seq) < orf_length):
				rand_seq = random.choice(sequences)
			seq = rand_seq[:orf_length]
			seq = remove_stop_codons(seq)
			#if mRNA size is long, add stop codons in other frames to make sure it CPAT doens't find a longer ORF in that frame
			if diff >= 9:
				seq_start_stop = 'ATG' + seq + 'TAA' + 'ATAAATAAA' + ('A' * (diff - 9))
			#if mRNA size is short, this doesn't really matter because CPAT won't find longer ORFs there anyway
			else:
				seq_start_stop = 'ATG' + seq + 'TAA' + ('A' * diff)
			record_id = "same_length_as_" + fields[0]
			my_seq = SeqRecord(Seq(seq_start_stop,unambiguous_dna), id = record_id)
			fasta_records.append(my_seq)
			print >>sys.stderr, "%d lines finished\r" % count,
	SeqIO.write(fasta_records, options.output, "fasta")
	orfs_file.close()

if __name__ == '__main__':
	main()
