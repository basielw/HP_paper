#!/usr/bin/python

#import modules
import os,sys
import subprocess
from optparse import OptionParser
from random import randrange
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna

class bcolors:
	HEADER = '\033[95m'
	FAIL = '\033[91m'
	WARNING = '\033[93m'
	OKGREEN = '\033[92m'
	ENDC = '\033[0m'

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage)
	parser.add_option("-c","--cgene",action="store",dest="coding_file",help="Protein coding transcripts (used to build logit model) should be mRNA sequences in FASTA format. The input FASTA file could be regular text file or compressed file (*.gz, *.bz2) or accessible url. NOTE: transcript ID should be unique.")
	parser.add_option("-n","--ngene",action="store",dest="noncoding_file",help="Non protein coding transcripts (used to build logit model) either in BED format or mRNA sequences in FASTA format: If this is BED format file, '-r' must be specified; if this is mRNA sequence file in FASTA format, ignore the '-r' option. The input BED or FASTA file could be regular text file or compressed file (*.gz, *.bz2) or accessible url.  NOTE: transcript ID should be unique.")
	parser.add_option("-o","--outdir",action="store",dest="out_dir",help="Directory where output should be saved. You are safe to remove everything but FINAL_OUTPUT if you wish to do so. FINAL_OUTPUT: Tab separated text file: geneID <tab> mRNA size <tab> ORF sequence <tab> ORF position <tab> ORF size <tab> Fickett Score <tab> Hexamer Score<tab>Coding Probability.")
	parser.add_option("-x","--hex",action="store",dest="hexamer_dat",help="Prebuilt hexamer frequency table (Human, Mouse, Fly, Zebrafish). Run 'make_hexamer_tab.py' to make this table out of your own training dataset.")
	parser.add_option("-s","--start",action="store",dest="start_codons",default='ATG',help="Start codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). default=%default")
	parser.add_option("-t","--stop",action="store",dest="stop_codons",default='TAG,TAA,TGA',help="Stop codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). Multiple stop codons should be separated by ','. default=%default")
	parser.add_option("-f","--folds",action="store",dest="folds",default=5,help="Amount of folds that you wish to do the cross validation with")
	parser.add_option("-l","--min_length",action="store",dest="min_length",default=25,help="Minimum length of ORF's to be found by CPAT. When set to a number, only ORFs that consist of this amount of nucleotides or more will be withheld. default=%default")
	parser.add_option("-k","--kind",action="store",dest="kind",default="ALL",help="Should be either ALL or LONGEST. User can define if they only want the longest ORF or every ORF to be analyzed. When set to LONGEST, minimum length '-l' can be ignored. default=%default")

	(options,args)=parser.parse_args()

	"""CHECKING ALL KINDS OF THINGS"""

	folds = int(options.folds)

	#check input files
	for file in ([options.coding_file,options.noncoding_file,options.hexamer_dat]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	#check if output directory was supplied
	if not options.out_dir:
		print >>sys.stderr, bcolors.FAIL + "\nError: no output directory defined\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	#check if output directory ends with / and if not, add it
	if options.out_dir[-1] != "/":
		options.out_dir = options.out_dir + "/"

	#check if kind was defined correctly
	if (str(options.kind) != "ALL") and (str(options.kind) != "LONGEST"):
		print >>sys.stderr, bcolors.FAIL + "\nError: '-k' should be ALL or LONGEST\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	#make LogitModel directory within output directory
	if not os.path.exists(str(options.out_dir) + "LogitModels"):
		os.makedirs(str(options.out_dir) + "LogitModels")

	#make CPAT_output directory within output directory
	if not os.path.exists(str(options.out_dir) + "CPAT_output"):
		os.makedirs(str(options.out_dir) + "CPAT_output")

	#make sets directory within output directory
	if not os.path.exists(str(options.out_dir) + "sets"):
		os.makedirs(str(options.out_dir) + "sets")

	"""MAGIC HAPPENS BELOW"""

	#making lists
	coding_list = []
	coding_sequences = SeqIO.parse(open(options.coding_file),'fasta')
	for fasta in coding_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		coding_list.append([randrange(1, folds+1), name, sequence])

	noncoding_list = []
	noncoding_sequences = SeqIO.parse(open(options.noncoding_file),'fasta')
	for fasta in noncoding_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		noncoding_list.append([randrange(1, folds+1), name, sequence])

	outfile = open(options.out_dir + "FINAL_OUTPUT", 'w')
	outfile.write("geneID\tmRNA_size\tORF_starting_position\tORF_size\tORF_sequence\tFickett_score\tHexamer_score\tCPAT_coding_prob\tactual\n")
	count = 0
	for i in range(1, folds+1):
		
		#MAKING FASTA FILES

		trainingset_coding = []
		trainingset_noncoding = []
		testset = []
		for entry in coding_list:
			if entry[0] == i:
				my_seq = SeqRecord(Seq(entry[2],unambiguous_dna), id = "CODING___" + entry[1])
				testset.append(my_seq)
			else:
				my_seq = SeqRecord(Seq(entry[2],unambiguous_dna), id = entry[1])
				trainingset_coding.append(my_seq)
		for entry in noncoding_list:
			if entry[0] == i:
				my_seq = SeqRecord(Seq(entry[2],unambiguous_dna), id = "NONCODING___" + entry[1])
				testset.append(my_seq)
			else:
				my_seq = SeqRecord(Seq(entry[2],unambiguous_dna), id = entry[1])
				trainingset_noncoding.append(my_seq)

		trainingset_coding_dir = str(options.out_dir) + "sets/trainingset_coding_" + str(i) + ".fasta"
		trainingset_coding_file = open(trainingset_coding_dir, 'w')
		SeqIO.write(trainingset_coding, trainingset_coding_file, "fasta")
		trainingset_coding_file.close()

		trainingset_noncoding_dir = str(options.out_dir) + "sets/trainingset_noncoding_" + str(i) + ".fasta"
		trainingset_noncoding_file = open(trainingset_noncoding_dir, 'w')
		SeqIO.write(trainingset_noncoding, trainingset_noncoding_file, "fasta")
		trainingset_noncoding_file.close()

		testset_dir = str(options.out_dir) + "sets/testset_" + str(i) + ".fasta"
		testset_file = open(testset_dir, 'w')
		SeqIO.write(testset, testset_file, "fasta")
		testset_file.close()

		#CARRYING OUT CPAT

		print >>sys.stderr, bcolors.HEADER + "\n\nMaking logit model number " + str(i) + "\n" + bcolors.ENDC
		logit_model_dir = str(options.out_dir) + "LogitModels/logitmodel_" + str(i)
		subprocess.call(["/home/basielw/programs/CPAT/bin/make_logitModel_modified.py", "-c", trainingset_coding_dir, "-n", trainingset_noncoding_dir, "-o", logit_model_dir, "-x", options.hexamer_dat, "-s", options.start_codons, "-t", options.stop_codons, "-l", str(options.min_length), "-k", options.kind])

		print >>sys.stderr, bcolors.HEADER + "\nRunning CPAT number " + str(i) + "\n" + bcolors.ENDC
		logit_model_dir_rdata = logit_model_dir + ".logit.RData"
		cpat_output_dir = str(options.out_dir) + "CPAT_output/output_" + str(i)
		subprocess.call(["/home/basielw/programs/CPAT/bin/cpat_modified.py", "-g", testset_dir, "-o", cpat_output_dir, "-x", options.hexamer_dat, "-d", logit_model_dir_rdata, "-s", options.start_codons, "-t", options.stop_codons, "-l", str(options.min_length), "-k", options.kind])

		#WRITING ALL TO FINAL_OUTPUT

		with open(options.out_dir + "CPAT_output/output_" + str(i)) as infile:
			for i, line in enumerate(infile):
				if i == 0:
					continue
				else:
					count += 1
					fields = line.split('\t')
					names = fields[0].split('___')
					coding_prob = fields[-1].rstrip("\n")
					outlist = [names[1]]
					outlist.extend(fields[1:-1])
					outlist.append(coding_prob)
					outlist.append(names[0])
					outline = "\t".join(outlist) + "\n"
					outfile.write(outline)

	print >>sys.stderr, "\n\nWrote %s lines to FINAL_OUTPUT" % (str(count))
	outfile.close()
	print >>sys.stderr, bcolors.OKGREEN + "Done!\n\n" + bcolors.ENDC,

if __name__ == '__main__':
	main()
