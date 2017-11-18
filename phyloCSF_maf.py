#!/usr/bin/python

'''-----------------------------------
Tool to parse BED files into MAF files
-----------------------------------'''

#import modules
import os,sys
import subprocess
from optparse import OptionParser
from tempfile import mkstemp
from shutil import move

class bcolors:
    FAIL = '\033[91m'
    WARNING = '\033[93m'
    OKGREEN = '\033[92m'
    ENDC = '\033[0m'

def main():
	usage = "\n%prog  [options]"
	desc = "\nThis script takes one BED file and outputs a FASTA file for every gene in the BED file. These FASTA files compare the same gene in different organisms (29mammals/46way). It retrieves this information from MAF-files and will thus create a MAF file for every gene first, which will then be transformed into FASTA. These FASTA files is what PhyloCSF needs as an input. The process takes about 2 hours for 24 chromosomes. The time it takes is independent of the size of your BED-file, because it has to read the whole MAF-files anyway. You will need to provide MAF-files for every chromosome you would like to analyse. MAF-files for hg19 can be downloaded from here: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz46way/\n"
	parser = OptionParser(usage,description=desc)
	parser.add_option("-b","--BED_file",action="store",dest="bed_file",help="BED file that needs to be processed")
	parser.add_option("-m","--MAF_file_directory",action="store",dest="maf_files",help="Directory where all the corresponding MAF files are located.")
	parser.add_option("-o","--output_directory",action="store",dest="output_dir",default="./",help="Directory where everything should be created. Default = current directory")
	parser.add_option("-i","--fasta_id",action="store",dest="change_fasta_id",help="Modify the fasta-IDs to match a .nh file. This will only work if you downloaded the MAF-files here: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz46way/ . If not, see the '-r' argument. Included .nh files: 29mammals (PhyloCSF), 46way. When using 29mammals, a lot of the animals in the FASTA-files will be skipped because they are not a part of 29mammals. This is completely normal. If you would like to specify your own fasta-ID's, you will have to add a custom file with the ID's by adding the -r argument followed by a path to the file. You have to pick either the -i argument or the -r argument.")
	parser.add_option("-r","--replacement_file",action="store",dest="replacement_dir",help="Add the path to the custom replacements file. The file should contain one line per organism, with each line having this format: 'code organism' for example 'hg19 human' (no parentheses). The code cannot contain spaces, however the organism name can if you wish. When there are animals in the FASTA-files that are not in your replacement file, these animals will be skipped. You have to pick either the -i argument or the -r argument.")

	(options,args)=parser.parse_args()

	"""DEFINE DEFAULT REPLACEMENT DICTS"""
	fourtysixway = {
		'hg19':'Human',
		'GRCh37':'Human',
		'panTro2':'Chimp',
		'gorGor1':'Gorilla',
		'ponAbe2':'Orangutan',
		'rheMac2':'Rhesus',
		'papHam1':'Baboon',
		'calJac1':'Marmoset',
		'tarSyr1':'Tarsier',
		'micMur1':'Mouse_lemur',
		'otoGar1':'Bushbaby',
		'tupBel1':'TreeShrew',
		'mm9':'Mouse',
		'rn4':'Rat',
		'dipOrd1':'Kangaroo_rat',	
		'cavPor3':'Guinea_Pig',
		'speTri1':'Squirrel',
		'oryCun2':'Rabbit',
		'ochPri2':'Pika',
		'vicPac1':'Alpaca',
		'turTru1':'Dolphin',
		'bosTau4':'Cow',
		'equCab2':'Horse',
		'felCat3':'Cat',
		'canFam2':'Dog',
		'myoLuc1':'Microbat',
		'pteVam1':'Megabat',
		'eriEur1':'Hedgehog',
		'sorAra1':'Shrew',
		'loxAfr3':'Elephant',
		'proCap1':'Rock_hyrax',
		'echTel1':'Tenrec',
		'dasNov2':'Armadillo',
		'choHof1':'Sloth',
		'macEug1':'Wallaby',
		'monDom5':'Opossum',
		'ornAna1':'Platypus',
		'galGal3':'Chicken',
		'taeGut1':'Zebra_finch',
		'anoCar1':'Lizard',
		'xenTro2':'X_tropicalis',
		'tetNig2':'Tetraodon',
		'fr2':'Fugu',
		'gasAcu1':'Stickleback',
		'oryLat2':'Medaka',
		'danRer6':'Zebrafish',
		'petMar1':'Lamprey',
		}

	twentyninemammals = {
		'hg19':'Human',
		'GRCh37':'Human',
		'panTro2':'Chimp',
		'rheMac2':'Rhesus',
		'tarSyr1':'Tarsier',
		'micMur1':'Mouse_lemur',
		'otoGar1':'Bushbaby',
		'tupBel1':'TreeShrew',
		'mm9':'Mouse',
		'rn4':'Rat',
		'dipOrd1':'Kangaroo_rat',	
		'cavPor3':'Guinea_Pig',
		'speTri1':'Squirrel',
		'oryCun2':'Rabbit',
		'ochPri2':'Pika',
		'vicPac1':'Alpaca',
		'turTru1':'Dolphin',
		'bosTau4':'Cow',
		'equCab2':'Horse',
		'felCat3':'Cat',
		'canFam2':'Dog',
		'myoLuc1':'Microbat',
		'pteVam1':'Megabat',
		'eriEur1':'Hedgehog',
		'sorAra1':'Shrew',
		'loxAfr3':'Elephant',
		'proCap1':'Rock_hyrax',
		'echTel1':'Tenrec',
		'dasNov2':'Armadillo',
		'choHof1':'Sloth',
		}

	"""A LOT OF CHECKING HAPPENING BELOW"""
	#check if bed_file was given
	if not options.bed_file:
		print >>sys.stderr, bcolors.FAIL + "\nError: no BED file defined\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	#check if maf_files was supplied
	if not options.maf_files:
		print >>sys.stderr, bcolors.FAIL + "\nError: no MAF files directory defined\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	#check if maf_files directory ends with / and if not, add it
	if options.maf_files[-1] != "/":
		options.maf_files = options.maf_files + "/"

	#check if output_directory was given and if not, use current directory
	if not options.output_dir:
		options.output_dir = "./"

	#check if output directory ends with / and if not, add it
	if options.output_dir[-1] != "/":
		options.output_dir = options.output_dir + "/"

	#check if output directory exists and if not, make it
	if not os.path.exists(options.output_dir):
		os.makedirs(options.output_dir)

	#if both -i and -r are supplied, an error should be given
	if options.change_fasta_id and options.replacement_dir:
		print >>sys.stderr, bcolors.FAIL + "\nError: both the '-i' argument and the '-r' argument were supplied. You have to pick one.\n" + bcolors.ENDC
		parser.print_help()
		sys.exit(0)

	#check if custom replacements file was supplied and if so, make a dict out of it
	if options.replacement_dir:
		replacements = {}
		replacement_file = open(options.replacement_dir)
		for line in replacement_file:
			fields = line.split()
			replacements[fields[0]] = ' '.join(fields[1:])
	#check if one of the default dicts should be used
	if options.change_fasta_id:
		if options.change_fasta_id == '46way':
			replacements = fourtysixway
		elif options.change_fasta_id == '29mammals':
			replacements = twentyninemammals
		else:
			print >>sys.stderr, bcolors.FAIL + "\nError: the '-i' argument has to be either '29mammals' or '46way'.\n" + bcolors.ENDC
			parser.print_help()
			sys.exit(0)

	#make MAF-files directory within output directory
	maf_output_dir = "%sMAF_files/" % options.output_dir
	if not os.path.exists(maf_output_dir):
		os.makedirs(maf_output_dir)

	#make FASTA-files directory within output directory
	fasta_output_dir = "%sFASTA_files/" % options.output_dir
	if not os.path.exists(fasta_output_dir):
		os.makedirs(fasta_output_dir)

	"""SPLITTING BED PER CHROMOSOME"""
	#get name of the BED file
	bed_name = os.path.basename(options.bed_file)
	bed_name = bed_name[:-4]

	#this is where the magic happens
	count = 0
	file_list = {}
	print >>sys.stderr, "\nSplitting BED...\n",
	bed_file = open(options.bed_file)
	for line in bed_file:
		count +=1
		fields = line.split()
		file_name = str(options.output_dir + str(bed_name) + '_' + str(fields[0]))

		#check if chromosome file is already opened and if not, open it
		#when the file doesn't exist yet, it will be created
		#when the file already exists, it will be cleared
		if file_name not in file_list:
			file = open('%s.bed' % file_name, 'w')
			file_list[file_name] = file
		
		#write lines to the correct chromosome file
		file_list[file_name].write(line)

		#show bed file count to not fall asleep
		print >>sys.stderr, "%d genes in the BED file splitted per chromosome\r" % count,
	print >>sys.stderr, bcolors.OKGREEN + "\nDone!\n" + bcolors.ENDC,

	#close all the files
	for key, value in file_list.iteritems():
		value.close()
	bed_file.close()

	"""GENERATING MAF FILES FROM BED"""
	count = 0
	gene_dict = {}
	print >>sys.stderr, "\nGenerating MAF files...\n",
	for bed_file in os.listdir(options.output_dir):
		if bed_file[-4:] == '.bed':
			#get chromosome name for MAF file
			chromosome = bed_file[-9:-4]
			if chromosome[0] != 'c':
				chromosome = chromosome[1:]
			if not os.path.isfile(options.maf_files + chromosome + '.maf'):
				print >>sys.stderr, bcolors.WARNING + "Warning: %s was skipped because no MAF file was found for this chromosome." % bed_file + bcolors.ENDC
				continue #if there is no MAF file for the chromosome, skip this BED file
			if chromosome[0:3] != 'chr':
				print >>sys.stderr, bcolors.FAIL + "\nError: BED files have incorrect names\n" + bcolors.ENDC + "\nThey should have this format: name_chr..\nFor example: whatever_chr1\n"
				sys.exit(0)

			count += 1
			
			#make folder for each chromosome
			maf_output_dir_chr = maf_output_dir + chromosome + '/'
			if not os.path.exists(maf_output_dir_chr):
				os.makedirs(maf_output_dir_chr)

			#get the full bed file directory
			bed_file_dir = options.output_dir + bed_file

			#generate dict with strand per gene
			bed_file_open = open(bed_file_dir)
			for line in bed_file_open:
				fields = line.split('\t')
				gene_dict[fields[3]] = fields[5]
			bed_file_open.close()

			#get the corresponding full maf file directory
			maf_file_dir = options.maf_files + chromosome + '.maf'

			#use maf_parse on the bed files
			subprocess.call(["maf_parse", "-g", bed_file_dir, "-P", "id", maf_file_dir], cwd=maf_output_dir_chr)

			#show bed file count to not fall asleep
			print >>sys.stderr, "%d BED files per chromosome finished\r" % count,
	print >>sys.stderr, bcolors.OKGREEN + "\nDone!\n" + bcolors.ENDC,

	"""GENERATING FASTA FILES FROM MAF"""
	count = 0
	print >>sys.stderr, "\nGenerating FASTA files...\n",
	warning_set = set()
	#loop through chromosome folders within MAF output directory
	for subdir, dirs, files in os.walk(maf_output_dir):
		for directory in dirs:
			maf_output_dir_chr = maf_output_dir + directory + '/'

			#actually make FASTA files from MAF files
			for maf_file in os.listdir(maf_output_dir_chr):
				#skip files that are not in gene_dict
				#this script won't create files that are not in gene_dict, obviously
				#however, the directory where it searches for MAF files could already contain files from before
				if maf_file[10:-4] not in gene_dict:
					continue

				count += 1
				maf_file_dir = maf_output_dir_chr + maf_file
				fasta_file_dir = fasta_output_dir + maf_file[10:-4] + '.fa'
				
				#for positive strands
				if gene_dict[maf_file[10:-4]] == '+':
					some_tmp = open(fasta_file_dir, 'w')
					subprocess.call(["msa_view", "-i", "MAF", "-o", "FASTA", "-m", "-G", "1", maf_file_dir], cwd=fasta_output_dir, stdout=some_tmp)
					some_tmp.close()
					if options.change_fasta_id or options.replacement_dir:
						#create temp file
						fh, abs_path = mkstemp()
						with open(abs_path,'w') as new_file:
							with open(fasta_file_dir) as old_file:
								for line in old_file:
									if line[0] == '>':
										if line[2:-1] in replacements:
											new_file.write('> ' + replacements[line[2:-1]] + '\n')
											flag = 1
										else:
											flag = 0
											warning_set.add(line[2:-1]) #make a set of items that were not found in 'replacements'
									elif flag == 1:
										new_file.write(line)
						os.close(fh)
						#remove original file
						os.remove(fasta_file_dir)
						#move new file
						move(abs_path, fasta_file_dir)
					
				#for negative strands
				elif gene_dict[maf_file[10:-4]] == '-':
					some_tmp = open(fasta_file_dir, 'w')
					subprocess.call(["msa_view", "-i", "MAF", "-o", "FASTA", "-m", "-V", "-G", "1", maf_file_dir], cwd=fasta_output_dir, stdout=some_tmp)
					some_tmp.close()
					if options.change_fasta_id or options.replacement_dir:
						#create temp file
						fh, abs_path = mkstemp()
						with open(abs_path,'w') as new_file:
							with open(fasta_file_dir) as old_file:
								for line in old_file:
									if line[0] == '>':
										if line[2:-1] in replacements:
											new_file.write('> ' + replacements[line[2:-1]] + '\n')
											flag = 1
										else:
											flag = 0
											warning_set.add(line[2:-1]) #make a set of items that were not found in 'replacements'
									elif flag == 1:
										new_file.write(line)
						os.close(fh)
						#remove original file
						os.remove(fasta_file_dir)
						#move new file
						move(abs_path, fasta_file_dir)
				
				#show maf file count to not fall asleep
				print >>sys.stderr, "%d MAF files finished\r" % count,
	print >>sys.stderr, bcolors.OKGREEN + "\nDone!\n" + bcolors.ENDC,

	#if there are items in the warning set, print a warning
	if len(warning_set) != 0:
		warning_list = ', '.join(warning_set)
		if options.replacement_dir:
			name = os.path.basename(options.replacement_dir)
		if options.change_fasta_id:
			name = options.change_fasta_id
		if len(warning_set) == 1:
			isorare = 'is'
			organism = 'this organism'
		if len(warning_set) > 1:
			isorare = 'are'
			organism = 'these organisms'
		print >>sys.stderr, bcolors.WARNING + "Warning: %s %s not in %s, so all the lines from %s in the FASTA-files were skipped" % (warning_list, isorare, name, organism) + bcolors.ENDC
		print >>sys.stderr, "\nThis is a list of all known ID's and replacements found in %s: " % name
		for keys, values in replacements.iteritems():
			print >>sys.stderr, "	" + keys + " : " + values
		print >>sys.stderr, ""

if __name__ == '__main__':
	main()
