#!/usr/bin/python

#import modules used here
import sys
import requests
import pymysql.cursors
from time import sleep

class bcolors:
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ENDC = '\033[0m'

infile = open("/Users/basielw/Downloads/results.txt")
outfile = open("/Users/basielw/Downloads/sorfs_org_converted_2", 'w')
outfile.write("name\tchrom\tstrand\tchromstart\tchromend\tORFstart\tORFstop\tORFsize\tORFseq\tORFscore\tcovF1prop\tpep_seq\tinLNCipedia\n")
select_seq_region_start_37 = "SELECT seq_region_start FROM homo_sapiens_core_75_37.transcript WHERE stable_id = %(ensemblid)s"
select_seq_region_end_37 = "SELECT seq_region_end FROM homo_sapiens_core_75_37.transcript WHERE stable_id = %(ensemblid)s"
select_seq_region_start_38 = "SELECT seq_region_start FROM homo_sapiens_core_83_38.transcript WHERE stable_id = %(ensemblid)s"

db_37 = pymysql.connect(host="ensembldb.ensembl.org",
					user="anonymous",
					db="homo_sapiens_core_75_37")
db_38 = pymysql.connect(host="ensembldb.ensembl.org",
					user="anonymous",
					db="homo_sapiens_core_83_38")
cur_37 = db_37.cursor()
cur_38 = db_38.cursor()

count = 0
unknown_count = 0
error_list = []
for i, line in enumerate(infile):
	if i == 0:
		continue
	else:
		error_37 = False
		count += 1
		line = line.rstrip()
		fields = line.split('\t')
		ensemblid = fields[0]
		if ensemblid == "":
			unknown_count += 1
			continue
		#get new chromstart and chromend (hg19)
		cur_37.execute(select_seq_region_start_37, {'ensemblid': ensemblid})
		row = cur_37.fetchone()
		if row == None:
			error_37 = True
		else:
			chromstart_37 = str(row[0])
		cur_37.execute(select_seq_region_end_37, {'ensemblid': ensemblid})
		row = cur_37.fetchone()
		if row == None:
			error_37 = True
		else:
			chromend_37 = str(row[0])
		#get old chromstart to modify ORFstart and ORFstop to match hg19
		cur_38.execute(select_seq_region_start_38, {'ensemblid': ensemblid})
		row = cur_38.fetchone()
		if row == None:
			print >>sys.stderr, bcolors.WARNING + "Warning: line %d: no such transcript found in Ensembl 83: %s, line skipped\n" % (count,ensemblid) + bcolors.ENDC
			continue
		else:
			chromstart_38 = str(row[0])
		#get lncipedia ID if available
		url = "http://lncipedia.org/api/search?id={}".format(ensemblid)
		while True:
			try:
				data = requests.get(url) #403
			except requests.exceptions.ConnectionError:
				print >>sys.stderr, bcolors.WARNING + "Warning: connection error occured, trying again  in 30 seconds\n" + bcolors.ENDC
				sleep(30)
				continue
			break
		#if transcript not in Ensembl 75, get chromstart and end from LNCipedia if possible
		if error_37 == True:
			if int(data.json()['count']) == 0:
				print >>sys.stderr, bcolors.WARNING + "Warning: line %d: no such transcript found in Ensembl 75: %s, line skipped\n" % (count,ensemblid) + bcolors.ENDC
				continue
			elif int(data.json()['count']) == 1:
				chromstart_37 = data.json()['transcripts'][0]['start']
				chromend_37 = data.json()['transcripts'][0]['end']
			elif int(data.json()['count']) > 1:
				error_list.append(ensemblid)
				continue
		#define all other variables
		chrom = 'chr' + fields[3]
		if fields[4] == '1':
			strand = '+'
		elif fields[4] == '-1':
			strand = '-'
		ORFstart = str(int(fields[1]) - int(chromstart_38) + int(chromstart_37))
		ORFstop = str(int(fields[2]) - int(chromstart_38) + int(chromstart_37))
		ORFsize = str(int(fields[2]) - int(fields[1]) + 1)
		ORFseq = fields[7]
		ORFscore = fields[5]
		covF1prop = fields[6]
		pep_seq = fields[8]
		#check if in LNCipedia or not
		if int(data.json()['count']) == 0:
			name = ensemblid
			inLNCipedia = "0"
			outfile.write(name + "\t" + chrom + "\t" + strand + "\t" + chromstart_37 + "\t" + chromend_37 + "\t" + ORFstart + "\t" + ORFstop + "\t" + ORFsize + "\t" + ORFseq + "\t" + ORFscore+ "\t" + covF1prop+ "\t" + pep_seq + "\t" + inLNCipedia + "\n")
		elif int(data.json()['count']) == 1:
			name = data.json()['transcripts'][0]['lncipediaTranscriptID']
			inLNCipedia = "1"
			outfile.write(name + "\t" + chrom + "\t" + strand + "\t" + chromstart_37 + "\t" + chromend_37 + "\t" + ORFstart + "\t" + ORFstop + "\t" + ORFsize + "\t" + ORFseq + "\t" + ORFscore+ "\t" + covF1prop+ "\t" + pep_seq + "\t" + inLNCipedia + "\n")
		elif int(data.json()['count']) > 1:
			error_list.append(ensemblid)
		print >>sys.stderr, "%d lines finished\r" % count,
print >>sys.stderr, "\n"

outfile.close()
infile.close()
cur_37.close()
cur_38.close()
db_37.close()
db_38.close()

if len(error_list) != 0:
	print >>sys.stderr, bcolors.WARNING + "Warning: these Ensembl ID's had multiple LNCpedia references. Please check them manually.\n" + bcolors.ENDC
	for item in error_list:
		print >>sys.stderr, '\t-\t' + item
if unknown_count != 0:
	print >>sys.stderr, bcolors.WARNING + "Warning: there were %d lines without transcript ID. These were skipped automatically.\n" % unknown_count + bcolors.ENDC
print >>sys.stderr, bcolors.OKGREEN + "DONE!\n" + bcolors.ENDC
