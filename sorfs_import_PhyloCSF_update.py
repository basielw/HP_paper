#!/usr/bin/python

# import modules used here
import os,sys
import argparse
import MySQLdb
import getpass

# INSERT query
add = ("INSERT INTO PhyloCSF "
		"(ORFs_idORFs, ORFs_Transcripts_idTranscripts, PhyloCSF_score, peptide_seq) "
		"VALUES (%(idORFs)s, %(idTranscripts)s, %(PhyloCSF_score)s, %(peptide_seq)s)")
select_trans = "SELECT idTranscripts FROM small_orfs.Transcripts WHERE LNCipediaID = %(name)s"
select_orf = "SELECT idORFs, ORFsequence FROM small_orfs.ORFs WHERE Transcripts_idTranscripts = %(idTranscripts)s"
select_PhyloCSF = "SELECT idPhyloCSF FROM small_orfs.PhyloCSF WHERE ORFs_idORFs = %(idORFs)s"

username = raw_input("Enter your small_orfs username:")
password = getpass.getpass("Enter your small_orfs password:")

# Gather our code in a main() function
def main():
 parser = argparse.ArgumentParser(description='Some very cool algoritm.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

 parser.add_argument('input_file', metavar='INPUT', type=str, help='Database file')

 args = parser.parse_args()
 #print(args.accumulate(args.integers))

 # connect to database
 db = MySQLdb.connect(host="localhost",     # your host, usually localhost
					user=username,         # your username
					passwd=password,    # your password
					db="small_orfs")        # name of the data base

 # you must create a Cursor object. It will let
 # you execute all the queries you need
 cur = db.cursor()

 codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
    'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }

 handle = open(args.input_file, "rU")
 count = 0
 for i, line in enumerate(handle):
 	if i == 0:
 		continue
	else:
		line = line.rstrip()
		fields = line.split('\t')
		count += 1
		name = str(fields[0])
		cur.execute(select_trans, {'name': name})
		row_trans = cur.fetchone()
		trans_id = row_trans[0]
		cur.execute(select_orf, {'idTranscripts': trans_id})
		numrows = cur.rowcount
		orf_id = None
		for x in xrange(0,numrows):
			row_orf = cur.fetchone()
			cds = str(row_orf[1])
			proteinsequence = ''
			for n in range(0,len(cds),3):
				try:
					proteinsequence += codontable[cds[n:n+3]]
				except KeyError:
					orf_id = row_orf[0]
					print >>sys.stderr, "There seems to be a mistake in the sequence of the ORF with ID %s\n" % orf_id
			if proteinsequence == str(fields[1]):
				orf_id = row_orf[0]
				break
			else: 
				continue
		if orf_id == None:
			print >>sys.stderr, "There was no ORF of transcript %s found that matched the PhyloCSF protein sequence\n" % name
			continue
		else:
			cur.execute(select_PhyloCSF, {'idORFs': int(orf_id)})
			row = cur.fetchone()
			if row == None:
				data = {
					'idORFs': int(orf_id),
					'idTranscripts': int(trans_id),
					'PhyloCSF_score': fields[2],
					'peptide_seq': str(fields[1])
					 }
				cur.execute(add, data)
			else:
				continue
		print >>sys.stderr, "%d lines finished\r" % count,
		if i % 1000 == 0:
			db.commit()
 handle.close()

 db.commit()
 cur.close()
 db.close()

# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
 main()
