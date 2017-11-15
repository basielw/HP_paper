#!/usr/bin/python

# import modules used here
import os,sys
import argparse
import MySQLdb
import getpass

# INSERT query
add = ("INSERT INTO PhyloCSF "
		"(ORFs_idORFs, ORFs_Transcripts_idTranscripts, PhyloCSF_score) "
		"VALUES (%(idORFs)s, %(idTranscripts)s, %(PhyloCSF_score)s)")
select_trans = "SELECT idTranscripts FROM pan_cancer.Transcripts WHERE Name = %(name)s"
select_orf = "SELECT idORFs FROM pan_cancer.ORFs WHERE Transcripts_idTranscripts = %(idTranscripts)s AND ORFstart = %(ORFstart)s"
select_score = "SELECT idPhyloCSF FROM pan_cancer.PhyloCSF WHERE ORFs_idORFs = %(idORFs)s"
add_orf = ("INSERT INTO ORFs "
		"(ORFstart, ORFsize, ORFsequence, Transcripts_idTranscripts, peptideSeq) "
		"VALUES (%(ORFstart)s, %(ORFsize)s, %(ORFsequence)s, %(idTranscripts)s, %(peptideSeq)s)")

username = raw_input("Enter your pan_cancer username:")
password = getpass.getpass("Enter your pan_cancer password:")

# Gather our code in a main() function
def main():
 parser = argparse.ArgumentParser(description='Some very cool algoritm.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

 parser.add_argument('input_file', metavar='INPUT', type=str, help='Database file')
 parser.add_argument('transcript_error_file', metavar='OUTPUT', type=str, help='Output file for errors in transcript names (no matching transcript found)')

 args = parser.parse_args()
 #print(args.accumulate(args.integers))

 # connect to database
 db = MySQLdb.connect(host="localhost",     # your host, usually localhost
					user=username,         # your username
					passwd=password,    # your password
					db="pan_cancer")        # name of the data base

 # you must create a Cursor object. It will let
 # you execute all the queries you need
 cur = db.cursor()

 outfile_transcript = open(args.transcript_error_file, 'w')

 handle = open(args.input_file, "rU")
 count = 0
 for i, line in enumerate(handle):
	line = line.rstrip()
	fields = line.split('\t')
	count += 1
	name = os.path.basename(str(fields[0]))[:-3]
	try:
		PhyloCSF_score = float(fields[2])
	except ValueError:
		continue
	ORFstart = int(fields[3]) + 1
	ORFsize = int(fields[4]) - ORFstart + 4
	peptideSeq = str(fields[5])
	cur.execute(select_trans, {'name': name})
	row_trans = cur.fetchone()
	if row_trans == None:
		outfile_transcript.write(name + "\n")
		continue
	trans_id = row_trans[0]
	cur.execute(select_orf, {'idTranscripts': trans_id, 'ORFstart': ORFstart})
	row_orf = cur.fetchone()
	if row_orf == None:
		data_orf = {
			'ORFstart': ORFstart,
			'ORFsize': ORFsize,
			'ORFsequence': None,
			'idTranscripts': trans_id,
			'peptideSeq': peptideSeq
		}
		cur.execute(add_orf,data_orf)
		orf_id = cur.lastrowid
	else:
		orf_id = int(row_orf[0])
	cur.execute(select_score, {'idORFs': orf_id})
	row_score = cur.fetchone()
	if row_score == None:
		data = {
			'idORFs': orf_id,
			'idTranscripts': trans_id,
			'PhyloCSF_score': PhyloCSF_score,
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
 outfile_transcript.close()

# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
 main()
