#!/usr/bin/python

# import modules used here
import sys
import argparse
import MySQLdb
import getpass

# INSERT query
add_CPAT = ("INSERT INTO CPAT "
		"(ORFs_idORFs, ORFs_Transcripts_idTranscripts, Fickett_score, Hexamer_score, CPAT_score, CPAT_score_newlogit) "
		"VALUES (%(idORFs)s, %(idTranscripts)s, %(Fickett_score)s, %(Hexamer_score)s, %(CPAT_score)s, %(CPAT_score_newlogit)s)")

username = raw_input("Enter your small_orfs username:")
password = getpass.getpass("Enter your small_orfs password:")

db = MySQLdb.connect(host="localhost",
					user=username,
					passwd=password,
					db="small_orfs")

cur = db.cursor()

CPAT_old_file = open("/home/basielw/Database/sorfs_output_old")
CPAT_new_file = open("/home/basielw/Database/sorfs_output_new")

count = 0
for i, (line_old, line_new) in enumerate(zip(CPAT_old_file,CPAT_new_file)):
	if i == 0:
		continue
	else:
		line_old = line_old.rstrip()
		fields_old = line_old.split('\t')
		line_new = line_new.rstrip()
		fields_new = line_new.split('\t')
		if fields_old[0] != fields_new[0]:
			raise ValueError("CPAT outputs don't align")
		count += 1
		idORFs = int((((fields_old[0].split('___'))[0]).split('='))[1])
		idTranscripts = int((((fields_old[0].split('___'))[1]).split('='))[1])
		data_CPAT = {
			'idORFs': idORFs,
			'idTranscripts': idTranscripts,
			'Fickett_score': fields_old[5],
			'Hexamer_score': fields_old[6],
			'CPAT_score': fields_old[7],
			'CPAT_score_newlogit': fields_new[7]
			}
		cur.execute(add_CPAT, data_CPAT)
		if (i % 1000) == 0:
			db.commit()
		print >>sys.stderr, "%d lines finished\r" % count,
db.commit()
CPAT_old_file.close()
CPAT_new_file.close()
cur.close()
db.close()
