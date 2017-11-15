#!/usr/bin/python

# import modules used here
import sys
import argparse
import MySQLdb
import getpass

# INSERT query
add_orf = ("INSERT INTO ORFs "
		"(ORFstart, ORFsize, ORFsequence, Transcripts_idTranscripts) "
		"VALUES (%(ORFstart)s, %(ORFsize)s, %(ORFsequence)s, %(idTranscripts)s)")
add_CPAT = ("INSERT INTO CPAT "
		"(ORFs_idORFs, ORFs_Transcripts_idTranscripts, Fickett_score, Hexamer_score, CPAT_score) "
		"VALUES (%(idORFs)s, %(idTranscripts)s, %(Fickett_score)s, %(Hexamer_score)s, %(CPAT_score)s)")
select_trans = "SELECT idTranscripts FROM pan_cancer.Transcripts WHERE Name = %(name)s"

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
					db="pan_cancer")        # name of the data base

	# you must create a Cursor object. It will let
	# you execute all the queries you need
	cur = db.cursor()

	handle = open(args.input_file, "rU")
	count = 0
	last_name = ('empty', 'empty')
	for i, line in enumerate(handle):
		if i == 0:
			continue
		else:
			line = line.rstrip()
			fields = line.split('\t')
			count += 1
			#as this script will use the same idTranscripts multiple times consecutively,
			#it would be a waste of time to search for it every time
			#thats why we introduce last_idTranscripts
			if last_name[0] == fields[0]:
				item = last_name[1]
			else:
				cur.execute(select_trans, {'name': fields[0]})
				row = cur.fetchone()
				item = row[0]
				last_name = (fields[0], item)
			data_orf = {
				'ORFstart': fields[2],
				'ORFsize': fields[3],
				'ORFsequence': fields[4],
				'idTranscripts': int(item)
					}
			cur.execute(add_orf, data_orf)
			data_CPAT = {
				'idORFs': cur.lastrowid,
				'idTranscripts': int(item),
				'Fickett_score': fields[5],
				'Hexamer_score': fields[6],
				'CPAT_score': fields[7],
				}
			cur.execute(add_CPAT, data_CPAT)
			if (i % 1000) == 0:
				db.commit()
			print >>sys.stderr, "%d lines finished\r" % count,
	handle.close()

	db.commit()
	cur.close()
	db.close()


# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
	main()
