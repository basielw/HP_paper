#!/usr/bin/python

# import modules used here
import sys
import argparse
import MySQLdb
import getpass

# INSERT query
select_trans = "SELECT idTranscripts FROM small_orfs.Transcripts WHERE Name = %(name)s"
select_orfs = "SELECT idORFs, ORFstart_mRNA FROM small_orfs.ORFs WHERE Transcripts_idTranscripts = %(idTranscripts)s"
update_orf = "UPDATE ORFs SET ORFstart_transcript = %(ORFstart)s WHERE idORFs = %(idORFs)s"

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

	handle = open(args.input_file, "rU")
	count = 0
	for i, line in enumerate(handle):
		line = line.rstrip()
		fields = line.split('\t')
		exonsizes = fields[10].rstrip(',\n').split(',')
		exonstarts = fields[11].rstrip( ',\n' ).split( ',' )
		cur.execute(select_trans, {'name': fields[3]})
		row = cur.fetchone()
		idTranscripts = row[0]
		cur.execute(select_orfs, {'idTranscripts': idTranscripts})
		for row in cur:
			idORFs = int(row[0])
			old_ORFstart = int(row[1])
			previoussizes = []
			print >>sys.stderr, "%d rows in small_orfs.ORFs finished\r" % count,
			for size, start in zip(exonsizes, exonstarts):
				previoussizes.append(int(size))
				if sum(previoussizes) > old_ORFstart:
					new_ORFstart = int(fields[1]) + int(start) + old_ORFstart - sum(previoussizes[:-1])
					break
			cur.execute(update_orf, {'ORFstart': new_ORFstart, 'idORFs': idORFs})
		if (i % 1000) == 0:
			db.commit()

	handle.close()
	db.commit()
	cur.close()
	db.close()


# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
	main()
