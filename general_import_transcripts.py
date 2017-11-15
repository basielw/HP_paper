#!/usr/bin/python

# import modules used here
import sys
import argparse
import MySQLdb
import getpass

# INSERT query
add = ("INSERT INTO Transcripts "
		"(Name, Chromosome, chromStart, chromEnd, Strand, mRNA_size) "
		"VALUES (%(Name)s, %(Chromosome)s, %(chromStart)s, %(chromEnd)s, %(Strand)s, %(mRNA_size)s)")

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
 for line in handle:
	count += 1
	line = line.rstrip()
	fields = line.split('\t')
	size = fields[10].split(',')
	size_int = []
	for item in size:
		if item == '':
			continue
		else:
			size_int.append(int(item))
	data = {
		'Name': fields[3],
		'Chromosome': fields[0],
		'chromStart': int(fields[1]) + 1,
		'chromEnd': fields[2],
		'Strand': fields[5],
		'mRNA_size': sum(size_int)
		 }
	cur.execute(add, data)
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
