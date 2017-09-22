#!/usr/bin/env python3
#
# extract_rbbh.py
#
# Script to extract RBBH matches from a local SQLite3 database, when
# provided with a protein locus tag identifier, and the path to a file
# containing sequences for all proteins in the database
#
# NOTE: This script is one-shot, and brittle to the presence of correct 
# input files

import logging
import os
import sqlite3
import sys

from argparse import ArgumentParser

from Bio import SeqIO

# Process command-line arguments
def parse_cmdline(args):
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="extract_rbbh.py")
    parser.add_argument("--db", dest="database",
                        action="store", default=None,
                        help="SQLite3 database of RBBH")
    parser.add_argument("--locus_tag", dest="locus_tag",
                        action="store", default=None,
                        help="Locus tag of protein query")
    parser.add_argument("--seqfile", dest="seqfile",
                        action="store", default=None,
                        help="FASTA file of all protein sequences")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    return parser.parse_args()




# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('extract_rbbh.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(args)
    logger.info("Command-line: %s" % ' '.join(sys.argv))

    # Have we got a database, locus tag, and sequence file
    if args.database is None:
        logger.error("No input database (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % args.database)
    if args.locus_tag is None:
        logger.error("No locus_tag provided (exiting)")
        sys.exit(1)
    logger.info("Locus tag for query: %s" % args.locus_tag)
    if args.seqfile is None:
        logger.error("No sequence file (exiting)")
        sys.exit(1)
    logger.info("Sequence file: %s" % args.seqfile)

    # Check for spaces in input: there should be none
    #arglist = [args.database, args.locus_tag, args.seqfile]
    #for fname in arglist:
    #    if ' ' in  os.path.abspath(fname):
    #        logger.error("File or directory '%s' contains whitespace" % fname)
    #        logger.error("(exiting)")
    #        sys.exit(1)

    # Connect to database and get cursor
    with sqlite3.connect(args.database) as conn:
        cursor = conn.cursor()
    logger.info("Opened database %s" % args.database)

    # Load in sequence data, keyed by locus tag
    with open(args.seqfile) as fh:
        seqdict = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    logger.info("Opened %s as dictionary; %d records" % (args.seqfile,
                                                         len(seqdict)))

    # Construct SQL query from locus_tag and execute
    sql = 'SELECT * FROM rbbh WHERE locus_tag_1=?'
    logger.info("Executing query: %s" % sql.replace('?', args.locus_tag))
    cursor.execute(sql, (args.locus_tag,))

    # Collect results
    results = cursor.fetchall()
    logger.info("Retrieved %d rows" % len(results))

    # Report results
    with sys.stdout as outfh:
        print("%s makes RBBH with..." % args.locus_tag, file=outfh)
        for r in results:
            print("\t%s" % r[2], file=outfh)

    # Write sequences to file
    outfname = "%s.rbbh.fasta" % args.locus_tag
    logger.info("Writing RBBH sequences to %s" % outfname)
    SeqIO.write([seqdict[args.locus_tag]] + [seqdict[r[2]] for r in results],
                outfname, 'fasta')

    # Say we're done
    logger.info("Done.")
            
            

