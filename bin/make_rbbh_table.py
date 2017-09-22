#!/usr/bin/env python3
# pylint: disable=C0103
# make_rbbh_table.py
# Script to generate a table of RBBH matches from a local SQLite3 database,
# when provided with a protein locus tag identifier.
#
# NOTE: This script is one-shot, and brittle to the presence of correct
# input files


"""
This script produces a table for RBBH
"""
import logging
import sqlite3
import sys

from argparse import ArgumentParser


# Process command-line arguments
def parse_cmdline(args):
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="make_rbbh_table.py")
    parser.add_argument(dest="locus_tag", type=str,
                        action="store", default=None,
                        help="Locus tag of protein query")
    parser.add_argument(dest="database", type=str,
                        action="store", default=None,
                        help="SQLite3 database of RBBH")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        type=str,
                        action="store", default=None,
                        help="Logfile location")
    return parser.parse_args()




# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('make_rbbh_table.py')
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

    # Connect to database and get cursor
    with sqlite3.connect(args.database) as conn:
        cursor = conn.cursor()
    logger.info("Opened database %s" % args.database)

    # Construct SQL query from locus_tag and execute
    sql = 'SELECT * FROM rbbh ' +\
          'INNER JOIN rbbh_rbbh_data ON ' +\
          'rbbh.rbbh_id=rbbh_rbbh_data.rbbh_id ' +\
          'INNER JOIN rbbh_data ON ' +\
          'rbbh_rbbh_data.data_id=rbbh_data.data_id ' +\
          'WHERE locus_tag_1=?'
    logger.info("Executing query: %s" % sql.replace('?', args.locus_tag))
    cursor.execute(sql, (args.locus_tag,))

    # Collect results
    results = cursor.fetchall()
    logger.info("Retrieved %d rows" % len(results))

    # Write table to file
    outfname = "%s.tab" % args.locus_tag
    logger.info("Writing RBBH summary table to %s" % outfname)
    header = '\t'.join(['locus_tag', 'rbbh', 'bitscore', 'pid', 'cov'])
    with open(outfname, 'w') as ofh:
        print(header, file=ofh)
        with sys.stdout as stdout:
            print("%s makes RBBH with..." % args.locus_tag, file=stdout)
            for r in results:
                print("\t%s" % r[2], file=stdout)
                print("%s" % '\t'.join([str(e) for e in
                                        [r[1], r[2], r[7], r[8], r[9]]]),
                      file=ofh)

    # Say we're done
    logger.info("Done.")
