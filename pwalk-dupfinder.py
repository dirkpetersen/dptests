#! /usr/bin/env python3

"""
pwalk-dupfinder aggregates pwalk output files and
finds dupicate files in different paths that have 
the same filename, modification time and size
"""
# internal modules
import sys, os, argparse, csv, platform, textwrap, inspect
if sys.platform.startswith('linux'):
    import getpass, pwd, grp
# stuff from pypi
import duckdb

__app__ = 'pwalk dupfinder'
__version__ = '0.0.1'

def main():

    #remove trailing slash
    if args.csvpath.endswith('/'):
        args.csvpath = args.csvpath[:-1]

    if args.subcmd in ['report', 'rep']:
    
        # Initialize DuckDB connection
        cores = 16
        conn = duckdb.connect(':memory:')
        conn.execute(f'PRAGMA threads={cores};')
                
        if os.path.isdir(args.csvpath):
            #print('Using:', args.csvpath)
            # if dir, create a union view over all the CSV files in the directory
            print('find csv files ....', flush=True)
            csv_files = conn.execute(f"SELECT * FROM glob('{args.csvpath}/*.csv')").fetchall()
            #print('csv_files', csv_files)
            query_parts = [f"SELECT * FROM read_csv_auto('{csv_file[0]}')" for csv_file in csv_files]
            union_query = " UNION ALL ".join(query_parts)
        elif os.path.isfile(args.csvpath):
            union_query = f"SELECT * FROM read_csv_auto('{args.csvpath}')"
        else:
            print(f"Path {args.csvpath} is not a file or directory")
            return False
        print("Execute query:", union_query)
        conn.execute(f"CREATE VIEW combined_csvs AS {union_query}")
        # Now you can query the combined data from all CSV files directly
        


def parse_arguments():
    """
    Gather command-line arguments.
    """       
    parser = argparse.ArgumentParser(prog='pwalk-summary ',
        description='A (mostly) automated tool for archiving large scale data ' + \
                    'after finding folders in the file system that are worth archiving.')
    parser.add_argument( '--debug', '-d', dest='debug', action='store_true', default=False,
        help="verbose output for all commands")
    parser.add_argument('--cores', '-c', dest='cores', action='store', default='4', 
        help='Number of cores to be allocated for the machine. (default=4)')
    parser.add_argument('--version', '-v', dest='version', action='store_true', default=False, 
        help='print Froster and Python version info')
    
    subparsers = parser.add_subparsers(dest="subcmd", help='sub-command help')
    # ***
    parser_config = subparsers.add_parser('report', aliases=['rep'], 
        help=textwrap.dedent(f'''
            Print a pwalk report             
        '''), formatter_class=argparse.RawTextHelpFormatter)
    parser_config.add_argument( '--monitor', '-m', dest='monitor', action='store', default='',
        metavar='<email@address.org>', help='setup as a monitoring cronjob ' +
        'on a machine and notify an email address')
    parser_config.add_argument('csvpath', action='store', default="", nargs='?',
        help='csv path can be a file or a folder' +
                '(default=~ home directory)  ')
    
    # ***

    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)               

    return parser.parse_args()

if __name__ == "__main__":
    try:
        args = parse_arguments()
        if main():
            sys.exit(0)
        else:
            sys.exit(1)
    except KeyboardInterrupt:
        print('\nExit !')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)                