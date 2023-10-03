#! /usr/bin/env python3

"""
Froster automates much of the challening tasks when 
archiving many Terabytes of data on large (HPC) systems.
"""
# internal modules
import sys, os, argparse, json, configparser, csv, platform, asyncio
import urllib3, datetime, tarfile, zipfile, textwrap, tarfile, time
import concurrent.futures, hashlib, fnmatch, io, math, signal, shlex
import shutil, tempfile, glob,  subprocess, itertools, socket, inspect
if sys.platform.startswith('linux'):
    import getpass, pwd, grp
# stuff from pypi
#import requests, duckdb, boto3, botocore, psutil
import duckdb

__app__ = 'pwalk query'
__version__ = '0.0.1'

def main():

    if args.subcmd in ['report', 'rep']:
    
        # Initialize DuckDB connection
        conn = duckdb.connect()

        # Create a union view over all the CSV files in the directory
        print('csv_files ....', flush=True)        
        csv_files = conn.execute(f"SELECT * FROM glob('{args.csvpath}/*.csv')").fetchall()
        print('query_parts ....', flush=True)        
        query_parts = [f"SELECT * FROM read_csv_auto('{csv_file[0]}')" for csv_file in csv_files]
        print('union_query ....', flush=True)        
        union_query = " UNION ALL ".join(query_parts)
        conn.execute(f"CREATE VIEW combined_csvs AS {union_query}")

        # Now you can query the combined data from all CSV files directly
        print('fetch_all ....', flush=True)        
        result = conn.execute("SELECT * FROM combined_csvs WHERE ...").fetchall()

        print(result)


class Reporter:
    def __init__(self, args):
        self.args = args

    def query(self, pwalkfolder):

        # Connect to an in-memory DuckDB instance
        con = duckdb.connect(':memory:')
        #con.execute('PRAGMA experimental_parallel_csv=TRUE;') # now standard 
        con.execute(f'PRAGMA threads={self.args.cores};')


        pwalkcsv = self.args.pwalkcsv
        with tempfile.NamedTemporaryFile() as tmpfile:
            with tempfile.NamedTemporaryFile() as tmpfile3:
                # removing all files from pwalk output, keep only folders
                mycmd = f'grep -v ",-1,0$" "{pwalkcsv}" > {tmpfile3.name}'
                self.cfg.printdbg(f' Running {mycmd} ...', flush=True)
                result = subprocess.run(mycmd, shell=True)
                if result.returncode != 0:
                    print(f"Folder extraction failed: {mycmd}")
                    return False
                # Temp hack: e.g. Revista_Española_de_Quimioterapia in Spellman
                # Converting file from ISO-8859-1 to utf-8 to avoid DuckDB import error
                # pwalk does already output UTF-8, weird, probably duckdb error 
                mycmd = f'iconv -f ISO-8859-1 -t UTF-8 {tmpfile3.name} > {tmpfile.name}'
                self.cfg.printdbg(f' Running {mycmd} ...', flush=True)  
                result = subprocess.run(mycmd, shell=True)
                if result.returncode != 0:
                    print(f"File conversion failed: {mycmd}")
                    return False
            
        sql_query = f"""SELECT UID as User,
                        st_atime as AccD, st_mtime as ModD,
                        pw_dirsum/1073741824 as GiB, 
                        pw_dirsum/1048576/pw_fcount as MiBAvg,                            
                        filename as Folder, GID as Group,
                        pw_dirsum/1099511627776 as TiB,
                        pw_fcount as FileCount, pw_dirsum as DirSize
                    FROM read_csv_auto('{tmpfile.name}', 
                            ignore_errors=1)
                    WHERE pw_fcount > -1 AND pw_dirsum > 0
                    ORDER BY pw_dirsum Desc
                """  # pw_dirsum > 1073741824
        self.cfg.printdbg(f' Running SQL query on CSV file {tmpfile.name} ...', flush=True)
        rows = con.execute(sql_query).fetchall()
        # also query 'parent-inode' as pi,
        
        # Get the column names
        header = con.execute(sql_query).description

        totalbytes=0
        agedbytes=[]
        for i in daysaged:
            agedbytes.append(0)
        numhotspots=0

        mycsv = self._get_hotspots_path(pwalkfolder)
        self.cfg.printdbg(f' Running filter and write results to CSV file {mycsv} ...', flush=True)

        with tempfile.NamedTemporaryFile() as tmpcsv:
            with open(tmpcsv.name, 'w') as f:
                writer = csv.writer(f, dialect='excel')
                writer.writerow([col[0] for col in header])
                # 0:Usr,1:AccD,2:ModD,3:GiB,4:MiBAvg,5:Folder,6:Grp,7:TiB,8:FileCount,9:DirSize
                for r in rows:
                    row = list(r)
                    if row[3] >= self.thresholdGB and row[4] >= self.thresholdMB:
                        row[0]=self.uid2user(row[0])                        
                        row[1]=self.daysago(self._get_newest_file_atime(row[5],row[1]))
                        row[2]=self.daysago(row[2])
                        row[3]=int(row[3])
                        row[4]=int(row[4])
                        row[6]=self.gid2group(row[6])
                        row[7]=int(row[7])
                        writer.writerow(row)
                        numhotspots+=1
                        totalbytes+=row[9]
                    for i in range(0,len(daysaged)):
                        if row[1] > daysaged[i]:
                            agedbytes[i]+=row[9]
            if numhotspots > 0:
                shutil.copyfile(tmpcsv.name,mycsv)


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
        TABLECSV = '' # CSV string for DataTable
        SELECTEDFILE = '' # CSV filename to open in hotspots 
        MAXHOTSPOTS = 0    
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