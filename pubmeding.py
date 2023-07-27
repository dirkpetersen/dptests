#! /usr/bin/env python3

"""
Pubmed find uses pubmed_sdk to
search for Pubmed Articles using CLI
"""

import sys, os, argparse, textwrap, json
from pubmed_sdk import PubMed
#json, configparser, csv, platform

__app__ = 'Pubmed Find'
__version__ = '0.1'

def main():
    
    if args.debug:
        pass

    if len(sys.argv) == 1:        
        print(textwrap.dedent(f'''\n
            For example, use one of these commands:
              pubmed-find --search "COVID-19"
              pubmed-find id1 id2 id3
            '''))

    if args.version:
        print(f'Pubmed-Find version: {__version__}')
        print(f'Python version:\n{sys.version}')
        return True
    
    pubmed = PubMed()

    if args.search:
        if args.year:
            args.search = f'{args.search} AND ("{args.year}"[Date - Publication] : "3000"[Date - Publication])'
        results = pubmed.search(args.search,term=f'retmax={args.retmax}')
        args.pmids = results['id_list']

    if args.pmids:
        details = pubmed.fetch_details(args.pmids).get('PubmedArticle')
        #This method accepts a list of IDs and returns a list of dictionaries containing the details of each article.

        for detail in details:
            article = detail['MedlineCitation']['Article']
            pretty=json.dumps(article, indent=2)
            print(pretty,'\n\n')
            #print(article['ArticleTitle'])
            #print(article['Abstract']['AbstractText'], '\n\n')

    print('List of PMIDs:')
    print(" ".join(args.pmids))
   

def parse_arguments():
    """
    Gather command-line arguments.
    """       
    parser = argparse.ArgumentParser(prog='pubmed  ',
        description='A (mostly) automated tool for archiving large scale data ' + \
                    'after finding folders in the file system that are worth archiving.')
    parser.add_argument( '--debug', '-d', dest='debug', action='store_true', default=False,
        help="verbose output for all commands")
    parser.add_argument('--version', '-v', dest='version', action='store_true', default=False, 
        help='print Froster and Python version info')      
    parser.add_argument('--search', '-s', dest='search', action='store', default='', 
        help='PubMED search string ')
    parser.add_argument('--retmax', '-r', type=int, dest='retmax', action='store', default=1000, 
        help='Maximum number of of results to return')    
    parser.add_argument('--year', '-y', dest='year', action='store', default='', 
        help='Limit publications to year newer than this one')
    parser.add_argument('pmids', action='store', default=[],  nargs='*',
        help='List of PMIDs you would like to get details on (separated by space), ' +
                ' ')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)               

    return parser.parse_args()

if __name__ == "__main__":
    if not sys.platform.startswith('linux'):
        print('This software currently only runs on Linux x64')
        sys.exit(1)
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
