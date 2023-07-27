#! /usr/bin/env python3 

# get citations  
# https://gist.github.com/mcfrank/c1ec74df1427278cbe53

# This script executes a simple pubmed search
# with publication year constraint  
#   Example:
# pubmed "<search1> <search2>"
# pubmed "<search1> <search2>" <./mydir/myfile>
# pubmed "<search1> <search2>" 2015 <./mydir/myfile>

import sys, os, operator 
import xml.etree.ElementTree as ET

#### Name Convention ###
NC="FL"
SINCE_YEAR='2020'
# FmL=Joel T Nigg
# LFm=Nigg, Joel T
# FL=Joel Nigg
# LF=Nigg, Joel
myxml=''
lines=[]

def myhelp():
    print("\nThis tool searches Pubmed publications for metadata")
    print('%s "<term1> <term2>"' % sys.argv[0])
    print("\noptionally add the cached pubmed-articles.xml file to prevent downloading from Pubmed API")
    print('%s "<term1> <term2>" pubmed-articles.xml' % sys.argv[0])

if len(sys.argv) < 2:
    myhelp()
    sys.exit()

SEARCH1=sys.argv[1]

if len(sys.argv) >= 3:
    filename=''
    if os.path.exists(sys.argv[2]):
        filename=sys.argv[2]
    else:
        SINCE_YEAR=sys.argv[2]
    if filename:
        try:
            with open(filename, "r") as fh:
                if fh.read(5) == '<?xml':
                    myxml='<?xml' + fh.read()
                else:
                    lines=fh.readlines()
        except:
            print('*** Could not open %s ! *** ' % filename)
            myhelp()
            sys.exit()

def main(myxml,lines):
    newxml=False
    if not myxml:
        if not lines:
            SSTR='"%s" AND ("%s"[Date - Publication] : "3000"[Date - Publication])' % (SEARCH1, SINCE_YEAR)
            lines=searchArticles(SSTR)
        if not lines:
            print("No articles found. Please change your search.")
            return 1
        myxml=getArticles(lines)
        newxml=True
        if not myxml:
            print("No articles could be retrieved. Please use a valid XML file")
            return 1
    root = ET.fromstring(myxml)
    SEARCH1DICT={}
    SEARCH1DICTFL={}
    AFFDICT={}
    for article in root.findall("PubmedArticle"):
        pmid = article.find("MedlineCitation/PMID").text
        year = article.find("MedlineCitation/Article/Journal/JournalIssue/PubDate/Year")
        title = article.find("MedlineCitation/Article/ArticleTitle").text
        if year is None: year = 'NA'
        else: year = year.text
        aulist = article.findall("MedlineCitation/Article/AuthorList/Author")
        #print(pmid, year, len(aulist), affiliation, aulist, ET.dump(root))
        SEARCH1L=''
        iauthors=0
        firstauthor = ''
        SEARCH1firstauthor = ''
        SEARCH1lastauthor = ''
        for author in aulist:
            iauthors+=1
            atext = ET.tostring(author).decode()
            affiliations = author.findall('AffiliationInfo/Affiliation')
            lastname = author.find('LastName') 
            forename = author.find('ForeName')
            if lastname == None:
                continue
            if forename == None:
                fullname = lastname.text
            else:
                firstname = forename.text.split()[0]
                if len(firstname) <= 1:
                    firstname = forename.text.split()[-1]
                if len(firstname) <= 1:
                    firstname = forename.text
                if NC=='FmL':
                    fullname = forename.text + ' ' + lastname.text
                elif NC=='LFm':
                    fullname = lastname.text + ', ' + forename.text
                elif NC=='FL':
                    fullname = firstname + ' ' + lastname.text
                elif NC=='LF':
                    fullname = lastname.text + ', ' + firstname
                
            if not firstauthor:
                firstauthor = fullname
            lastauthor = fullname
            if SEARCH1 in atext:
                SEARCH1L+=fullname+'\n'
                SEARCH1DICT[fullname] = SEARCH1DICT.get(fullname,0) + 1
                if not SEARCH1firstauthor:
                    SEARCH1firstauthor = fullname
                SEARCH1lastauthor = fullname
                for aff in affiliations:
                    AFFDICT[aff.text] = AFFDICT.get(aff.text,0) + 1                            

        if SEARCH1firstauthor != firstauthor:
           SEARCH1firstauthor = ''
        if SEARCH1lastauthor != lastauthor:
           SEARCH1lastauthor = ''

        if SEARCH1firstauthor:
             SEARCH1DICTFL[SEARCH1firstauthor]=SEARCH1DICTFL.get(SEARCH1firstauthor,0) + 1
        if SEARCH1lastauthor:
             SEARCH1DICTFL[SEARCH1lastauthor]=SEARCH1DICTFL.get(SEARCH1lastauthor,0) + 1

        if SEARCH1L:
            print('\n\n***** Article: %s / Year: %s / Total Authors: %s *****' % (pmid,year,iauthors))
            print('*** Title: %s\n' % title)
            
            print ("* "+SEARCH1+':\n'+SEARCH1L)
             
            print ('First Author: '+firstauthor)
            print ('Last Author: '+lastauthor)
            
    SEARCH1DICT = dict( sorted(SEARCH1DICT.items(), key=operator.itemgetter(1),reverse=True))
    SEARCH1DICTFL = dict( sorted(SEARCH1DICTFL.items(), key=operator.itemgetter(1),reverse=True))
    AFFDICT = dict( sorted(AFFDICT.items(), key=operator.itemgetter(1),reverse=True))

    dlen=25
    l=list(AFFDICT)
    if len(AFFDICT) < dlen:
        dlen=len(AFFDICT)
    print ('\n\n****** Top %s affiliations by frequency *******' % dlen)

    for i in range(0,dlen-1):
        print(str('(%s) %s' % (AFFDICT[l[i]],l[i]), "utf-8"))
        
    print ('\n\n****** All collaborators by frequency *******')

    print ("\n**************** "+SEARCH1+':')
    print (SEARCH1DICT)

    print ('\n\n****** First and last authors by frequency *******')
    print ('****** (inaccurate when shared first/last authorship) *******')

    print ("\n**************** "+SEARCH1+':')
    print (SEARCH1DICTFL)

    mincollab=2
    print ('\n\n****** eMail lists of authors with %s or more collaborations *******' % mincollab)
    print ('****** paste these in Outlook and hit "Check Names"')
    print ('****** (might need to change variable NC)')

    print ("\n**************** "+SEARCH1+':')
    elist=""
    for k,v in SEARCH1DICT.items():
        if v>=mincollab:
            elist+=k+'; '
    print (elist)
   
    if newxml:
       print ('\nDownloaded articles to pubmed-articles.xml. You can now run this command without querying Pubmed')
       print ('%s "%s" "%s" ./pubmed-articles.xml' % (sys.argv[0],sys.argv[1],sys.argv[2]))
     
class redirector(object):
    def __init__(self, filename="pubmed-articles.xml"):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        #self.terminal.write(message)
        self.log.write(message)


def searchArticles(searchStr):
    try:
        import entrezpy.esearch.esearcher
    except:
        print('"entrezpy" not installed. Please run:')
        print('pip3 install --user --upgrade entrezpy')
        return ""

    e = entrezpy.esearch.esearcher.Esearcher("entrezpy",
                                             "test@testorg.io",
                                             apikey=None,
                                             apikey_var=None,
                                             threads=None,
                                             qid=None)
    analyzer = e.inquire({'db' : 'pubmed',
                          'term' : searchStr,
                          'retmax' : '100000',
                          'rettype' : 'uilist'})
    return analyzer.result.uids

def getArticles(pmids):
    try:
        import entrezpy.efetch.efetcher
    except:
        print('"entrezpy" not installed. Please run:')
        print('pip3 install --user --upgrade entrezpy')
        return ""

    e = entrezpy.efetch.efetcher.Efetcher('entrezpy',
                                          'test@testorg.io',
                                          apikey=None,
                                          apikey_var=None,
                                          threads=None,
                                          qid=None)

    oldout = sys.stdout
    sys.stdout = redirector()

    analyzer = e.inquire({'db' : 'pubmed',
                          'id' : pmids,
                          'retmax' : '100000',
                          'retmode' : 'xml',
                          'rettype' : 'abstract'})
    sys.stdout = oldout

    with open("pubmed-articles.xml", "r") as fh:
        myxml=fh.read()

    return myxml


if __name__ == "__main__":
    sys.exit(main(myxml,lines))


