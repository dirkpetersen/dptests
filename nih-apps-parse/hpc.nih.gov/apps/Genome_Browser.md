

document.querySelector('title').textContent = 'UCSC Executables';
UCSC Executables
The [UCSC Genome Browser](http://genome.ucsc.edu) is maintained using a variety of executables that perform
functions ranging from sequence analysis and format conversion, to basic number crunching and statistics, to complex
database generation and manipulation.


These executables have been compiled and made available on the Helix Systems.


How to use
----------


The easiest way is to load the [ucsc environment module](/apps/modules.html):



```
$ module load ucsc
```

The **module load** statement can be placed in your [startup files](/docs/startup_files.html) for permanency.


Then, just call the executable of choice on the commandline:



```
$ faSplit byname scaffolds.fa outRoot/
```

List of current executables
---------------------------




```

aarToAxt		      addAveMedScoreToPsls	    addCols			  affyPairsToSample
agpAllToFaFile		      agpCloneCheck		    agpCloneList		  agpMergeChromScaf
agpToFa			      agpToGl			    agxToBed			  agxToIntronBeds
agxToTxg		      allenCollectSeq		    altAnalysis			  altPaths
altSplice		      altSummary		    aNotB			  apacheMonitor
assessLibs		      autoDtd			    autoSql			  autoXml
ave			      aveCols			    averagExp			  averageZoomLevels
avgTranscriptomeExps	      axtAndBed			    axtBest			  axtCalcMatrix
axtChain		      axtDropOverlap		    axtDropSelf			  axtFilter
axtForEst		      axtIndex			    axtPretty			  axtQueryCount
axtRecipBest		      axtRescore		    axtSort			  axtSplitByTarget
axtSwap			      axtToBed			    axtToChain			  axtToMaf
axtToPsl		      bamToPsl			    barChartMaxLimit		  bedClip
bedCommonRegions	      bedCons			    bedCoverage			  bedDown
bedExtendRanges		      bedGeneParts		    bedGraphPack		  bedGraphToBigWig
bedInGraph		      bedIntersect		    bedItemOverlapCount		  bedJoinTabOffset
bedJoinTabOffset.py	      bedMergeAdjacent		    bedOrBlocks			  bedPartition
bedPileUps		      bedRemoveOverlap		    bedRestrictToPositions	  bedSort
bedSplitOnChrom		      bedToBigBed		    bedToExons			  bedToFrames
bedToGenePred		      bedToPsl			    bedToTxEdges		  bedUp
bedWeedOverlapping	      bigBedInfo		    bigBedNamedItems		  bigBedSummary
bigBedToBed		      bigChainBreaks		    bigGenePredToGenePred	  bigGuessDb
bigHeat			      bigMafToMaf		    bigPslToPsl			  bigWigAverageOverBed
bigWigCat		      bigWigCluster		    bigWigCorrelate		  bigWigInfo
bigWigMerge		      bigWigSummary		    bigWigToBedGraph		  bigWigToWig
binFromRange		      blastToPsl		    blastXmlToPsl		  blat
blatServersCheck	      blatz			    blatzClient			  blatzServer
borfBig			      buildReleaseLog		    calc			  calcGap
catDir			      catUncomment		    ccCp			  chainAntiRepeat
chainBridge		      chainCleaner		    chainDbToFile		  chainFilter
chainMergeSort		      chainNet			    chainPreNet			  chainScore
chainSort		      chainSplit		    chainStitchId		  chainSwap
chainToAxt		      chainToPsl		    chainToPslBasic		  checkableBorf
checkAgpAndFa		      checkBigDbSnp		    checkCardinality		  checkChain
checkCoverageGaps	      checkHgFindSpec		    checkSgdSync		  checkTableCoords
checkUrlsInTable	      chopFaLines		    chromGraphFromBin		  chromGraphToBin
chromToUcsc		      clusterGenes		    clusterMatrixToBarChartBed	  clusterPsl
clusterRna		      colTransform		    consForBed			  convolve
countChars		      countNib			    cpg_lh			  createSageSummary
crTreeIndexBed		      crTreeSearchBed		    ctgFaToFa			  ctgToChromFa
dbDbToHubTxt		      dbFindFieldsWith		    dbSnoop			  dbSnpJsonToTab
dbTrash			      detab			    dnaMotifFind		  dnsInfo
dumpNib			      edwAddAssembly		    edwAddQaContamTarget	  edwAddQaEnrichTarget
edwAddSubscriber	      edwChangeFormat		    edwComparePeaks		  edwCorrectFileTags
edwCreateUser		      edwFakeManifestFromSubmit	    edwJob			  edwMakeContaminationQa
edwMakeEnrichments	      edwMakePairedEndQa	    edwMakeRepeatQa		  edwMakeReplicateQa
edwMakeValidFile	      edwQaEvaluate		    edwReallyRemoveFiles	  edwReplicatedPeaks
edwRetryJob		      edwRsyncEdwUser		    edwRunDaemon		  edwRunOnIds
edwSamPairedEndStats	      edwSolexaToSangerFastq	    edwSubmit			  eisenInput
emblMatrixToMotif	      embossToPsl		    endsInLf			  estLibStats
estOrient		      exonerateGffDoctor	    expMatrixToBarchartBed	  expToRna
faAlign			      faCmp			    faCount			  faFilter
faFilterN		      faFlyBaseToUcsc		    faFrag			  faGapLocs
faGapSizes		      fakeFinContigs		    fakeOut			  faNcbiToUcsc
faNoise			      faOneRecord		    faPolyASizes		  faRandomize
faRc			      faRenameRecords		    faSimplify			  faSize
faSomeRecords		      faSplit			    fastqStatsAndSubsample	  fastqToFa
faToFastq		      faToNib			    faToTab			  faToTwoBit
faToVcf			      faTrans			    faTrimPolyA			  faTrimRead
featureBits		      fetchChromSizes		    ffaToFa			  findMotif
findStanAlignments	      fishClones		    fixCr			  fixHarbisonMotifs
fixStepToBedGraph.pl	      fqToQa			    fqToQac			  fragPart
gapSplit		      gapToLift			    gbGetEntries		  gbOneAcc
gbSeqCheck		      gbToFaRa			    gcForBed			  gencodeVersionForGenes
genePredCheck		      genePredFilter		    genePredHisto		  genePredSingleCover
genePredToBed		      genePredToBigGenePred	    genePredToFakePsl		  genePredToGtf
genePredToMafFrames	      genePredToProt		    genePredToPsl		  geoMirrorNode
getChroms		      getFeatDna		    getRna			  getRnaPred
gfClient		      gff3ToGenePred		    gff3ToPsl			  gffPeek
gfPcr			      gfServer			    git-reports			  gmtime
gpcrParser		      gpStats			    gpToGtf			  groupSamples
gsBig			      gtfToGenePred		    hapmapPhaseIIISummary	  headRest
hgAddLiftOverChain	      hgAvidShortBed		    hgBbiDbLink			  hgBioCyc
hgCeOrfToGene		      hgCGAP			    hgChroms			  hgClonePos
hgClusterGenes		      hgCountAlign		    hgCtgPos			  hgDeleteChrom
hgDropSplitTable	      hgEmblProtLinks		    hgExonerate			  hgExperiment
hgExtFileCheck		      hgFakeAgp			    hgFiberglass		  hgFindSpec
hgFlyBase		      hgGcPercent		    hgGeneBands			  hgGenericMicroarray
hgGetAnn		      hgGnfMicroarray		    hgGoAssociation		  hgGoldGapGl
hgGtex			      hgGtexAse			    hgGtexExonBed		  hgGtexGeneBed
hgJaxQtl		      hgKegg			    hgKegg2			  hgKegg3
hgKgGetText		      hgKgMrna			    hgKnownGeneList		  hgKnownMore
hgKnownToSuper		      hgLoadBed			    hgLoadBlastTab		  hgLoadChain
hgLoadChromGraph	      hgLoadEranModules		    hgLoadGap			  hgLoadGenePred
hgLoadItemAttr		      hgLoadMaf			    hgLoadMafFrames		  hgLoadMafSummary
hgLoadNet		      hgLoadNetDist		    hgLoadOut			  hgLoadOutJoined
hgLoadPsl		      hgLoadRnaFold		    hgLoadSample		  hgLoadSeq
hgLoadSqlTab		      hgLoadWiggle		    hgLsSnpPdbLoad		  hgMapMicroarray
hgMapToGene		      hgMapViaSwissProt		    hgMaxExp			  hgMedianMicroarray
hgMrnaRefseq		      hgNearTest		    hgNetDist			  hgNibSeq
hgPar			      hgPepPred			    hgPhMouse			  hgProtIdToGenePred
hgRatioMicroarray	      hgRenameSplitTable	    hgRnaGenes			  hgSanger20
hgSanger22		      hgSelect			    hgSgdGff3			  hgSgdGfp
hgSgdPep		      hgSoftberryHom		    hgSoftPromoter		  hgSpeciesRna
hgsql			      hgsqladmin		    hgsqldump			  hgsqlimport
hgsqlSwapTables		      hgsqlTableDate		    hgStanfordMicroarray	  hgStsAlias
hgStsMarkers		      hgSuperfam		    hgTablesTest		  hgTomRough
hgTpf			      hgTraceInfo		    hgTrackDb			  hgTracksRandom
hgvsToVcf		      hgWaba			    hgWiggle			  hgWormLinks
hgYeastRegCode		      hicInfo			    hprdP2p			  htmlCheck
htmlPics		      hubCheck			    hubClone			  hubCrawl
hubPublicCheck		      intronEnds		    iriToControlTable		  iriToDnaMotif
isPcr			      jkUniq			    joinableFields		  joinerCheck
joinerRoute		      knownToHprd		    knownToVisiGene		  knownVsBlat
kvsSummary		      lavToAxt			    lavToPsl			  ldHgGene
lfsOverlap		      libScan			    liftAcross			  liftAgp
liftFrags		      liftOver			    liftOverMerge		  liftPromoHits
liftUp			      lineCount			    linesToRa			  localtime
mafAddIRows		      mafAddQRows		    mafCoverage			  mafFetch
mafFilter		      mafFrag			    mafFrags			  mafGene
mafMeFirst		      mafNoAlign		    mafOrder			  mafRanges
mafsInRegion		      mafSpeciesList		    mafSpeciesSubset		  mafSplit
mafSplitPos		      mafToAxt			    mafToBigMaf			  mafToPsl
mafToSnpBed		      makeTableDescriptions	    makeTableList		  makeTrackIndex
maskOutFa		      matrixClusterColumns	    matrixMarketToTsv		  matrixNormalize
matrixToBarChartBed	      maxTranscriptomeExps	    mdbPrint			  mdbUpdate
mdToNcbiLift		      mgcFastaForBed		    mktime			  motifLogo
motifSig		      mousePoster		    mrnaToGene			  mysqlSecurityCheck
netChainSubset		      netClass			    netFilter			  netSplit
netStats		      netSyntenic		    netToAxt			  netToBed
netToBedWithId		      newProg			    newPythonProg		  nibbImageProbes
nibbNameFix		      nibbParseImageDir		    nibbPrepImages		  nibFrag
nibSize			      normalizeSampleFile	    nt4Frag			  oligoMatch
orf			      orfStats			    orthoEvaluate		  orthologBySynteny
orthoMap		      orthoPickIntron		    orthoSplice			  overlapSelect
paraFetch		      paraSync			    patCount			  pepPredToFa
phToPsl			      phyloRenameAndPrune	    polyInfo			  positionalTblCheck
promoSeqFromCluster	      pslCat			    pslCDnaFilter		  pslCheck
pslCoverage		      pslDiff			    pslDropOverlap		  pslFilter
pslFilterPrimers	      pslFixCdsJoinGap		    pslGlue			  pslHisto
pslHitPercent		      pslIntronsOnly		    pslLiftSubrangeBlat		  pslMap
pslMapPostChain		      pslMismatchGapToBed	    pslMrnaCover		  pslPairs
pslPartition		      pslPosTarget		    pslPretty			  pslProtToRnaCoords
pslQuickFilter		      pslRc			    pslRecalcMatch		  pslRemoveFrameShifts
pslReps			      pslScore			    pslSelect			  pslSimp
pslSomeRecords		      pslSort			    pslSortAcc			  pslSplitOnTarget
pslStats		      pslSwap			    pslToBed			  pslToBigPsl
pslToChain		      pslToPslx			    pslToXa			  pslUniq
pslUnpile		      pslxToFa			    qacAgpLift			  qacToQa
qacToWig		      qaToQac			    randomLines			  raSqlQuery
raToCds			      raToLines			    raToTab			  refreshNamedSessionCustomTrac
refSeqGet		      regionPicker		    relPairs			  reviewIndexes
rikenBestInCluster	      rmFaDups			    rmskAlignToPsl		  rnaFoldBig
rowsToCols		      safePush			    samHit			  sanger22gtf
scaffoldFaToAgp		      scaleSampleFiles		    scanRa			  scrambleFa
semiNorm		      seqCheck			    sequenceForBed		  sim4big
simpleChain		      sizeof			    snpException		  snpMaskAddInsertions
snpMaskCutDeletions	      snpMaskSingle		    snpNcbiToUcsc		  snpValid
spacedToTab		      spideyToPsl		    splitFa			  splitFaIntoContigs
splitFile		      splitFileByColumn		    splitSim			  sqlToXml
stageMultiz		      stanToBedAndExpRecs	    strexCalc			  stringify
subChar			      subColumn			    subs			  subsetAxt
subsetTraces		      tableSum			    tailLines			  tdbQuery
tdbRename		      tdbSort			    testSearch			  textHist2
textHistogram		      tfbsConsLoc		    tfbsConsSort		  tickToDate
timePosTable		      toDev64			    toLower			  toUpper
trackDbIndexBb		      trackDbPatch		    trackDbRaFormat		  trackDbToTxt
trackOverlap		      transMapPslToGenePred	    trfBig			  twinOrf
twinOrf2		      twinOrf3			    twinOrfStats		  twinOrfStats2
twinOrfStats3		      twoBitDup			    twoBitInfo			  twoBitMask
twoBitToFa		      txAbFragFind		    txBedToGraph		  txCdsBadBed
txCdsCluster		      txCdsEvFromBed		    txCdsEvFromBorf		  txCdsEvFromProtein
txCdsEvFromRna		      txCdsGoodBed		    txCdsOrfInfo		  txCdsOrtho
txCdsPick		      txCdsPredict		    txCdsRaExceptions		  txCdsRefBestEvOnly
txCdsRepick		      txCdsSuspect		    txCdsSvmInput		  txCdsToGene
txCdsWeed		      txgAddEvidence		    txgAnalyze			  txGeneAccession
txGeneAlias		      txGeneAltProt		    txGeneCanonical		  txGeneCdsMap
txGeneColor		      txGeneExplainUpdate1	    txGeneFromBed		  txGeneProtAndRna
txGeneSeparateNoncoding	      txGeneXref		    txgGoodEdges		  txgToAgx
txgToXml		      txgTrim			    txInfoAssemble		  txOrtho
txPslFilter		      txPslToBed		    txReadRa			  txWalk
ucscApiClient		      udcCleanup		    undupFa			  uniqSize
updateStsInfo		      upper			    utrFa			  validateFiles
validateManifest	      varStepToBedGraph.pl	    vcfFilter			  vcfRenameAndPrune
vcfToBed		      vcfToHgvs			    vegaBuildInfo		  venn
verticalSplitSqlTable	      webSync			    weedLines			  whyConserved
wigBedToStep		      wigCorrelate		    wigEncode			  wigToBigWig
wordLine		      xmlCat			    xmlToSql

```

Documentation
-------------


For most of the executables, issuing the command without any arguments will give a brief description of the 
executable.


For example:



```

$ twoBitToFa 
twoBitToFa - Convert all or part of .2bit file to fasta
usage:
   twoBitToFa input.2bit output.fa
options:
   -seq=name - restrict this to just one sequence
   -start=X  - start at given position in sequence (zero-based)
   -end=X - end at given position in sequence (non-inclusive)
   -seqList=file - file containing list of the desired sequence names 
                    in the format seqSpec[:start-end], e.g. chr1 or chr1:0-189
                    where coordinates are half-open zero-based, i.e. [start,end)
   -noMask - convert sequence to all upper case
   -bpt=index.bpt - use bpt index instead of built in one
   -bed=input.bed - grab sequences specified by input.bed. Will exclude introns
   -bedPos        - with -bed, to use chrom:start-end as the fasta ID in output.fa

Sequence and range may also be specified as part of the input
file name using the syntax:
      /path/input.2bit:name
   or
      /path/input.2bit:name
   or
      /path/input.2bit:name:start-end

```

Many of the executables are poorly documented. For further help, you can submit questions to the main
UCSC discussion list. See <http://genome.ucsc.edu/contacts.html>.




