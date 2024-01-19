

document.querySelector('title').textContent = 'ROSE: Rank Ordering of Super-Enhancers ';
**ROSE: Rank Ordering of Super-Enhancers** 


|  |
| --- |
| 
Quick Links
[Documentation](#doc)
[Notes](#notes)
[Interactive job](#int) 
[Batch job](#sbatch) 
[Swarm of jobs](#swarm) 
 |


  


ROSE (Rank Ordering of Super-Enhancers) is a tool for   

(1) creating stitched enhancers, and   

(2) separating super-enhancers from typical enhancers.   

given sequencing data (.bam) and a file of previously identified constituent enhancers (.gff)



  

### References:


* Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young   

*Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes*  

[Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003929) 2013, **153** : 307-319   
* Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vako
c, James E. Bradner, Tong Ihn Lee, and Richard A. Young   

*Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers*  

[Cell](https://www.sciencedirect.com/science/article/pii/S0092867413003930) 2013, **153**: 320-334


Documentation
* [ROSE on Github](https://github.com/stjude/ROSE)


Important Notes
* Module Name: ROSE (see [the modules page](https://hpc.nih.gov/apps/modules.html) for more information)
* Unusual environment variables set
	+ **ROSE\_HOME**  ROSE installation directory
	+ **ROSE\_BIN**       ROSE bin directory
	+ **ROSE\_SRC**       ROSE source directory
	+ **ROSE\_DATA**    ROSE sample data directory



Interactive job
[Interactive jobs](/docs/userguide.html#int) should be used for debugging, graphics, or applications that cannot be run as batch jobs.
Allocate an [interactive session](/docs/userguide.html#int) and run the program. Sample session:



```

[user@biowulf]$ **sinteractive** 
[user@cn3107 ~]$ **module load ROSE**
[+] Loading gcc  7.2.0  ...
[+] Loading GSL 2.4 for GCC 7.2.0 ...
[+] Loading openmpi 3.0.0  for GCC 7.2.0
[+] Loading pandoc  2.10  on cn2385
[+] Loading R 3.4.4
[+] Loading samtools 0.1.19  ...
[+] Loading ROSE  20200707

```

Copy sample data to the current folder:

```

[user@cn3107 ~]$ **cp -rP $ROSE\_DATA/\* .** 

```

The ROSE usage is as follows:

```

[user@cn3107 ~]$ **ROSE\_main.py -h** 
Usage: ROSE_main.py [options] -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]

Options:
  -h, --help            show this help message and exit
  -i INPUT, --i=INPUT   Enter a .gff or .bed file of binding sites used to
                        make enhancers
  -r RANKBY, --rankby=RANKBY
                        bamfile to rank enhancer by
  -o OUT, --out=OUT     Enter an output folder
  -g GENOME, --genome=GENOME
                        Enter the genome build (MM9,MM8,HG18,HG19,MM10,HG38)
  -b BAMS, --bams=BAMS  Enter a comma separated list of additional bam files
                        to map to
  -c CONTROL, --control=CONTROL
                        bamfile to rank enhancer by
  -s STITCH, --stitch=STITCH
                        Enter a max linking distance for stitching
  -t TSS, --tss=TSS     Enter a distance from TSS to exclude. 0 = no TSS
                        exclusion

```

Run ROSE on the sample data:

```

[user@cn3107 ~]$ **ROSE\_main.py -g hg18 -i HG18\_MM1S\_MED1.gff -r MM1S\_MED1.hg18.bwt.sorted.bam -o my\_output\_dir**
USING HG18_MM1S_MED1.gff AS THE INPUT GFF
USING hg18_refseq.ucsc AS THE GENOME
('upper(genome)=', 'HG18_REFSEQ.UCSC')
[user@cn2380 ROSE]$ vi ROSE_main.py
[usera@cn2380 ROSE]$ python ROSE_main.py -g hg18 -i HG18_MM1S_MED1.gff -r MM1S_MED1.hg18.bwt.sorted.bam -o my_output_dir
USING HG18_MM1S_MED1.gff AS THE INPUT GFF
USING hg18 AS THE GENOME
('upper(genome)=', 'HG18')
MAKING START DICT
LOADING IN GFF REGIONS
CHECKING INPUT TO MAKE SURE EACH REGION HAS A UNIQUE IDENTIFIER
REFERENCE COLLECTION PASSES QC
STITCHING REGIONS TOGETHER
PERFORMING REGION STITCHING
MAKING GFF FROM STITCHED COLLECTION
WRITING STITCHED GFF TO DISK AS my_output_dir/gff/HG18_MM1S_MED1_12KB_STITCHED.gff
OUTPUT WILL BE WRITTEN TO  my_output_dir/HG18_MM1S_MED1_12KB_STITCHED_ENHANCER_REGION_MAP.txt
python ROSE_bamToGFF.py -f 1 -e 200 -r -m 1 -b MM1S_MED1.hg18.bwt.sorted.bam -i my_output_dir/gff/HG18_MM1S_MED1_12KB_STITCHED.gff
-o my_output_dir/mappedGFF/HG18_MM1S_MED1_12KB_STITCHED_MM1S_MED1.hg18.bwt.sorted.bam_MAPPED.gff &
python ROSE_bamToGFF.py -f 1 -e 200 -r -m 1 -b MM1S_MED1.hg18.bwt.sorted.bam -i HG18_MM1S_MED1.gff -o my_output_dir/mappedGFF/HG18_
MM1S_MED1_MM1S_MED1.hg18.bwt.sorted.bam_MAPPED.gff &
PAUSING TO MAP
{'matrix': '1', 'extension': '200', 'floor': '1', 'sense': 'both', 'output': 'my_output_dir/mappedGFF/HG18_MM1S_MED1_12KB_STITCHED_
MM1S_MED1.hg18.bwt.sorted.bam_MAPPED.gff', 'bam': 'MM1S_MED1.hg18.bwt.sorted.bam', 'rpm': True, 'input': 'my_output_dir/gff/HG18_MM
1S_MED1_12KB_STITCHED.gff'}
[]
mapping to GFF and making a matrix with fixed bin number
{'matrix': '1', 'extension': '200', 'floor': '1', 'sense': 'both', 'output': 'my_output_dir/mappedGFF/HG18_MM1S_MED1_MM1S_MED1.hg18
.bwt.sorted.bam_MAPPED.gff', 'bam': 'MM1S_MED1.hg18.bwt.sorted.bam', 'rpm': True, 'input': 'HG18_MM1S_MED1.gff'}
[]
mapping to GFF and making a matrix with fixed bin number
WAITING FOR MAPPING TO COMPLETE. ELAPSED TIME (MIN):
0
using a MMR value of 17.4141
using a MMR value of 17.4141
has chr
has chr
Number lines processed
0
Number lines processed
0
100
100
200
200
...
19300
19400
MAPPING TOOK 30 MINUTES
BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS
FORMATTING TABLE
1000
2000
...
13000
14000
GETTING MAPPED DATA
GETTING MAPPING DATA FOR  MM1S_MED1.hg18.bwt.sorted.bam
OPENING my_output_dir/mappedGFF/HG18_MM1S_MED1_12KB_STITCHED_MM1S_MED1.hg18.bwt.sorted.bam_MAPPED.gff
MAKING SIGNAL DICT FOR MM1S_MED1.hg18.bwt.sorted.bam
CALLING AND PLOTTING SUPER-ENHANCERS
R --no-save my_output_dir/ my_output_dir/HG18_MM1S_MED1_12KB_STITCHED_ENHANCER_REGION_MAP.txt HG18_MM1S_MED1 NONE < ROSE_callSuper.
R
...
 version 3.4.4 (2018-03-15) -- "Someone to Lean On"
opyright (C) 2018 The R Foundation for Statistical Computing
latform: x86_64-pc-linux-gnu (64-bit)

 is free software and comes with ABSOLUTELY NO WARRANTY.
ou are welcome to redistribute it under certain conditions.
ype 'license()' or 'licence()' for distribution details.

 Natural language support but running in an English locale

 is a collaborative project with many contributors.
ype 'contributors()' for more information and
citation()' on how to cite R or R packages in publications.

ype 'demo()' for some demos, 'help()' for on-line help, or
help.start()' for an HTML browser interface to help.
ype 'q()' to quit R.

 #============================================================================
 #==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
 #============================================================================
>
>  #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
>  calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
+       inputVector <- sort(inputVector)
+       inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
+       slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is t
he diagonal.
+       xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find
 the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
+       y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
+
+       if(drawPlot){  #if TRUE, draw the plot
+               plot(1:length(inputVector), inputVector,type="l",...)
+               b <- y_cutoff-(slope* xPt)
+               abline(v= xPt,h= y_cutoff,lty=2,col=8)
+               points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
+               abline(coef=c(b,slope),col=2)
+               title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFol
d over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
+               axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
+       }
+       return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
+ }
>
> #this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
> numPts_below_line <- function(myVector,slope,x){
+       yPt <- myVector[x]
+       b <- yPt-(slope*x)
+       xPts <- 1:length(myVector)
+       return(sum(myVector<=(xPts*slope+b)))
+ }
>
>
> convert_stitched_to_bed <- function(inputStitched,trackName,trackDescription,outputFile,splitSuper=TRUE,score=c(),superRows=c(),b
aseColor="0,0,0",superColor="255,0,0"){
+       outMatrix <- matrix(data="",ncol=4+ifelse(length(score)==nrow(inputStitched),1,0),nrow=nrow(inputStitched))
+
+       outMatrix[,1] <- as.character(inputStitched$CHROM)
+       outMatrix[,2] <- as.character(inputStitched$START)
+       outMatrix[,3] <- as.character(inputStitched$STOP)
+       outMatrix[,4] <- as.character(inputStitched$REGION_ID)
+       if(length(score)==nrow(inputStitched)){
+               score <- rank(score,ties.method="first")
+               score <- length(score)-score+1  #Stupid rank only does smallest to largest.
+               outMatrix[,5] <- as.character(score)
+       }
+       trackDescription <- paste(trackDescription,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
+       trackDescription <- gsub("\n","\t", trackDescription)
+       tName <- gsub(" ","_",trackName)
+       cat('track name="', tName,'" description="', trackDescription,'" itemRGB=On color=',baseColor,"\n",sep="",file=outputFile)
+       write.table(file= outputFile,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
+       if(splitSuper==TRUE){
+               cat("\ntrack name=\"Super_", tName,'" description="Super ', trackDescription,'" itemRGB=On color=', superColor,"\n"
,sep="",file=outputFile,append=TRUE)
+               write.table(file= outputFile,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE
)
+       }
+ }
>
>
>
> convert_stitched_to_gateway_bed <- function(inputStitched,outputFileRoot,splitSuper=TRUE,score=c(),superRows=c()){
+       outMatrix <- matrix(data="",ncol=6,nrow=nrow(inputStitched))
+
+       outMatrix[,1] <- as.character(inputStitched$CHROM)
+       outMatrix[,2] <- as.character(inputStitched$START)
+       outMatrix[,3] <- as.character(inputStitched$STOP)
+       outMatrix[,4] <- as.character(inputStitched$REGION_ID)
+
+       if(length(score)==nrow(inputStitched)){
+               score <- rank(score,ties.method="first")
+               score <- length(score)-score+1  #Stupid rank only does smallest to largest.
+               outMatrix[,5] <- as.character(score)
+       }
+
+       outMatrix[,6] <- as.character(rep('.',nrow(outMatrix)))
+
+
+       outputFile1 = paste(outputFileRoot,'_Gateway_Enhancers.bed',sep='')
+       write.table(file= outputFile1,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
+       if(splitSuper==TRUE){
+               outputFile2 = paste(outputFileRoot,'_Gateway_SuperEnhancers.bed',sep='')
+
+               write.table(file= outputFile2,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRU
E)
+       }
+ }
>
>
>
>
> writeSuperEnhancer_table <- function(superEnhancer,description,outputFile,additionalData=NA){
+       description <- paste("#",description,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
+       description <- gsub("\n","\n#",description)
+       cat(description,"\n",file=outputFile)
+       if(is.matrix(additionalData)){
+               if(nrow(additionalData)!=nrow(superEnhancer)){
+                       warning("Additional data does not have the same number of rows as the number of super enhancers.\n--->>> AD
DITIONAL DATA NOT INCLUDED <<<---\n")
+               }else{
+                       superEnhancer <- cbind(superEnhancer,additionalData)
+                       superEnhancer = superEnhancer[order(superEnhancer$enhancerRank),]
+
+               }
+       }
+       write.table(file=outputFile,superEnhancer,sep="\t",quote=FALSE,row.names=FALSE,append=TRUE)
+ }
>
>
>
> #============================================================================
> #===================SUPER-ENHANCER CALLING AND PLOTTING======================
> #============================================================================
>
>
>
> args <- commandArgs()
>
> print('THESE ARE THE ARGUMENTS')
[1] "THESE ARE THE ARGUMENTS"
> print(args)
[1] "/usr/local/apps/R/3.4/3.4.4/lib64/R/bin/exec/R"
[2] "--no-save"
[3] "my_output_dir/"
[4] "my_output_dir/HG18_MM1S_MED1_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
[5] "HG18_MM1S_MED1"
[6] "NONE"
>
> #ARGS
> outFolder = args[3]
> enhancerFile = args[4]
> enhancerName = args[5]
> wceName = args[6]
>
>
>
> #Read enhancer regions with closestGene columns
> stitched_regions <- read.delim(file= enhancerFile,sep="\t")
>
>
> #perform WCE subtraction. Using pipeline table to match samples to proper background.
> rankBy_factor = colnames(stitched_regions)[7]
>
> prefix = unlist(strsplit(rankBy_factor,'_'))[1]
>
> if(wceName == 'NONE'){
+
+
+       rankBy_vector = as.numeric(stitched_regions[,7])
+
+ }else{
+       wceName = colnames(stitched_regions)[8]
+       print('HERE IS THE WCE NAME')
+       print(wceName)
+
+       rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])
+ }
>
>
> #SETTING NEGATIVE VALUES IN THE rankBy_vector to 0
>
> rankBy_vector[rankBy_vector < 0] <- 0
>
>
> #FIGURING OUT THE CUTOFF
>
> cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankByFactor,'
 Signal','- ',wceName),lwd=2,col=4)
>
>
> #These are the super-enhancers
> superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
> typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)
> enhancerDescription <- paste(enhancerName," Enhancers\nCreated from ", enhancerFile,"\nRanked by ",rankBy_factor,"\nUsing cutoff
of ",cutoff_options$absolute," for Super-Enhancers",sep="",collapse="")
>
>
> #MAKING HOCKEY STICK PLOT
> plotFileName = paste(outFolder,enhancerName,'_Plot_points.png',sep='')
> png(filename=plotFileName,height=600,width=600)
> signalOrder = order(rankBy_vector,decreasing=TRUE)
> if(wceName == 'NONE'){
+       plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy
_factor,' Signal'),pch=19,cex=2)
+
+ }else{
+       plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy
_factor,' Signal','- ',wceName),pch=19,cex=2)
+ }
> abline(h=cutoff_options$absolute,col='grey',lty=2)
> abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
> lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')
> text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnha
ncerRows)),pos=4)
>
> dev.off()
null device
          1
>
>
>
>
> #Writing a bed file
> bedFileName = paste(outFolder,enhancerName,'_Enhancers_withSuper.bed',sep='')
> convert_stitched_to_bed(stitched_regions,paste(rankBy_factor,"Enhancers"), enhancerDescription,bedFileName,score=rankBy_vector,sp
litSuper=TRUE,superRows= superEnhancerRows,baseColor="0,0,0",superColor="255,0,0")
>
>
> bedFileRoot = paste(outFolder,enhancerName,sep='')
>
> convert_stitched_to_gateway_bed(stitched_regions,bedFileRoot,splitSuper=TRUE,score=rankBy_vector,superRows = superEnhancerRows)
>
>
> #This matrix is just the super_enhancers
> true_super_enhancers <- stitched_regions[superEnhancerRows,]
>
> additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
> colnames(additionalTableData) <- c("enhancerRank","isSuper")
> additionalTableData[,1] <- nrow(stitched_regions)-rank(rankBy_vector,ties.method="first")+1
> additionalTableData[,2] <- 0
> additionalTableData[superEnhancerRows,2] <- 1
>
>
> #Writing enhancer and super-enhancer tables with enhancers ranked and super status annotated
> enhancerTableFile = paste(outFolder,enhancerName,'_AllEnhancers.table.txt',sep='')
> writeSuperEnhancer_table(stitched_regions, enhancerDescription,enhancerTableFile, additionalData= additionalTableData)
Warning message:
In write.table(file = outputFile, superEnhancer, sep = "\t", quote = FALSE,  :
  appending column names to file
>
> superTableFile = paste(outFolder,enhancerName,'_SuperEnhancers.table.txt',sep='')
> writeSuperEnhancer_table(true_super_enhancers, enhancerDescription,superTableFile, additionalData= additionalTableData[superEnhan
cerRows,])
Warning message:
In write.table(file = outputFile, superEnhancer, sep = "\t", quote = FALSE,  :
  appending column names to file

```

End the interactive session:

```

[user@cn3107 ~]$ **exit**
salloc.exe: Relinquishing job allocation 46116226
[user@biowulf ~]$

```





