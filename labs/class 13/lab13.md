Lab 13: RNA-seq analysis
================
Torrey Rhyne (A14397504)

Setup:

``` r
library(DESeq2)
```

## (3) Import countData and colData

Read csv:

``` r
# Complete the missing code
counts <- read.csv("../airway_scaledcounts.csv", row.names=1)
metadata <-  read.csv("../airway_metadata.csv")
```

View:

``` r
head(counts)
```

                    SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
    ENSG00000000003        723        486        904        445       1170
    ENSG00000000005          0          0          0          0          0
    ENSG00000000419        467        523        616        371        582
    ENSG00000000457        347        258        364        237        318
    ENSG00000000460         96         81         73         66        118
    ENSG00000000938          0          0          1          0          2
                    SRR1039517 SRR1039520 SRR1039521
    ENSG00000000003       1097        806        604
    ENSG00000000005          0          0          0
    ENSG00000000419        781        417        509
    ENSG00000000457        447        330        324
    ENSG00000000460         94        102         74
    ENSG00000000938          0          0          0

``` r
head(metadata)
```

              id     dex celltype     geo_id
    1 SRR1039508 control   N61311 GSM1275862
    2 SRR1039509 treated   N61311 GSM1275863
    3 SRR1039512 control  N052611 GSM1275866
    4 SRR1039513 treated  N052611 GSM1275867
    5 SRR1039516 control  N080611 GSM1275870
    6 SRR1039517 treated  N080611 GSM1275871

How many genes are in this dataset?

``` r
nrow(counts)
```

    [1] 38694

How many ‘control’ cell lines do we have?

``` r
table(metadata$dex)
```


    control treated 
          4       4 

## (4) Toy differential gene expression

Lets perform some exploratory differential gene expression analysis.
Note: this analysis is for demonstration only. NEVER do differential
expression analysis this way!

Step 1: Identify and extract “control” columns. Step 2: Calculate the
mean value per gene for all these “control” columns. Step 3: Do the same
for “treated”. Step 4: Compare the “control.mean” and “treated.mean”
values.

Step 1:

``` r
control <- metadata[metadata[,"dex"]=="control",]
```

Step 2:

``` r
control.counts <- counts[ ,control$id]
control.mean <- rowSums( control.counts )/4 
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

Alternative with dplyr:

``` r
suppressMessages(library(dplyr))
control <- metadata %>% filter(dex=="control")
control.counts <- counts %>% select(control$id) 
control.mean <- rowSums(control.counts)/4
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

Step 3:

``` r
suppressMessages(library(dplyr))

treated <- metadata %>% filter(dex=="treated")
treated.counts <- counts %>% select(treated$id) 
treated.mean <- rowMeans(treated.counts) # more robust
head(treated.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             658.00            0.00          546.00          316.50           78.75 
    ENSG00000000938 
               0.00 

``` r
# Barry's way (combine)
treated.mean <- rowMeans(counts[,metadata$dex == "treated"])
head(treated.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             658.00            0.00          546.00          316.50           78.75 
    ENSG00000000938 
               0.00 

Step 4 (combine into one dataframe and plot):

``` r
mean.counts <- data.frame(control.mean, treated.mean)
plot(mean.counts)
```

![](lab13_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
suppressMessages(library(ggplot2))

ggplot(mean.counts) + 
  aes(control.mean, treated.mean) +
  geom_point(alpha = 0.2)
```

![](lab13_files/figure-commonmark/unnamed-chunk-11-1.png)

Skewed data. Do a log transformation:

``` r
plot(mean.counts, log = "xy")
```

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0 omitted
    from logarithmic plot

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0 omitted
    from logarithmic plot

![](lab13_files/figure-commonmark/unnamed-chunk-12-1.png)

Logs are super useful when we have such skewed data and/or are
interested in big changes (orders of magnitude).

``` r
# treated / control
10/10
```

    [1] 1

``` r
log2(10/10) # useful because log1 = 0 ("no change" - like in volcano plots)
```

    [1] 0

``` r
log2(40/10)
```

    [1] 2

Add log2 fold-change column to mean.counts data:

``` r
suppressMessages(library(dplyr))

mean.counts <- mean.counts %>% mutate("log2fc" = log2(treated.mean/control.mean))
head(mean.counts)
```

                    control.mean treated.mean      log2fc
    ENSG00000000003       900.75       658.00 -0.45303916
    ENSG00000000005         0.00         0.00         NaN
    ENSG00000000419       520.50       546.00  0.06900279
    ENSG00000000457       339.75       316.50 -0.10226805
    ENSG00000000460        97.25        78.75 -0.30441833
    ENSG00000000938         0.75         0.00        -Inf

I need to exclude any genes with zero counts. We can’t say anything
about them anyways from this experiment and it causes me math pain.

``` r
# what values in the first two cols = 0?
# remember T = 1, F = 0 (so we want to keep F + f = 0)
keep.inds <- rowSums(mean.counts[,1:2] == 0) == 0
mycounts <- mean.counts[keep.inds,]

# Barry uses which in the notebook
which(T,F,T)
```

    [1] 1

``` r
head(which(mean.counts[,1:2] == 0))
```

    [1]   2  65  70  73 121 123

How many genes do I have left?

``` r
nrow(mycounts)
```

    [1] 21817

How many are up-regulated using a threshold of log2fc = 2?

``` r
sum(mycounts$log2fc > +2)
```

    [1] 250

How many are down-regulated?

``` r
sum(mycounts$log2fc < -2)
```

    [1] 367

Are we satisfied? No. No measure of significance. Use DESeq2. You could
get this much variation by chance within the control samples, for
example.

## (5) settin up DESeq

Like many bioconductor analysis packages, DESeq wants it’s input in a
very particular way.

``` r
dds <- DESeqDataSetFromMatrix(countData = counts, 
                       colData = metadata,
                       design = ~dex)
```

    converting counts to integer mode

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

## (6) PCA

Can use PCA as sort of a QC (should see separation between control and
treated).

## (7) DESeq analysis

To run DESeq analysis we call the main function from the package
`DESeq(dds)`:

``` r
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

To get the results out of the `dds` object we can use the DESeq
`results()` function:

``` r
res <- results(dds)
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 6 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj
                    <numeric>
    ENSG00000000003  0.163035
    ENSG00000000005        NA
    ENSG00000000419  0.176032
    ENSG00000000457  0.961694
    ENSG00000000460  0.815849
    ENSG00000000938        NA

## (9) Data visualization

A common summary visualization is called a volcano plot.

``` r
plot(res$log2FoldChange, res$padj) 
```

![](lab13_files/figure-commonmark/unnamed-chunk-22-1.png)

Flip the y-axis, do a transformation, and add labels:

``` r
plot(res$log2FoldChange, -log(res$padj),
     xlab = "log2 fold-change",
     ylab = "-log agjusted p-value")
abline(v=c(-2,2), col="red")
abline(h=-log(0.05), col="blue")
```

![](lab13_files/figure-commonmark/unnamed-chunk-23-1.png)

Much better! Now add color:

``` r
mycols <- rep("grey",nrow(res))
mycols[res$log2FoldChange > 2 | res$log2FoldChange < -2] <- "black"
mycols[res$padj > 0.05] <- "grey"


plot(res$log2FoldChange, -log(res$padj), col=mycols,
     xlab = "log2 fold-change",
     ylab = "-log agjusted p-value")
abline(v=c(-2,2), col="red")
abline(h=-log(0.05), col="blue")
```

![](lab13_files/figure-commonmark/unnamed-chunk-24-1.png)

Summarize and save results so far:

``` r
res05 <- results(dds, alpha=0.05)
summary(res05)
```


    out of 25258 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 1236, 4.9%
    LFC < 0 (down)     : 933, 3.7%
    outliers [1]       : 142, 0.56%
    low counts [2]     : 9033, 36%
    (mean count < 6)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

## (8) Adding annotation data.

We need to translate or “map” our ensemble IDs into more understandable
gene names and the identifiers that other useful databases use.

Packages:

``` r
library(AnnotationDbi)
```


    Attaching package: 'AnnotationDbi'

    The following object is masked from 'package:dplyr':

        select

``` r
library("org.Hs.eg.db") # data package
```

``` r
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

``` r
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),     # Our gene names
                     keytype="ENSEMBL",       # The format of our gene names
                     column="SYMBOL",         # The new format we want to add
                     multiVals="first")       # Choosing most abundant here
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 7 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj      symbol
                    <numeric> <character>
    ENSG00000000003  0.163035      TSPAN6
    ENSG00000000005        NA        TNMD
    ENSG00000000419  0.176032        DPM1
    ENSG00000000457  0.961694       SCYL3
    ENSG00000000460  0.815849       FIRRM
    ENSG00000000938        NA         FGR

Run the `mapIds()` function three more times to add the Entrez ID,
UniProt accession, and Genename.

``` r
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),     # Our gene names
                     keytype="ENSEMBL",       # The format of our gene names
                     column="ENTREZID",         # The new format we want to add
                     multiVals="first")       # Choosing most abundant here
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$uniprot <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),     # Our gene names
                     keytype="ENSEMBL",       # The format of our gene names
                     column="UNIPROT",         # The new format we want to add
                     multiVals="first")       # Choosing most abundant here
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$genename <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),     # Our gene names
                     keytype="ENSEMBL",       # The format of our gene names
                     column="GENENAME",         # The new format we want to add
                     multiVals="first")       # Choosing most abundant here
```

    'select()' returned 1:many mapping between keys and columns

## Pathway analysis

Packages:

``` r
library(pathview)
```

    ##############################################################################
    Pathview is an open source software package distributed under GNU General
    Public License version 3 (GPLv3). Details of GPLv3 is available at
    http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    formally cite the original Pathview paper (not just mention it) in publications
    or products. For details, do citation("pathview") within R.

    The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ##############################################################################

``` r
library(gage)
```

``` r
library(gageData)
```

Examine the first 2 pathways in this kegg set for humans:

``` r
data(kegg.sets.hs)
head(kegg.sets.hs, 2)
```

    $`hsa00232 Caffeine metabolism`
    [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   

    $`hsa00983 Drug metabolism - other enzymes`
     [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"   "1551"  
     [9] "1553"   "1576"   "1577"   "1806"   "1807"   "1890"   "221223" "2990"  
    [17] "3251"   "3614"   "3615"   "3704"   "51733"  "54490"  "54575"  "54576" 
    [25] "54577"  "54578"  "54579"  "54600"  "54657"  "54658"  "54659"  "54963" 
    [33] "574537" "64816"  "7083"   "7084"   "7172"   "7363"   "7364"   "7365"  
    [41] "7366"   "7367"   "7371"   "7372"   "7378"   "7498"   "79799"  "83549" 
    [49] "8824"   "8833"   "9"      "978"   

Note - higher chance of hitting pathway with more genes involved.

Run gage (need vector of foldchanges first, not whole df):

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

           7105       64102        8813       57147       55732        2268 
    -0.35070302          NA  0.20610777  0.02452695 -0.14714205 -1.73228897 

``` r
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

Examine:

``` r
attributes(keggres)
```

    $names
    [1] "greater" "less"    "stats"  

``` r
head(keggres$less, 3)
```

                                          p.geomean stat.mean        p.val
    hsa05332 Graft-versus-host disease 0.0004250461 -3.473346 0.0004250461
    hsa04940 Type I diabetes mellitus  0.0017820293 -3.002352 0.0017820293
    hsa05310 Asthma                    0.0020045888 -3.009050 0.0020045888
                                            q.val set.size         exp1
    hsa05332 Graft-versus-host disease 0.09053483       40 0.0004250461
    hsa04940 Type I diabetes mellitus  0.14232581       42 0.0017820293
    hsa05310 Asthma                    0.14232581       29 0.0020045888

Lets have a look at one of these pathways

``` r
pathview(gene.data=foldchanges, pathway.id="hsa05310")
```

    Info: Downloading xml files for hsa05310, 1/1 pathways..

    Info: Downloading png files for hsa05310, 1/1 pathways..

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/TorreyRhyne/Desktop/BioSci PhD/classes/BGGN213_bioinformatics/R/labs

    Info: Writing image file hsa05310.pathview.png

![](hsa05310.pathview.png)
