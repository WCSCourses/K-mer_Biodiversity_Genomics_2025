# 3. What is in my metagenome? Metagenomic profiling via FracMinHash

[HackMD version](https://hackmd.io/@bluegenes/HybXTfFzlg)

**This follows from the [Introduction to FracMinHash sketching for sequence comparison](https://hackmd.io/@bluegenes/H1ItsC_fxx) and [Comparing genomes metagenomes using FracMinHash sketches](https://hackmd.io/@bluegenes/Byh_JfKGgg) modules. Please refer to the installation instructions in the introduction for installing sourmash and related dependencies.**

With metagenomic sequencing, a first step is often to determine what organisms are present in the sequenced community.

Here, we demonstrate comprehensive metagenome compositional analysis using [sourmash](sourmash.readthedocs.io) FracMinHash sketching, and aggregation to taxonomic groupings. This analysis is entirely reference-dependent. The quality of your results will depend in part on the quality and completeness of the reference database.

## Before we begin:
Do you have sourmash installed?

`sourmash --version`

> If not, run `conda activate sourmash` (if using conda and you've previously created the environment) or go back to the installation instructions in the first part.

## Genome-based compositional analysis of a metagenome

We're going to use the same metagenome used in part 2. Instead of searching with a single genome, we're going to compare the metagenome against a comprehensive list of reference genomes using the "minimum metagenome cover approach". This approach uses two steps internally: first, search all reference genomes for k-mer overlaps, and second, use a greedy method to build the smallest list of reference genomes that contain all shared (known) k-mer content.

We will start by searching only bacterial and archaeal reference genomes, as the database is much smaller. If we have time at the end, we will also search a database the includes eukaryotic genomes, since we can provide that database in a shared location without requiring donwload.

In sourmash, we can use the `gather` command to run this analysis.

Datasets we will use:
  - **Metagenome sample**: [SRR29654720](https://www.ncbi.nlm.nih.gov/sra/?term=SRR29654720) from the [BioProject PRJNA1127303](https://www.ncbi.nlm.nih.gov/bioproject/1127303).
  - **Associated paper**: [Risk assessment of antibiotic resistance genes in rhizosphere soil during tomato growth under bio-control bacterial inoculation](https://www.sciencedirect.com/science/article/pii/S0959652625002616)
    > The study aims to investigate the effects of different fertilization practices during various growth stages of tomatoes on the abundance and dynamics of antibiotic resistance genes (ARGs) in the soil. Understanding the mechanisms behind these impacts is critical, as ARGs pose a significant threat to public health by promoting the development and spread of antibiotic-resistant bacteria. By examining the influence of organic, inorganic, and combined fertilization methods on ARG profiles, as well as considering factors such as microbial community composition and soil physicochemical properties, this research seeks to elucidate how fertilization strategies might mitigate or exacerbate the propagation of ARGs in agricultural environments. The findings from this study will provide valuable insights for developing sustainable farming practices that minimize the risk of ARG proliferation, thereby contributing to the broader effort of combating antibiotic resistance.


### Find your metagenome signature

You should have this file available locally, as we downloaded it in part 2.
```
ls SRR29654720.sig.zip
```
> If you don't see this file, you can re-download by uncommenting and running this command:
> ```
> # curl -JLO https://github.com/WCSCourses/K-mer_Biodiversity_Genomics_2025/raw/refs/heads/main/course_data_2025/module8_sketching/SRR29654720.sig.zip
> ```

Examine this file with `sourmash`:
```
sourmash sig describe SRR29654720.sig.zip
```

You should see:
```
signature filename: SRR29654720.sig.zip
signature: SRR29654720
source file: SRR29654720_2.fastq.gz
md5: 57e826e874fe03b8e4438b19d64b8652
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=1
size: 2714833
sum hashes: 8461460
signature license: CC0

---
signature filename: SRR29654720.sig.zip
signature: SRR29654720
source file: SRR29654720_2.fastq.gz
md5: 4781783025b1e0bf76b06171631c897f
k=51 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=1
size: 2561376
sum hashes: 7030076
signature license: CC0

loaded 2 signatures total, from 1 files
```
> This file contains k=31, k=51 sourmash sketches of SRA run SRR29654720

### Find and examine the database

We will be searching against Genome Taxonomy Database (GTDB) for refseq release 226 (Sept 2024). This release contains 732k bacterial and archael genomes from 143k species.

![image](https://hackmd.io/_uploads/HyLZcXYGle.png)

**This sourmash database is constructed at k=31, scaled=10,000**, meaning we keep roughly 1/10000 k-mers per genome. We're using this low resolution to keep the database small and the analysis fast during the workshop. Higher resolution databases are available for follow-up analyses.

This database can take a little time to download, so we started the download in the lecture portion of the course.

You should now see a compressed version of the database available:
```
ls gtdb-rs226.k31-sc10k.sig.rocksdb.tar.gz
```
> The download command was: 
> ```
> #curl -JLO https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226.k31-sc10k.sig.rocksdb.tar.gz
> ```

Now uncompress and examine the database:
```
tar xvf gtdb-rs226.k31-sc10k.sig.rocksdb.tar.gz
sourmash sig summarize gtdb-rs226.k31-sc10k.sig.rocksdb
```
You should see:

```
** loading from 'gtdb-rs226.k31-sc10k.sig.rocksdb'
path filetype: DiskRevIndex
location: gtdb-rs226.k31-sc10k.sig.rocksdb
is database? yes
has manifest? yes
num signatures: 732034
** examining manifest...
total hashes: 256986256
summary of sketches:
   732034 sketches with DNA, k=31, scaled=10000       256986256 total hashes
```

This database is a 'rocksdb' index, a sourmash format that maps k-mers (hashes) as keys to their datasets (values) and is optimized for fast on-disk search using a single thread. This is the same type of database that supports fast search across SRA metagenomes on [branchwater.jgi.doe.gov](https://branchwater.jgi.doe.gov).

:::spoiler **If you can't download this database, follow these instructions**
::: warning
If you're unable to download the full database due to space or internet connection, you can use a tiny database where I've extracted only signatures that matched our metagenome sample. **If you're on the GitPod, this database should be pre-loaded.**

Here are the commands for that:

download:
```
curl -JLO https://github.com/WCSCourses/K-mer_Biodiversity_Genomics_2025/raw/refs/heads/main/course_data_2025/module8_sketching/gtdb-rs226-SRR29654720subset.k31-sc10k.rocksdb.tar.gz
```
uncompress:
```
tar xvf gtdb-rs226-SRR29654720subset.k31-sc10k.rocksdb.tar.gz
```

examine:
```
sourmash sig summarize gtdb-rs226-SRR29654720subset.k31-sc10k.rocksdb.tar.gz
```

run gather:
```
sourmash gather SRR29654720.sig.zip gtdb-rs226-SRR29654720subset.k31-sc10k.rocksdb.tar.gz \
                -k 31 --scaled 10_000 --threshold-bp 3_000_000 \
                -o SRR29654720-x-gtdb-rs226.gather.csv
```

The results should be the same as if using the full database. Here we output to the same filename as in the main code so we can use the same visulization code.
:::


### Run the `sourmash gather` analysis 


```
sourmash gather SRR29654720.sig.zip gtdb-rs226.k31-sc10k.sig.rocksdb \
                -k 31 --scaled 10_000 --threshold-bp 3_000_000 \
                -o SRR29654720-x-gtdb-rs226.gather.csv
```
Command options: 
> - `--scaled 10_000` The database is scaled at 10k, while our query signatures are scaled=1000. We need to downsample the query to 10k to be able to use this search database. This is done on the fly during the search.
> - `-k 31` Search with k=31. Since the database is k=31, this is the only k-mer that will work with this database.
> - `--threshold-bp` This is the total amount of unique sequence overlap required to report a match. Here, I've set `--threshold-bp 3_000_000`, which is 3 million base pairs (3Mb). At `scaled=10_000`, this corresponds to 300 k-mers matched (`threshold-bp/scaled`). _This high threshold reduces the overall runtime and (more importantly for now) makes the results easier to display and visualize below. If time remains at the end, try rerunning this analysis with a lower threshold (e.g. 1Mb or 100kb) to see additional reference matches with <3Mb of unique sequence overlap._

You should see:
```
selecting specified query k=31
loaded query: SRR29654720 tomato rhizosphere... (k=31, DNA)
downsampling query from scaled=1000 to 10000
--
loaded 732034 total signatures from 1 locations.
after selecting signatures compatible with search, 732034 remain.

Starting prefetch sweep across databases.
Prefetch found 368229 signatures with overlap >= 3.0 Mbp.
Doing gather to generate minimum metagenome cover.

overlap     p_query p_match avg_abund
---------   ------- ------- ---------
6.4 Mbp        1.7%   93.9%      22.9    GCF_003951285.1 Variovorax beijingensis
6.3 Mbp        0.5%   97.8%       6.8    GCF_001976025.1 Rhodococcus sp. D-1
6.1 Mbp        0.4%   79.5%       5.3    GCF_001461685.1 Sinorhizobium sp. GL28
5.5 Mbp       36.5%   95.7%     568.1    GCA_013375735.1 Ralstonia solanacearum
5.2 Mbp        0.4%   89.7%       6.1    GCA_900472815.1 Hyphomicrobiales bacterium
5.0 Mbp        0.2%   68.7%       3.3    GCF_008824165.1 Pseudomonas umsongensis
5.1 Mbp        0.4%   77.3%       7.0    GCF_035610345.1 Sinorhizobium meliloti
4.7 Mbp        0.5%   88.5%       8.4    GCF_900013505.1 Agrobacterium deltaense RV3
5.0 Mbp        0.1%   61.1%       2.5    GCA_900472255.1 Hyphomicrobiales bacterium
4.7 Mbp        0.4%   75.5%       7.8    GCF_004343065.1 Neorhizobium sp. S3-V5DH
4.6 Mbp        0.1%   64.1%       2.3    GCA_035565435.1 Dyadobacter sp.
4.6 Mbp        0.4%   74.6%       7.4    GCF_022347545.1 Mesorhizobium retamae
4.6 Mbp        0.3%   77.7%       5.7    GCA_003400185.1 Bacillus sp. RC
5.1 Mbp        0.2%   64.7%       3.0    GCA_013421725.1 Ensifer adhaerens
4.5 Mbp        0.2%   75.6%       3.4    GCF_008807375.1 Pseudomonas lalkuanensis
4.4 Mbp        0.3%   76.1%       6.1    GCF_018447135.1 Enterobacter hormaechei subsp. xiangfangensis
4.4 Mbp        0.2%   61.5%       3.7    GCF_029841185.1 Achromobacter spanius
4.5 Mbp        0.2%   69.6%       3.0    GCA_900473555.1 Hyphomicrobiales bacterium
4.3 Mbp        0.2%   66.6%       3.4    GCF_031822455.1 Pseudomonas sp. HTZ1
4.1 Mbp        0.1%   67.0%       3.1    GCF_025989245.1 Rhodococcus pyridinivorans
4.1 Mbp        0.1%   56.7%       3.2    GCF_010994755.1 Shinella lacus
3.9 Mbp        0.2%   76.6%       3.8    GCF_022763525.1 Acidovorax sp.
3.9 Mbp        0.1%   37.3%       2.1    GCF_040446465.1 Streptomyces sp. RG80
3.9 Mbp        0.1%   60.5%       2.9    GCF_013283895.1 Pseudomonas sp. MPDS
3.8 Mbp        0.2%   78.8%       5.6    GCA_031995645.1 Leclercia adecarboxylata
4.0 Mbp        0.2%   76.4%       3.9    GCF_023751425.1 Enterobacter kobei
3.5 Mbp        0.1%   69.2%       2.7    GCF_037027445.1 Lysobacter firmicutimachus
3.6 Mbp        0.1%   55.3%       2.6    GCF_009601395.1 Sinorhizobium saheli
4.0 Mbp        0.3%   70.4%       8.9    GCF_000482285.1 Agrobacterium sp. UNC420CL41Cvi
3.4 Mbp        0.1%   64.9%       2.6    GCF_019976975.1 Arthrobacter sp. NtRootA1
3.7 Mbp        0.2%   59.3%       4.0    GCF_916855635.1 Acidovorax sp. Leaf73
3.2 Mbp        0.2%   81.5%       5.8    GCF_022662495.1 Lysobacter sp. M2-1
3.3 Mbp        0.1%   65.4%       3.4    GCF_037402075.1 Lysobacter sp. CCNWLW52
3.6 Mbp        0.1%   54.6%       2.8    GCF_014076515.1 Pseudomonas putida
3.1 Mbp        0.1%   43.0%       1.8    GCF_927798175.1 Paenibacillus sp. JJ-223
3.3 Mbp        0.1%   65.9%       3.6    GCA_016892445.1 Enterobacteriaceae bacterium RIT 814
3.7 Mbp        0.1%   60.7%       3.3    GCA_905331205.2 Enterobacter cloacae
found less than 3.0 Mbp in common. => exiting

found 37 matches total;
the recovered matches hit 45.6% of the abundance-weighted query.
the recovered matches hit 5.8% of the query k-mers (unweighted).
```
    

#### Key elements seen above:
Columns:
- `overlap` - estimated (unique) matching sequence
- `p_query` - percent of the query metagenome matched by this reference genome (containment of metagenome)
- `p_match` - percent of the reference genome matched by the query metagenome (containment of reference genome)
- `avg_abund` - the average abundance of k-mers matching this reference genome. e.g. 22.9 means we have an average of 23x coverage over this genome.
- **Note: the containment values are weighted by the k-mer abundance in the metagenome.**

Final reporting:
- `% of abundance-weighted query` -  fraction of the metagenome matched when accounting for k-mer abundance

**This information (and a lot more) is also in the output `csv` file, `SRR29654720-x-gtdb-rs226.gather.csv`. Descriptions of all columns are available in the sourmash documentation, [here](https://sourmash.readthedocs.io/en/latest/classifying-signatures.html#appendix-d-gather-csv-output-columns).**

## 2. Taxonomic profiling: summarize genome-level information by integrating taxonomy

`sourmash gather` uses a "minimum metagenome cover" approach to obtain the smallest list of reference genomes that best explain the query metagenome sequence. While this list of genomes is useful for follow-up via e.g. mapping approaches, we can also generate a 'taxonomic profile' of the metagenome, a summary of the taxonomic groupings present in the sequenced community.

To do this, we apply taxonomic lineage information from each genome, and optionally summarize to higher taxonomic rank. We can visualize the taxonomic profile using a sankey diagram.


### Obtain the taxonomic 'lineages' file

Sourmash provides a 'lineages' file for each database that links the genome accession to taxonomic information from domain/superkindom down to species (and strain, if provided by the reference database).

Download this file:
```
curl -JLO https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv
```

### Use `sourmash tax` to annotate the genomic results with their taxonomic lineages

```
sourmash tax annotate --gather-csv SRR29654720-x-gtdb-rs226.gather.csv \
                       --taxonomy gtdb-rs226.lineages.csv 
```

You should see:

```
Starting annotation on 'SRR29654720-x-gtdb-rs226.gather.csv'. Using ID column: 'name'
saving 'annotate' output to 'SRR29654720-x-gtdb-rs226.gather.with-lineages.csv'.
Annotated 37 of 37 total rows from 'SRR29654720-x-gtdb-rs226.gather.csv'.
```


### Visualize the taxonomic profile

To run this step, you'll need the [sourmash_plugin_betterplot](https://github.com/sourmash-bio/sourmash_plugin_betterplot) plugin installed. We installed this in module 1, but you can obtain this via: `pip install sourmash_plugin_betterplot`.

First, check that you can see the betterplot `sankey` command
```
sourmash info -v
```

You should see the following amongst the listed commands:
```the following plugins are installed:

plugin type          from python module             v     entry point name
-------------------- ------------------------------ ----- --------------------
sourmash.cli_script  sourmash_plugin_betterplot     0.5.1 sankey
```

Now, run build a `sankey` plot from the results:

```
sourmash scripts sankey --annotate-csv SRR29654720-x-gtdb-rs226.gather.with-lineages.csv -o SRR29654720-x-gtdb-rs226.sankey.png
```

![SRR29654720-x-gtdb-rs226.k31-sc10k-t3Mb.sankey](https://hackmd.io/_uploads/rkX1xEtGxl.jpg)


### Summarize at desired taxonomic rank:

There can be quite a lot of genomes (and quite a lot of species) represented in the metagenome. To help us make sense of the results, we can summarize at higher taxonomic ranks. Here, we'll look at summarizing at the `order` level, for simplicity. 

The `human` output format produces human readable output - there are a number of other formats that you may find useful. See documentation [here](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-metagenome-summarize-metagenome-content-from-gather-results) or run `sourmash tax metagenome --help`. You can produce multiple output formats with a single command, and can specify the output basename via `-o`.

**Build taxonomic summary using `tax metagenome`:**
```
sourmash tax metagenome --gather-csv SRR29654720-x-gtdb-rs226.gather.csv \
                        --taxonomy gtdb-rs226.lineages.csv \
                        --output-format human --rank order
```


You should see:

```
loaded 1 gather results from 'SRR29654720-x-gtdb-rs226.gather.csv'.
loaded results for 1 queries from 1 gather CSVs
sample name    proportion   cANI   lineage
-----------    ----------   ----   -------
SRR29654720 tomato rhizosphere   54.4%     -      unclassified
SRR29654720 tomato rhizosphere   38.7%     85.8%  d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales
SRR29654720 tomato rhizosphere    3.5%     88.2%  d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales
SRR29654720 tomato rhizosphere    1.0%     85.0%  d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales
SRR29654720 tomato rhizosphere    0.8%     85.4%  d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales
SRR29654720 tomato rhizosphere    0.6%     83.6%  d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales
SRR29654720 tomato rhizosphere    0.5%     83.4%  d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Xanthomonadales
SRR29654720 tomato rhizosphere    0.3%     81.4%  d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales
SRR29654720 tomato rhizosphere    0.1%     81.4%  d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Cytophagales
SRR29654720 tomato rhizosphere    0.1%     80.6%  d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Actinomycetales
SRR29654720 tomato rhizosphere    0.1%     80.9%  d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Streptomycetales
SRR29654720 tomato rhizosphere    0.1%     80.3%  d__Bacteria;p__Bacillota;c__Bacilli;o__Paenibacillales
```
This order-level describes the 11 orders that we can see by looking at the order (`o__`) level of the sankey diagram. This summarization can be done at any rank, or at all ranks by omitting the `--rank` parameter.


Output the full summary to a csv file by running this slightly modified version:

```
sourmash tax metagenome --gather-csv SRR29654720-x-gtdb-rs226.gather.csv \
--taxonomy gtdb-rs226.lineages.csv \
--output-format csv_summary -o SRR29654720-x-gtdb-rs226.gather
```

## 3. Searching eukaryotic genomes

This metagenome sample is a tomato rhizosphere (root environment), so we might expect to see some tomato sequence also present in our sequence. While we've limited this module to searching microbial genomes, FracMinHash sketching also works well for eukaryotic genomes. The challenge is that since eukaryotic genomes are bigger, the databases are also correspondingly larger! To get around that for this session, I've built a subset database just including genomes that matched to this sample. We can prepare that here:

download the database and the full lineages file:
```
curl -JLO https://github.com/WCSCourses/K-mer_Biodiversity_Genomics_2025/raw/refs/heads/main/course_data_2025/module8_sketching/gbentire-SRR29654720subset.k51-sc10k.rocksdb.tar.gz
curl -JLO https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/entire-2025-01-21/entire-2025-01-21.lineages.csv
```

untar
```
tax xzf gbentire-SRR29654720subset.k51-sc10k.rocksdb.tar.gz
```

Run `gather`
```
sourmash gather SRR29654720.sig.zip gbentire-SRR29654720subset.k51-sc10k.rocksdb \
                --scaled 10_000 -k 51 --threshold-bp 3_000_000 \
                -o SRR29654720-x-gbentire.k51.gather.csv
```

Command line output:

```
selecting specified query k=51
loaded query: SRR29654720... (k=51, DNA)
downsampling query from scaled=1000 to 10000
--
loaded 1296 total signatures from 1 locations.
after selecting signatures compatible with search, 1296 remain.

Starting prefetch sweep across databases.
Prefetch found 1296 signatures with overlap >= 3.0 Mbp.
Doing gather to generate minimum metagenome cover.

overlap     p_query p_match avg_abund
---------   ------- ------- ---------
26.5 Mbp       0.6%   19.2%       1.6    GCA_036172945.1 southern root-knot nematode (Meloidogyne incognita)
23.6 Mbp       0.5%    3.3%       1.4    GCF_036512215.1 tomato (Solanum lycopersicum)
6.6 Mbp        0.5%   94.9%       5.3    GCF_001976025.1 Rhodococcus sp. D-1 strain=D-1, ASM197602v1
6.2 Mbp        1.3%   90.3%      14.9    GCF_003951285.1 Variovorax sp. 502 strain=502, ASM395128v1
6.2 Mbp        0.1%   16.7%       1.3    GCF_023168085.1 Purpureocillium lilacinum
6.2 Mbp        0.3%    1.8%       3.8    GCA_907164665.1 Labrador sulphur (Colias nastes)
5.5 Mbp        0.3%   76.6%       3.7    GCA_900472695.1 Rhizobiales bacterium, AFS086625
5.5 Mbp       35.1%   94.5%     445.8    GCA_013306015.1 Ralstonia solanacearum strain=FJAT91.F50 ASM1330601v1
4.9 Mbp        0.3%   87.9%       4.8    GCA_900472815.1 Rhizobiales bacterium, AFS089153
4.5 Mbp        0.4%   89.1%       5.7    GCF_002008205.1 Agrobacterium sp. YIC 4121 strain=YIC 4121, ASM200820v1
4.5 Mbp        0.3%   72.8%       4.8    GCF_002565155.1 Bacillus megaterium strain=AFS075626, ASM256515v1
4.5 Mbp        0.3%   64.8%       5.0    GCF_002807095.1 Sinorhizobium meliloti strain=AC50a, ASM280709v1
4.4 Mbp        0.3%   73.2%       5.5    GCF_004368875.1 Neorhizobium sp. R1-B strain=R1-B, ASM436887v1
4.5 Mbp        0.1%   59.1%       2.3    GCA_013421725.1 Ensifer adhaerens strain=HP1, ASM1342172v1
4.2 Mbp        0.2%   70.1%       3.1    GCF_020422905.1 Pseudomonas lalkuanensis strain=ACYW.190, ASM2042290v1
4.2 Mbp        0.3%   64.9%       5.3    GCF_022347545.1 Mesorhizobium sp. IRAMC:0171 strain=IRAMC:0171, ASM2234754v1
4.0 Mbp        0.1%   53.3%       2.6    GCF_008824165.1 Pseudomonas umsongensis strain=GO16, ASM882416v1
4.1 Mbp        0.1%   49.2%       2.1    GCF_024380055.1 Ensifer adhaerens strain=NER9, ASM2438005v1
3.9 Mbp        0.2%   74.1%       4.4    GCF_002588255.1 Leclercia adecarboxylata strain=FDAARGOS_404, ASM258825v1
3.6 Mbp        0.1%   59.3%       2.7    GCF_024494485.1 Pseudomonas asiatica strain=MD9, ASM2449448v1
3.6 Mbp        0.1%   50.4%       2.9    GCF_029841185.1 Achromobacter spanius strain=GD03843 ASM2984118v1
3.5 Mbp        0.1%   48.5%       2.2    GCF_010994755.1 Shinella sp. CPCC 100929 strain=CPCC 100929, ASM1099475v2
3.7 Mbp        0.1%   73.8%       2.8    GCF_001022355.1 Enterobacter kobei strain=GN02366, ASM102235v1
3.4 Mbp        0.1%   53.2%       2.0    GCF_900156295.1 Pseudomonas sp. A214 strain=A214, IMG-taxon 2681812807 annotated assembly
3.4 Mbp        0.1%   56.9%       2.5    GCA_900473555.1 Rhizobiales bacterium, R129_B
3.5 Mbp        0.3%   67.6%       6.7    GCF_000482285.1 Agrobacterium sp. UNC420CL41Cvi strain=UNC420CL41Cvi, ASM48228v1
3.4 Mbp        0.1%   49.8%       2.2    GCF_024520415.1 Pseudomonas putida strain=T815, ASM2452041v1
3.1 Mbp        0.1%   55.0%       2.2    GCF_000336385.2 Lysobacter antibioticus HS124 strain=HS124, HS124_Split10plusN
3.1 Mbp        0.1%   53.9%       2.2    GCF_021026015.1 Rhodococcus rhodochrous strain=IEGM 1360, ASM2102601v1
3.0 Mbp        0.1%   59.6%       2.4    GCF_019976975.1 Arthrobacter sp. NtRootA1 strain=NtRootA1, ASM1997697v1
found less than 3.0 Mbp in common. => exiting

found 30 matches total;
the recovered matches hit 43.1% of the abundance-weighted query.
the recovered matches hit 6.7% of the query k-mers (unweighted).
```

We use the full lineages file to annotate the gather results:
```
sourmash tax annotate -g SRR29654720-x-entire-2025-01-21.k51.gather.csv \
                      -t entire-2025-01-21.lineages.csv
```

### Visualize the output:

...and then use the `sankey` command to visualize:
```
sourmash scripts sankey --annotate-csv SRR29654720-x-entire-2025-01-21.k51.gather.with-lineages.csv -o SRR29654720-x-entire-2025-01-21.k51.sankey.png
```

![SRR29654720-x-gtdb-rs226.k31-sc10k-t3Mb.sankey](https://hackmd.io/_uploads/H1nAGEFfgl.jpg)

### Why are we using k=51 instead of k=31?

Eukaryotic genomes on GenBank seem to share a lot more k-mers than the GTDB microbial database genomes do. This may partly be because eukaryotes as a group are more closely related than bacteria and archaea, and longer ksizes help to distinguish closely related genomes. However, there are also lots of contamination issues! For example, we've seen large sections of human sequence - Mb of exact sequence matches - within draft genomes of fish. Keep in mind that reference-based approaches can only be as good as their reference databases!

## Appendix

### A. Using other sourmash databases

You can download a prepared databases of sketches from the [prepared databases](https://sourmash.readthedocs.io/en/latest/databases.html) page. Running `sourmash gather` on sourmash `zip` files will be slightly slower than using `rocksdb` databases. Please look into `fastgather` and `fastmultigather` in the [branchwater plugin](https://github.com/sourmash-bio/sourmash_plugin_branchwater) for faster implementations.


### Downloading the GenBank 'entire' database for use (27G)
The `entire` database is not yet available on the list of sourmash prepared databases, partly because it is so large (27G!). However, you can download and use it like so:

Download rocksdb:
```
wget -c --recursive --no-parent -nH --cut-dirs=3 --reject "index.html*"  https://farm.cse.ucdavis.edu/\~ctbrown/sourmash-db/entire-2025-01-21/entire-2025-01-21.k51.rocksdb
```
Download lineages(same as above):
```
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/entire-2025-01-21/entire-2025-01-21.lineages.csv
```

`gather`:

```
sourmash gather SRR29654720.sig.zip entire-2025-01-21.k51.rocksdb \
                --scaled 10_000 -k 51 --threshold-bp 3_000_000 \
                -o SRR29654720-x-gbentire.k51.gather.csv
```
You can reduce or even remove the `--threshold-bp 3_000_000` to get more matches, but it will take a bit longer to run. When I remove it entirely, using the default threshold of 50 kb, the final output is updated to:

```
found 1296 matches total;
the recovered matches hit 54.4% of the abundance-weighted query.
the recovered matches hit 20.6% of the query k-mers (unweighted).
```
This shows that we were able to match more of the original metagenome when we included those smaller-overlap matches. If you do taxonomic summarization, you can see that many of the additional matched genomes represent strain variants of prior matches.

 
### B. How does sourmash gather analysis compare with alignment approaches?

For the _R. nicotianae_ (_R. pseudosolanacearum_)** genome and our metagenome, we will map metagenome reads via `minimap2` and compare the fraction found via FracMinHash vs mapping.

**_There's debate over this naming, as the proposed "nicotianae" species is a member of the _R. pseudosolanacearum_ species based on ANI value. Keep in mind taxonomic debates are common and names are subject to change over time!_


#### Install `minimap2` and `samtools`

```
conda install -c conda-forge -c bioconda minimap2 samtools
```

#### Index the genome for minimap2

```
minimap2 -d GCF_013306015.mmi GCF_013306015.1_ASM1330601v1_genomic.fna.gz
```

You should see:
```
[M::mm_idx_gen::0.120*1.00] collected minimizers
[M::mm_idx_gen::0.155*1.40] sorted minimizers
[M::main::0.223*1.21] loaded/built the index for 2 target sequence(s)
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 2
[M::mm_idx_stat::0.231*1.20] distinct minimizers: 972790 (91.92% are singletons); average occurrences: 1.126; average spacing: 5.362; total length: 5874107
[M::main] Version: 2.29-r1283
[M::main] CMD: minimap2 -d GCF_013306015.mmi GCF_013306015.1_ASM1330601v1_genomic.fna.gz
[M::main] Real time: 0.235 sec; CPU: 0.281 sec; Peak RSS: 0.097 GB
```


#### Map reads to the genome
This step takes a little while, since we are mapping all the reads onto the genome.
```
minimap2 -ax sr GCF_013306015.mmi SRR29654720_1.fastq.gz SRR29654720_2.fastq.gz | samtools view -bS - | samtools sort -o SRR29654720-x-GCF_013306015.sorted.bam
```

Note that `minimap2` uses minimizers internally to identify initial exact matches (to anchor reads and speed up final alignment).


#### Index the BAM

```
samtools index SRR29654720-x-GCF_013306015.sorted.bam
```

#### Count average depth over this genome

```
samtools depth  SRR29654720-x-GCF_013306015.sorted.bam >  SRR29654720-x-GCF_013306015.depth.txt
```

Take a quick look at the file: each line is the depth at a particular position
```
head SRR29654720-x-GCF_013306015.depth.txt
```

Sum column 3 and divide by the total count to get the average depth
```
awk '{sum+=$3; count++} END {print "Average depth:", sum/count}' SRR29654720-x-GCF_013306015.depth.txt
```

    Average depth: 731.563
    
> This is a bit higher than we saw with `gather`. The reason is likely two-fold: First, scaled=10,000 is a pretty lossy scaled value for microbial genomes. We recommend using scaled=1000 for most microbial searches. Second, exact k-mer matching will miss regions of the genome with SNPs, while mapping (with mismatches) will be able to cover these.


#### How much of the genome was covered?


```
awk '$3 > 0 {covered++} {total++} END {print "Covered:", covered, "Total:", total, "Fraction covered:", covered/total}' SRR29654720-x-GCF_013306015.depth.txt
```

    Covered: 5744411 Total: 5744564 Fraction covered: 0.999973


Here we see 99% of the genome was covered.
