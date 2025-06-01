# 3. What is in my metagenome? Metagenomic profiling via FracMinHash


**This follows from the [Introduction to FracMinHash sketching for sequence comparison](https://hackmd.io/@bluegenes/H1ItsC_fxx) and [Comparing genomes metagenomes using FracMinHash sketches](https://hackmd.io/@bluegenes/Byh_JfKGgg) modules. Please refer to the installation instructions in the introduction for installing sourmash and other dependencies.**

With metagenomic sequencing, a first step is often to determine what organisms are present in the sequenced community.

Here, we demonstrate comprehensive metagenome compositional analysis using [sourmash](sourmash.readthedocs.io) FracMinHash sketching, and aggregation to taxonomic groupings. This analysis is entirely reference-dependent. The quality of your results will depend in part on the quality and completeness of the reference database.

> Checks:
>
> Do you have sourmash installed?
> `sourmash --version`

## Genome-based compositional analysis of a metagenome

We're going to use the same metagenome used in part 2. Instead of searching with a single genome, we're going to compare the metagenome against a comprehensive list of reference genomes using the "minimum metagenome cover approach". This approach uses two steps internally: first, search all reference genomes for k-mer overlaps, and second, use a greedy method to build the smallest list of reference genomes that contain all shared (known) k-mer content.

We will start by searching only bacterial and archaeal reference genomes, as the database is much smaller. If working on the course instances, we will also search a database the includes eukaryotic genomes, since we can provide that database in a shared location without requiring donwload.

In sourmash, we can use the `gather` command to run this analysis.

Datasets we will use:
  - **Metagenome sample**: [SRR29654720](https://www.ncbi.nlm.nih.gov/sra/?term=SRR29654720) from the [BioProject PRJNA1127303](https://www.ncbi.nlm.nih.gov/bioproject/1127303).
  - **Associated paper**: [Risk assessment of antibiotic resistance genes in rhizosphere soil during tomato growth under bio-control bacterial inoculation](https://www.sciencedirect.com/science/article/pii/S0959652625002616)
    > The study aims to investigate the effects of different fertilization practices during various growth stages of tomatoes on the abundance and dynamics of antibiotic resistance genes (ARGs) in the soil. Understanding the mechanisms behind these impacts is critical, as ARGs pose a significant threat to public health by promoting the development and spread of antibiotic-resistant bacteria. By examining the influence of organic, inorganic, and combined fertilization methods on ARG profiles, as well as considering factors such as microbial community composition and soil physicochemical properties, this research seeks to elucidate how fertilization strategies might mitigate or exacerbate the propagation of ARGs in agricultural environments. The findings from this study will provide valuable insights for developing sustainable farming practices that minimize the risk of ARG proliferation, thereby contributing to the broader effort of combating antibiotic resistance.

If you don't have the file locally available (e.g. you are following along at a later date), you can download this file here:
```
# (generated in part 2)
!curl -L -o SRR29654720.sig.zip https://osf.io/download/hrdsw
!ls SRR29654720.sig.zip
```

Examine this file with `sourmash`:
```
sourmash sig describe SRR29654720.sig.zip
```

You should see:
```
signature filename: SRR29654720.sig.zip
signature: SRR29654720 tomato rhizosphere
source file: SRR29654720_2.fastq.gz
md5: 57e826e874fe03b8e4438b19d64b8652
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=1
size: 2714833
sum hashes: 8461460
signature license: CC0
```
> This file contains a sourmash sketch of SRA run SRR29654720

### Obtain a database and examine it

We are downloading sketches of the Genome Taxonomy Database's (GTDB) "species representatives" for refseq release 226 (Sept 2024). This release contains 732k bacterial and archael genomes, from 143k species.

![image](https://hackmd.io/_uploads/HyLZcXYGle.png)


**This sourmash database is constructed at k=31, scaled = 10,000**, meaning we keep roughly 1/10000 k-mers per genome. We're using this low resolution to keep the database small and the analysis fast during the workshop. Higher resolution databases are available for follow-up analyses.

The database should be pre-loaded onto your instance. You should see it available via:

```
ls gtdb-rs226.k31-sc10k.sig.rocksdb
```

If you need to re-download at a later date:

```
# Download the database
#curl -JLO https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226-reps.k31-sc10k.sig.rocksdb.tar.gz 
#curl -JLO https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226.k31-sc10k.sig.rocksdb.tar.gz 
and untar like so:
#tar xzf gtdb-rs226-reps.k31-sc10k.sig.rocksdb.tar.gz
#tar xzf gtdb-rs226.k31-sc10k.sig.rocksdb.tar.gz
```
> The `gtdb-rs226-reps.k31-sc10k.sig.rocksdb` database is a smaller version of the full GTDB database, containing only the "representative" genomes for each species. The `gtdb-rs226.k31-sc10k.sig.rocksdb` database contains all genomes in the GTDB release. You can use either for this analysis, but downloading the representative database will be faster and use less disk space. The commands below specifically use the `gtdb-rs226.k31-sc10k.sig.rocksdb` database (update if needed).


### Examine the database with sourmash
This database is a 'rocksdb' index, a sourmash format that maps k-mers (hashes) as keys to their datasets (values) and is optimized for fast on-disk search using a single thread.


<!---
reps version
```
#sourmash sig summarize gtdb-rs226-reps.k31-sc10k.sig.rocksdb
```
--->

```
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

### Run the analysis 
Command options: 
> - `--scaled 10_000` The database is scaled at 10k, while our query signatures are scaled=1000. We need to downsample the query to 10k to be able to use this search database. This is done on the fly during the search.
> - `-k 31` Search with k=31. Since the database is k=31, this is the only k-mer that will work with this database.
> - `--threshold-bp` This is the total amount of unique sequence overlap required to report a match. Here, I've set `--threshold-bp 3_000_000`, which is 3 million base pairs (3Mb). At `scaled=10_000`, this corresponds to 300 k-mers matched (`threshold-bp/scaled`). _This high threshold reduces the overall runtime and (more importantly for now) makes the results easier to display below. If time remains at the end, rerun this analysis with a lower threshold (e.g. 1Mb or 100kb) to see additional reference matches with <3Mb of unique sequence overlap._

<!---
reps version
```
#sourmash gather SRR29654720.sig.zip gtdb-rs226-reps.k31-sc10k.sig.rocksdb \
                -k 31 --scaled 10_000 --threshold-bp 3_000_000 \
                -o SRR29654720-x-gtdb-rs226-reps.gather.csv
```
--->

```
sourmash gather SRR29654720.sig.zip gtdb-rs226.k31-sc10k.sig.rocksdb \
                -k 31 --scaled 10_000 --threshold-bp 3_000_000 \
                -o SRR29654720-x-gtdb-rs226.gather.csv
```

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

**This information (and a lot more) is also in the output `csv` file, `SRR29654720-x-gtdb-rs226.gather.csv`.**

## 2. Taxonomic profiling: summarize genome-level information by integrating taxonomy

`sourmash gather` uses a "minimum metagenome cover" approach to obtain the smallest list of reference genomes that best explain the query metagenome sequence. While this list of genomes is useful for follow-up via e.g. mapping approaches, we can also generate a 'taxonomic profile' of the metagenome, a summary of the taxonomic groupings present in the sequenced community.

To do this, we apply taxonomic lineage information from each genome, and optionally summarize to higher taxonomic rank. We can visualize the taxonomic profile using a sankey diagram.

### Obtain the taxonomic 'lineages' file

Sourmash provides a 'lineages' file for each database that links the genome accession to taxonomic information from domain/superkindom down to species (and strain, if provided by the reference database).

You should have this file locally. You can download later using:
```
# download and check this file
#curl -JLO https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv
wc -l gtdb-rs226.lineages.csv
```

### Use `sourmash tax` to annotate the genomic results with their taxonomic lineages

<!---
```
sourmash tax annotate --gather-csv SRR29654720-x-gtdb-rs226-reps.gather.csv \
                       --taxonomy gtdb-rs226.lineages.csv 
```
--->
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
sourmash.cli_script  sourmash_plugin_betterplot     0.5.0 sankey
```

Now, run build a `sankey` plot from the results:

<!---
```
sourmash scripts sankey --annotate-csv SRR29654720-x-gtdb-rs226-reps.gather.with-lineages.csv
```
--->
```
sourmash scripts sankey --annotate-csv SRR29654720-x-gtdb-rs226.gather.with-lineages.csv -o SRR29654720-x-gtdb-rs226.k31-sc10k-t3Mb.sankey.png
```

![SRR29654720-x-gtdb-rs226.k31-sc10k-t3Mb.sankey](https://hackmd.io/_uploads/rkX1xEtGxl.jpg)


### Summarize at desired taxonomic rank:

Here, we'll look at summarizing at the `order` level, for simplicity. Default output produces the summary at each rank.

The `human` output format produces human readable output - there are a number of other formats that you may find useful. See documentation [here](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-tax-metagenome-summarize-metagenome-content-from-gather-results) or run `sourmash tax metagenome --help`. You can produce multiple output formats with a single command, and can specify the output basename via `-o`.

<!---
```
### Build taxonomic summary using `tax metagenome`:
!sourmash tax metagenome --gather-csv SRR29654720-x-gtdb-rs226-reps.gather.csv \
                         --taxonomy gtdb-rs226.lineages.csv \
                         --output-format human --rank order
```
--->
```
### Build taxonomic summary using `tax metagenome`:
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

This sample is a tomato rhizosphere (root environment), so we might expect to see some tomato sequence also present in our sequence. Sourmash also provides search across all GenBank genomes, including eukaryotes. While search can be done very quickly with `rocksdb` indexes, the database is too large to download during this session.

The database should be installed on your instance. If you're interested in running later, you can download the database and lineages file at a later date -- the code will be provided at the very bottom.

Run the following command:
```
sourmash gather SRR29654720.sig.zip entire-2025-01-21.k31.rocksdb \
                -o SRR29654720-x-entire-2025-01-21.k31-sc10k-t3Mb.gather.csv --scaled 10_000 \
                -k 31 --threshold-bp 3_000_000
```

Command line output:

```
selecting specified query k=31
loaded query: SRR29654720 tomato rhizosphere... (k=31, DNA)
downsampling query from scaled=1000 to 10000
loading from '/group/ctbrowngrp5/sourmash-db/entire-2025-01-21/entire-2025-01-21.k31.rocksdb'...
loaded 729299 total signatures from 1 locations.
after selecting signatures compatible with search, 729299 remain.
Starting prefetch sweep across databases.
Prefetch found 351580 signatures with overlap >= 3.0 Mbp.
Doing gather to generate minimum metagenome cover.


overlap     p_query p_match avg_abund
---------   ------- ------- ---------
30.6 Mbp       0.6%   24.0%       1.7    GCA_036172945.1 southern root-knot n...
27.4 Mbp       0.5%    4.3%       1.5    GCF_036512215.1 tomato (Solanum lyco...
8.1 Mbp        0.1%   21.1%       1.5    GCF_023168085.1 Purpureocillium lila...
6.5 Mbp        0.4%    2.0%       5.4    GCA_907164665.1 Labrador sulphur (Co...
6.4 Mbp        1.7%   93.9%      22.9    GCF_003951285.1 Variovorax sp. 502 s...
6.3 Mbp        0.5%   97.8%       6.8    GCF_001976025.1 Rhodococcus sp. D-1 ...
6.1 Mbp        0.4%   79.5%       5.3    GCF_001461685.1 Sinorhizobium sp. GL...
5.5 Mbp       36.5%   95.8%     568.1    GCA_013306015.1 Ralstonia solanacear...
5.2 Mbp        0.4%   89.7%       6.1    GCA_900472815.1 Rhizobiales bacteriu...
5.0 Mbp        0.2%   68.7%       3.3    GCF_008824165.1 Pseudomonas umsongen...
4.9 Mbp        0.2%    0.3%       3.8    GCA_964133935.1 divine cactus (Lopho...
4.7 Mbp        0.5%   88.3%       8.4    GCF_900013505.1 Agrobacterium deltae...
5.0 Mbp        0.1%   61.1%       2.6    GCA_900472255.1 Rhizobiales bacteriu...
4.7 Mbp        0.4%   75.5%       7.8    GCF_004343065.1 Neorhizobium sp. S3-...
4.6 Mbp        0.4%   74.7%       7.5    GCF_022347545.1 Mesorhizobium sp. IR...
4.6 Mbp        0.3%   77.5%       5.1    GCA_003400185.1 Bacillus sp. RC stra...
5.1 Mbp        0.2%   64.7%       3.0    GCA_013421725.1 Ensifer adhaerens st...
4.6 Mbp        0.4%   72.7%       6.9    GCA_029201645.1 Sinorhizobium sp. K1...
4.5 Mbp        0.2%   74.9%       3.4    GCF_008807375.1 Pseudomonas sp. PE08...
4.4 Mbp        0.2%   61.3%       3.6    GCF_029841185.1 Achromobacter spaniu...
4.5 Mbp        0.2%   69.4%       3.0    GCA_900473555.1 Rhizobiales bacteriu...
4.3 Mbp        0.2%   63.5%       3.3    GCF_001306495.1 Pseudomonas putida s...
4.1 Mbp        0.1%   56.7%       3.2    GCF_010994755.1 Shinella sp. CPCC 10...
4.1 Mbp        0.1%   67.0%       3.1    GCF_025989245.1 Rhodococcus pyridini...
3.9 Mbp        0.2%   76.6%       3.8    GCF_022763525.1 Acidovorax sp. ASM22...
3.9 Mbp        0.1%   60.3%       2.9    GCF_013283895.1 Pseudomonas fluoresc...
3.8 Mbp        0.2%   78.3%       5.6    GCF_029892245.1 Leclercia adecarboxy...
3.5 Mbp        0.1%   68.2%       2.8    GCF_000336385.2 Lysobacter antibioti...
3.6 Mbp        0.1%   55.3%       2.6    GCF_009601395.1 Sinorhizobium saheli...
4.0 Mbp        0.3%   70.4%       8.9    GCF_000482285.1 Agrobacterium sp. UN...
4.0 Mbp        0.1%   71.2%       3.6    GCF_023751425.1 Enterobacter kobei s...
3.4 Mbp        0.1%   64.9%       2.6    GCF_019976975.1 Arthrobacter sp. NtR...
3.7 Mbp        0.2%   59.3%       4.0    GCF_916855635.1 Acidovorax sp. Leaf7...
3.2 Mbp        0.2%   81.0%       5.8    GCF_022662495.1 Lysobacter sp. M2-1 ...
3.1 Mbp        0.0%    7.0%       1.4    GCF_001653235.2 Pochonia chlamydospo...
3.1 Mbp        0.1%   43.0%       1.8    GCF_927798175.1 Paenibacillus sp. JJ...
3.6 Mbp        0.1%   54.1%       2.7    GCF_014076515.1 Pseudomonas putida s...
3.3 Mbp        0.1%   65.7%       3.5    GCA_016892445.1 Enterobacteriaceae b...
found less than 3.0 Mbp in common. => exiting
```

Use the lineages file to annotate the gather results:
```
sourmash tax annotate -g SRR29654720-x-entire-2025-01-21.k31-sc10k-t3Mb.gather.csv \
                      -t entire-2025-01-21.lineages.csv
```

### Visualize the output:

We can again use the `sankey` command to visualize:
```
sourmash scripts sankey --annotate-csv SRR29654720-x-entire-2025-01-21.k31-sc10k-t3Mb.gather.with-lineages.csv -o SRR29654720-x-entire-2025-01-21.k31-sc10k-t3Mb.sankey.png
```

![SRR29654720-x-gtdb-rs226.k31-sc10k-t3Mb.sankey](https://hackmd.io/_uploads/H1nAGEFfgl.jpg)


## Appendix

### Using different databases

You can download a prepared databases of sketches from the [prepared databases](https://sourmash.readthedocs.io/en/latest/databases.html) page. Running `sourmash gather` on sourmash `zip` files will be slightly slower than using `rocksdb` databases. Please look into `fastgather` and `fastmultigather` in the [branchwater plugin](https://github.com/sourmash-bio/sourmash_plugin_branchwater) for faster implementations.


### Downloading the 'entire' database for use (17G)
The `entire` database is not yet available on the list of sourmash prepared databases. However, you can download it like so:

Download rocksdb:
```
wget -c --recursive --no-parent -nH --cut-dirs=3 --reject "index.html*"  https://farm.cse.ucdavis.edu/\~ctbrown/sourmash-db/entire-2025-01-21/entire-2025-01-21.k31.rocksdb
```
Download lineages:
```
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/entire-2025-01-21/entire-2025-01-21.lineages.csv
```

Use the `gather` commands above to run the analysis.

 
 ### How does sourmash gather analysis compare with alignment approaches?

For the _R. nicotianae_/_R. solanacearum_ genome and our metagenome, we will map metagenome reads via `minimap2` and compare the fraction found via FracMinHash vs mapping.

Note that `minimap2` uses minimizers internally to identify initial exact matches (to anchor reads and speed up final alignment).

#### First, let's index the genome for minimap2


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
Note that `minimap2` is using minimizers for faster initial seeding of reads to the genome!

#### Map reads to the genome
This step takes a little while, since we are mapping all the reads onto the genome.
```
minimap2 -ax sr GCF_013306015.mmi SRR29654720_1.fastq.gz SRR29654720_2.fastq.gz | samtools view -bS - | samtools sort -o SRR29654720-x-GCF_013306015.sorted.bam
```

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


#### How much of the genome was covered?


```
awk '$3 > 0 {covered++} {total++} END {print "Covered:", covered, "Total:", total, "Fraction covered:", covered/total}' SRR29654720-x-GCF_013306015.depth.txt
```

    Covered: 5744411 Total: 5744564 Fraction covered: 0.999973


Here we see 99% of the genome was covered.
