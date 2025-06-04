# 2. Comparing genomes and metagenomes using FracMinHash sketches

[HackMD version](https://hackmd.io/@bluegenes/Byh_JfKGgg)

**This follows from the [Introduction to FracMinHash sketching for sequence comparison](https://hackmd.io/@bluegenes/H1ItsC_fxx) module. Please see the installation instructions there for installing sourmash and other dependencies.**

A common application for sketching is detecting genomes of interest within metagenomes to determine presence/absence of organisms of interest, such as detecting pathogenic microbes in wastewater metagenomic sequencing.

[sourmash](sourmash.readthedocs.io) FracMinHash sketching is well-suited to conducting these comparisons at scale. We first look at genome -x- metagenome comparisons and then move to finding genomes(s) in metagenomes at scale.

## Genome x Metagenome comparisons
In this section, we will compare a query genome to a metagenome sample. The goal is to determine if the query genome is present in the metagenome and, if so, how much of it is present.

Datasets we will use:
- **Query microbial genome**: [GCA_013306015 "FJAT91.F50"](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013306015.1) _Ralstonia nicotianae_ / _Ralstonia psuedosolanacearum_.A
- **Metagenome**: Tomato rhizosphere metagenome, collected from a field in China, treated with organic fertilizer and inoculated with _Ralstonia pseudosolanacearum_ to examine disease progression. It was collected during the fruiting stage of tomato growth.
  - **Metagenome sample**: [SRR29654720](https://www.ncbi.nlm.nih.gov/sra/?term=SRR29654720) from the [BioProject PRJNA1127303](https://www.ncbi.nlm.nih.gov/bioproject/1127303).
  - **Associated paper**: [Risk assessment of antibiotic resistance genes in rhizosphere soil during tomato growth under bio-control bacterial inoculation](https://www.sciencedirect.com/science/article/pii/S0959652625002616)
    > The study aims to investigate the effects of different fertilization practices during various growth stages of tomatoes on the abundance and dynamics of antibiotic resistance genes (ARGs) in the soil. Understanding the mechanisms behind these impacts is critical, as ARGs pose a significant threat to public health by promoting the development and spread of antibiotic-resistant bacteria. By examining the influence of organic, inorganic, and combined fertilization methods on ARG profiles, as well as considering factors such as microbial community composition and soil physicochemical properties, this research seeks to elucidate how fertilization strategies might mitigate or exacerbate the propagation of ARGs in agricultural environments. The findings from this study will provide valuable insights for developing sustainable farming practices that minimize the risk of ARG proliferation, thereby contributing to the broader effort of combating antibiotic resistance.


### Download and sketch the query genome
Here we follow the same process as in the previous section, with a new genome

Download
```
curl -JLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/306/015/GCF_013306015.1_ASM1330601v1/GCF_013306015.1_ASM1330601v1_genomic.fna.gz
```

Sketch
```
sourmash sketch dna -p k=31,scaled=1000,abund GCF_013306015.1_ASM1330601v1_genomic.fna.gz \
                     --name "GCF_013306015.1 Ralstonia nicotianae FJAT91" -o GCF_013306015.sig

```
> We will add the 'abund' parameter to the sketch parameters (after `-p`). This will allow us to compare k-mer abundances between the two sketches.

### Download a sketch of the metagenome sample and examine it

To save time and space on our devices, we'll download a pre-built sketch of the metagenome sample. For commands to rebuild from reads, see the Appendix section.

```
curl -JLO https://github.com/WCSCourses/K-mer_Biodiversity_Genomics_2025/raw/refs/heads/main/course_data_2025/module8_sketching/SRR29654720.sig.zip

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
> There are two signatures included in this file, one for each k-mer size.


You will notice that this file has a `sig.zip` extension, which means this is a compressed archive of the signature file. While this format is particularly useful for large datasets and for storing multiple signatures in a single file (e.g. databases), we recommend using it for any number of signatures.



## Compare the query genome to the metagenome

We will use the `sourmash sig overlap` command to compare the query genome sketch to the metagenome sketch. This will give us an idea of how many k-mers from the query genome are present in the metagenome.
```
sourmash sig overlap GCF_013306015.sig SRR29654720.sig.zip -k 31
```

Here is the relevant comparison output:

```
--- Similarity measures ---
jaccard similarity:          0.00200
first contained in second:   0.96858 (cANI: 0.99897)
second contained in first:   0.00200 (cANI: 0.81832)
average containment ANI:     0.90865


--- Hash overlap summary ---
number of hashes in first:   5601
number of hashes in second:  2714833

number of hashes in common:  5425
only in first:               176
only in second:              2709408
total (union):               2715009


--- Abundance-weighted similarity: ---
angular similarity:          0.60318
first contained in second (weighted): 0.96862
second contained in first (weighted): 0.36035

number of hashes in first (weighted): 5736
number of hashes in second (weighted): 8461460
```
Here we see that the Jaccard similarity is quite small, which makes sense given the size difference between the query genome and the metagenome. However, the containment values are more informative: the query genome is contained in the metagenome at a level of 97%, meaning that 97% of the unique k-mers in the genome are found in the metagenome. This is a high level of containment, indicating that the query genome is likely present in the metagenome.

Even more interesting are the Abundance-weighted similarity values. The "second contained in first (weighted)" value of 36% indicates that nearly 36% of the abundance-weighted k-mers in the metagenome are k-mers from this genome. This roughly corresponds to the proportion of the metagenome that is _Ralstonia_: aka, ~36% of all content in the metagenome is _Ralstonia_ sequence. This suggests that the query genome was present in the community at really high abundance, which is consistent with the fact that this sample was innoculated with _Ralstonia pseudosolanacearum_.


## Finding genomes in metagenomes at scale

In this trivial example, we knew we would find this genome in this metagenome. However, we can use the same strategy to search for a genome of interest in arbitrarily large collections of metagenomes. We will explore this via web-based search of over 1.1 million metagenomes in the Sequence Read Archive at [branchwater.jgi.doe.gov](https://branchwater.jgi.doe.gov/)

1. **Go to the branchwater website, and click on the 'Examples' tab.**

![image](https://hackmd.io/_uploads/HkfeGXKMxg.png)

2. **Select 'Ralstonia solanacearum'**

![image](https://hackmd.io/_uploads/rkgtGXtGeg.png)

3. **Click "Submit".**

4. **Wait for Results (<30 seconds)**

![image](https://hackmd.io/_uploads/HkpHX7tGex.png)

> If you're waiting a minute or more for this example search, you may need to refresh your page and start over.

5. **Try filtering for a minimum of 98% cANI (0.98)**

The default parameters may return results with genus-level similarity (k=21). To narrow the the results, you can select on the reported cANI value to ensure high similarity with your query. Type `0.98` in the cANI 'Min' box and hit 'Enter'/'Return' on your keyboard. You should see that now the list of results is reduced to just two SRA metagenomes, one from the United States, and a second within a "synthetic metagenome" from Switzerland.

You can also try subsetting in other categories. For example, type 'soil' in the 'organism' box and hit "Enter"/"Return". This will limit to just metagenomes with "soil" in the SRA metadata column. Alternatively, you can download the full CSV and do any subselection on the results file afterwards.

### Notes and Limitations

Branchwater web search has a few limitations at the moment. First, it uses a 'scaled' value of 1000, meaning that it will not be useful for finding small queries < 10kb, and may even struggle with queries <50kb. If you want to search for small queries (<1kb) or want to search other (non-metagenome) datasets, check out [Logan Search](https://logan-search.org/), which searches in pre-built assemblies of >27 million SRA runs.

Second, the query sketch is limited to 5Mb, which means you can't query with e.g. metagenome samples (most genomes should be fine; extremely large eukaryotic genomes may present issues). 

The branchwater SRA search is web-based because the underlying metagenome database is over 1Tb in size, meaning it is impractical to provide the database for download. However, you can build sketch databases with your own local metagenomes and search them quickly using the sourmash [branchwater](https://github.com/sourmash-bio/sourmash_plugin_branchwater) plugin, e.g. via `manysearch`. For very large databases, we recommend indexing the files using sourmash's `rocksdb` index type. We will use this in the next hands-on section.


## Appendix

If you want to download the full metagenome reads and sketch at a later date, uncomment the code below:

Download:
```
#curl -JLO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/020/SRR29654720/SRR29654720_1.fastq.gz
#curl -JLO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/020/SRR29654720/SRR29654720_2.fastq.gz
```

sketch:
```
#sourmash sketch dna -p k=31,k=51,scaled=1000,abund SRR29654720_1.fastq.gz SRR29654720_2.fastq.gz \
                    --name SRR29654720 -o SRR29654720.sig.zip
```
> -p k=31,k=51,scaled=1000,abund are the sketch parameters. This means we are sketching at both k=31 and k=51, using a scaled value of 1000 (keep ~ 1/1000 k-mers, and keeping track of k-mer abundance)
> *Note, if you do this in the same directory as an existing `SRR29654720.sig.zip` file, the sketches will be automatically added to the same zipfile.*


If you want to sketch a bit faster, we recommend using the [branchwater](https://github.com/sourmash-bio/sourmash_plugin_branchwater) plugin. Here is the equivalent command using `sourmash scripts singlesketch`. Again, you need to uncomment to run.


```
#conda install -c conda-forge -c bioconda sourmash_plugin_branchwater
#sourmash scripts singlesketch -p k=31,k=51,scaled=1000,abund SRR29654720_1.fastq.gz SRR29654720_2.fastq.gz \
                              --name SRR29654720 -o SRR29654720.sig.zip
```

This plugin also contains the `sourmash scripts manysketch` command, which is the fastest way to sketch many files at once. 

Run the search with `manysearch`. In this example I'm using `--scaled 1000`, but you can use `--scaled 10_000` to make this faster and mimic the `overlap` analysis above.

```
sourmash scripts manysearch GCF_013306015.sig SRR29654720.sig.zip -k 31 --scaled 1000 --output GCF_013306015-x-SRR29654720.k31-sc1000.manysearch.csv
```
