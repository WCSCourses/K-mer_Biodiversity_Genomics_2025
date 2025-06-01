# 2. Comparing genomes metagenomes using FracMinHash sketches

**This follows from the [Introduction to FracMinHash sketching for sequence comparison](https://hackmd.io/@bluegenes/H1ItsC_fxx) module. Please see the installation instructions there for installing sourmash and other dependencies.**

A common application for sketching is detecting genomes of interest within metagenomes to determine presence/absence of organisms of interest, such as detecting pathogenic microbes in wastewater metagenomic sequencing.

[sourmash](sourmash.readthedocs.io) FracMinHash sketching is well-suited to conducting these comparisons at scale. This module demonstrates that and includes a comparison of FracMinHash results vs alignment (minimap2). Finally, we discuss finding genomes(s) in metagenomes at scale.

## Genome x Metagenome comparisons
In this section, we will compare a query genome to a metagenome sample. The goal is to determine if the query genome is present in the metagenome and, if so, how much of it is present.

Datasets we will use:
- **Query microbial genome**: [GCA_013306015 "FJAT91.F50"](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013306015.1) _Ralstonia nicotianae_ / _Ralstonia psuedosolanacearum_.A
- **Metagenome**: Tomato rhizosphere metagenome, collected from a field in China, treated with organic fertilizer and inoculated with _Ralstonia pseudosolanacearum_ to examine disease progression. It was collected during the fruiting stage of tomato growth.
  - **Metagenome sample**: [SRR29654720](https://www.ncbi.nlm.nih.gov/sra/?term=SRR29654720) from the [BioProject PRJNA1127303](https://www.ncbi.nlm.nih.gov/bioproject/1127303).
  - **Associated paper**: [Risk assessment of antibiotic resistance genes in rhizosphere soil during tomato growth under bio-control bacterial inoculation](https://www.sciencedirect.com/science/article/pii/S0959652625002616)
    > The study aims to investigate the effects of different fertilization practices during various growth stages of tomatoes on the abundance and dynamics of antibiotic resistance genes (ARGs) in the soil. Understanding the mechanisms behind these impacts is critical, as ARGs pose a significant threat to public health by promoting the development and spread of antibiotic-resistant bacteria. By examining the influence of organic, inorganic, and combined fertilization methods on ARG profiles, as well as considering factors such as microbial community composition and soil physicochemical properties, this research seeks to elucidate how fertilization strategies might mitigate or exacerbate the propagation of ARGs in agricultural environments. The findings from this study will provide valuable insights for developing sustainable farming practices that minimize the risk of ARG proliferation, thereby contributing to the broader effort of combating antibiotic resistance.

These files should already be downloaded to your working directory. If you want to rerun this at a later date, you can download by uncommenting and running the following:
```
#curl -JLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/306/015/GCF_013306015.1_ASM1330601v1/GCF_013306015.1_ASM1330601v1_genomic.fna.gz
#curl -JLO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/020/SRR29654720/SRR29654720_1.fastq.gz
#curl -JLO ftp://ftp.sra.ebi.ac.uk/vol1/fAastq/SRR296/020/SRR29654720/SRR29654720_2.fastq.gz

```

## Sketch the microbial query genome

We will sketch this genome with the 'abund' parameter. This allows us to compare k-mer abundances between the two sketches.
```
sourmash sketch dna -p k=31,scaled=1000,abund GCF_013306015.1_ASM1330601v1_genomic.fna.gz \
                     --name "GCF_013306015.1 Ralstonia nicotianae FJAT91" -o GCF_013306015.sig
```
<!---
This is the genome gather finds with the full gtdb database:
```
sourmash sketch dna -p k=31,scaled=1000,abund GCF_013375735.1_ASM1337573v1_genomic.fna.gz  --name "GCF_013375735 Ralstonia nicotianae FJAT91-F8" -o GCF_013375735.sig
```
--->


Now sketch the metagenome sample:
```
sourmash sketch dna -p k=31,scaled=1000,abund SRR29654720_1.fastq.gz SRR29654720_2.fastq.gz \
                    --name SRR29654720 -o SRR29654720.sig.zip
```
> Here we are using the `sig.zip` extension, which is a compressed archive of the signature file. This format is particularly useful for large datasets and for storing multiple signatures in a single file (e.g. databases). Note we are also including the `abund` option, which will include k-mer abundances.


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
Here we see that the Jaccard similarity is quite small, which makes sense given the size difference between the query genome and the metagenome. However, the containment values are more informative: the query genome is contained in the metagenome at a level of 97%, meaning that 96% of the unique k-mers in the genome are found in the metagenome. This is a high level of containment, indicating that the query genome is likely present in the metagenome.

Even more interesting are the Abundance-weighted similarity values. The "second contained in first (weighted)" value of 36% indicates that nearly 36% of the abundance-weighted k-mers in the metagenome are k-mers from this genome. This suggests that the query genome was present in the community at really high abundance, which is consistent with the fact that this sample was innoculated with _Ralstonia pseudosolanacearum_.


## Finding genomes in metagenomes at scale

In this trivial example, we knew we would find this genome in this metagenome. However, sourmash can search for a genome of interest in arbitrarily large collections of metagenomes. We will explore this via web-based search of over 1.1 million metagenomes in the Sequence Read Archive at [branchwater.jgi.doe.gov](https://branchwater.jgi.doe.gov/)

1. **Go to the branchwater website, and click on the 'Examples' tab.**

![image](https://hackmd.io/_uploads/HkfeGXKMxg.png)

2. **Select 'Ralstonia solanacearum'**

![image](https://hackmd.io/_uploads/rkgtGXtGeg.png)

3. **Click "Submit".**

4. **Wait for Results**

![image](https://hackmd.io/_uploads/HkpHX7tGex.png)

5. **Try filtering for a minimum of 98% cANI (0.98)**

The default parameters may return results with genus-level similarity (k=21). To narrow the the results, you can select on the reported cANI value to ensure high similarity with your query. Type `0.98` in the cANI 'Min' box and hit 'Enter'/'Return' on your keyboard. You should see that now the list of results is reduced to just two SRA metagenomees, one from the United States, and a second within a "synthetic metagenome" from Switzerland.

You can also try subsetting in other categories. For example, type 'soil' in the 'organism' box and hit "Enter"/"Return". This will limit to just metagenomes with "soil" in the SRA metadata column. Alternatively, you can download the full CSV and do any subselection on the results file afterwards.

### Notes and Limitations

Branchwater search has a few limitations at the moment. First, it uses a 'scaled' value of 1000, meaning that it will not be useful for finding small queries < 10kb, and may even struggle with queries <50kb. Second, the query sketch is limited to 5Mb, which means you can't query with e.g. metagenome samples (though most genomes should be fine). 

The branchwater search is web-based because the underlying search database is over 1Tb in size, meaning it is impractical to provide the metagenome sketches for download. However, you can build sketch databases with your own local metagenomes and search them quickly using the sourmash [branchwater](https://github.com/sourmash-bio/sourmash_plugin_branchwater) plugin via `manysearch`. For very large databases, we recommend indexing the files using the `rocksdb` index type. We will use this in the next hands-on section.


