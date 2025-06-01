# 1. An introduction to FracMinHash sketching for sequence comparisons

## Prerequisites

This module assumes you have a basic understanding of k-mers. Here, we will use FracMinHash sketching to conduct fast yet accurate comparisons of genomic datasets. Specifically, we will explore the following comparisons with k-mers and FracMinHash sketching:
    - Jaccard
    - Containment
    - ANI estimation
    - Comparing many genomes

Software used in this session:
- [sourmash](sourmash.readthedocs.io) 

_Note: Here, we will exclusively use the sourmash command line interface. However, sourmash has a robust python API and an experimental Rust API that are also available for use. It also contains a series of software plugins for even faster and more scalable analysis; see sourmash.readthedocs.io for documentation._

## Installations

If you are running this on your own computer, you will need to install [sourmash](https://github.com/sourmash-bio/sourmash) via conda. There are instructions for setting up conda in the [README.md](README.md). We will also install the `sourmash_plugin_betterplot` plugin, which allows us to visualize some of these k-mer comparisons.

Install sourmash via conda
```
conda install -c bioconda sourmash
```

Install the sourmash `betterplot` plugin for enhanced plotting functions
```
pip install sourmash_plugin_betterplot
```

Now check the installations:
```
sourmash info -v
```

## Pairwise comparisons of genomes

First, we'll download a couple of genomes and sketch them using sourmash. This will allow us to compare them using various methods. We will use genomes from the _Ralstonia solanacearum_ species complex, which is a group of bacteria that cause 'bacterial wilt' and 'brown rot' disease in plants. Pandemic lineages of this species complex have caused significant agricultural losses worldwide. RSSC is the subject of extensive ongoing research to understand its evolution, host range, and spread. We'll use the following genomes:

pandemic lineage:
- [GCF_002251655.1 "UW551"](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002251655.1) is a member of the pandemic lineage (_Ralstonia solanacearum sequevar 1 subclade of phylotype IIB_).

comparison genomes:
- [GCF_000009125.1 "GMI1000"](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009125.1) is a member of the _Ralstonia pseudosolanacearum_ species, which is closely related to _Ralstonia solanacearum_.
- [GCF_001644795.1 "CIP120"](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001644795.1) is a member of the _Ralstonia solanacearum_ species, but part of Phylotype IIA (not part of the pandemic lineage).
- [GCF_001587155.1 "IBSBF1503"](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001587155.1) is another member of the _Ralstonia solanacearum_ Phylotype IIB (closely related strain).

These genomes can be downloaded from NCBI or ENA. For this tutorial, they have been included in the `course_data_2025` folder of the github repository, which you should have cloned to your local machine.


### Setup folder for this module and link in the course data
```
mkdir -p module8_sketching
ln -s ../course_data_2025/module8_sketching/* module8_sketching/
```

### If you do need to download the genomes:
```
curl -JLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/251/655/GCF_002251655.1_ASM225165v1/GCF_002251655.1_ASM225165v1_genomic.fna.gz
curl -JLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/125/GCF_000009125.1_ASM912v1/GCF_000009125.1_ASM912v1_genomic.fna.gz
curl -JLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/644/795/GCF_001644795.1_ASM164479v1/GCF_001644795.1_ASM164479v1_genomic.fna.gz
curl -JLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/587/155/GCF_001587155.1_ASM158715v1/GCF_001587155.1_ASM158715v1_genomic.fna.gz
```

### Sketch the genomes

### Sketch the UW551 pandemic lineage genome
```
sourmash sketch dna -p k=31,scaled=1000 GCF_002251655.1_ASM225165v1_genomic.fna.gz \
                     --name "GCF_002251655.1 Ralstonia solanacearum IIB UW551" -o GCF_002251655.sig
```
> Here, we use a k-mer size of 31 and a scaled value of 1000. The `--name` option allows us to assign a name to the signature, which will be useful for identification later.

You should see output like this:

```
== This is sourmash version 4.9.0. ==
== Please cite Irber et. al (2024), doi:10.21105/joss.06830. ==
 
Computing signatures for files: GCF_002251655.1_ASM225165v1_genomic.fna.gz
Computing a total of 1 signature(s) for each input.
... reading sequences from GCF_002251655.1_ASM225165v1_genomic.fna.gz
... GCF_002251655.1_ASM225165v1_genomic.fna.gz 2 sequences
calculated 1 signature for 2 sequences taken from 1 files
saved 1 signature(s) to 'GCF_002251655.sig'
```

### Examine the signature file with the `sourmash describe` command:
```
sourmash sig describe GCF_002251655.sig
```
This will show you the signature metadata, including the name, k-mer size, and other parameters used to create the signature.

```
signature filename: GCF_002251655.sig
signature: GCF_002251655.1 Ralstonia solanacearum IIB UW551
source file: GCF_002251655.1_ASM225165v1_genomic.fna.gz
md5: 17c4bbc60a07474e73275b998b060b2e
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=0
size: 5209
sum hashes: 5209
signature license: CC0
```
> Here, we did not track the abundance of k-mers, so the `track_abundance` parameter is set to 0. We will look at abundance-weighted comparisons later in this module.

You can also look at the signature file itself, which is a zip file containing the k-mer signatures. It is currently stored as `json` format and can be viewed with any text editor, or the `cat` command in the terminal.


### Now sketch the comparison genomes
We'll use the same command here, but with different input files and names for the signatures. Note, that there are multiple ways to automate this process, either via bash scripting or using the `sourmash sketch fromfile` command. For simplicity, we will sketch each genome individually here.

```
sourmash sketch dna -p k=31,scaled=1000 GCF_000009125.1_ASM912v1_genomic.fna.gz \
                     --name "GCF_000009125.1 Ralstonia pseudosolanacearum_I_GMI1000" -o GCF_000009125.sig
sourmash sketch dna -p k=31,scaled=1000 GCF_001644795.1_ASM164479v1_genomic.fna.gz \
                     --name "GCF_001644795.1 Ralstonia solanacearum_IIA_CIP120" -o GCF_001644795.sig
sourmash sketch dna -p k=31,scaled=1000 GCF_001587155.1_ASM158715v1_genomic.fna.gz \
                     --name "GCF_001587155.1 Ralstonia solanacearum_IIB_IBSBF1503" -o GCF_001587155.sig

```

### Comparing a pair of genomes

First, let's compare the pandemic lineage genome (UW551) with the _R. pseudosolanacearum_ genome (GMI1000). We will use the `sourmash sig overlap` command to get a comprehensive set of comparison information.

```
sourmash sig overlap GCF_002251655.sig GCF_000009125.sig -k 31
``` 
`overlap` first prints descriptive information for each signature, then it prints a summary of the comparison, including the Jaccard similarity, containment ratios, and average containment ANI. The output will look something like this:
```
--- Similarity measures ---
jaccard similarity:          0.06967
first contained in second:   0.13553 (cANI: 0.93757)
second contained in first:   0.12540 (cANI: 0.93522)
average containment ANI:     0.93639


--- Hash overlap summary ---
number of hashes in first:   5209
number of hashes in second:  5630

number of hashes in common:  706
only in first:               4503
only in second:              4924
total (union):               1013
```
> We can look first to the **Hash overlap summary** section, which tells us how many k-mers are present in each genome and how many are shared between them. The `number of hashes in first` and `number of hashes in second` values indicate the total number of unique k-mers in each genome. The `number of hashes in common` value indicates how many k-mers are shared between the two genomes, while the `only in first` and `only in second` values indicate how many k-mers are unique to each genome.

These hash values are used to calculate the Jaccard similarity and containment, above. Recall the `jaccard similarity` value is the ratio of the number of shared k-mers to the total number of unique k-mers in both genomes, which gives us a measure of how similar the two genomes are based on their k-mer content.

### Visualize this comparison
The `betterplot` plugin for sourmash allows us to visualize this comparison using a Venn diagram. We can use the `sourmash scripts venn` command to create this visualization.

```
sourmash scripts venn GCF_002251655.sig GCF_000009125.sig -k 31 \
                 --name1 "Rsol UW551" --name2 "Rpseudosol GMI1000" \
                 -o venn_GCF_002251655-x-GCF_000009125.png
```

![venn_GCF_002251655-x-GCF_000009125](https://hackmd.io/_uploads/B1yj5gtzee.png)

<!---
![Venn diagram of GCF_002251655 and GCF_000009125](module8_sketching/_static/venn_GCF_002251655_GCF_000009125.png)
-->

Since these are members of different species, we expect to see a low Jaccard similarity and containment. The Venn diagram shows that there are some shared k-mers, but the majority of k-mers are unique to each genome. We can use the containment to estimate the average ANI (average nucleotide identity) between the two genomes, which is a measure of how similar the genomes are at the nucleotide level. Here, we get ~93% ANI.

### How do we know the sketching isn't impacting the results?
Since sourmash allows for flexible sketching parameters, we can check how the scaled value impact the results. We can do this by running the same comparison with different k-mer sizes and scaled values. For example, we can run the same comparison with a scaled value of 10.

For this, we will re-sketch both genomes with a scaled value of 10, which will keep more k-mers in the sketch, and then compare them again.
```
sourmash sketch dna -p k=31,scaled=10 GCF_002251655.1_ASM225165v1_genomic.fna.gz \
                    --name "GCF_002251655.1 Ralstonia solanacearum IIB UW551" -o GCF_002251655-scaled10.sig
sourmash sketch dna -p k=31,scaled=10 GCF_000009125.1_ASM912v1_genomic.fna.gz \
                    --name "GCF_000009125.1 Ralstonia pseudosolanacearum_I_GMI1000" -o GCF_000009125-scaled10.sig
sourmash sig overlap GCF_002251655-scaled10.sig GCF_000009125-scaled10.sig -k 31
```
> Note: You'll see a warning as you're sketching that we recommend scaled values of >=100. Low scaled values such as 10 or even 1 (keep all k-mers in the "sketch") are completely valid, but sourmash is not optimized for large-scale comparisons with this strategy. Higher resolution will create large\(r) sketch files that are a bit slower to compare!
```
--- Similarity measures ---
jaccard similarity:          0.06802
first contained in second:   0.13265 (cANI: 0.93691)
second contained in first:   0.12252 (cANI: 0.93452)
average containment ANI:     0.93572


--- Hash overlap summary ---
number of hashes in first:   520788
number of hashes in second:  563839

number of hashes in common:  69081
only in first:               451707
only in second:              494758
total (union):               1015546
```
The largest difference here is in the total number of hashes we are counting and comparing; the Jaccard similarity and containment ratios are very similar to the previous comparison, indicating that the sketching parameters did not significantly impact the results. The average containment ANI is also very similar, at ~93.6%. This shows that the sketching minimally impacts our similarity comparisons.


### Comparing the pandemic lineage genome with a member of the same species
Now, let's compare the _R. solanacearum_ pandemic lineage genome (UW551) with the other _R. solanacearum_ genome (CIP120). We will use the same `sourmash sig overlap` command to get a comprehensive set of comparison information.

```
sourmash sig overlap GCF_002251655.sig GCF_001644795.sig -k 31
```

We see:
```
--- Similarity measures ---
jaccard similarity:          0.23171
first contained in second:   0.39086 (cANI: 0.97015)
second contained in first:   0.36266 (cANI: 0.96781)
average containment ANI:     0.96898


--- Hash overlap summary ---
number of hashes in first:   5209
number of hashes in second:  5614

number of hashes in common:  2036
only in first:               3173
only in second:              3578
total (union):               8787
```

We can see that the Jaccard similarity is higher than in the previous comparison, indicating that these two genomes are more similar. The containment ratios are also higher, with an average ANI of ~97%. We can clearly see these genomes are more closely related, as they share a larger proportion of k-mers.

We can also visualize this comparison using the `sourmash scripts venn` command:

```
sourmash scripts venn GCF_002251655.sig GCF_001644795.sig -k 31 \
                 --name1 "Rsol UW551" --name2 "Rsol CIP120" \
                 -o venn_GCF_002251655-x-GCF_001644795.png
```

![venn_GCF_002251655-x-GCF_001644795](https://hackmd.io/_uploads/rkeFxsgFMlx.png)

Finally, we can compare the pandemic lineage genome (UW551) with another member of the same strain (IBSBF1503). This will allow us to see how similar these two genomes are within the same strain.

```
sourmash sig overlap GCF_002251655.sig GCF_001587155.sig -k 31
```

We see:

```
--- Similarity measures ---
jaccard similarity:          0.33882
first contained in second:   0.51776 (cANI: 0.97899)
second contained in first:   0.49504 (cANI: 0.97757)
average containment ANI:     0.97828

--- Hash overlap summary ---
number of hashes in first:   5209
number of hashes in second:  5448

number of hashes in common:  2697
only in first:               2512
only in second:              2751
total (union):               7960
```
> We can see that the Jaccard similarity is even higher than in the previous comparison, indicating that these two genomes are very similar. The containment ratios are also higher, with an average ANI of ~98%. This is expected, as these two genomes are quite closely related. If we were to add a genome in the exact pandemic lineage, we would see the similarity jump to >99%.

We can visualize this comparison using the `sourmash scripts venn` command:

```
sourmash scripts venn GCF_002251655.sig GCF_001587155.sig -k 31 \
                 --name1 "Rsol UW551" --name2 "Rsol IBSBF1503" \
                 -o venn_GCF_002251655-x-GCF_001587155.png
```

![venn_GCF_002251655-x-GCF_001587155](https://hackmd.io/_uploads/rkZBoetzll.png)


### Comparing many genomes at once

`sourmash sig overlap` provides a wealth of information for interactive pairwise comparisons, but it can be cumbersome to use for many genomes. Instead, we can use the `sourmash sig compare` command to compare many genomes at once. This command will compute pairwise comparisons between all signatures in a given directory and output a matrix of pairwise Jaccard similarities, containment ratios, or average containment ANI values.

```
sourmash compare -k 31 GCF_002251655.sig GCF_000009125.sig GCF_001644795.sig GCF_001587155.sig -o RSSC.jaccard
```
This will create a file called `RSSC.np` that contains the pairwise Jaccard similarities between all signatures. You can also specify the `--containment` or `--ani` options to compute containment ratios or average containment ANI values, respectively.
The matrix will also print to the screen, which will look like this:
``` 
0-GCF_002251655.1...	[1.    0.07  0.232 0.339]
1-GCF_000009125.1...	[0.07  1.    0.067 0.068]
2-GCF_001644795.1...	[0.232 0.067 1.    0.225]
3-GCF_001587155.1...	[0.339 0.068 0.225 1.   ]
``` 
We can then visualize this matrix using the `sourmash plot` command, which will create a heatmap of the pairwise Jaccard similarities.
```
sourmash plot RSSC.jaccard
```
This will create a file called `RSSC.jaccard.matrix.png` that contains the heatmap of the pairwise Jaccard similarities. The heatmap will show the pairwise Jaccard similarities between all signatures, with darker colors indicating higher similarities.

![RSSC.jaccard.matrix](https://hackmd.io/_uploads/H1E8jxtzeg.png)


## Concluding thoughts

The comparisons here are primarily to get you thinking about how k-mer sketches can be used to compare sequence datasets. Since we're just using a small number of genomes, it would be simple to compare these with more precise measures. However, this sketching strategy enables comparisons between many different types of datasets, including finding genomes of interest in metagenomes and comprehensive metagenome characterization, which we'll explore in the next hands-on modules.


## Scaling up these comparisons

> Here we have used primarily base `sourmash` commands, but if you want to scale to millions of pairwise comparisons, it's best to install the [branchwater](https://github.com/sourmash-bio/sourmash_plugin_branchwater) plugin for sourmash, which allows you to compute comparisons in parallel. Here we would instead use the `sourmash scripts pairwise` command. The [betterplot](https://github.com/sourmash-bio/sourmash_plugin_betterplot) plugin also contains a few different visualizations you may want to explore.
