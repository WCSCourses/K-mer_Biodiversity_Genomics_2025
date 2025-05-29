## Prerequisites

The starting point of this module's practicals is understanding (at least kind of) what are k-mers. If you are puzzled by the concept, do ask any of the instructors to go over it again with you, it is rather important. We will build a k-mer spectrum and try to get an intuitive understanding of what it means as well as learn how to use GenomeScope to fit genome models.

Software used this session:
 - FastK
 - R
 - genomeScope

## k-mer counting

k-mer counting is nothing else then calculating the number of times all the k-mers occur in a dataset. It is a simple task, but there are many k-mer counters out there, why is that? Well, it is a task that needs to be done a lot, so having the counting optimise is rather important, but on top of that, k-mer counters come bundled in tools that allow to do a lot more than just count k-mers, they can do all sort of tranformations and operations. 

If you want to know details about speed of individual k-mer counters, check this [benchmarking blogpost](https://kamilsjaron.github.io/posts/2021/12/2021-12-14-testing-fastk-on-a-difficult-genome/
), but in general, FastK and KMC are very fast, but a bit harder to use. Jellyfish is in comparison to them a bit slow, but very easy to use. Both FastK and KMC are a lot more powerful and if you have more specific k-mer tasks, you might want to use one of the more advanced tools.

In this workshop, we will use FastK mostly because it is very fast as well as compatible with many common k-mer tools we will learn about in this workshop.

### building the first k-mer spectrum using FastK

Either download the dataset, or use the pre-downloaded dataset in the gitpod.

```
mkdir -p data/Scer && wget -P data/Scer https://sk13.cog.sanger.ac.uk/SRR3265401_{1,2}.fastq.gz
```

This is a short read dataset of a yeast. It should be ~800MB of data. You can check how big are the files and if the content of files looks alright...

```
ls -lh data/Scer/
zcat data/Scer/SRR3265401_1.fastq.gz | head -8
```

Alright, if all looks good. You can calculate a k-mer database, the parameters you need to specify are the minimal k-mer coverage to compute (`-t1`), the choice of k (`-k31`) resources that can be used by this operation (`-M120 -T4`) and the name of the calculated k-mer database (`-Ndata/Scer/FastK_Table`). All together looks like this

```
FastK -v -t1 -k31 -M120 -T4 data/Scer/SRR3265401_[12].fastq.gz -Ndata/Scer/FastK_Table
```

Look what files were generated

```
ls -lh data/Scer/
```

You should see `FastK_Table.hist`, `FastK_Table.ktab` and several dotted files (which store the actual k-mer database). You will need to use various commands from the FastK suite to work with it. One of the functionalities is to extract the k-mer histogram as follows. Notice the histogram is streamed to `stdout`, to get it into a file, you need to redirect the stream using `>` (which is not anything specific to k-mers).

```
Histex -G data/Scer/FastK_Table > data/Scer/SRR3265401_k31.hist
```

Congratulations, you made your first k-mer histogram. You can take a look how it looks like

```
head data/Scer/SRR3265401_k31.hist
```

You should see two column tab-delimited file. The first column is the coverage, and the second is the number of distinct k-mers in the dataset that have that coverage. This second column is often referred to as frequency.

## plot k-mer spectrum

Now that we know how to build histograms, let's try to interpret them and extract as much information about our genomes as we can from them Try to load the k-mer histogram and plot it in **R**.

```
kmer_hist <- read.table('data/Scer/SRR3265401_k21.hist', col.names = c('cov', 'freq'))
plot(kmer_hist$cov, kmer_hist$freq, type = 'l', xlab = 'Coverage', ylab = 'Frequency', main = 'Naive plot')
```

<details> 
  <summary>What do you see? </summary>

Most of sequencing errors generate k-mers with very small coverage and as a consequence, the coverages 1, 2, 3... are usually by far with the highest number of k-mers, therefore the y-axis is streched mostly because of the sequencing errors. The k-mer histogram reported by FastK contains coverages all the way till 1000x, which is stretching the x-axis to repetitive part of the gneome we are not interested to see. The part you are really interested is are those tiny peaks in the bottom right corner. 

</details>


<details> 
  <summary>Can you plot it in a more useful way? How many peaks do we see? </summary>

```
cov_range <- c(1:120) # there are some smart ways too, but you can also just look where the peaks approximately are 
ymax <- max(kmer_hist$freq[kmer_hist$cov > 10]) # here I say.

plot(kmer_hist$cov[cov_range], kmer_hist$freq[cov_range], ylim = c(0, ymax), type = 'l', xlab = 'Coverage', ylab = 'Frequency', main = 'Rescaled plot')
```

Ha, peaks! Three in fact (one is blended with others). These peaks correspond to k-mers in 1, 2 and 4 genomic copies in the genome. The first peak is centered around 1n coverage, and the other two are centered around 1n coverage multiplied by their respective copy numbers. **The peaks in k-mer spectrum are evenly spaced!**. The heights of the individual peaks give us an idea about the amount of k-mers in this copy number.

</details>

This genome is complicated, we well later on return to this spectrum to analyse it in more details. In the meantime you should appreciate:
- the most of the distinct k-mers in the dataset are sequencing errors
- seqencing errors are relatively easy to see and separate from the rest of the k-mers
- k-mer spectra is made of several peaks
- peaks are evenly spaced with centers of individual peaks around copy_number * 1n_coverage

## fit a simple GenomeScope model

We have an intuition what is a k-mer spectrum is - distribution of k-mer coverages. We also intuitively understand it carries information about genome size, heterozygosity and coverage. We can fit all these using genome profiling techniques, such as GenomeScope.

Let's start with a simple spectrum before analysing the one we visualised before. _Timema_monikensis_ is a completelly homozygous diploid organism - it's a parthenogenetic stick insect.

```bash
cd course_data_2025/histograms/homozygous_diploid
genomescope.R -i Timema_monikensis_k21.hist -o . -n timema -k 21
```

This is a simple fit, with default parameters the model converge to a genome size that is in the right ballpark (a bit over 1Gbp) and the species is predicted to be homozygous. This is what we excpect when things are easy.

## Fitting more models

The directory `course_data_2025/histograms` contains several k-mer spectra. Try to fit them one by one and make sense out of the fits. Many of the fits will require non-default parameters, you can look at the GenomeScope help (`GenomeScope.R --help`) and try to figure out how to do it, underneath there is a walkthrough with explanations.

The recommended order is begonia, strawberry, human, and toad.

### begonia

You can try also a more heterozygous organism, such as begonia (`begonia_simple_heterozygous_diploid/Begonia_luxurians_k21.hist`)

<details> 
  <summary>Does the model look right? What is the heterozygosity level? Do you think it is a high or low for a plant?</summary>

This one converges well with the default

```bash
cd course_data_2025/histograms/begonia_simple_heterozygous_diploid 
genomescope.R -i Begonia_luxurians_k21.hist -o . -n begonia -k 21
```

</details>

### strawberry

Try to fit a strawberry k-mer spectrum now (`course_data_2025/histograms/convergence_problem`). It is again a diploid organism.

<details> 
  <summary>What is the monoploid genome size of a strawberry (`Fragaria`), did the model converge to something sensible? How do we know if all the peaks are part of the genome or not? Is there a parameter of GenomeScope you could explpoit to get a nicer fit?</summary>

The default run converges to a heterozygous and by far too small genome. If you know about strawberries, both those will be already suspicious, but there is more to it - there is clearly one peak that is 1/2 of the coverage of the first peak included in the model. This is strongly suggesting it is part of the same genome and therefore should be part of the model. Seems the 1n coverage was estimated to be 2x as high as it should have been. You can give GenomeScope prior on coverage using parameter `-l`.

```bash
cd course_data_2025/histograms/convergence_problem
genomescope.R -i Fragaria_iinumae_k21.hist -k 21 -p 2 -o . -n Fragaria # does not converge well
genomescope.R -i Fragaria_iinumae_k21.hist -k 21 -p 2 -o . -n Fragaria_fixed -l 140 # yep, this is it
```

</details>

### Human

Another simple case is the Human genome. This is specficially one of the sequencing runs of H002, which is a male sample. Try to fit the model to it (`course_data_2025/histograms/human`).

```bash
cd course_data_2025/histograms/human
genomescope.R -i m84005_220919_232112_s2.hifi_reads.bam.21.kmc.nozero.hist -o . -n hsap -k 21
```

what do you see? Do you trust the genome size / heterozygosity estimates? What could possibly be throwing them off?

### Toad

The toad dataset 

```bash
cd course_data_2025/histograms/toad_very_large_genome
```
The two k-mer histograms in this directory are corresponding to a k-mer histogram build with def
ault parameters of kmc (`bombina_naive_kmc_k21.hist`), which don't count coverage over 10,000x and all k-mers with the coverage higher are simply reported as 10,000x k-mers. And a histogram with all the k-mers counted (`bombina_sp_k21.hist`).

<details> 
  <summary>What do you think would be the consequences of fitting models to one or the other for the gneome size estimation?</summary>

```
genomescope.R -i bombina_naive_kmc_k21.hist -o . -n naive
```

Hmmm, the model looks quite alright, but the genome is much smaller than what one would expect for a firebelly toad. Try to fit the same model to a histogram with all the k-mers counted (`bombina_sp_k21.hist`).

```
genomescope.R -i bombina_sp_k21.hist -o . -n kmc_full
```

</details>

## Polyploids

Let's return to our very first k-mers spectrum Saccharomyces cerevisiae (`Scer`). We can either visualise the spectrum again, or simply run genomescope. If you do that later, you must remember that if the model (the black line) mismatches the data (the blue histogram) you need to ingnore all the estimated numbers and focus only on the shape of the histogram itself. And it looks like this

![image](https://github.com/user-attachments/assets/29b4100b-70aa-4a45-a93b-6fdd108e2ccc)

To determine the ploidy of this species, it would be helpful to know if there are k-mers that are similar to each other - those that are coming from heterozygous sites. Smudgeplot is a tool that implement this approach. The first step is to chose a theshold for excluding low frequencing k-mers that will be considered errors with the aim of processing mostly real genomic k-mers. That choice is not too difficult to make by looking at the k-mer spectra. 

In this example, a meaningful error threshold would be 10. As a rule of thumb, no dataset should have this threshold much below <10 unless it's a very clean sequencing run (you can separate error and genomic peaks). Furthermore, it is not the end of the world if we lose a bit of the real genomic k-mers (as long as there is enough signal). However these are just some gudances, what is sensible really depends on each individual datasets!

Now we can finally find k-mer pairs by running `smudgeplot.py hetmers`. This command will interanlly call a C-kernel optimised for the searched designed by Gene Myers. We specify `-L 10` (the coverage threshold for genomic k-mers) and use 4 threads (`-t`). The parameter `-o` just specifies the pattern for the output files

```
smudgeplot hetmers data/Scer/FastK_Table -L 10 -t 4 -o data/Scer/smudgeplot_pairs
```

This algorithm found all k-mer pairs (`A` and `B`) that are exactly one nt away from each other and form a unique k-mer pair. What we generated now is sort of a 2D histogram of k-mer coverages. For each corresponding coverage combination (of the k-mer A and k-mer B) we have a respective frequency that is seen within the dataset.

```
head data/Scer/smudgeplot_pairs.smu
10      10      300
10      11      566
10      12      530
11      11      254
10      13      576
11      12      552
10      14      552
11      13      546
12      12      328
10      15      572
```

You could play with 2D visualisation or any other tranformation, or you can pass this 2D histogram to smudgeplot inference algorithm. That is called as follows

```
smudgeplot all data/Scer/smudgeplot_pairs.smu -o data/Scer/smudgeplot
```

<details>
<summary><b> What do you think about ploidy of this species? </b></summary>
Tetraploid, specifically of `AAAB` type. Notably, this constitution does not necesarily indicate one of the haplotypes is more divergent to others because we the B k-mers can be on different haplotype for each individual k-mer pair possibly making the 4 haplotypes equidistant. We can refute a hypothesis of two and two haplotypes that are closer to each genomes as subgenomes, as those would generate prominent AABB smudge, hence the genome is quite possibly autotetraploid.  
</details>
