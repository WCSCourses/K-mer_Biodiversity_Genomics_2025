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

## fit genomeScope model

We have an intuition what is a k-mer spectrum is - distribution of k-mer coverages. We also intuitively understand it carries information about genome size, heterozygosity and coverage. We can fit all these using genome profiling techniques, such as GenomeScope.

Let's start with several simpler spectra before analysing the spectrum above.

## Genome size estimation



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







FIt a model to a simple k-mer spectrum
Expectations over homozygous and heterozygous genomes

Appreciate what could and could not be fit
Estimates of genome sizes
Common pitfalls due to wrong convergence, or limited k-mer counting

See spectra that do NOT look like things went well
Signatures of different problems
Low coverage
Contamination (two or multiple sources)
Sequencing biases
Poor model convergence / wrong model
Genome Size estimate
Explain the magical last position
