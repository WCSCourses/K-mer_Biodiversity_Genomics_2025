## Prerequisites

The starting point of this module's practicals is familiarity with genome profiling using GenomeScope. It's also expected the user will know how to compute a k-mer spectrum, in this practicals we will use precomputed spectra to save time.

Software used this session:
 - R
 - GenomeScope

## Details of GenomeScope 

In the following section we will construct our own GenomeScope-like model in R. This will help us apprecite more what the model is exactly predicting and how it works. Subsequentially, we will use the newly acquired skill to create alternative models.

### heterozygosity

To construct out first model we use idealised dataset from californian stick insect.

```bash
cd course_data_2025/histograms/advanced/Timema/
```

All the subsequent code is R. Load the spectrum and visualise it.

```R
ideal_kmer_spec <- read.table('Timema_cristinae_kmer_k21_simplified.hist', col.names = c('cov', 'freq'))
plot(ideal_kmer_spec$cov, ideal_kmer_spec$freq, type = 'l', xlab = 'Coverage', ylab = 'Frequency', main = 'ideal')
```

You can see that there are no errors in the spectrum (flat right side of the spectrum), not there are any repetitive seqeunces (flat right side of the spectrum). There is no point in fitting more than two peaks, so the model we will fit will be something like

```
frequency ~ 1n peak +
            2n peak
```

We could choose from various distributions to model the two peaks, but it just happens, that negative binomial, fits sequencing coverage very well. Negative bionomial has usually two parameters stopping parameter and the success probability, however, that is not a convinient notion for us, we will therefore use a reparametrised version of negative binomial with mean 'mu', and 'size', the dispersion parameter. So in our model, we will fit mean, and the dispersion parameter and the model can look like this

```
Probability ~ α NB(  kcov, dispersion) +
              β NB(2*kcov, dispersion2)
```

Now, the α and β will be the contributions of the two peaks. Notice that we fit the mean of the first peak with the same parameter as the mean of the second peak s, just multiplied by a factor of 2. That is because the heterozygous loci have 1/2 of the coverage of the homozygous loci and therefore, we actually desire to co-estimate the haploid coverage jointly for both of the at the same time.

There are still a few adjustments we need to do to the model to finally make a reasonable fit. In the last model we parametrise the dispersion parameters of the two models independently, that is however unnecessary - with increased coverage, the variance in sequencing depth increases too, but the "overdispersal" scales the same for the whole sequencing run. We therefore use the same trick we used for the coverage and fit the overdispersion as a single parameter scaled by the mean

```
Probability ~ α NB(  kcov,   kcov / overdispersion) +
              β NB(2*kcov, 2*kcov / overdispersion)
```

We can express α and β using probability of a heterozygous nucleotide (r) and the k-mer length (k). Look at the cascade of following expressions

![Screenshot 2023-02-07 at 13 17 51](https://user-images.githubusercontent.com/8181573/217258225-a3927db7-e8e4-484a-9a03-b78f7be0110c.png)

With this we could estimate the heterozygosity using `k` as a constant, while assuming that the heteorzygosity (probability of a nucleotidy being heterozygous) is the same for each nucleotide across the genome. Note there is a one less parameter too, instead of α and β (or `a` and `b`), there is just `r`, but there is one more consideration - while our previous model features large α and β coeficients, this model expects them to be relative proportions of the two peaks (the probability of heterozygous vs homozygous k-mers), which means we need to multiple the joint distribution by the genome length - this will be a one more added parameter to the equation, leading to the same total number as it was in the last case.

```R
y <- ideal_kmer_spec$freq
x <- ideal_kmer_spec$cov
k <- 21

gm2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kcov   / overdispersion, mu = kcov    ) +
                ((1-r)^k)       * dnbinom(x, size = kcov*2 / overdispersion, mu = kcov * 2)) * length,
           start = list(r = 0, kcov = 25, overdispersion = 0.5, length = 1000e6),
           control = list(minFactor=1e-12, maxiter=40))
summary(gm2)
```

So, what do you think? The heterozygosity seems to be in the right ballpark (0.99% by GenomeScope vs 0.98% by our model), K-mer coverage is 27.6x in both models, and the genome size is also remarkably similar ~760Mbp.

Now, try to fit the same model, while changing the starting values (e.g. r = 0.02 or 0.2). Have the starting condition affected the model? 

<details> 
  <summary>What do you think singular gradient means? </summary>


That is when the least-square surface is flat around the initial values, which leads to impossibility to optimize the model, in other words, the initial guesstimates need to be in the right ballpark. 
</details>

Till now we had an idealised case - trimmed of all repetitions and errors. In truth, the unidealised version of the k-mer spectra still looks pretty good. This is how the GenomeScope of that spectrum looks like

![real_linear_plot](https://user-images.githubusercontent.com/8181573/217610706-8a9c1cc1-3b93-408c-8bb7-11efa70123af.png)

So what's the effect of the "unidelised version"? Is the fit still good? Which estimated values get affected the most?

### genome size again

The last time we evaluated the effect of counting repetitive k-mers. Here we would like to present an alternative to counting all the k-mers as well as discussing how variable ploidy levels among chromosomes will affect / change the estimate. 

#### Counting the sum of k-mers instead of the full histogram

The "peaky" section of the k-mer histogram is used to fit the 1n k-mer coverage as well as heterozygosity. However, the right side of the k-mer histogram is used for nothing else than the genome size, which is the sum of all coverages multiplied by their respective frequencies and divided by the 1n coverage and ploidy. Therefore, we don't need to know the full k-mer histogram as long as we keep a sum of coverages*frequencies for all the repetitive k-mers. FastK does this trick by default and for a high repetitive threshold (1000 by default), the reported frequency is representing a sum of all the coverages * frequencies divided by that threshold. 

To get an intuition, load 3 spectra of a toad (`course_data_2025/histograms/toad_very_large_genome`) in R - one calculated by KMC with the default parameters (`bombina_naive_kmc_k21.hist`), one calculated by KMC that includes all the repetitive sequences (`bombina_sp_k21.hist`) and one by FastK that includes this magical position at the position 1000. Use the genome size estimation equation to estimate the genome size for the toad.

<details> 
  <summary>Can you see a difference between FastK and KMC full histogram? And what about truncated (default/naive) histogram by KMC? </summary>

```R
kmc_default <- read.table("course_data_2025/histograms/toad_very_large_genome/bombina_naive_kmc_k21.hist")
fastk_default <- read.table("course_data_2025/histograms/toad_very_large_genome/bombina_sp_FastK_k21.hist") 
kmc_full <- read.table("course_data_2025/histograms/toad_very_large_genome/bombina_sp_k21.hist")   

# you could fit either your simplified genomescope model, or fit regular genomescope to get 1n coverage estimate, or use the fit from yesterday, which I will do here; the only thing you need is the 1n k-mer coverage
cov <- 67.08
ploidy <- 2

# I will add 1e9 to the denominator to make the genome estimate in Gbp (easier to read)
sum(kmc_default[, 1] * kmc_default[, 2]) / ((cov * ploidy) * 1e9)
sum(fastk_default[, 1] * fastk_default[, 2]) / ((cov * ploidy) * 1e9)
sum(kmc_full[, 1] * kmc_full[, 2]) / ((cov * ploidy) * 1e9)
# we see that kmc full and FastK give the same estimates, while kmc default has much smaller one

kmc_default[1000,]
# There are 78124 k-mers with coverage 1000
tail(kmc_default, 1)
# and 5628862 k-mers with coverage 10,000 or greater

fastk_default[1000, ]
# The frequency of coverage 1000 is 574997621; however, this is the magical position

repetitive_kmers <- kmc_full[1000:nrow(kmc_full), ]
sum(repetitive_kmers[, 2]) # there are 66078908 distinct k-mers with coverages >1000
sum(repetitive_kmers[, 2] * repetitive_kmers[, 1])  # and 574,997,620,805 total
# if we divide the total number of k-mers by the coverage of the last magical position in the FastK histogram (1000) ...
sum(repetitive_kmers[, 2] * repetitive_kmers[, 1]) / 1000 
# ... we get 574997621 - the number FastK reports to capture the totality of repetitive k-mers
```

</details>

<details> 
  <summary>Why do you think the GenomeScope and this genome size estimate differ a bit?  </summary>

  That's because GenomeScope also substract the left-side residuals of the model - i.e. what GenomeScope models as an error. If you would like a challenge, you can subtract the errors from the k-mer histograms and see if you can match the genome size estimates perfectly.

</details>

#### Heterogametic Sex genome size

The in-built command `nls` is using Gauss-Newton algorithm by detail. Experimentally, I found that for sex chromosomes it is easier to use Levenberg-Marquardt instead - it practically change very little for fitting the models, but we will use `nlsLM` function from `minipack.lm` package. In fact, it's not just the sex chromosomes, genomescope2 is using Levenberg-Marquardt too.

Now, we want to define a function that will fit a model given **heterozygosity** and `k` while estimating the proportion of the heterogametic genome. It is analogous to the functions we have developed above.

```R
nlsLM_2peak_heterogametic_model_fixed_r <- function(x, y, estKmercov, estLength, heterogameticSize, r, k){
        nlsLM(y ~ (((2*(1-(1-r)^k))  * dnbinom(x, size = kmercov   / overdispersal, mu = kmercov)) +
                   ((1-r)^k)         * dnbinom(x, size = kmercov*2 / overdispersal, mu = kmercov*2)) * (length - hetSize) +
                                       dnbinom(x, size = kmercov   / overdispersal, mu = kmercov) * hetSize,
                  start = list(kmercov = estKmercov, overdispersal = 0.5, length = estLength, hetSize = heterogameticSize),
                  control = list(minFactor=1e-12, maxiter=40))
}
```

Now that we got the model, let's load the data. Because human is well sequenced, we also save all the real values so we have something to compare our fit with...

```R
library('minpack.lm')

HG002 <- read.table('course_data_2025/histograms/human/m84005_220919_232112_s2.hifi_reads.bam.21.kmc.nozero.hist', col.names = c('cov', 'freq'))

# these are the real results for comparison
human_wo_Y <- 3.055e9 # from T2T paper which is without Y
X <- 154259566	
Y <- 62460029 # from complete Y paper
PAR <- 3.03e6
total_human <-	human_wo_Y + Y
male_disomic <- total_human - X - Y + PAR
male_monosomic <- X + Y - 2 * PAR
```

And now we would like to plot it. Here we will define ourself a neat little function so we can do it for any histogram afterwards.

```R
coverage_barplot <- function(bar_heights, bar_positions, xlim = c(0, 0), ylim = c(0, 0), font_size = 1, width = 0.5){

  if (ylim[2] == 0){
      ylim[2] <- max(bar_heights)
  }
  if (xlim[2] == 0){
      xlim = range(bar_positions)
  }
  
  plot(bar_heights, type="n", xlab="Coverage", ylab="Frequency",
       ylim = ylim, xlim= xlim,
       cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  for ( i in 1:length(bar_heights)){
    rect(bar_positions[i] - width, 0, bar_positions[i] + width, bar_heights[i], col = 'deepskyblue', border = F)
  }
}

# plot it
fitting_range <- 6:60 # this is a prior knowledge of where is the sensible peak, normally I would look and see and then specify this variable
coverage_barplot(HG002[1:max(fitting_range), 'freq'], HG002[1:max(fitting_range), 'cov'], ylim = c(0, max(HG002[fitting_range, 'freq']))) 
```

<details> 
  <summary>Alright, table `HG002` contains the histogram, we know how the histogram looks like and what is the sensible coverage range. Throw it in the modeling function we defined with some sensible starting values... Now, using this function, we can try to fit the size of the heterogametic part of the human genome. One consideration is what heterozygosity we will assume. The textbook choice for human would be `0.001` (1 in 1000), but if you would like to see something slightly different that what is in this tutorial, you can make a slightly different choice.
 </summary>

```R
# let's
x <- HG002$cov
y <- HG002$freq
k <- 21
estKmercov <- 15
estLength <- 3e9
estR <- 0.02
max_iterations <- 20
heterogameticSize <- male_monosomic # here we cheat a bit, we can take the real value

het_model2 <- nlsLM_2peak_heterogametic_model_fixed_r(x[fitting_range], y[fitting_range], estKmercov, estLength, heterogameticSize, 0.001, k)

summary(het_model2)

# Parameters:
#                Estimate Std. Error t value Pr(>|t|)    
# kmercov       1.544e+01  2.110e-02  731.63   <2e-16 ***
# overdispersal 3.464e-01  1.698e-02   20.40   <2e-16 ***
# length        2.183e+09  1.090e+07  200.38   <2e-16 ***
# hetSize       1.895e+08  7.751e+06   24.45   <2e-16 ***
```

It seems it converged (it would be good to plot it to make sure). Nevertheless, we have estimates for disomic portion of the genome (meaning autosomes and PAR fitted as the parameter `length`) and the size of the monosomic (heterogametic portions of X and Y fitted as the parameters `hetSize`). We see that the sum of the two is far far below the actual human genome size, that is because we did not include repetitive part of the spectrum in the fit (we fit only k-mers in one or two genomic copies).

</details>

We can plot the model. To do that, We would like to see the individual components of the model. The in-built predict function won't do the job, becuase it would select the whole model at the same time, so we need to define our own functins for this prediction 

```R
predict_disomic_portion_2peak <- function(x, k, kmercov, overdispersal, est_length, heterogameticSize, r){
    (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / overdispersal, mu = kmercov)) +
    ((1-r)^k)         * dnbinom(x, size = kmercov*2 / overdispersal, mu = kmercov*2)) * (est_length - heterogameticSize)
}

predict_monosomic_portion_2peak <- function(x, k, kmercov, overdispersal, est_length, heterogameticSize, r){
    dnbinom(x, size = kmercov   / overdispersal, mu = kmercov) * heterogameticSize
}
```

And then plotting of the whole model will be a combination of plotting the raw data with these predictios on top of it. Let's make a function for that.

```R
plot_heterogametic_model_fixed_r_model <- function(het_model){
       # variable from the model environment
       x <- 1:max(fitting_range)
       y <- HG002[x, 'freq']
       k <- het_model$m$getEnv()$k
       r <- het_model$m$getEnv()$r
      #  cov_xlim <- 1:max(fitting_range)
  	   predicted_by_model <- predict(het_model, response = T, newdata = list(x = x))
       
      coverage_barplot(y, x, ylim = c(0, max(HG002[fitting_range, 'freq']))) 
      lines(predicted_by_model ~ x, lwd = 3)

       estimates <- coef(het_model)
       # extracting all the fitted values
       kmercov <-  estimates['kmercov']
       overdispersal <-  estimates['overdispersal']
       est_length <- estimates['length']
  	   heterogameticSize <- estimates['hetSize']

       disomic_prediction <- predict_disomic_portion_2peak(x, k, kmercov, overdispersal, est_length, heterogameticSize, r)
       monosomic_prediction <- predict_monosomic_portion_2peak(x, k, kmercov, overdispersal, est_length, heterogameticSize, r)

       lines(disomic_prediction ~ x, lwd = 3, col = 'darkgoldenrod1')
       lines(monosomic_prediction ~ x, lwd = 3, col = 'red')

       legend('topright',
              c('kmer histogram','full model', 'autosomes', 'Y + X chromosomes'),
              col = c('deepskyblue','black', 'darkgoldenrod1', 'red'),
              lty = c(NA, 1, 1, 1), lwd = 3,
              pch = c(15, NA, NA, NA), bty = 'n')

       total_genome <- round(est_length / 1e6, 2)
       X_chrom_size_est <- round(heterogameticSize / 1e6, 2)
       het <- round(r * 100, 2)
       title(paste0("Heterozygosity: ", het , '% X + Y: ', X_chrom_size_est, ' Mbp out of ', total_genome, ' Mbp'))
}
```

Now we can plot the model we created by calling this function

```R
plot_heterogametic_model_fixed_r_model(het_model2)
```

### Refining genome size estimate

So, how can we refine the total genome size estimate to include repeats? In regular GenomeScope, if all sequenced had the same ploidy, you could simply sum all the k-mer frequencies and divide it by ploidy and 1n coverage

```R
sum(HG002$cov * as.numeric(HG002$freq)) / (2 * coef(het_model2)['kmercov'])
# 2964003697
```

If you see `NAs produced by integer overflow`, you need to convert intigers to type `numeric` using `as.numeric()`, which is capable of storing much larger numbers. 

Of course, GenomeScope actually substract all the error k-mers from the equation above. Errors in GenomeScope are defined as the left-side residual of the genome model.

```R

freq_wo_errors <- HG002[, 'freq']
errors <- (y[1:20] - predict(het_model2, newdata = list(x = 1:20)))
# errors I will pick all the values of 1:6 which is where most of the errors are
freq_wo_errors[1:6] <- freq_wo_errors[1:6] - errors[1:6]

sum(HG002$cov * as.numeric(freq_wo_errors)) / (2 * coef(het_model2)['kmercov'])
# Our estimate 2,895,567,201 bp
# Compare with genomescope: 2,870,876,412 bp
```

Alright, this is finally a true-GenomeScope-like genome size estimate, but still does not take into account the X and Y. So, in absence of better model, we will simply assume that the repetitiveness of X and Y is the same as of autosommes - yes, that's not exactly right, mostly because of Y, but we will still manage to make a more accurate prediction than ignoring the sex chromosomes completelly. The trick we will do is calculating average ploidy - meaning weighting 1 and 2 by the proportion of disomic and monosomic part of the genome.

```R
genome_size_est <- coef(het_model2)['length'] + coef(het_model2)['hetSize']
disomic_prop <- coef(het_model2)['length'] / genome_size_est
monosomic_prop <- coef(het_model2)['hetSize'] / genome_size_est

average_ploidy <- monosomic_prop + (2 * disomic_prop)

refined_genome_size_est <- sum(HG002$cov * as.numeric(freq_wo_errors)) / (average_ploidy * coef(het_model2)['kmercov'])
# 3,016,012,444
```

This estimate is a pretty damn close to the assembled genome size (`3.11 Gbp`). And of course, using the disomic and monosomic proportion we can refine also the estimates of the autosomes + par and X + Y

```R
refined_monosomic <- refined_genome_size_est * monosomic_prop
refined_disomic <- refined_genome_size_est * disomic_prop

male_disomic - refined_disomic
# Our estimate is 128,648,475 Mbp lower for autosomes + PAR
male_monosomic - refined_monosomic
# Our esitmate is 30,230,890 Mbp higher for X + Y
```

Which is quite interesting, but perhaps not all THAT suprising, given we just assumed heterozygosity of the sample. Perhaps if the real heterozygosity would be a lower, the estimates would fit the known truth better. Further refining the estimate is a problem for another day though...

<details> 
  <summary>(optional only if you have time, go to lichen first) You can fit a grit of heteroxygosities and find out which heterozygosity is closest to the known genome sizes...  </summary>

```R
heterozygosities <- seq(0.0001,0.003, by = 0.0001)

human_wrapper <- function(one_het){
    nlsLM_2peak_heterogametic_model_fixed_r(x[fitting_range], y[fitting_range], estKmercov, estLength, 200e6, one_het, k)
}

het_models <- lapply(heterozygosities, human_wrapper)

refine_estimates <- function(genome_model, x, y){
	genome_size_est <- coef(genome_model)['length'] + coef(genome_model)['hetSize']
	disomic_prop <- coef(genome_model)['length'] / genome_size_est
	monosomic_prop <- coef(genome_model)['hetSize'] / genome_size_est
	genome_model$average_ploidy <- monosomic_prop + (2 * disomic_prop)

	genome_model$refined_genome_size_est <- sum(x * y) / (genome_model$average_ploidy * coef(genome_model)['kmercov'])

	genome_model$refined_monosomic <- genome_model$refined_genome_size_est * monosomic_prop
	genome_model$refined_disomic <- genome_model$refined_genome_size_est * disomic_prop
	return(genome_model)
}

refined_models <- lapply(het_models, function(hetm){ refine_estimates(hetm, x, y) })

refined_disomic <- sapply(refined_models, function(x){x$refined_disomic})
refined_mono <- sapply(refined_models, function(x){x$refined_monosomic})
refined_total <- sapply(refined_models, function(x){x$refined_genome_size_est})

plot(heterozygosities ~ c(refined_mono / 1e6), pch = 20, xlab = 'Heterogametic (monosomic) genome size', ylab = 'Heterozygosity')

refined_mono + refined_disomic

which.min(abs(refined_mono - male_monosomic))
which.min(abs(total_human - refined_total))

heterozygosities[5]

plot_heterogametic_model_fixed_r_model(het_models[[5]])
```

</details>


## The stochiometry principle

Stochiometry as the main principle for diagnosis of the problems. If peaks are not in even ratios, there is somethng funny happening - there always must be a biological or technical reason for it. One good example is sequencing low-complexity metagenomes of only of a few co-bionts, for example lichens. The default genomescope of course does not work on these kind of spectra, so let's try to write our own models that could help us guessing something about our samples...

### Co-bionts dog lichen

The first up is `dog_lichen`. Let's load the spectrum and visualise it.

```R
dog_lichen <- read.table('course_data_2025/histograms/dog_lichen/glPelHori1.k31.hist.txt', col.names = c('cov', 'freq'))

fitting_range <- 40:950
coverage_barplot(dog_lichen[1:max(fitting_range), 'freq'], dog_lichen[1:max(fitting_range), 'cov'], ylim = c(0, max(dog_lichen[fitting_range, 'freq']))) 
```

Right, there are three things in there and none of them are in the even ratios, meaning all of them are probably different genomes. Let's fit a 3 genome model.

<details> 
  <summary>If you want to try to construct it yourself, it is very similar to the original diploid genomescope model, but you need only a single negative bionomial per genome, and no need for the heterozygosity terms. The three coverages and genome lengths need to be different parameters, but the overdispersal can be the same (it's more of a property of the sequencing run, than samples, but don't forget to scale it with coverage). </summary>

```R
x <- dog_lichen$cov[fitting_range]
y <- dog_lichen$freq[fitting_range]
lengthEst <- 1e7
kmerEst <- 220
kmerEst2 <- 300
kmerEst3 <- 840

# (x, y, kmerEst, lengthEst, hetEst = 0.6)
dog_lichen_model <- nlsLM(y ~ (dnbinom(x, mu = kmercov, size = kmercov/overdispersion) * length1) + 
                          (dnbinom(x, mu = kmercov2, size = kmercov2/overdispersion) * length2) +
                          (dnbinom(x, mu = kmercov3, size = kmercov3/overdispersion) * length3),
                        start = list(kmercov = kmerEst, kmercov2 = kmerEst2, kmercov3 =  kmerEst3,
                        overdispersion = 0.5, length1 = lengthEst, length2 = lengthEst, length3 = lengthEst),
                        control = list(minFactor = 1e-12, maxiter = 40))
```

Nice, the model converged, three genomes with 214x, 318x and 839x coverage respectivelly, and genome sizes 53.3Mbp, 41Mbp and 83Mbp respectively.

</details>

<details> 
  <summary>Plotting the model </summary>

```R
model_env <- dog_lichen_model$m$getEnv()

left_peak <- dnbinom(model_env$x, size = model_env$kmercov/model_env$overdispersion, mu = model_env$kmercov) * model_env$length1
middle_peak <- dnbinom(model_env$x, size = model_env$kmercov2/model_env$overdispersion, mu = model_env$kmercov2 ) * model_env$length2
right_peak <- dnbinom(model_env$x, size = model_env$kmercov3/model_env$overdispersion, mu = model_env$kmercov3 ) * model_env$length3

pal <- c("chocolate", "darkorchid4", "pink")

lines(left_peak ~ fitting_range, lwd = 3, col = pal[1])
lines(middle_peak ~ fitting_range, lwd = 3, col = pal[2])
lines(right_peak ~ fitting_range, lwd = 3, col = pal[3])
```

Which looks quite alrigt. It is hard to say at this point which genome is which, if there are two fungi and one photobiont or other way around, but at least we know what to expect in general. One could for example subset reads by k-mer coverages using Gene's profiles and assemble the individual components separately, but that carries a risk of excluding all the repetitive sequences which will fragment the assembly. ALternatively, one could use the knowledge of the expected coverage to process the assembly and separate the individual cobionts.

</details>

### Co-bionts - wolf lichen

```R
wolf_lichen <- read.table('course_data_2025/histograms/Letharia_vulpina_k21.hist', col.names = c('cov', 'freq'))

fitting_range <- 20:400 # I always plot it with a coverage range and redefine it, there would be a more sensible way to do this, but it never bothered me enough to actually implement it
coverage_barplot(wolf_lichen[1:max(fitting_range), 'freq'], wolf_lichen[1:max(fitting_range), 'cov'], ylim = c(0, max(wolf_lichen[fitting_range, 'freq']))) 
```

So, looking at this spectrum, we see two or three genomes. Potentially the first two peaks could be a single genome diploid gneome as their coverage ratio is close to 1:2, while the third peak (~190x) has no other peak it could be associated with so should be modelled as separated genome. 

<details> 
  <summary>So, let's fit the model as such. If you want to try to construct it yourself, it is combination of the original diploid genomescope model and the one we just made for the dog lichen. There will be one diploid and one haploid genome. The two coverages and genome lengths need to be different parameters, but again, the overdispersal can be the same. </summary>

```R
x <- wolf_lichen$cov[fitting_range]
y <- wolf_lichen$freq[fitting_range]
lengthEst <- 1e7
kmerEst <- 60
kmerEst2 <- 190
rEst <- 0.5
# (x, y, kmerEst, lengthEst, hetEst = 0.6)

lichen_model <- nlsLM(y ~ ((2*(1-(1-r)^k))  *  dnbinom(x, size = kmercov/overdispersal, mu = kmercov ) +
                           ((1-r)^k) *         dnbinom(x, size = 2*kmercov/overdispersal , mu = 2*kmercov )) * length1 + (dnbinom(x, size = kmercov2/overdispersal, mu = kmercov2) * length2), 
                        start = list(kmercov = kmerEst, kmercov2 = kmerEst2, 
                        overdispersal = 0.5, length1 = lengthEst, length2 = lengthEst, r = rEst), 
                        control = list(minFactor = 1e-12, maxiter = 40))
 
```

Ha, converged somewhere... The diploid genome seems to be 15.65 Mbp long, the haploid genome is 26.81Mbp. The heterozygosity of the diploid sample is ~2%. This of course counting only the non-repetitive parts of the genome and also assuming the model we decided on. If for example the diploid genome is not really a diploid genome, they are just two different organisms, this model would give very wrong numbers. Usually it is good to create a set of possible realities using these models and then combine those with the assembly to figure out what is going on.

You should also plot the model to visually confirm that the model worked.

</details>



