# Part B — Using **Logan Search** to Locate Sequences Across the SRA

> “Query the planet in minutes.”  

Logan Search lets you search up to **1 kb** of DNA/RNA and instantly identifies every public *Sequence Read Archive* (SRA) dataset that is likely to contain it. It does so by finding shared k-mers between the query sequence and each SRA accession. It’s built on a petabyte-scale Bloom filter index covering >27 million runs assembled by **Logan**.

---

## 0 . Prerequisites

| Tool | Why you need it | Install |
|------|-----------------|---------|
| Web browser & email | Submit queries + receive result links | Any modern browser |
| *Optional* `aws` CLI | Fetch assemblies for the hits | See Part A or `aws help` |
| *Optional* `parallel`, `seqkit` | Batch downloads & sequence fiddling | `GitHub`, `brew`, `apt`, etc. |

---


## 1 . Prepare a Query Sequence

* **Limit:** ≤ 1 000 bp (service cap)  
* **Formats:** FASTA
* **Tip:** pick a region *specific* to your organism/gene—avoid conserved rRNA or adapters.

```bash
# Example: extract 800 bp from the E.coli assembly from PartA
seqkit subseq -r 1001:1800 SRR2584403.contigs.fa --quiet | seqkit seq -m 800 -w 0  | head -n 2 > query.fa
```

Trivia: what does the command above really do?

But actually, let us not query that particular sequence, because E. coli is everywhere, there will be too many hits and Logan Search takes more time.

Let us instead query this mysterious sequence:

```
TATGTCCTGTAACCATGAATTAAACTGCACCATATGATTGACTGATGCGGATGAACAGTATTGATCAAGTCCAGCTACAATACATTCTCTGCGACGTAGGCCTATTTACTAGGCCGTAGCTGAGTCACAGCAACTTCCCCCCCCAAAATTAAAAGATACCAAATGCATGACACCACAGTTATGTGTGATCATGGGCTGATTAGTCATTAGTCCATCGAGATGTCCCATTTGCAGGATCAGTAGGATTGAGTATGCGATGGCATGGCCCTGAAGTAAGAACCAGATGCCAGGTATAGTTCCAGTATAGAAACCCCCACATTTTATGGGCCCGGAGCGAGAAGAGGTACACTTCGCGCAAGGATTGATGGTTTCTCGAGGCTTGGTGATTAAGCTCGTGATCTAGGTGTGTTGCACCTATGTAAGATGAACCATAATATGTACATGCTTATATGCATGGGGCAAACCACCAATGCACGACGTACATAGGGTAGTTTAAAGGTTATATATGGGGTGGTGCAGTGGGTATTTGATATGTAGGTGTTGGGGTGGG
```


---

## 2 . Open the Interface

1. Navigate to **<https://logan-search.org>**.  
2. Click **“Query the planet”** (upper‑right).  
3. Under **Input**, choose **Text** (paste) *or* **File** (upload) and provide your sequence.  
4. *(Optional)* adjust the sidebar settings (see table below).  
5. Hit **Submit**. Wait about 5 minutes.

Pro tip: search `All no viral` instead of `all` groups for faster search. (Because search scales with number of SRA runs, and there are a lot of Covid19 runs.)

### Sidebar settings

| Setting | What it does |
|---------|--------------|
| **Groups** | Search a subset of SRA (*Genomic*, *Single‑cell*, etc). By default, use *all* |
| **Minimum coverage** | Fraction of 31‑mers that must match (default = 0.5). Raising it reduces sensitivity but increases specificity. |
| **Email** | Necessary, and results links are sent there. |

> **Service limits**  
> • One sequence per submission  
> • Top 20 000 hits returned  
> • Results retained for **1 month**—download them!

---

## 3 . Look At The Results

Once Logan Search has finished searching (takes about 5 minutes), browse through the results webpage, and take a look at the map and plots!

The first tab, "Table", show the list of SRA accessions where the query sequence is present. 

![image](https://github.com/user-attachments/assets/e4f6d588-1e7c-4580-ad4a-3581c1839fc6)

Key columns are:

* `ID`: SRA accession 
* `kmer_coverage`: how many k-mers of the query are also present in the accession
*  all the other columns are from the SRA metadata

The next tab, "Map", shows the location of SRA samples according to ther geolocalization metadata provided by the data submitters.

![image](https://github.com/user-attachments/assets/08e56e6f-8b4c-4fdc-9968-8a719d045366)


The color of the circles correspond to `kmer_coverage`, and the size of the circle is proportional to the number of SRA samples at that location.

The next tab, "Plots", shows various plots for the metadata!

![image](https://github.com/user-attachments/assets/056088e6-4705-4b09-a826-6b8124fe7b90)


Be sure to check the bottom-right button:

![image](https://github.com/user-attachments/assets/16a3b4d5-0851-4050-a41c-a1da77933434)

To see many other plot templates, i.e. for sequencing centers, tissues.

**Q1.** What do you think is the organism?

**Q2.** Where is it mostly found?

**Q3.** Is the sequence specific to that organism?



## 4. Get a BLAST-like alignment

Logan Search main results table does not return alignments, instead, it just returns "coverage" of each accession for your query sequence. Think of it as "how likely is the sequence present in that accession", or, "how similar is the query sequence to the sequence present in the accession", although both statements are not strictly true.

You may verify that the sequence is indeed present by manually asking Logan Search to run an alignment of the query sequence to the unitigs or contigs of a particular accession. To do this, click on the "BLAST-like alignment" tab:

![image](https://github.com/user-attachments/assets/64ff2f9d-d5b8-4582-b91e-ac7260374df0)

Then select accession `SRR12518690`, select contigs, and press search.

**Q4.** How many search results?

**Q5.** Is the query sequence fully contained in the contig?

**Q6.** Is the query sequence exactly present in the contigs, or are there mismatches?

In some cases we're not so lucky to have the whole query sequence present in a single contig. Try for example to align to `SRR26176344`, *unitigs* this time.

## 5 . Retrieve Your Results

The confirmation e‑mail contains **two links**:

| Link URL  | Purpose |
|-------------|---------|
| `/dashboard` | Raw output: `.tsv`, `session.json`, your query. |
| `/api/download` | Interactive *kmviz* dashboard (taxonomy barplot, hit table, coverage slider). |

```bash
# Unzip & extract the accession column
unzip kmviz-<session>.zip
cut -f1 <queryname>.tsv | tail -n +2 > hits.acc   # skip header
```

---

## 6 . Download Assemblies for the Hits

```bash
# Fetch the top 25 hits using 4 parallel threads
head -n 25 hits.acc | parallel -j4 '
  aws s3 cp s3://logan-pub/c/{}/{}.contigs.fa.zst . --no-sign-request &&
  zstd -d {}.contigs.fa.zst
'
```

You now have contigs ready for alignment, variant calling, or pangenome analysis.

---


## 7 . Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| “Too many concurrent queries” | High server load | Wait a few minutes, then resubmit. |
| Empty result set | Sequence too short or repetitive | Lengthen the query or choose a unique region. |
| No e‑mail received | Spam filter or typo | Check spam; resubmit with correct address. |
| “Top 20 k cutoff” notice | Query matched *too many* datasets | Raise **Minimum coverage** or pick a more specific sequence. |

---

## 8 . Next Steps

* **Variant scanning:** align all retrieved contigs to a reference.  
* **Phylogenetics:** build a tree from sequences extracted from the matched runs. 
* **Epidemiology:** combine `.tsv` Logan Search results with sampling dates/locations.  
* **Iterative discovery:** pull novel haplotypes → search again → expand the network.

---

### 🏁 Quick Recap

1. Craft ≤ 1 kb sequence.  
2. Submit on Logan Search, tweak filters, add your e‑mail.  
3. Explore *kmviz*, download the ZIP.  
4. Script‑fetch assemblies, and dive into analysis.

You’ve just queried **every** public sequencing run on Earth in minutes—welcome to planetary‑scale genomics!



Now, let's move on to [Part C](PartC.md).
