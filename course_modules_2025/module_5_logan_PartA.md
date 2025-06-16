## Part A: Accessing Logan data

### 0. Preparation

You need the AWS command-line interface (CLI) installed either on your computer or in a virtual environment (if one is provided by the tutorial).

To check if it is installed, run the command in a terminal.

```
   aws help
```

If it is not installed, you may download it from https://aws.amazon.com/cli/

Or, one handy command line provided by tester Sam Ebdon, to install AWS CLI and everything else in the practical, is:

    mamba install awscli seqkit minimap2 samtools blast bioconda::bandage_ng bioconda::sra-tools

Other prerequisites are Linux tools `shuf` and `parallel`.

---

### 1. Exploring an SRA accession

The Short Read Archive is the central repository of all public sequencing data. Access it at this location: https://www.ncbi.nlm.nih.gov/sra

![image](https://github.com/user-attachments/assets/a0a73965-5764-4dff-acba-71645dde48fc)

We will focus on a E. coli long-term evolution run, identifier *SRR2584403*, because it is small enough to test.

![image](https://github.com/user-attachments/assets/a4e931d7-35f5-420c-b2b4-1ef138ca9900)

You may see some more information about the run (SRR2584403) by clicking on it. Head to the `Data access` tab to see how to download that run.

![image](https://github.com/user-attachments/assets/d9e6c218-cb9b-4cf3-88e7-dc39a7f83cff)

Now, how to download a FASTQ file using the command line? It is not immediate. One way to get a FASTQ file is to first get a SRA file and then convert that SRA file to FASTQ.

There is an easier way, but for now we will go this route, because it will teach us important things about how to get data from the cloud.

To get the SRA file, look at the `SRA archive data` section and download any of the files.

![image](https://github.com/user-attachments/assets/0567c726-37f1-41bb-ade7-64c59df6eddc)

You may download using the first link by typing:

    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR2584403/SRR2584403
    mv SRR2584403 SRR2584403.sra

But, do not do this, see instead if you can download any of the other links.

Then once you have done the download, what you have is a file in the SRA (or SRAlite) format that needs to be converted to FASTQ. To do that, install the SRA toolkit https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
and type:

    fasterq-dump SRR2584403.sra

If all went well, you should see files `SRR2584403_1.fastq` and `SRR2584403_2.fastq` in your folder (these are paired-end reads), around 700 MB each.

You may attempt to download the data in the `Original format`, at location:

    s3://sra-pub-src-8/SRR2584403/Ara-1_500gen_762B_R1.fastq.gz

by typing:

    aws s3 cp s3://sra-pub-src-8/SRR2584403/Ara-1_500gen_762B_R1.fastq.gz . --no-sign-request

Although this is a S3 URL, try to download it, you'll likely get an error 403: Forbidden. Turns out the way to access this data is through `Cloud Data Delivery` which requires NCBI to retrieve it from cold storage, so it is not immediate.

Now that you have `fastq` files, you may do any processing you like on those reads, e.g. assemble them or align them. But we won't do that in this tutorial, we will instead show you how to access
this same accession using Logan instead.

You may wonder if there was an easier way to get FASTQ files from SRA directly, as alluded above. There is. The `fasterq-dump` program also takes a run identifier as argument, so that if you type:

    fasterq-dump SRR2584403

It will find the location of the SRA file automatically, download it, convert it to FASTQ, and delete it, so that you end up with only the two FASTQ files. We did all the steps above to see what happened behind the scenes.

Also, if you are bothered with `fasterq-dump` not showing any output until it has finished, you may add the `--progress` command line argument.


### 1. Accessing SRA data through Logan

[Logan](https://github.com/IndexThePlanet/Logan/) has assembled all SRA data until late 2023, so this run `SRR2584403` is part of Logan. You may download the Logan assembly using this command:

    aws s3 cp s3://logan-pub/c/SRR2584403/SRR2584403.contigs.fa.zst . --no-sign-request

Now you have a FASTA file of the assembly of that run, compressed in Zstandard format.

See if you have the `zstd` and `zstdcat` commands installed on your system by just typing:

    zstdcat

If this doesn't work, install Zstd here: https://github.com/facebook/zstd/releases

Then extract the assembly by typing:

    zstd -d SRR2584403.contigs.fa.zst

Now you have download your (first?) Logan assembly!

You may quickly check some basic statistics about the assembly using Seqkit: https://github.com/shenwei356/seqkit/releases

    seqkit stats SRR2584403.contigs.fa -a

Should display:

    file  format  type  num_seqs    sum_len  min_len  avg_len  max_len  Q1  Q2     Q3  sum_gap     N50  Q20(%)  Q30(%)  GC(%)
    -     FASTA   DNA      1,297  4,556,063       31  3,512.8   81,999  32  57  1,015        0  22,171       0       0  50.72


### 2. Quick note about bacterial data

Above we downloaded a small bacterial assembly. But, when it comes to bacterial isolates, it is actually not recommended to use Logan. A better resource is AllTheBacteria: https://allthebacteria.readthedocs.io/en/latest/

Because it has bacterial assemblies that are of higher contiguity. Though, note that AllTheBacteria uses SRA sample identifiers (SAMxxxxxx), not run identifiers (SRRxxxxxx). Not to worry, you can find the sample identifier in the SRA page: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2584403&display=metadata

![image](https://github.com/user-attachments/assets/8fbecc76-2b98-4e63-8cb1-68524cf9edbb)

Then you can download the AllTheBacteria assembly like this:

    aws s3 cp --no-sign-request s3://allthebacteria-assemblies/SAMN04096045.fa.gz .
    gunzip SAMN04096045.fa.gz

And notice that it has improved assembly stats:

    seqkit stats -a  SAMN04096045.fa

gives:

    file             format  type  num_seqs    sum_len  min_len   avg_len  max_len     Q1      Q2      Q3  sum_gap      N50  Q20(%)  Q30(%)  GC(%)
    SAMN04096045.fa  FASTA   DNA        118  4,544,381      202  38,511.7  348,273  1,034  14,303  46,440        0  110,753       0       0  50.72


Notice how the N50 is greater, though the total lengths are about the same.

### 2. Finding many SRA datasets

Now, downloading a single SRA accession is neat, but what about many many runs? Let us download 100 E.coli assemblies. But how to find the corresponding runs?

One way is to browse the SRA website, search for "ecoli", and painstakingly copy the run names of the 100 first search results. That would take too much time.

A more clever approach is to export a list of search results:

![image](https://github.com/user-attachments/assets/b62787f6-1818-4a7b-b751-fca523b70181)

And then parse that CSV to get just the accession names:

    head SraAccList.csv # to inspect
    grep -v "sra" SraAccList.csv | head -n 100 > ecoli.acc.txt

Now you have in `ecoli.acc.txt` a list of 100 ecoli accessions. Admittedly, they aren't in a random order, so it may not reflect a broad diversity of E. coli genomes.
This can be mitigated by typing instead:

    grep -v "sra" SraAccList.csv | shuf | head -n 100 > ecoli.acc.txt

A more advanced version of this section is available as an independent Logan tutorial: https://github.com/IndexThePlanet/Logan/blob/main/SRA_list.md

### 3. Downloading many SRA datasets

Now you have a truly random sample of E.coli accessions from the SRA. You could download them one by one:

    first_acc=$(head -n 1 ecoli.acc.txt)
    echo $first_acc
    fasterq-dump $first_acc

and repeat for all the other accessions in the file, e.g. with (but do not run that command):

    for acc in $(cat ecoli.acc.txt); do echo $acc; fasterq-dump $acc; done

But doing this is outside the scope of this tutorial, so we won't do that. 

Instead, to download all the Logan assemblies for this list of accessions, type:

    for acc in $(cat ecoli.acc.txt); do echo $acc; aws s3 cp s3://logan-pub/c/$acc/$acc.contigs.fa.zst . --no-sign-request; done

You do not need to wait for this command to complete, because, how long do you think it will take? Around a few seconds per accession, and we have 100 accessions, so it will take a dozen minutes or so. I suggest that you interrupt it with Control+C. 

Also, you will notice that some accessions cannot be found. This is because Logan has currently assembled only data prior to 2024, so anything uploaded later is not there.

The point was to show how to download many accessions. 

You may instead of doing a `for` loop, use the `parallel` program (https://ftp.gnu.org/gnu/parallel/). It will download in parallel instead of one after the other. You may add the command line argument `-j 4` to limit it to 4 threads at the same time (by default it uses the number of CPU cores).

Try it, and stop it too with Control+C.


### 4. Recap

What have we seen?


| 📚 Step                         | 📝 What you did                                                             |
|  ------------------------------- | --------------------------------------------------------------------------- |
| 🛠️ **AWS CLI check**            | `aws help` to confirm CLI, install if missing                               | 
| 🔍 **SRA deep-dive**            | Visited SRA page → located run **SRR2584403** → peeked at “Data access” tab |
| ⬇️ **Grab raw SRA file**        | `wget`/`aws s3 cp` (learned about cold storage 403)                         |
| 🔄 **Convert to FASTQ**         | `fasterq-dump SRR2584403.sra` (or run without pre-download)                 | 
| 🛰️ **Logan data**              | `aws s3 cp s3://logan-pub/...contigs.fa.zst` Logan ships assembled data      | 
| 🎈 **Decompress + QC**          | `zstd -d` ➜ `seqkit stats -a`                                               | 
| 🦠 **Better for bacteria?**     | Checked **AllTheBacteria**, higher contiguity                               |
| 🔎 **Harvest 100 E. coli runs** | Export CSV from SRA → `grep` / `shuf` → `ecoli.acc.txt`                      | 
| 📥 **Mass download**            | `for acc … aws s3 cp …` or `parallel -j 4`                                  |


Now, let's move on to [Part B](PartB.md).
