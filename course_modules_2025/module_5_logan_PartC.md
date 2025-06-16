# Part C — Assembly‑Graph Exploration with **Bandage**

In Parts A and B you learned how to grab Logan contigs and pinpoint SRA runs that contain your favourite 1 kb sequence.  
Now you’ll **see** where that sequence sits in the assembly graph, follow alternative paths, and pull out neighbouring genes — all without writing a custom parser.

We’ll combine three tools:

| Tool | Purpose | Install |
|------|---------|---------|
| [`minimap2`](https://github.com/lh3/minimap2) | Fast alignment of the query against contigs | `brew install minimap2` / `conda install -c bioconda minimap2` |
| [`BandageNG`](https://github.com/asl/BandageNG) | Interactive GFA viewer with built‑in BLAST | Download binaries for Win / macOS / Linux |

We’ll work with **SRR2584403** (E. coli LTEE) because it’s tiny and quick. Adapt the run ID for your own data.

---

## 0 Download contigs & graph

```bash
# 0.1 Fetch the contigs (≈12 MB compressed)
aws s3 cp s3://logan-pub/c/SRR2584403/SRR2584403.contigs.fa.zst . --no-sign-request
zstd -d SRR2584403.contigs.fa.zst

# --- reconstruct the GFA (follows https://github.com/IndexThePlanet/Logan/blob/main/Unitigs.md#assembly-graph)
wget https://raw.githubusercontent.com/GATB/bcalm/refs/heads/master/scripts/convertToGFA.py
sed -i 's/>.*_/>/' SRR2584403.contigs.fa
# for Mac, apparently the command is:
# sed -i '' 's/>.*_/>/' SRR2584403.contigs.fa
python convertToGFA.py SRR2584403.contigs.fa SRR2584403.contigs.gfa 31
```

> **Why the extra step?** Logan’s public bucket stores contigs in FASTA format, but a the assembly graph is easy to reconstruct from it.

Inspect the GFA file:

    head SRR2584403.contigs.gfa -c 100

should return:

    H       VN:Z:1.0        ks:i:31
    S       0       GATAAGACGCGCCAGCGTCGCATCAGGCAAGACCGTATTAATTCGGCGTATAGCCCATTACCACGCCACTTAAGCCA

which is the beginning of a valid GFA file. To learn more about the GFA format, look here: https://gfa-spec.github.io/GFA-spec/GFA1.html

---

## 1 Open the graph in **Bandage**

1. Launch Bandage → **File › Load Graph** → `SRR2584403.gfa`.  
2. The spaghetti resolves into an E. coli‑sized blob (~5 Mb).  


## 2 Investigate parts of the graph

You should see a weird "balloon" that stick out of the hairball:

![image](https://github.com/user-attachments/assets/83a1d21c-d978-454f-88c2-cb9237214cda)

What is it? Let us investigate.

We will select one node that is long enough:

![image](https://github.com/user-attachments/assets/14e597a5-ec2d-444d-ad08-550297ff9828)

and search it on BLAST:

![image](https://github.com/user-attachments/assets/538ebbc1-161d-40c3-a3b1-c558e249e14f)

or alternatively, you could have also copied the sequence to the clipboard:

![image](https://github.com/user-attachments/assets/16a9e97b-e11c-425c-add1-41576c6dc23c)

But then, you'll notice that actually this sequence is not worth being searched:

![image](https://github.com/user-attachments/assets/e5c8ed99-cad4-4ea7-ab5d-5cb3cc755465)

because it is low-complexity. What's going on? To investigate further, select _all_ the sequences in that "balloon":

![image](https://github.com/user-attachments/assets/858cefed-e34f-408e-9044-67b5cb0103d3)

and copy those sequences to the clipboard. It should contain:

    AAAAAAAAAAAAAAAAAAAAAAAAAAAAACA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAGA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAATA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAATT
    ...

They are all low-complexity, so, not worth investigating further, but at least now we know. Why were they sequenced? This is what Chatgpt observes when given that list of sequences:

| Observation                                         | Why it screams “artefact”                                                                                                 |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| **>90 % A/T, frequently pure A or T for 30–40 bp**  | *E. coli* rarely has homopolymers longer than \~9 nt—its AT-rich intergenic terminators top out around 8–10 T’s.          |
| **Step-wise ladders of ‘A…C’, ‘A…G’, ‘A…T’**        | Exactly what you get when a de Bruijn graph assembler walks a homopolymer and the k-mer window shifts one base at a time. |
| **Motifs like `AATGAT`, `ATCATT`, `CGTATCATT`**     | Sub-strings of Illumina P5/P7 adapters and Nextera/Tn5 mosaic ends (“AATGATAC…”, “CTGTCTCTTATACAC…”).                     |
| **Blocks of `TTTT…AATGAT(A/C/T/G)`**                | Classic read-through into adapter after the insert ends—library QC/trim step missed them.                                 |
| **Appear as dozens of tiny unitigs, each 31–80 bp** | The k-mer (k=31) graph can’t merge pure homopolymers with adapter seed—gives many short dead-end nodes.                   |


## 2 Map the 16S with **minimap2**

Let's search for the 16S gene in that assembly. First we'll use minimap2. 

To find the sequence of 16S, let us enter "rRNA gene ecoli ncbi" on Google. First hit is: https://www.ncbi.nlm.nih.gov/nuccore/J01859

![image](https://github.com/user-attachments/assets/9b1b47e0-29b3-4049-b4e8-5f3622e1f33f)

Go there, export the sequence using the "Send to" button:

![image](https://github.com/user-attachments/assets/0cb2bb5f-b20b-441a-885c-bd7cf99e4510)

then type:

```bash
minimap2 -x asm5 -t4 SRR2584403.contigs.fa sequence.fa > query_vs_SRR2584403.paf
column -t query_vs_SRR2584403.paf | head
```

* **`asm5` preset** tolerates 5 % divergence; good for Illumina assemblies.
* You may instead use the `-a` file to ask for a SAM alignment instead of default PAF alignment format (https://github.com/lh3/miniasm/blob/master/PAF.md). 
* * The PAF tells you which contig(s) contain the sequence and at what coordinates.

You should see:

    J01859.1  1541  1  1539  +  360  1732  147  1686  1506  1539  60  tp:A:P  cm:i:148  s1:i:1506  s2:i:0  dv:f:0.0004  rl:i:0

It is is an alignment in PAF format (https://github.com/lh3/miniasm/blob/master/PAF.md).  Try again with the `-a` argument to see:

    J01859.1        0       360     147     60      88M1D1453M      *       0       0       AAATTGAAGAGTTT[...]  *       NM:i:3  ms:i:1459       AS:i:1459       nn:i:0  tp:A:P  cm:i:148        s1:i:1506   s2:i:0  de:f:0.0019     rl:i:0

Which I find easier to read, thank to the CIGAR string and the NM flags, indicating that the sequence is essentially found with little mutations.

---
### 2.1 Run *Bandage BLAST*

You will need to have BLAST installed on your system for this feature to be activated in Bandage.

1. In BandageNG, go to: Create/view graph search -> Build Blast DB -> Load from fasta file -> Select sequence.fasta
4. Click **Run BLAST search**
5. Color the nodes by double clicking on Annotations - Blast hits:

![image](https://github.com/user-attachments/assets/68402ecb-4db6-4b7b-bcf8-3fd33dfd87b7)
 
You should see a small portion of the graph colored (in rainbow color if you selected it so).

5. Extract just the portion of the graph containing the alignment:

![image](https://github.com/user-attachments/assets/dd53d94d-cf31-4486-820e-dfb41cda28f9)

Then what do you notice? The 16S sequence is flanked by variations and repeats! This is because the assembly is not fully repeat-resolved, due to the short reads sequencing.

![image](https://github.com/user-attachments/assets/39336156-f88c-4334-bfa7-dbb6dc1d2531)

But at least, it is present, and you can visualize the genomic neighborhood. For instance, you may ask, what's following it? Select the node right after (or before) the bubble and Web blast it:

![image](https://github.com/user-attachments/assets/aa3b385d-3c86-4979-a190-9d22732c9f21)

In the BLAST results, you can then click on "Graphics" in one of the hits:

![image](https://github.com/user-attachments/assets/6e4aadb3-f9b4-41e4-bb3e-ec216404a4bd)

And then you should see a genome browser:

![image](https://github.com/user-attachments/assets/4b2b37d8-f2a5-42d3-b003-c207be36b174)

With the annotation track displayed! (For some of the BLAST hits, it won't display the annotation track.)

So, this indicates that our query is within the 23S gene. This is apparently expected, as 23s follows 16s in E. coli:

![image](https://github.com/user-attachments/assets/9abfd6b4-d92e-4401-82ff-de6a876986bb)

(Image from https://www.biorxiv.org/content/10.1101/2020.08.03.235192v1.full)

That's it! To recap:

* We've searched for a known sequence in the assembly graph (16S)
* Looked at the graph neighborhood
* Reverse searched for an unknown sequence in the graph near the known one, and found it was a 23S sequence

This allows to do genomic exploration in longer ranges than what is contained in a single Logan contig.

---


