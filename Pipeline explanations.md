# Sequence processing pipeline
Here is a summary of the pipeline used to process the Ion Torrent sequences. The pipeline is based on the ["alternative vsearch pipeline"](https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline) proposed by Torbjørn Rognes, and based on the approach that usearch uses. In short, after quality filtering, OTUs with 97% similarity are produced from all sequences except the singleton sequences. Then, these OTUs are checked for chimeras. Finally, all original sequences (including singletons) are mapped to these created OTUs using the 97% similarity threshold.

In steps 1-4, the samples are being processed separately, so we use a for loop that iterates over all files for these steps. Afterwards, samples are concatenated and analysed together.
A good habit to take is to check after each step how many reads are being removed, as some surprises may appear (more about this in the decontamination protocol).

The following software/packages are needed to run the pipeline:
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [vsearch](https://github.com/torognes/vsearch)
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
- [LULU](https://github.com/tobiasgf/lulu)
- [CREST4](https://github.com/xapple/crest4)
- You will also need a copy of the SILVA 138 nr database if you want to run the chimera check against it.

### 1. Removal of the forward primers.
When we receive our sequences back, they are already demultiplexed, meaning that each sample has its own fastq file, and the adaptors and barcodes are already removed. However, the forward primer is still present. We want to remove it, as it does not contain any phylogenic information. We use the cutadapt tool. 

```bash
for f in *.fastq; do
    s=$(awk -F'.' '{print $1}' <<< $f)
    cutadapt -j 10 -g CAGCMGCCGCGGTAA --no-indels --discard-untrimmed -o $s.noprimers.fastq $f
done
```

-  `-j`: Run the command over 10 CPU.
- `-g `: The primer sequence to be removed. This primer is the 519F, but for some of the files the sequenced used was the 515F.
- `--no-indels`:  We do not allow for insertion or deletion in the primer sequence, ie. we want to find exactly the good sequence.
- `--discard-untrimmed`: Throw away the whole sequence if the primer is not found in it.
- `-o`: Name of the output.
- `$f`:  This is the input file.

## Trim the sequences at 220bp.  
This can be shortened or lengthened based on what you saw in FastQC. 220 is a good length usually. 

# Trim the sequences at 220bp.
for f in *noprimers.fastq; do
    s=$(awk -F'.' '{print $1 "." $2}' <<< $f)
    $VSEARCH --threads $THREADS \
        --fastq_filter $f \
        --fastq_maxns 0 \
        --fastq_trunclen 220 \
        --fastqout $s.trimmed.fastq
done

# Amount reads in after trimming.
for f in *trimmed.fastq; do
     awk '{s++}END{print s/4}' $f
done

# Quality filtering at maxee = 2 and relabel reads within files.
for f in *trimmed.fastq; do
    s=$(awk -F'.' '{print $1}' <<< $f)
    $VSEARCH --threads $THREADS \
        --fastq_filter $f \
        --relabel $s. \
        --fastq_maxee 2 \
        --fastaout $s.fa \
        --fasta_width 0
done

# Amount reads in fasta file. Little change of strategy here as we are just looking at the amount of ">" in each file.
for f in *.fa; do
    grep -c "^>" $f
done

# Dereplicate at sample level and relabel with sample_n
for f in *.fa; do
    s=$(awk -F'.' '{print $1}' <<< $f)
    $VSEARCH --threads $THREADS \
        --derep_fulllength $f \
        --strand plus \
        --sizeout \
        --relabel $s. \
        --fasta_width 0 \
        --output $s.derep.dfa
done 

# Amount reads in fasta file. Little change of strategy here as we are just looking at the amount of ">" in each file.
for f in *.dfa; do
    grep -c "^>" $f
done



# Concatenate all dfa files
cat *.dfa > AllSeqs.fa

# Remove no more used files, but keep the original files.
rm *noprimers*.fastq *.dfa *.fasta


# Amount reads in fasta file. Little change of strategy here as we are just looking at the amount of ">" in each file.
for f in *.fa; do
    grep -c "^>" $f
done


# Dereplication of the pooled fasta file, removing singletons
$VSEARCH --derep_fulllength AllSeqs.fa \
    --threads $THREADS \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc AllSeqs.derep.uc \
    --output AllSeqs.derep.fasta


grep -c "^>" AllSeqs.derep.fasta

# Clustering at 97% similarity.
$VSEARCH --cluster_size AllSeqs.derep.fasta \
    --threads $THREADS \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --centroids Centroids.fasta


# Sort centroids and remove singletons
$VSEARCH --sortbysize Centroids.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --output Centroids.sorted.fasta

# denovo chimera detection.
$VSEARCH --uchime_denovo Centroids.sorted.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --nonchimeras Centroids.denovo.nonchimeras.fasta \


# Reference chimera detection against SILVA138.
$VSEARCH --uchime_ref Centroids.denovo.nonchimeras.fasta \
    --threads $THREADS \
    --db /home/svenlemoinebauer/Desktop/Work/Projects/202106_AnnaSequencing/SequencingProcessing/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --dbmask none \
    --nonchimeras Centroids.nonchimeras.fasta

# Relabel OTUs.
$VSEARCH --fastx_filter Centroids.nonchimeras.fasta \
    --threads $THREADS \
    --sizein \
    --fasta_width 0 \
    --relabel OTU_ \
    --fastaout Otus.fasta

# Map sequences to OTUs.
$VSEARCH --usearch_global AllSeqs.fa \
    --threads $THREADS \
    --db Otus.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --dbmask none \
    --otutabout Otutab.txt

 # Sort OTU table numerically.
sort -k1.5n Otutab.txt > Otutab.sorted.txt

# To run Lulu, we first need to make a database out of the OTU file.
makeblastdb -in Otus.fasta -dbtype nucl
blastn -db Otus.fasta -outfmt '6 qseqid sseqid pident' -out LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query Otus.fasta -num_threads 12

Run the LULU.R script

# Remove sequences from fasta file that got removed using LULU.
python3 filter_seq_by_OTUtable.py Otus.fasta Otutab_curated.tsv > Otus_curated.fasta

# Sort OTU table numerically.
sort -k1.5n Otutab_curated.tsv > Otutab_curated.sorted.tsv

# Amount reads in fasta file. 
grep -c "^>" Otus_curated.fasta

# Run Crest against the moded SILVA138 database
crest4 -f Otus_curated.fasta -t 12
