## NAME

iBWA - Iterative Burrows-Wheeler Alignment

## SYNOPSIS

<pre class='terminal'>ibwa index -a bwtsw primary_ref.fa
ibwa index -a bwtsw alt_ref.fa

ibwa aln -b1 primary_ref.fa reads.bam > primary_ref_r1.sai
ibwa aln -b2 primary_ref.fa reads.bam > primary_ref_r2.sai
ibwa aln -b1 alt_ref.fa reads.bam > alt_ref_r1.sai
ibwa aln -b2 alt_ref.fa reads.bam > alt_ref_r2.sai

ibwa sampe -R primary_ref.fa primary_ref_r1.sai primary_ref_r2.sai \
    reads.bam reads.bam \
    alt_ref.fa alt_ref_r1.sai alt_ref_r2.sai > output.sam
</pre>

## NOTES

iBWA is a fork of [Heng Li's BWA aligner][bwa] with support for iteratively adding alternate haplotypes, reference patches, and variant hypotheses.

For additional information about the original BWA please see Heng Li's [BWA @ SourceForge][bwa] or [BWA's manual page][bwaman].

## COMMANDS AND OPTIONS

iBWA adds the `-R` option to `sampe` to enable compound sequence remapping. Since iBWA is based on BWA v0.5.9 it also supports those commands and options found in that version of BWA.
See [BWA's manual page][bwaman] or run the command for additional details.

## FILE FORMATS

iBWA accepts reads in FASTQ and BAM format, and it can align to references in FASTA format. Remap files are required when using the `-R` option to enable compound sequence remapping.

For instance, if a user has a reference FASTA named `primary_ref.fa`, and has a second reference named `alternates.fa` that contains alternate sequences for `primary_ref.fa`, a `alternates.fa.remap` file must also exist next to `patches.fa`.

The remap file format is as follows:

<pre class='terminal'>&gt;seqid-chrom|start|stop
cigar
</pre>

* `seqid` is the sequence ID from alternates.fa
* `chrom`, `start`, and `stop` explains the region of the primary reference where the alternate sequence is going to remap (note that this can be much longer than the alternate sequence in the case of sequences with long deletions from the reference)
* `cigar` is a cigar string that maps the alternate sequence back onto the primary reference; the cigar string has a format like `5M10D5M10I5M`, where `M` indicates that reads line up (they could match or mismatch), `I` indicates an insertion (the alternate sequence has extra bases, not the other way around), and `D` indicates a deletion.

The sequence ids in the `.fa` and the `.fa.remap` files must match and must be in the same order, and `seqid` should not contain dashes.

# DOWNLOADING PRE-MADE REFERENCE SEQUENCES

All URLs below are accessible via http and also ftp.

The Genome Institute hosts a copy of the "lite" build 37 (hg19) human reference:
<pre class='terminal'>
wget http://genome.wustl.edu/pub/software/ibwa/hs37lite.fa.gz
gunzip *.gz
</pre>

TGI also hosts the latest patches from The Genome Reference Consortium pre-built to work with iBWA.  This fasta includes all alternate haplotypes as of the current patch level, in addition to all "fix patches" correcting errors in the reference.
<pre class='terminal'>
wget http://genome.wustl.edu/pub/software/ibwa/hs37patch10.fa.gz
wget http://genome.wustl.edu/pub/software/ibwa/hs37patch10.fa.remap.gz
gunzip *.gz
</pre>

TGI also hosts an extension the reference made from dbSNP variants where the variants have a global minor allele frequencey (GMAF) of at least 5%.

<pre class='terminal'>
wget http://genome.wustl.edu/pub/software/ibwa/37dbsnp137cutoff5.fa.gz
wget http://genome.wustl.edu/pub/software/ibwa/37dbsnp137cutoff5.fa.remap.gz
gunzip *.gz
</pre>


# EXAMPLE PIPELINE

---

## STEP 1
Create an aligner index for each reference fasta:

<pre class='terminal'>
ibwa index -a bwtsw hs37lite.fa
ibwa index -a bwtsw hs37patch10.fa
ibwa index -a bwtsw hs37dbsnp135cutoff5.fa

&#35; If small_ref.fa is &lt;= 10 MB use `-a is`:
ibwa index -a is small_ref.fa
</pre>

This will create index files next to the corresponding reference FASTAs.

There are two algorithms available. Specify `-a bwtsw` to use BWT-SW, and use `-a is` to use the IS algorithm. If the reference FASTA is less than 11,000,000 bytes, use the IS algorithm.

From the [BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml#3):
<dl>
<dt markdown='1'>`is`</dt>
<dd>IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.</dd>
<dt markdown='1'>`bwtsw`</dt>
<dd>Algorithm implemented in BWT-SW. This method works with the whole human genome, but it does not work with database smaller than 10MB and it is usually slower than IS.</dd>
</dl>

## STEP 2
Align each input file (FASTQ or BAM) to each reference (FASTA):

<pre class='terminal'>
ibwa aln hs37lite.fa reads1.fq > hs37lite_reads1.sai
ibwa aln hs37lite.fa reads2.fq > hs37lite_reads2.sai

ibwa aln hs37patch10.fa reads1.fq > hs37patch10_reads1.sai
ibwa aln hs37patch10.fa reads2.fq > hs37patch10_reads2.sai

ibwa aln hs37dbsnp135cutoff5.fa reads1.fq > hs37dbsnp135cutoff5_reads1.sai
ibwa aln hs37dbsnp135cutoff5.fa reads2.fq > hs37dbsnp135cutoff5_reads2.sai
</pre>

The `-b` option is required to use BAM files as an input. The `-b1` and `-b2` options can be used to select paired end reads in a BAM file that contains both strands of paired end data. For a BAM file that contains both single and paired end reads, use `-b0` to specify that only single-end reads should be mapped.

For example:

<pre class='terminal'>
ibwa aln hs37lite.fa -b1 allreads.bam > hs37lite_reads1.sai
ibwa aln hs37lite.fa -b2 allreads.bam > hs37lite_reads2.sai
</pre>

## STEP 3
Run `samse` (for single-ended data) or `sampe` (for paired-end data) to generate a SAM file:

<pre class='terminal'>
ibwa sampe -R \
    hs37lite.fa hs37lite_reads1.sai hs37lite_reads2.sai reads1.fq reads2.fq \
    hs37patch9.fa hs37patch10_reads1.sai hs37patch10_reads2.sai \
    hs37dbsnp135cutoff5.fa hs37dbsnp135cutoff5_reads1.sai \
    hs37dbsnp135cutoff5_reads2.sai \
    > aln.sam
</pre>

The `-R` option enables iBWA's remapping mode. Each alternate reference requires an `alt_ref.fa`, followed by the `alt_ref_reads.sai` files generated for that reference.

For example:

<pre class='terminal'>
ibwa sampe -R \
     pri_ref.fa pri_ref_reads1.sai pri_ref_reads2.sai reads1.fq reads2.fq \
     alt_ref1.fa alt_ref1_reads1.sai alt_ref1_reads2.sai \
     alt_ref2.fa alt_ref2_reads1.sai alt_ref2_reads2.sai \
     ...
     alt_refN.fa alt_refN_reads1.sai alt_refN_reads2.sai \
     > aln.sam
</pre>

## INTERPRETATION

All alignments to contigs in references with a .remap file will be converted into primary reference space.  This is a translation of the alignment position and cigar from the alternate, not a re-alignment.  The original refseq name, position and cigar are preserved in a ZR tag.

The resulting SAM/BAM file can be used with standard variant detectors, and does not necessarily require futher special processing.  When reads hit an alternate sequence or a patch, the "variants" which would describe that alternate or patch will appear with a stronger signal in the output.  New variants in the vicinity of the alternate are more likely to be detectable.  False positives commonly occurring around complex variations will be reduced or removed.

There are two things which may be useful to down-stream analysis:

Down-stream analysis and annotation will possibly want to remove variant calls which are related to fix-patches in the reference.  Removing these variants would be of benefit even with regular BWA on the "lite" reference, but those variants would be more difficult to even detect.  The normalization effect created by the alternate sequence will help make these more detectable, preventing things like the shortening of a repeat from appearing in inconsistent ways in the read alignments.

Variants aligning to alternate haplotypes can be interpreted at the whole-haplotype level instead of on a per-variant basis.  The later is detecting variations by counting occurrences of the ZR tag for a given alternate, versus other alternates or versus the primary reference.   A simple execution of samtools view, and piping through grep, awk, sort and uniq will summarize the known haplotypes present in the sample.  Counting occurrences provides an initial perspective on which haplotypes are present.  It is possible that re-aligning against a narrower list of haplotypes (say selecting no more than two in any region of a diploid genome) after examining initial alignments could produce an even clearer picture of the genome.

## ISSUES

The BWA alignment algorithm intentionally scores alignments poorly where certainty of position in the genome is ambiguous.  The iBWA modification mitigates this loss of confidence explicitly where competing alignments re-map to the exact same location.  

In cases where the cigar string associated with a patch lacks nuance, iBWA will not be able to protect the alignment from a scoring discount.  Similarly, the re-location of a piece of the genome may only receive clear scoring improvements where paired-end data is able to recognize that the relocated region is in play versus the original.

[bwa]: http://bio-bwa.sourceforge.net
[bwaman]: http://bio-bwa.sourceforge.net/bwa.shtml
