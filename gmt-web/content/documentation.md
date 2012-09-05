## NAME
iBWA - Iterative Burrows-Wheeler Alignment Tool

## SYNOPSIS
<pre class='terminal'>ibwa index -a bwtsw ref.fa

ibwa aln -b1 ref0.fa reads.bam > ref0_r1.sai
ibwa aln -b2 ref0.fa reads.bam > ref0_r2.sai
ibwa aln -b1 ref1.fa reads.bam > ref1_r1.sai
ibwa aln -b2 ref1.fa reads.bam > ref1_r2.sai

ibwa sampe -R ref0.fa ref0_r1.sai ref0_r2.sai reads.bam reads.bam
    ref1.fa ref1_r1.sai ref1_r2.sai > output.sam</pre>

## NOTES
iBWA is a a fork of [Heng Li's BWA aligner](http://bio-bwa.sourceforge.net), with support for iteratively adding alternate haplotypes, patches, and variant hypotheses. For information about other BWA features, please see the [latest BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml).

## COMMANDS AND OPTIONS
<dl>
<dt markdown='1'>`sampe -R`</dt>
<dd>enable compound sequence remapping</dd>
</dl>

See the latest [BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml) for a description of other commands and options. 

## FILE FORMATS

iBWA accepts reads in FASTQ and BAM format, and it can align to references in FASTA format. Remap files are required when using the `-R` option to enable compound sequence remapping.

For instance, if a user has a reference FASTA named `base_ref.fa`, and has a second reference named `alternates.fa` that contains alternate sequences for `base_ref.fa`, a `alternates.fa.remap` file must also exist next to `patches.fa`.

The remap file format is as follows:

<pre class='terminal'>>seqid-chrom|start|stop
cigar</pre>

* `seqid` is the sequence ID from alternates.fa
* `chrom`, `start`, and `stop` explains the region of the primary reference where the alternate sequence is going to remap (note that this can be much longer than the alternate sequence in the case of sequences with long deletions from the reference)
* `cigar` is a cigar string that maps the alternate sequence back onto the primary reference; the cigar string has a format like `5M10D5M10I5M`, where `M` indicates that reads line up (they could match or mismatch), `I` indicates an insertion (the alternate sequence has extra bases, not the other way around), and `D` indicates a deletion.

The sequence ids in the `.fa` and the `.fa.remap` files must match and must be in the same order, and `seqid` should not contain dashes.

# EXAMPLE PIPELINE

---

## STEP 1
Create aligner index for each reference fasta.

<pre class='terminal'>
ibwa index -a bwtsw hs37lite.fa
ibwa index -a bwtsw hs37patch9.fa
ibwa index -a bwtsw hs37dbsnp135cutoff5.fa
</pre>

This will create several index files next to the reference FASTA.

There are two algorithms available. Specify `-a bwtsw` to use BWT-SW, and use `-a is` to use the IS algorithm. If the reference fasta is less than 11,000,000 bytes, use the IS algorithm.

From the [BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml#3):
    > *is* IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.
    > *bwtsw* Algorithm implemented in BWT-SW. This method works with the whole human genome, but it does not work with database smaller than 10MB and it is usually slower than IS. 

<pre class='terminal'>
ibwa index -a is small_ref.fa
</pre>

## STEP 2
Align each input file to each references. FASTQ and BAM files are valid inputs.

<pre class='terminal'>
ibwa aln hs37lite.fa reads1.fq > hs37lite_reads1.sai
ibwa aln hs37lite.fa reads2.fq > hs37lite_reads2.sai

ibwa aln hs37patch9.fa reads1.fq > hs37patch9_reads1.sai
ibwa aln hs37patch9.fa reads2.fq > hs37patch9_reads2.sai

ibwa aln hs37dbsnp135cutoff5.fa reads1.fq > hs37dbsnp135cutoff5_reads1.sai
ibwa aln hs37dbsnp135cutoff5.fa reads2.fq > hs37dbsnp135cutoff5_reads2.sai
</pre>

The `-b` option is required to use BAM files as an input. The `-b1` and `-b2` options can be used to select paired end reads in a BAM file that contains both strands of paired end data. For a BAM file that contains both single and paired end reads, use `-b0` to specify that only single-end reads should be mapped. For instance:

<pre class='terminal>
ibwa aln hs37lite.fa -b1 allreads.bam > hs37lite_reads1.sai
ibwa aln hs37lite.fa -b2 allreads.bam > hs37lite_reads2.sai
</pre>

## STEP 3
Run `samse` (for single-ended data) or `sampe` (for paired-end data) to generate a SAM file.

<pre class='terminal'>
ibwa sampe -R
    hs37lite.fa hs37lite_reads1.sai hs37lite_reads2.sai reads1.fq reads2.fq
    hs37patch9.fa hs37patch9_reads1.sai hs37patch9_reads2.sai
    hs37dbsnp135cutoff5.fa hs37dbsnp135cutoff5_reads1.sai
    hs37dbsnp135cutoff5_reads2.sai
    > aln.sam
</pre>

The `-R` option enables iBWA's remapping mode. Each additional reference requires the `reference.fa`, followed by the `reads.sai` files generated for that reference. More generally:

<pre class='terminal'>
ibwa sampe -R
     baseref.fa baseref_reads1.sai baseref_reads2.sai reads1.fq reads2.fq
     ref1.fa ref1_reads1.sai ref1_reads2.sai
     ref2.fa ref2_reads1.sai ref2_reads2.sai
     ...
     refN.fa refN_reads1.sai refN_reads2.sai
     > aln.sam
</pre>

