## NAME
iBWA - Iterative Burrows-Wheeler Alignment Tool

## SYNOPSIS
<pre class='terminal'>ibwa index -a bwtsw ref.fa

ibwa aln -b0 ref0.fa reads.bam > ref0_r1.sai
ibwa aln -b1 ref0.fa reads.bam > ref0_r2.sai
ibwa aln -b0 ref1.fa reads.bam > ref1_r1.sai
ibwa aln -b1 ref1.fa reads.bam > ref1_r2.sai

ibwa sampe -R ref0.fa ref0_r1.sai ref0_r2.sai reads.bam read.bam
    ref1.fa ref1_r1.sai ref1_r2.sai > output.sam</pre>

## NOTES
iBWA is a a fork of [Heng Li's BWA aligner](http://bio-bwa.sourceforge.net), with support for iteratively adding alternate haplotypes, patches, and variant hypotheses. For information about other BWA features, please see the [latest BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml).

## COMMANDS AND OPTIONS
<dl>
<dt markdown='1'>`sampe -R`</dt>
<dd>enable compound sequence remapping</dd>
</dl>

See the latest [BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml) for a description of other commands and options. 

# EXAMPLE PIPELINE

---

## STEP 1
Create aligner index for each reference fasta.

<pre class='terminal'>
ibwa index -a bwtsw hs37lite.fa
ibwa index -a bwtsw hs37patch9.fa
ibwa index -a bwtsw hs37dbsnp135cutoff5.fa
</pre>

If the reference fasta is less than 11,000,000 bytes, use the IS algorithm instead of BWT-SW algorithm:

<pre class='terminal'>
ibwa index -a is small_ref.fa
</pre>

## STEP 2
Align each input file to each references.

<pre class='terminal'>
ibwa aln hs37lite.fa reads1.fq > hs37lite_reads1.sai
ibwa aln hs37lite.fa reads2.fq > hs37lite_reads2.sai

ibwa aln hs37patch9.fa reads1.fq > hs37patch9_reads1.sai
ibwa aln hs37patch9.fa reads2.fq > hs37patch9_reads2.sai

ibwa aln hs37dbsnp135cutoff5.fa reads1.fq > hs37dbsnp135cutoff5_reads1.sai
ibwa aln hs37dbsnp135cutoff5.fa reads2.fq > hs37dbsnp135cutoff5_reads2.sai
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

