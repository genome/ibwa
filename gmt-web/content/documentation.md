# Running iBWA

***

To run ibwa, first use ibwa index to create an index for your base reference, and each reference to remap against. Next, use ibwa aln to align each read file against each reference index created by ibwa index. Finally, use ibwa samse (single-end) or ibwa sampe (paired-end) to create a .sam output. 

# EXAMPLE PIPELINE

***

Please see below for detailed information about the ibwa commands.

## STEP 1
Create aligner index for each reference fasta.

<p class='terminal' markdown='1'>
ibwa index -a bwtsw hs37lite.fa
ibwa index -a bwtsw hs37patch9.fa
ibwa index -a bwtsw hs37dbsnp135cutoff5.fa
</p>

If the reference fasta is less than 11,000,000 bytes, use the IS algorithm instead of BWT-SW algorithm:

<p class='terminal' markdown='1'>
ibwa index -a is small_ref.fa
</p>

## STEP 2
Align each input file agaisnt each references.

<p class='terminal' markdown='1'>
ibwa aln hs37lite.fa reads1.fq > hs37lite_reads1.sai
ibwa aln hs37lite.fa reads2.fq > hs37lite_reads2.sai
ibwa aln hs37patch9.fa reads1.fq > hs37patch9_reads1.sai
ibwa aln hs37patch9.fa reads2.fq > hs37patch9_reads2.sai
ibwa aln hs37dbsnp135cutoff5.fa reads1.fq > hs37dbsnp135cutoff5_reads1.sai
ibwa aln hs37dbsnp135cutoff5.fa reads2.fq > hs37dbsnp135cutoff5_reads2.sai
</p>

## STEP 3
Run samse (for single-ended data) or sampe (for paired-end data) to generate a sam file .

<p class='terminal' markdown='1'>
ibwa sampe -R hs37lite.fa hs37lite_reads1.sai hs37lite_reads2.sai reads1.fq reads2.fq hs37patch9.fa hs37patch9_reads1.sai hs37patch9_reads2.sai hs37dbsnp135cutoff5.fa hs37dbsnp135cutoff5_reads1.sai hs37dbsnp135cutoff5_reads2.sai > aln.sam
</p>

The -R option enables ibwa's remapping mode. Each additional reference requires the reference .fa, followed by the .sai files generated for that reference:

<p class='terminal' markdown='1'>
ibwa sampe -R ref0.fa ref0_reads1.sai ref0_reads2.sai reads1.fq reads2.fq ref1.fa ref1_reads1.sai ref1_reads2.sai ... refN.fa refN_reads1.sai refN_reads2.sai > aln.sam
</p>

## NAME
ibwa - Iterative Burrows-Wheeler Alignment Tool

## SYNOPSIS
    ibwa index -a bwtsw ref.fa

    ibwa aln ref.fa reads.fq > ref_reads.sai

    ibwa samse ref.fa ref_reads.sai > aln.sam

    ibwa sampe -R ref0.fa ref0_reads1.sai ref0_reads2.sai reads1.fq reads2.fq ref1.fa ref1_reads1.sai ref1_reads2.sai > aln.sam

## NOTES
A fork of [Heng Li's BWA aligner](http://bio-bwa.sourceforge.net), with support for iteratively adding alternate haplotypes, patches, and variant hypotheses.

## OPTIONS
<dl>
<dt>-R</dt>
<dd>Enables remapping against multiple references.</dd>
</dl>

