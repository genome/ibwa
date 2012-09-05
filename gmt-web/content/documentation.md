# Running iBWA

***

To run iBWA, first use iBWA index to create an index for your base reference, and each reference to remap against. Next, use iBWA aln to align each read file against each reference index created by iBWA index. Finally, use iBWA samse (single-end) or iBWA sampe (paired-end) to create a .sam output. 

# EXAMPLE PIPELINE

***

Please see below for detailed information about the iBWA commands.

## STEP 1
Create aligner index for each reference fasta.

<p class='terminal' markdown='1'>
ibwa index -a bwtsw hs37lite.fa<br/>
ibwa index -a bwtsw hs37patch9.fa<br/>
ibwa index -a bwtsw hs37dbsnp135cutoff5.fa
</p>

If the reference fasta is less than 11,000,000 bytes, use the IS algorithm instead of BWT-SW algorithm:

<p class='terminal' markdown='1'>
ibwa index -a is small_ref.fa
</p>

## STEP 2
Align each input file agaisnt each references.

<p class='terminal' markdown='1'>
ibwa aln hs37lite.fa reads1.fq > hs37lite_reads1.sai<br/>
ibwa aln hs37lite.fa reads2.fq > hs37lite_reads2.sai
<p class='terminal' markdown='1'>
</p>
ibwa aln hs37patch9.fa reads1.fq > hs37patch9_reads1.sai<br/>
ibwa aln hs37patch9.fa reads2.fq > hs37patch9_reads2.sai
<p class='terminal' markdown='1'>
</p>
ibwa aln hs37dbsnp135cutoff5.fa reads1.fq > hs37dbsnp135cutoff5_reads1.sai<br/>
ibwa aln hs37dbsnp135cutoff5.fa reads2.fq > hs37dbsnp135cutoff5_reads2.sai
</p>

## STEP 3
Run samse (for single-ended data) or sampe (for paired-end data) to generate a sam file .

<p class='terminal' markdown='1'>
ibwa sampe -R hs37lite.fa hs37lite_reads1.sai hs37lite_reads2.sai reads1.fq reads2.fq<br/>
    &nbsp;&nbsp;&nbsp; hs37patch9.fa hs37patch9_reads1.sai hs37patch9_reads2.sai<br/>
    &nbsp;&nbsp;&nbsp; hs37dbsnp135cutoff5.fa hs37dbsnp135cutoff5_reads1.sai hs37dbsnp135cutoff5_reads2.sai<br/>
    &nbsp;&nbsp;&nbsp; > aln.sam
</p>

The -R option enables iBWA's remapping mode. Each additional reference requires the reference .fa, followed by the .sai files generated for that reference. More generally:

<p class='terminal' markdown='1'>
ibwa sampe -R baseref.fa baseref_reads1.sai baseref_reads2.sai reads1.fq reads2.fq<br/>
    &nbsp;&nbsp;&nbsp; ref1.fa ref1_reads1.sai ref1_reads2.sai<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ...<br/>
    &nbsp;&nbsp;&nbsp; refN.fa refN_reads1.sai refN_reads2.sai<br/>
    &nbsp;&nbsp;&nbsp; > aln.sam
</p>

## NAME
iBWA - Iterative Burrows-Wheeler Alignment Tool

## SYNOPSIS
    ibwa index -a bwtsw ref.fa

    ibwa aln ref.fa reads.fq > ref_reads.sai

    ibwa sampe -R baseref.fa baseref_reads1.sai baseref_reads2.sai reads1.fq reads2.fq<br/>&nbsp;&nbsp;&nbsp;&nbsp; ref1.fa ref1_reads1.sai ref1_reads2.sai > aln.sam

## NOTES
iBWA is a a fork of [Heng Li's BWA aligner](http://bio-bwa.sourceforge.net), with support for iteratively adding alternate haplotypes, patches, and variant hypotheses. For information about other BWA features, please see the [latest BWA documentation](http://bio-bwa.sourceforge.net/bwa.shtml).

## COMMANDS AND OPTIONS
<dl>
<dt>-R</dt>
<dd>Enables remapping against multiple references.</dd>
</dl>

