# bamcmp
<div style="text-align: left">
bamcmp is a tools for deconvolving host and graft reads. The algorithm makes use of full-length alignments and their scores and can use values
extracted from the CIGAR string or mapping quality scores when
alignment scores are unavailable. With paired end data, conflicts are resolved at
the individual read level, not the fragment level, allowing more
data to be retained. This approach allows a weak but significant
match to one organism to be ignored in favor of a stronger match
to the other.
</div>

## Quick Install

bamcmp requires htslib (http://www.htslib.org).  

Please set the environment variable $HTSLIBDIR to the proper folder
before running make

``` bash
export HTSLIBDIR=/nonstandard/htslib
cd bamcmp/
make
```

## Usage

``` bash
bamcmp -n -1 ABC_human.bam -2 ABC_mouse.bam -a ABC_humanOnly.bam  
-A ABC_humanBetter.bam -b ABC_mouseOnly.bam -B ABC_mouseBetter.bam –C
ABC_humanLoss.bam -D ABC_mouseLoss.bam -s [as/match/mapq/balwayswins].
```


## Citation

Garima Khandelwal, Maria Girotti, Christopher Smowton, Sam Taylor, Chris Wirth, Marek Dynowski, Kris Frese, Ged Brady, Deborah Burt, Richard Marais, Crispin Miller.  <a href="http://mcr.aacrjournals.org/content/15/8/1012.long">Next-Gen Sequencing Analysis and Algorithms for PDX and CDX Models.</a> Molecular Cancer Research. 2017, 15:8, PMID: 28442585 DOI: 10.1158/1541-7786.MCR-16-0431

## License
This project is licensed under <a href="https://opensource.org/licenses/GPL-3.0">GPL-3.0</a>.
