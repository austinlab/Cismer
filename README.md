Cismer is a discriminative motif scanner and mapper.  It will scan a collection of PSSMs in JASPAR format against a foreground set
of promoter sequences and compare it against an expected background distribution built by randomly sampling the background
set to generate 1000 sets of equal size to the foreground set.  It will then calculate a z-score statistic to determine the 
statistical significance of said motif in the foreground.  Use of cismer without a background set (with the -m switch) will map 
PSSMs against a collection of promoters, returning the hit position, strand and functional depth of the match in BED format.

cismer v2.3

usage:
cismer -p _pssms_ -g _genomebackground.fasta_ -f _cluster.fasta_ [-hvai][-r cutoff][-d cutoff][-s seed][-z cutoff]

	-	- read the pssm output from stdin

	-p	- a file listing formatted PSSMs

	-g	- the genome background in fasta format

	-f	- the cluster file of target genes in fasta format

	-b	- the number of bootstraps to perform (default:1000)

	-z	- the z-score cutoff to use for significance (default: none

	-d	- the functional depth cutoff to use (default: 0)

	-o	- return log-odds PWMs instead of PSSMs

	-m	- just map the pssms against the fasta file (-f)

	-h	- print this usage message


Cite:  Ryan S Austin, Shu Hiu, Jamie Waese, Matthew Ierullo, Asher Pasha, Ting Ting Wang, Jim Fan, Curtis Foong, Robert Breit, 
Darrell Desveaux, Alan Moses, Nicholas J Provart.  **New BAR tools for mining expression data and exploring Cis‐elements in Arabidopsis
 thaliana.** (2016) _The Plant Journal_ 88:3 490-504.
