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
