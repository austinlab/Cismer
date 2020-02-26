/* fastaSeq.c - three models for the randomization of protein strings and associated helper functions and
 * a function for pulling sequence information from a fasta file stream */ 

#include "seq.h"

/**********************************************************************************************
 * getSeq - pull annotation and protein sequence information from a fasta file stream, get seq 
 * returns no '\n' characters or other special characters that may occur in the protein sequence.
 * requires  - a file stream to a fasta file 
 * returns   - a pointer to a Seq structure with pointers to the next annotation and sequence in the file 
 ***********************************************************************************************/

struct Seq *getSeq(FILE *fasta_stream) 
	{
	int seq_mem = CHUNK_SIZE;
	int c = 0;
	int i = 0;
	int cr = 0;		/* truth value for carriage return after annotation encountered */
	long block = 0;		/* memory used thus far */
	
	char *whole_sequence, *seq_p, *tmp_p;
 	static char *annotation,*sequence;
	static struct Seq seq, *struct_p;

	/* initialize all the pointers to NULL for cleanliness */
	seq_p = tmp_p = whole_sequence = annotation = sequence = NULL;
	struct_p = NULL;

	/* allocate memory for sequence to grab */
	if ((whole_sequence = calloc(sizeof(char),CHUNK_SIZE)) == NULL)
		{
		fprintf(stderr, "error: memory allocation failure\n");
		exit(1);
		}
	seq_p = whole_sequence;	  /* set the placement pointer */

	/* check the stream  */
	if (fasta_stream == NULL)
		{
		fprintf(stderr, "error: received null file stream\n");
		exit(1);
		}			
	c = getc(fasta_stream);
	if (feof(fasta_stream))
		{
		return(NULL);
		}

	/* check we have fasta */
	if (c != '>') 
		{
		fprintf(stderr, "error: failed to find fasta formatted entry -%c-\n",c);
		exit(1);
		}
	*seq_p++ = c; 	 /* put the '>' to our sequence string */
	seq_mem--;		/* adjust the memory tracking */
	block++;
	cr = 0;
	
	/* grab the next annotation and sequence pair */
	for (;;)
		{		
		if (((c = getc(fasta_stream)) == '>') && (cr == 1))
			{
			break;		/* quit at end of sequence entry */ 
			}
		if ((feof(fasta_stream)) || (ferror(fasta_stream))) 
			{
			break;		 /* quit if end of file */
			}		
		if (!seq_mem)	/* allocate more memory as needed */
			{			
			seq_mem = CHUNK_SIZE;
			if ((tmp_p = realloc(whole_sequence,((seq_mem + block + 1) * sizeof(char)))) == NULL)
				{
				fprintf(stderr, "error: memory allocation failure\n");
				exit(1);
				}			
			whole_sequence = seq_p = tmp_p;	 /* reset pointers that may have changed in reallocation */ 
			seq_p += block;
			}
	
		if (c == '\n')		/* only take the line break after the annotation */
			{						/* so we remove all other breaks from the protein string */	
			if(!cr)
				{
				*seq_p++ = c;
				seq_mem--;
				block++;
				cr = 1;
				}
			} 
		else 					/* take everything else for now */
			{
			*seq_p++ = c;
			seq_mem--;
			block++;
			}
		}

	*seq_p = '\0';
	ungetc(c, fasta_stream);	/* put the last '>' back for next pass over stream */
	
	/* allocate memory for splitting the sequence into */
	if ((annotation = malloc(sizeof(char) * (1 + strlen(whole_sequence)))) == NULL)
		{
		fprintf(stderr, "error: memory allocation failure\n");
		exit(1);
		}
	if ((sequence = malloc(sizeof(char) * (strlen(whole_sequence)))) == NULL)
		{
		fprintf(stderr, "error: memory allocation failure\n");
		exit(1);
		}
	struct_p = &seq;

	/* grab the annotation and point the struct to it */
	tmp_p = annotation;
	for (i=0 ; whole_sequence[i] != '\n'; *tmp_p++ = whole_sequence[i++])
		;
	*tmp_p = '\0';
	i++;	/* skip the \n */
	struct_p->annotation =  annotation;
	
	/* grab the sequence and point the struct to it */
	for (seq_p = sequence; whole_sequence[i] != '\0'; i++)
		{
		if(isalpha(whole_sequence[i]))			/* only take valid residues */
			{
			*seq_p++ = whole_sequence[i];
			}
		}
	*seq_p = '\0';
	struct_p->sequence = sequence;

	/* clean up & return the struct */
	free(whole_sequence);	
	whole_sequence = NULL;	
	return(struct_p);
	}

/*return the number of records in the file */
unsigned int count_seqs(char *fasta_file)
	{
	struct Seq *next_seq;
	FILE *fasta;
	unsigned int total_genes = 0;
	
	if ((fasta = fopen(fasta_file,"r")) == NULL)
		{
		fprintf(stderr,"error: unable to open fasta file\n");
		exit(1);
		}

	while ((next_seq = getSeq(fasta)) != NULL)
		{
		if ((next_seq->annotation) != NULL)
			{
			if (strlen(next_seq->sequence) != 0)
				{
				total_genes++;
				}
			}
		}
	return total_genes;
	}

/* reverse complement a sequence */
char *rc(char *seq, int length)
    {
      int x,idx;
      int array_size = length + 1;  // include the terminating '\0'
      char *rc;
      if ((rc = calloc(array_size, sizeof(char))) == NULL)
	{
	  fprintf(stderr, "error: failed to allocate revcomp string memory\n");
	}
      for (x = 0; x < length; x++)
	{
	  idx = length - x - 1;   // get the revcomp offset, building backwards
	  if (seq[x] == 'A')
	    {
	      rc[idx] = 'T';
	    }
	  if (seq[x] == 'C')
	    {
	      rc[idx] = 'G';
	    }
	  if (seq[x] == 'G')
	    {
	      rc[idx] = 'C';
	    }
	  if (seq[x] == 'T')
	    {
	      rc[idx] = 'A';
	    }
	}
      rc[length] = '\0';  // finally, terminate the new seq and return
      return(rc);
    }
