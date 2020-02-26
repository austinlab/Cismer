/* cismer 
 * usage: ./cismer <pssm_output> <cluster.fasta> <genome.fasta>
 *
 * Given a list of PWMs, a fasta file of origin for those PWMs and a fasta background 
 * cismer will compare via non-parametric bootstrapping a distribution for the source and
 * background distributions using a simple measure of the number of sd the mean of a markov
 * sampled version of your cluster lies away from the mean of a markov sampled cluster of randomly pulled genes 
 * 
 * Each random sampling of clusters of equal size to the user provided cluster is mapped with the PWM and 
 * assessed using the simple question of how many of these elements exist in the sampled cluster
 * 
 * Author: Ryan Austin, University of Toronto
 * Department Cell & Systems Biology 2004-2011
 * Present Address : Agriculture & Agri-Food Canada
 * ryan.austin@agr.gc.ca
 *
 * -------------------------------------
 * TODO:
 * - error check on fasta files, pssm file input 
 * - error check on the size of the background and foreground distributions - mismatch in seq lengths can f*up the Z-score
 * - fix memory leak 
 *
 * Revision 2.3  - Nov 2016
 * - changes in menu options and defaults 
 * - increased scan_pwm chunk size to 50,000 from 5,000 to fix whole genome leafy scan bug 
 * - various other things, warnings fix, etc
 * 
 * Revision 2.0.5 
 * - added funcDepth cutoff if statement to mapping output 
 * - modifed output to provide BED format mappings
 * - fixed annotation on genes bug , data being freed in fasta_read 
 * - improved verbosity 
 * - grabbing annotation info 
 * - returning mapping sequence annotation properly 
 * 
 * Revision 2.0.4 - Apr 15
 * - ACGT on rowlabels now optional in read_pssm 
 * - changed -c option to more logical -f 
 * - changed name of program to cismer
 * - mapping just prints first field with no preceeding white space
 * - fixed a bug with BOOSTRAP number on mapping i.e n=1 causing very small realloc's - set to 5000 default 
 * 
 * Revision 1.96 - Apr 7, 2014
 * - implemented mapping function against fasta  cluster 
 * - allowed the use of '>' as an annotation marker for PSSMs
 * - adjusted verbosity behaviour slightly 
 * - added FD calc, carry over of strandedness back to main and reporting on mapping info option
 *
 * Revision 1.95 - Jul 6, 2010
 * rewrote sections of read_pssm to be more osx friendly
 * replaced strcpy/strdup with strncpy/strndup and reworked annotation string concatenation 
 * 
 * Revision 1.94 - Jul 6, 2010
 * changed MAXFIELDS from 50 to 250
 * rare bug for 1000 up data set segfaulting on malloc in scan_pwm
 *
 * Revision 1.93  - Nov 11, 2009
 * changed malloc in log odds function from num-pfms (dumb) to length x width x double
 * miscellaneous comment adds and removals
 * in scan_pwm 	    adjusted - if (num_hits + (seqlen * 2) + 1 > chunk_size * round) - hopefully fixes longstanding segfault issue
 *
 * Revision 1.92 - Sept 16, 2009
 * ok, totally toned down the constants
 * now bs issues
 * maxline 500, maxfields 50  drop from like 5000 ridiculous mis-set
 *
 * Revision 1.91 - Aug 10, 2009
 * -- 11pm - more adjustments , i think these may work chunk 2048 and MAXLINE MAXFIELDS MAXPWMS = 250000
 * - persistent problems with segfaults on reading AlignAce files and segfaults on scan_pwm memory allocation
 * - adjusted callocs to mallocs in scan_pwm, adjusted constants from read_pwm down to moderate levels 
 * - adjusted chunk size for reallocation in scan_pwm  
 * - this is all trial and error at this point but I don't know what else to do with this curse of an issue
 *
 * Revision 1.9 - Sept 2, 2008
 *  - added PSSM reading verbose commenting
 *  - expanded MAXFIELDS to 16020 - sept 4 - changed to 6408
 * 
 * Revision 1.8 - Jun 2, 2008
 * - moved MAXPWMs up to 200,000
 * - MAXFIELDS also adjusted to 1648 at a prior date
 *
 * Revision 1.7 - Mar 12,2008
 * - added switch and code to search revcomp sequence
 *
 * Revision 1.5 - Mar 6, 2008
 * - added callocs to read_pssm to fix segfaults in meme, alignace ? 
 * - extended MAXFIELDS from 30 to 50 - this seems to have worked
 * - MAXLINE to 10,000, was 1000 i think
 *
 * Revision: 1.4 - Mar 4, 2008
 * - extended MAXLINES to fix one ace bug but still another
 * - beginnings of FD return but copped out
 * 
 * Revision: 1.3 - Jan 10, 2008
 * - added protection from division by zero in z-score calc
 * - fixed bug in annotation print (switched a malloc to a calloc)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "seq.h"

#define MAXPWMS  550000    // want to go over 100,000 with 8mer4mm so set to 200K - 215K and 1750 works well but segs on n=9
#define MAXLINE 1500           
#define MAXBOOTSTRAPS  10000
#define MAXSEQS  50000000				 // set to 50Mbases apr15/2014
#define BUFFER_SIZE 2056       
#define MAXFIELDS 250         

/* Methods */
struct fasta_set *read_fasta(char *fasta_file);
struct fasta_set *sample_fasta(struct fasta_set *fasta, int seed, int sample_number);
struct pfm **read_pssm(FILE *promomer);
int pssm_to_pwm(struct pfm *pssms[]);
struct distro *scan_pwm(struct pfm *pwm, struct fasta_set *fasta);
double cmp_distro(struct distro *cluster_distro, struct distro *bkgrnd_distro, double ic_max, double ic_min);
double **log_odds(int **matrix, float GC_content, int length);
unsigned int get_offset(char letter);
double std(int *array, int length);
double ic(double **matrix, int minmax,int length);
void usage();
void help();

/* Typedefs */
struct distro
  {
    int num_hits;
    int *hit_seqs;
    int *hit_sites;    //  3 columns: seq offset information_content
		int *strand;
    float *information_scores;
  };
struct fasta_set     // note the struct contains the whole genome
  {
    char **genes;
    char **annotation;   
    int num_genes;
  };
struct pfm           // this is used as an array of pfm structs for dynamic access
  {
    char *annotation;
    int **matrix;
    int length;        // the length of the motif
    double info_content;
    double **log_matrix;
  };

/* Globals */
static char *version = "2.3.2";   // VERSION NUMBER 

unsigned int bootstraps = 1000;
unsigned int seed = 0;
unsigned int num_pfms = 0;
unsigned int alphabet_size = 4;
unsigned int size_of_coreg_cluster = 0;
unsigned int use_func_depth = 1;
unsigned int mapping = 0;
unsigned int revcomp = 1;
float func_depth_cutoff = 0;  // no FD cutoff by default  
float cutoff = -100;          // print all Z by default 
float GC_content = 0.32;   // GC content of the promoter regions of TAIR7_upstream_500
unsigned int verbose = 0;
unsigned int logodds = 0;  // return a PSSM or a log-odds PWM


/* MAIN */
int main(int argc, char **argv)
	{
	  int c;
	  int promomer_in = 0;   // logical on whether to read pssms on stdin           
	  char *promomer_out = "";
	  char *genome_fasta = "";
	  char *cluster_fasta = "";
    fprintf(stderr, "cismer v%s\n",version);   // display the working version so there's never any confusion
	  while ( ( c = getopt(argc,argv, "p:g:f:h-z:s:b:iwd:vol:m")) != EOF)
	    {
	    switch(c)
	      {
	      case 'p': 
		promomer_out = optarg;
		break;
	      case '-':
		promomer_in = 1;
		break;
	      case 'g':
		genome_fasta = optarg;
		break;
	      case 'f':
		cluster_fasta = optarg;
		break;
	      case 'z':
		cutoff = atof(optarg);
	        break; 
	      case 's':
		seed = atoi(optarg);
		break;
	      case 'b':
		bootstraps = atoi(optarg);
		break;
	      case 'i':
		use_func_depth = 0;
		break;
	      case 'w':
		revcomp = 0;
		break;
	      case 'd':
		func_depth_cutoff = atof(optarg);
		break;
	      case 'v':
		verbose = 2;    
	        break;
	      case 'l':
		GC_content = atof(optarg);
		break;
	      case 'o':
		logodds = 1;
		break;
				case 'm':
					mapping = 1;
					break;
	      case 'h':
		help();
		break;
	      case '?':
	      default:
		usage();
		exit(1);
	      }
	    }
	if (strcmp(promomer_out,"") == 0 && promomer_in == 0) usage();
	if (strcmp(cluster_fasta,"") == 0) usage();
	if (strcmp(genome_fasta,"") == 0 && mapping == 0) usage();

	/* check bootstrapping boundaries */
	if (bootstraps > MAXBOOTSTRAPS)
	  {
	  fprintf(stderr, "error: bootstraps exceed MAX_BOOTSTAPS allowed, decrease.\n");
	  exit(1);
	  }
	if (mapping)
		{
		bootstraps = 1;
		}

	/* load the cluster and background genome sequences  */
	struct fasta_set *cluster = read_fasta(cluster_fasta);
	struct fasta_set *genome;
	if (mapping == 0)  genome = read_fasta(genome_fasta); // don't read genome if just mapping
	
	/* sample gene clusters from fasta */
	unsigned int sample_size = bootstraps * cluster->num_genes;
	struct fasta_set *cluster_sample = sample_fasta(cluster,seed,sample_size);
	struct fasta_set *genome_sample;
	if (mapping == 0) genome_sample = sample_fasta(genome,seed,sample_size);
	size_of_coreg_cluster = cluster->num_genes;
	if (mapping == 0)  // free genome once sampled if bootstrapping 
		{
		free(genome->genes);
		free(genome->annotation);
		free(genome);
		}

	/*  read in PSSMs and convert to log-odds PWM */
	struct pfm **PSSM;
	FILE *file;

	if (verbose) fprintf(stderr,"Beginning PSSM read  . . .\n");
	if(promomer_in)
	  {
	    PSSM = read_pssm(stdin);
	  }
	else
	  {
	    if ((file = fopen(promomer_out,"r")) == NULL)
	      {
				perror("fopen");
				exit(1);
	      }
	    PSSM = read_pssm(file);
	  }
	if (verbose) fprintf(stderr,"Pssm reading complete . . .\n");

	/* convert the frequency matrices to log-odds */
	pssm_to_pwm(PSSM);
	struct pfm **PWM = PSSM;
	if (verbose) fprintf(stderr,"Got log odds . . .\n");

	/*  scan the PWMs one at a time against sampled sequences */
	if (mapping == 0)
		{
		int x;
		for (x = 0; x < num_pfms; ++x)
		  {
		  struct pfm *next_pwm;
		  if ((next_pwm = PWM[x]) == NULL) break;
		  struct distro *bkgrnd_distro;
		  struct distro *cluster_distro;
		  if (verbose == 2) fprintf(stderr,"Scanning cluster distro %i . . .\n",x);
		  cluster_distro = scan_pwm(next_pwm,cluster_sample);
		  if (verbose == 2) fprintf(stderr,"Scanning bkg distro %i . . .\n",x);
			bkgrnd_distro = scan_pwm(next_pwm,genome_sample);
		  
		  /* k, got the distros, do the stats and report for this PWM if appropriate */
		  double imax = ic(next_pwm->log_matrix,next_pwm->length,1);
		  double imin = ic(next_pwm->log_matrix,next_pwm->length,0);
			if (verbose) fprintf(stderr,"%s  ic  max: %f  min: %f\n",next_pwm->annotation,imax,imin);
			double z = cmp_distro(cluster_distro, bkgrnd_distro,imax,imin);
		  if (z > cutoff)
		    {
		      int j,k;
		      printf(">%s Z:%3.2f\n",PWM[x]->annotation,z);
		      for (j = 0; j < alphabet_size; ++j)
						{
			  		for (k = 0; k < PWM[x]->length; ++k)
			    		{
			      	if (logodds)
								{
				  			printf("%2.2f\t",PWM[x]->log_matrix[j][k]);
								}
			      	else
								{
				  			printf("%i\t",PWM[x]->matrix[j][k]);
								}
			    		}
			  		printf("\n");
						}
		      printf("\n");
		    }
	
		  /*FREE the distribution memory for the next PWM */
		  free(cluster_distro->hit_seqs);
		  free(cluster_distro->hit_sites);
		  free(cluster_distro->information_scores);
		  free(cluster_distro);
			free(bkgrnd_distro->hit_seqs);
 		 	free(bkgrnd_distro->hit_sites);
 	 		free(bkgrnd_distro->information_scores);
 	 		free(bkgrnd_distro); 
	 	 }

		/* free the sampling distributions */
		free(cluster_sample->genes);
		free(cluster_sample);
		free(genome_sample->genes);
		free(genome_sample);
		}
	  /* Report on just mapping info  */
	else if (mapping == 1) 
	  {
			int x = 0;
			for (x = 0; x < num_pfms; ++x)
			  {
		 		if (verbose) fprintf(stderr,"Mapping %i . . .\n",x);
			  struct pfm *next_pwm;
			  if ((next_pwm = PWM[x]) == NULL)   break;
				struct distro *map = scan_pwm(next_pwm,cluster);
		  	double imax = ic(next_pwm->log_matrix,next_pwm->length,1);
		 		double imin = ic(next_pwm->log_matrix,next_pwm->length,0);
				double fdepth;
				if (verbose) fprintf(stderr,"ic  max: %f  min: %f\n",imax,imin);
				if (verbose) fprintf(stderr,"Mapping hits found:%i\n",map->num_hits);
				int t = 0;
	      for(t = 0; t < map->num_hits; ++t)
					{
					if (imax == imin)  // do this if the motif is a consensus pattern with no FD variation, avoid illegal /0
						{
						fdepth = 1;
						}
					else
						{
						fdepth = (map->information_scores[t] - imin) / (imax - imin);
						}

					char *a = strdup(PWM[x]->annotation);
					while (isspace(*a)) a++;    // move to character dropping whitespace between annotation marker and motif id
    			char *qp = strtok(a, " \t");		// just take the motif name, leave meta data 
					int pos = map->hit_sites[t] + 1;
					int posStop = map->hit_sites[t] + next_pwm->length;
					int seqIndex = map->hit_seqs[t] - 1;
					char *an = strdup(cluster->annotation[seqIndex]);  // just the first field of the annotation
					char *ann = strtok(an," \t");
					if (fdepth >= func_depth_cutoff) 
						{
						if (map->strand[t] == 0)
							{
							// BED format output:  chrom, start, end, name, score, strand
							printf("%s\t%i\t%i\t%s\t%3.2f\t+\n",ann,pos,posStop,qp,fdepth);
							}
						else
							{
							printf("%s\t%i\t%i\t%s\t%3.2f\t-\n",ann,pos,posStop,qp,fdepth);
							}
						}
					}
	    }
		}
		else
			{
			fprintf(stderr,"Error: invalid mapping value\n"); // should never hit this
			}	

	return(0);
	}


  /***********************/
 /*      METHODS        */
/***********************/

/* READ_FASTA  
 * TAKES: a fasta file
 * DEPENDS: Seq struct and methods from seq.h; MAX_SEQS global definition
 * RETURNS: an array af fasta sequence strings
 */
struct fasta_set *read_fasta(char *fasta_file)
   {
   FILE *fasta;
   char *sequence;
	 char *annot;
   struct Seq *next_seq;
   unsigned int total_genes = 0;
   char **sequences = calloc(MAXSEQS, sizeof(char *));
   char **annots = calloc(MAXSEQS, sizeof(char *));
   struct fasta_set *result = calloc(10,sizeof(struct fasta_set *)); // got me but calloc with size_t of 10 fixes bug on simple run with no first annotation   

   if ((fasta = fopen(fasta_file,"r")) == NULL)
      {
      fprintf(stderr, "error: unable to open fasta file\n");
      exit(1);
      }
   while ((next_seq = getSeq(fasta)) != NULL)
      {
      if ((next_seq->annotation) != NULL)
         {
				 annot = next_seq->annotation; 
				 ++annot;  // move past fasta identifier '>'
				 while (isspace(*annot)) annot++; // remove whitespace between '>' and gene id
				 annots[total_genes] = annot;
         if (strlen(next_seq->sequence) != 0)
            {
            sequence = next_seq->sequence; 
	   			  sequences[total_genes] = sequence;
            total_genes++;
            }
         }
      }
   fclose(fasta);
   result->num_genes = total_genes;
   result->genes = sequences;
	 result->annotation = annots;
   return(result);
   }


/* SAMPLE_FASTA  
 * TAKES: random seed number and the number of sequences to randomly sample
 * RETURNS:  a character array of fasta sequences randomly pulled from fasta with resampling 
 */
struct fasta_set *sample_fasta(struct fasta_set *fasta, int seed, int sample_number)
   {
   unsigned int x;
   struct fasta_set *result = malloc(sizeof(struct fasta_set *));
   char **random_fasta = calloc(MAXSEQS, sizeof(char *));
   int *random_array = (int *)malloc(sizeof(int) * (sample_number) + 2);

   /* check sample sizes */
   if (sample_number > MAXSEQS)
     {
       fprintf(stderr, "error: sample number exceeds MAXSEQS, try decreasing fasta size or bootstrap number\n");
       exit(1);
     }
   
   /* get an array of random numbers to use for sampling */
   srand(time(NULL));
   for (x = 0; x < sample_number; ++x)
      random_array[x] = ( rand() % fasta->num_genes );

   /* pull the random sets from the fasta set */
   if ((result = malloc(sizeof(struct fasta_set))) == NULL)
     {
       fprintf(stderr,"error: memory allocation failure in sample_fasta\n");
       exit(0);
     }
   result->num_genes = sample_number;

   char **sequences = fasta->genes;
   for (x = 0; x < sample_number; ++x)
      {
       int index = random_array[x];
       *(random_fasta + x) = strdup(sequences[index]);
      }
   result->genes = random_fasta;
   free(random_array);
   return(result);
   }


/* delspace - remove spaces from a string 
 *  failed attempt to fix a nasty bug in some bioprospector output - leading spaces on matrix values 
 * lifted from http://www.daniweb.com/forums/thread21071.html# 
 */
void delspace (char *Str)
	{
	int Pntr = 0;
	int Dest = 0;
	while (Str[Pntr])
		{
		if (Str[Pntr] != ' ')
			{
			Str[Dest++] = Str[Pntr];
			}
		Pntr++;	
		}
	Str[Dest] = 0;
	}

/* READ_PSSM
 * TAKES:  the name of a file containing PSSMs ala promomer2 output format
 * RETURNS: an array of PFM structs containing the PSSM matrices
 * GLOBALS: uses none but sets the global num_pfms for the rest of the routines to use, saving memory
 */
struct pfm **read_pssm(FILE *promomer)
   {
   char buf[BUFFER_SIZE], *p, *q;
   char *annotation_marker = "//";      // what does the first field of an annotation line look like
	 char *annotation_marker2 = ">";      // optional annotation marker 
   num_pfms = 0;               // start at -1 so the first assignment is 0 after autoincrementer on GLOBAL
   unsigned int row = -1;
   unsigned int annotation_logical = 0;
   unsigned int i;
	 unsigned int rowlabels = 0;
   struct pfm **result = malloc(MAXPWMS * sizeof(struct pfm *));
   
   /* parse the PSSM file for the PSSM outputs 
    * need an array of pfms, each with a matrix and information content */
while (1)
  {
  fgets(buf, sizeof buf, promomer);
  if (feof(promomer))  break;
  for (p = strtok(buf, "\n"); p; p = strtok((char *)NULL, "\n"))
    {
    for (q = strtok(buf, "\t"), i=0; q; q = strtok((char *)NULL, "\t"),i++)
    {
    if (i == 0)  // FIRST FIELD OF LINE BLOCK - see what we have
      {
	 	  if ((strcmp(q,annotation_marker) == 0) || strcmp(q,annotation_marker2) == 0)             // set annotation flag and increment pfm count
		     {
		     annotation_logical = 1;   // ANNOTATION BLOCK
		    	row = -1; 
		     /* DYNAMIC allocation of next PWM struct */
		     int **mtx_pointer;                               // pwm placeholder
		     char *annot_pointer;
		     if ((result[num_pfms] = malloc(sizeof(struct pfm))) == NULL)   // the pwm struct
		       {
					 fprintf(stderr,"error: memory allocation failure in read_psm1\n");
					 exit(1);
		       }
		     if ((mtx_pointer = calloc(alphabet_size, sizeof(int *))) == NULL)
		       {
					 fprintf(stderr, "error: memory allocation failure in read_psm2\n");
					 exit(1);
		       }
		     if ((annot_pointer = calloc(MAXLINE, sizeof(char))) == NULL)
		       {
					 fprintf(stderr, "error: memory allocation failure in read_psm3\n");
					 exit(1);
		       }
		     int v;
		     for (v = 0; v < alphabet_size; ++v)
		       {
					 mtx_pointer[v] = calloc(MAXFIELDS, sizeof(int));
		       }
		     result[num_pfms]->matrix = mtx_pointer;       // set the data initial values
		     result[num_pfms]->length = (int)NULL;
		     result[num_pfms]->annotation = annot_pointer;
		     num_pfms++;       // increment our pfm count here ( note we'll be 1 up from present pwm now)
		     }
		 //  NOT ANNOTATION - HAVE ROW LABELS
		 else if (strcmp(q,"A") == 0 || strcmp(q,"C") == 0 || strcmp(q,"G") == 0 || strcmp(q,"T") == 0)                      // start of a matrix, reset the annotation flag
		   {
		 	 if (strcmp(q,"A") == 0) // this block is when there are row labels in matrix
				{ 
				row = 0;  // probably not necessary, maybe assert? 
				}
			else
			 	{
				row++;
				}
			rowlabels = 1;
		 	annotation_logical = 0;
		  }
		 else         // THIS IS WHEN THERE ARE NO ROW LABELS                                   
		  {
		 	 row++;
			 annotation_logical = 0;
  		 int **mtx = result[num_pfms - 1]->matrix;
			 mtx[row][i] = atoi(q);
		  }
		 }
   else   // ALL OTHER FIELDS IN THE LINE BLOCK 
     {
		 if (annotation_logical)   
		     {
			   char *r = result[num_pfms - 1]->annotation;
		     if (r == NULL)
				 	{
				 	r = strncpy(r,q,(size_t)MAXLINE -1);
				 	}
		     else
			 		{
			   	r = strncat(r,"\t",(size_t)2);
			   	r = strncat(r,q,(size_t)500);
			 	  result[num_pfms - 1]->annotation = r;
					if (verbose) { fprintf(stderr,"Read PSSM %s . . .\n",r); }
			 		}
		     }
		 else        // store the field value in the matrix 
		     {
				 int pos;
		  	 int **mtx = result[num_pfms -1]->matrix;
				 if (rowlabels) { pos = i-1; } else { pos = i; }
		     mtx[row][pos] = atoi(q);
	  	   }
		 }
   }
	 int len;
	 if (rowlabels) { len = i-1; } else { len = i; }
	 result[num_pfms - 1]->length = len;
	}
 }
fclose(promomer);
return(result);
}


/* PSSM_TO_PWM
 * TAKES: an array of PFMs in the form of an array of pfm structs
 * RETURNS: the same array of PFMs, with the log_odds matrix added 
 * GLOBALS: uses num_pfms for better memory allocation and GC_content
 */
int pssm_to_pwm(struct pfm *pssms[])
  {
    int x;
    for ( x = 0; x < num_pfms; ++x)
      {
	int **mtx = pssms[x]->matrix;
	int length = pssms[x]->length;
	pssms[x]->log_matrix = log_odds(mtx,GC_content,length);
      }
    return(0);
  }

/* LOG_ODDS  -  a helper function of PSSM_TO_PWM
 * TAKES: a 2D matrix of ints
 * RETURNS: a 2D matrix of floats of the int matrix converted to log-odds ratios
 */
double **log_odds(int **matrix, float GC, int length)
  {
    int a,b,x,y;
    int N = 0;        // the total number of measurements
    float G = GC / 2;
    float AT_content = (double) 1 - GC;
    float A = AT_content / 2;
  
    /* allocate space for the log-odds matrix */
    if (verbose == 1) fprintf(stderr,"Allocating for log-odds...\n");
    double **lo_matrix = malloc(sizeof(double *) * alphabet_size * length);    // changed from num_pfms - nov11-09
    for (a = 0; a < alphabet_size; ++a)
	lo_matrix[a] = malloc(sizeof(double) * length + 1);

    /* calculate N by adding the first column*/
    for (b = 0; b < alphabet_size; ++b)
	    N += matrix[b][0];

    /* column wise log odds scoring 
     * this is based on log2( Fi,j / Pi )  where Pi is the base prob.(GC)
     * and Fi,j is the base frequency (count at position i,j / N )
     * Stormo 2000
     */
    for (x = 0; x < alphabet_size; ++x)
      {
	for (y = 0; y < length; ++y)
	  { 
	    float base_frequency = (float)matrix[x][y] / N;	    
	    if (((x+1) % alphabet_size) < 2)  // being clever here on the ACGT selection 1 and 4 (A and T) mod 4 < 2
	      {
		lo_matrix[x][y] = log(base_frequency/A)/log(2);
    	      }
	    else
	      {
		lo_matrix[x][y] = log(base_frequency/G)/log(2);
	      }
 	  }
      }
    if (verbose) fprintf(stderr,"Got log-odds...\n");
    return(lo_matrix);
  }

 /* SCAN_PWM 	
 *  TAKES: a pfm struct and a fasta_set struct
 *  RETURNS: a distribution struct detailing the occurances of all pwms within the fasta set exceeding a significance score
 *  GLOBALS:  num_pfm
 */
struct distro *scan_pwm(struct pfm *pwm, struct fasta_set *fasta)
{
	int x,y;
	unsigned int mer_length = pwm->length;
	unsigned int gene_num = fasta->num_genes;
	double **logodds = pwm->log_matrix;
	char **seqs = fasta->genes;
	char *next_seq;
	char *rc_seq = "";      // as the calloc is in the rc method, initialize to null to avoid warning
	char letter, next_letter;
	char *motif = pwm->annotation;
	int seqlen, idx, index, offset;
	float info_content;
	int round = 1;           // how many memory expansions of the dynamic arrays we've done

	/* DISTRIBUTION struct MEM-ALLOCATION */
	int chunk_size = 20000 * 2056; //  5000 * 2056;   /// the number of array elements to expand the array by when needed.  from 2048 -nov11-09
	int *hseqs;                 // changed this from 10000 to 50000 Apr5-09 to fix lfy pseudoscan segfault - 20k aug6-09
	int *hsites;
	int *strand;
	float *iscores;
	struct distro *distribution;
	if ((distribution = malloc(sizeof (struct distro))) == NULL)
	  {
	    fprintf(stderr, "error: memory allocation failure in scan_pwm\n");
	    exit(1);
	  }
	if((hseqs = malloc(chunk_size * sizeof(int))) == NULL)    // 3 column array of sequence:motif_offset:information_content
	  {
	    fprintf(stderr,"error:memory allocation failure in scan_pwm:hit_seqs\n");
	    exit(1);
	  }
	if ((hsites = malloc(chunk_size * sizeof(int))) == NULL)
	  {
	    fprintf(stderr,"error:memory allocation failure in scan_pwm:hit_sites\n");
	    exit(1);
	  }
	if ((strand = malloc(chunk_size * sizeof(int))) == NULL)
		{
	    fprintf(stderr,"error:memory allocation failure in scan_pwm:hit_sites\n");
	    exit(1);
		}
	if((iscores = malloc(chunk_size * sizeof(float))) == NULL)
	  {
	    fprintf(stderr,"error:memory allocation failure in scan_pwm:info_scores\n");
	    exit(1);
	  }

	// scanning each sequence one at time
	int num_hits = 0;
	for (x = 0; x < gene_num; ++x)
	  {
	    if ((next_seq = seqs[x]) == NULL)
	      {
				fprintf(stderr,"error: failed to retrieve next sequence\n");
				exit(0);
	      }
	    seqlen = strlen(next_seq);
	    if (seqlen < mer_length)  break;

	    /* grab the reverse complement sequence if flag set */
	    if (revcomp) rc_seq = rc(next_seq,seqlen);

	    /* realloc more memory to the arrays if needed */
	    if (num_hits + (seqlen * 2) + 1 > chunk_size * round)
	      {
				round++;    // increment the round
				if((hseqs = realloc(hseqs,chunk_size * round * sizeof(int))) == NULL)
				  {
		 		   fprintf(stderr,"error: memory reallacation error in scan_pwm\n");
			    exit(1);
				  }
				if((hsites = realloc(hsites,chunk_size * round * sizeof(int))) == NULL)
				  {
			    fprintf(stderr, "error: memory reallocation error in scan_pwm\n");
			    exit(1);
				  }
				if ((strand = realloc(strand,chunk_size * round * sizeof(int))) == NULL)
					{
		 		   fprintf(stderr, "error: memory reallocation error in scan_pwm\n");
		 		   exit(1);
					}
				if((iscores = realloc(iscores,chunk_size * round * sizeof(float))) == NULL)
				  {
		    	fprintf(stderr, "error: memory reallocation error in scan_pwm\n");
		   	  exit(1);
				  }
	      }

	    /* look at every position along the seq length */
	    for (y = 0; y < seqlen - mer_length + 1; ++y)
	      {
				letter = next_seq[y];
				offset = 0;
				info_content = -1000;
			
				if ((idx = get_offset(letter)) == -1) continue;

			// if the first letter is in the matrix try the rest of the motif, short circuiting if we hit a -inf
	  		if ((info_content = logodds[idx][0]) > -100 )     // officially log-odds yields -Inf but -100 close enough for niceness
				 {
	 		  index = y + 1;
	 		  for (offset = 1; offset < mer_length; ++offset,++index)  // run through the rest of the matrix columns
	       {
			 		next_letter = next_seq[index];
		 			if ((idx = get_offset(next_letter)) == -1) break;
					if (logodds[idx][offset] < -100) break;
					info_content += logodds[idx][offset];
					if (offset == mer_length - 1)
					  {
				    if(verbose == 2) printf("%i\t%s\tsequence:%i\tposition:%i\tIC:%f\tstrand:+\n",num_hits,motif,x + 1,y,info_content);
				    hseqs[num_hits] = x + 1;
				    hsites[num_hits] = y;
						strand[num_hits] = 0;
				    iscores[num_hits] = info_content;
				   	num_hits++;
				  	}
			    }
			   }
	    	}
	    /* now do the reverse complement strand if needed */
	    if (revcomp)
	      {
			for (y = 0; y < seqlen - mer_length + 1; ++y)
			  {
		    letter = rc_seq[y];
		    offset = 0;
		    info_content = -1000;
		    if ((idx = get_offset(letter)) == -1) continue;

		    // if the first letter is in the matrix try the rest of the motif, short circuiting if we hit a -inf
		    if ((info_content = logodds[idx][0]) > -100 )     // officially log-odds yields -Inf but -100 close enough for niceness
		      {
					index = y + 1;
					for (offset = 1; offset < mer_length; ++offset,++index)  // run through the rest of the matrix columns
					  {
			    	next_letter = rc_seq[index];
			  	  if ((idx = get_offset(next_letter)) == -1) break;
			    	if (logodds[idx][offset] < -100) break;
			 	    info_content += logodds[idx][offset];
			   	  if (offset == mer_length - 1)
			     	 {
			 			if(verbose == 2) printf("%i\t%s\tsequence:%i\tposition:%i\tIC:%f\tstrand:-\n",num_hits,motif,x + 1,y,info_content);
					 		hseqs[num_hits] = x + 1;
							hsites[num_hits] = seqlen - y - mer_length;  // account for the rc 
							strand[num_hits] = 1;
							iscores[num_hits] = info_content;
							num_hits++;
			       }
			  	 	}
		     	}
		  	}
				free(rc_seq);
	     }
	  }
	/* load our result object */
	distribution->hit_seqs = hseqs;
	distribution->hit_sites = hsites;
	distribution->information_scores = iscores;
	distribution->num_hits = num_hits;
	distribution->strand = strand;
	return(distribution);
	}

// Take a character ACGT bp and convert to 0 1 2 3 offsets respectively - helper to scan_pwm 
unsigned int get_offset(char letter)
  {
    unsigned int idx;
  if (letter == 'A')         // translate the residue into matrix row offset
    idx = 0;
  else if (letter == 'T')  // promoters are AT rich so put them first for optimization
    idx = 3;
  else if (letter == 'G')
    idx = 2;
  else if (letter == 'C')
    idx = 1;
  else
    return(-1);       // Thus note, X's or N's in promoter sequence won't match a motif
  return(idx);
  }

/* CMP_DISTRO_STDEV  
 * cmp_distro_stdev	
 * TAKES:  two distribution structs, one for the sampled cluster and one for sample background
 * RETURNS: a distro_stat struct detailing the statistical significance of the comparison
 * GLOBALS: size_of_coreg_cluster which is not reflected in the distro sample
 */
double cmp_distro(struct distro *cluster_distro, struct distro *bkgrnd_distro, double ic_max, double ic_min)
	{
	  int cluster_counts[bootstraps];
	  int bkgrnd_counts[bootstraps];

	  int clusterhitcount = cluster_distro->num_hits;	
	  int bkgrndhitcount = bkgrnd_distro->num_hits;

	  int cluster_size = size_of_coreg_cluster;  // global variable
	  int *cluster_hits = cluster_distro->hit_seqs;
	  int *bkgrnd_hits = bkgrnd_distro->hit_seqs;
	  
	  float *cluster_ic = cluster_distro->information_scores;
	  float *genome_ic = bkgrnd_distro->information_scores;

	  int cluster_total, bkgrnd_total;
	  double bkgrnd_mean, cluster_mean, bkgrnd_stdev;
	  double z_score = 0.0;
	  double funcdepth = 0;
	  
	  int fd = use_func_depth;
	  int x;
	  for (x = 0; x < bootstraps; ++x)  {   cluster_counts[x] = 0;  bkgrnd_counts[x] = 0;   }
	  if (ic_max == ic_min)  fd = 0;  // turn off the funcdepth cutoff if the max and min scores are equal, i.e. no degeneracy in the matrix

	  /*
	    For each bootstrap, there is a cluster of number genes equal in size to the number of genes submitted by the user
	     Scores are tallied by gene number of bootstrap * cluster_size samplings.  We need to run through the hits, see if the 
	     hit_number falls within the limits of that cluster (i.e. cluster 12 will have seq values between 12*cluster_size and 13 * cluster_size -1)
	     and tally the counts for that cluster and add it to the counts array 
	  */
	  int cluster, hit_idx;
	  cluster_total = 0;
	  for (cluster = 1, hit_idx = 0; cluster <= bootstraps; ++cluster)
	    {
	      int c_tally = 0;
	      while (cluster_hits[hit_idx] <= cluster_size * cluster && hit_idx < clusterhitcount)
					{
				  if (verbose == 2) fprintf(stderr,"hit: %i\n",cluster_hits[hit_idx]);
				  if (fd)
		  		  {
		 			   funcdepth = (cluster_ic[hit_idx] - ic_min)/(ic_max - ic_min);
		 			   if (funcdepth < func_depth_cutoff)
		    		  {
							hit_idx++;
							continue;
		    		  }
		   		 }
		   		 c_tally++;
		   		 hit_idx++;
					}
	      cluster_counts[cluster] = c_tally;
	     	if (verbose == 2) fprintf(stderr,"cluster: %i counts %i\n",cluster,c_tally);
	    }
	  for (x = 0; x < bootstraps; ++x) cluster_total += cluster_counts[x];
	  cluster_mean = (double) cluster_total / (double) bootstraps;

	  /* now the background totals */
	  for (cluster = 1, hit_idx = 0; cluster <= bootstraps; ++cluster)
	    {
	      int b_tally = 0;
	      while (bkgrnd_hits[hit_idx] <= cluster_size * cluster && hit_idx < bkgrndhitcount)
					{
				  if (verbose == 2) fprintf(stderr,"bhit: %i\n",bkgrnd_hits[hit_idx]);
				  if (fd)
				    {
		   		   funcdepth = (genome_ic[hit_idx] - ic_min)/(ic_max - ic_min);
		   	  	 if (funcdepth < func_depth_cutoff)
							{
						  hit_idx++;
			 				continue;
							}
		    		}
		  			b_tally++;
		  			hit_idx++;
					}
	      bkgrnd_counts[cluster] = b_tally;
	     	if (verbose == 2) fprintf(stderr,"b_cluster: %i counts %i\n",cluster,b_tally);
	    }
		  bkgrnd_total = 0;
	 		for (x = 0; x < bootstraps; ++x)
	      bkgrnd_total += bkgrnd_counts[x];
	  	bkgrnd_mean = (double) bkgrnd_total / (double) bootstraps;
	  	bkgrnd_stdev = std(bkgrnd_counts,bootstraps);

	  /* calculate the z-enrichment score */
		if (bkgrnd_stdev == 0)  // protect against illegal division by zero
			{
			z_score = 0.0;
			}
		else
			{
		 	 z_score = (cluster_mean - bkgrnd_mean) / bkgrnd_stdev;
			}
	    if (verbose == 2) fprintf(stderr,"Cmean:%f Bmean:%f Bsd:%f Z:%f\n",cluster_mean,bkgrnd_mean,bkgrnd_stdev,z_score);
	  return(z_score);
	}

/* return the standard deviation of an array */
double std(int *array, int length)
{
  double mu = 0;
  double sd = 0;
  double residuals = 0;
  
  // get the mean of the array first
  int y;
  for (y = 0; y < length; ++y)
    {
      mu += (double) array[y];   // get the sum of the array
    }
  mu /= (double)length;  // divide by N for mean

  // now calculate the residuals 
  for (y = 0; y < length; ++y)
    {
      double dif = (double) array[y] - mu;  // observed - expected
      dif *= dif;                            // squared residual
      residuals += dif;                     // summed
    }
  sd = (double) residuals / length;       // averaged to get variance
  sd = sqrt(sd);                         // sqrt to get stdev
  return(sd);
}


/* return the maximum or minimum information content value for a pwm */
double ic(double **pwm, int length, int minmax)
{
  double ic_minmax;
  int x,y;

  if (minmax)  // determine max IC score
    {
      ic_minmax = 0;
      for (x = 0; x < length; ++x)
	{
	  double best = -100;   // set this to our standard really low number (again the matrix will be -inf)
	  for (y = 0; y < alphabet_size; ++y)
	    {
	      if (pwm[y][x] > -100 && pwm[y][x] > best)
		best = pwm[y][x];
	    }
	  ic_minmax += best;
	}
      return(ic_minmax);
    }
  else     // determine min IC score
    {
      ic_minmax = 0;
      for (x = 0; x < length; ++x)
	{
	  double best = 100;
	  for (y = 0; y < alphabet_size; ++y)
	    {
	      if (pwm[y][x] > -100 && pwm[y][x] < best)
		best = pwm[y][x];
	    }
	  ic_minmax += best;
	}
      return(ic_minmax);
    }
}

/* print a usage message */
void usage()
{
  printf("cismer -p <pssms> -g <genome.fasta> -f <target.fasta> [-hvaim][-r cutoff][-d cutoff][-s seed][-z cutoff]\n");
  exit(1);
}


/* print a detailed usage and help message */
void help()
{
  printf("cismer -p <pssms> -g <genome_background.fasta> -f <cluster.fasta> [-hvmo][-b bootstraps][-d cutoff][-z cutoff]\n");
  printf("\t-\t- read the pssm output from stdin\n");
  printf("\t-p\t- a file listing formatted PSSMs\n");
  printf("\t-g\t- the genome background in fasta format\n");
  printf("\t-f\t- the cluster file of target genes in fasta format\n");
  printf("\t-b\t- the number of bootstraps to perform (default:1000)\n");
  printf("\t-z\t- the z-score cutoff to use for significance (default: none\n");
  printf("\t-d\t- the functional depth cutoff to use (default: 0)\n");
 // printf("\t-s\t- a random seed to use (default systime)\n");
 // printf("\t-i\t- turn off information content functional depth limit (default: on)\n");
 // printf("\t-w\t- only search the Watson strand (default: both)\n");
  printf("\t-o\t- return log-odds PWMs instead of PSSMs\n");
	printf("\t-m\t- just map the pssms against the fasta file (-f)\n");
  printf("\t-h\t- print this usage message\n");
  exit(1);
}
