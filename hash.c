/* functions for deploying hash lookup table */

#include "hash.h"

/* insert a string into the hash table */
int install(char *name, int i)
	{
	struct nlist *np;
	unsigned int hashval;
	int hashsize = HASHSIZE;

	// hash entry doesn't exist, create it
	if ((np = lookup(name)) == NULL)
		{
		  np = malloc(sizeof(struct nlist));  // this was changed from cast to (struct nlist *)
		  int k;
		  for (k = 0; k < MAX_SEQS; k++)  // set the hits array to 0
		    np->hits[k] = 0;
		  if ((np->name = strdup(name)) == NULL)
		    {
		      fprintf(stderr,"error: NULL value for strdup in hash.c\n");
		      return(1);
		    }
		  hashval = hash(name,hashsize);
		  np->next = hashtab[hashval];		// point next to hashtable entry (NULL or last struct)
		  np->hits[i]++;
		  hashtab[hashval] = np;				// point hashtable entry to new struct
		}
	// hash entry exists, increment the hit counter
	else
	  {
	    np->hits[i]++;
	  }
	return(0);
	}

/* find a particular string in the hash table */
struct nlist *lookup(char *s)
	{
	struct nlist *np;
	int hashsize = HASHSIZE;
	for (np = hashtab[hash(s,hashsize)]; np != NULL; np = np->next)
	  {
	    if (strcmp(s, np->name) == 0)
	      {
		return np;
	      }
	  }
	return NULL;					
	}

/* the hash key - D.J. Bernstein hash  */
unsigned int hash(char* str, int hash_size)
	{
	unsigned int hash = 5381;
	unsigned int i    = 0;
	unsigned int len = strlen(str);
	for(i = 0; i < len; str++, i++)
		{
		hash = ((hash << 5) + hash) + (*str);
		}
	hash =  (hash & 0x7FFFFFFF);
	return hash % hash_size;
	}



/* the hash key - K&R 2nd edition, p144 
unsigned int hash(char *s, int hash_size)
	{
	unsigned int hashval;
	for (hashval = 0; *s != '\0'; s++)
		{
		hashval = *s + 31 * hashval;
		}
	return hashval % hash_size;
	}
*/

/* return the keys to the hash 
char **keys()
	{
	char *keys[HASHSIZE * 2];
	int count = 0;
	struct nlist *np;
	
	for (int x = 0; x < HASHSIZE; ++x)
		{
		for(np = hashtab[x]; np != NULL; np = np->next)
			{
			keys[count++] = np->name;		
			}
		}
	return keys;
	}
*/
