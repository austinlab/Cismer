/* hash table functions header file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define HASHSIZE 89000
#define MAX_SEQS 201

struct nlist {
	struct nlist *next;
	char *name;
	int hits[MAX_SEQS];
	};

struct nlist *hashtab[HASHSIZE];

int install(char *name, int i);
struct nlist *lookup(char *s);
unsigned int hash(char* str, int hashsize);
//char *strdup(char *);
char **keys();

