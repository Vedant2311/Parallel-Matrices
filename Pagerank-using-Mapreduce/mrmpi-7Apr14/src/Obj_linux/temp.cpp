/* ----------------------------------------------------------------------
   MR-MPI = MapReduce-MPI library
   http://www.cs.sandia.gov/~sjplimp/mapreduce.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the modified Berkeley Software Distribution (BSD) License.

   See the README file in the top-level MapReduce directory.
------------------------------------------------------------------------- */

// MapReduce word frequency example in C++
// Syntax: wordfreq file1 dir1 file2 dir2 ...
// (1) read all files and files in dirs
// (2) parse into words separated by whitespace
// (3) count occurrence of each word in all files
// (4) print top 10 words

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include <bits/stdc++.h>

using namespace MAPREDUCE_NS;
using namespace std;

int num_webpages;
vector<vector<int> > outgoing_links(100000);
vector<double> pageranks(100000, 0.0f);
vector<double> pageranks_up(100000, 0.0f);

// void fileread(int, char *, KeyValue *, void *);
void fileread(int itask, KeyValue *kv, void *ptr);
void update(uint64_t, char *, int, char *, int, KeyValue *, void *);
void sum(char *, int, char *, int, int *, KeyValue *, void *);
int ncompare(char *, int, char *, int);
// void output(uint64_t, char *, int, char *, int, KeyValue *, void *);
/*
struct pagerank {
  int n,limit,flag;
};*/

/* ---------------------------------------------------------------------- */

int main(int narg, char **args)
{
  MPI_Init(&narg,&args);

  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  double alpha = 0.85;
	double conv = 0.00001;
	int a,b,i;
	ifstream fopen;
	fopen.close();
	fopen.open(args[1]);
	while(!fopen.eof()){
		fopen>>a>>b;
		num_webpages = max(num_webpages,max(a,b));
    outgoing_links[a].push_back(b);
	}
	fopen.close();
	num_webpages++;

  num_webpages = 20000;

  double pgr = double(1.0f/num_webpages);
	for(int i=0; i<num_webpages; i++){
		pageranks[i] = pgr;
	}
  int tt = 0;

  double tstart = MPI_Wtime();
  MapReduce *mr;
  while (true) {
    double dp =0.0f;
    for(i=0; i<num_webpages; i++){
      if(outgoing_links[i].size()==0){
    		dp += pageranks[i];	// sumlist.push_back(value);
    	}
    }
    dp = double(dp/num_webpages);

    MPI_Barrier(MPI_COMM_WORLD);
    mr = new MapReduce(MPI_COMM_WORLD);
    // mr->verbosity = 2;
    // mr->timer = 1;
    //mr->memsize = 1;
    //mr->outofcore = 1;

    int nwords = mr->map(nprocs,fileread,NULL);
    int nfiles = mr->mapfilecount;

    // cout << "gather -> " <<  << endl;
    mr->gather(1);
    mr->convert();

    int nunique = mr->reduce(sum,NULL);
    mr->broadcast(0);
    // mr->sort_keys(1);

    MPI_Barrier(MPI_COMM_WORLD);
    mr->map(mr,update,NULL);



    delete mr;
    for(i=0; i<num_webpages; i++){
			pageranks_up[i] = (double)alpha*pageranks_up[i] + alpha*dp + (double)(1-alpha)/num_webpages;
		}
    bool converging = true;
		for(i=0; i<num_webpages; i++){
			if(pageranks[i]-pageranks_up[i]>conv)
				converging=false;
        pageranks[i] = pageranks_up[i];
  			// pageranks_up[i] = 0;
		}
    // double sumer = 0.0;
  	// for(int i=0; i<num_webpages; i++){
  	// 	// cout<<i<<" = "<<pageranks[i]<<endl;
  	// 	sumer += pageranks[i];
  	// }
  	// cout<<"sum "<<sumer << endl;
    // if(me == 0){
    //   std::cout << "pageranks[0] = " << dp << '\n';
    // }
		if(converging && (me == 0))
			break;

    tt++;
  }

  std::cout << "Number of loop -> " << tt << '\n';

  double sumer = 0.0;
	for(int i=0; i<1000; i++){
		cout << i << " = " << pageranks[i] << endl;
		sumer += pageranks[i];
	}

  // MPI_Finalize();

	return 0;


}

void fileread(int itask, KeyValue *kv, void *ptr)
{
  // filesize = # of bytes in file

  int world_rank, num_procs;
  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // cout << " In World rank -> " << ((world_rank*20000)/num_procs) << " " << (((world_rank+1)*20000)/num_procs) << endl;

  for (int i = ((world_rank*20000)/num_procs); i < (((world_rank+1)*20000)/num_procs); i++) {
    // cout << i << endl;
    double pgi = 0;
    kv->add((char *) &i,sizeof(int),(char *) &pgi,sizeof(double));
    if(outgoing_links[i].size() == 0){
      continue;
    }
    for(int j =0; j<outgoing_links[i].size(); j++){
      // cout << j << " of " << outgoing_links[i].size() << " of webpage " << i << " of " << (world_rank+1)*(20000/num_procs) << " (" << world_rank << ") " << endl;
      double pg = (double)(pageranks[i]/outgoing_links[i].size());
      // std::cout << "OK " << pg << '\n';
      kv->add((char *) &outgoing_links[i][j],sizeof(int),(char *) &pg,sizeof(double));
		}
  }

  // std::cout << "____________" << '\n';

  // test = world_rank;


}

void update(uint64_t itask, char *key, int keybytes, char *value,
	    int valuebytes, KeyValue *kv, void *ptr)
{
  int keyint = *(int *) key;
  double pgrank = *(double *) value;
  pageranks_up[keyint] = pgrank;
}

/* ----------------------------------------------------------------------
   read a file
   for each word in file, emit key = word, value = NULL
------------------------------------------------------------------------- */

/*
void fileread(int itask, char *fname, KeyValue *kv, void *ptr)
{
  // filesize = # of bytes in file

  cout << "______" << itask  << "\n";

  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  int filesize = stbuf.st_size;

  FILE *fp = fopen(fname,"r");
  char *text = new char[filesize+1];
  int nchar = fread(text,1,filesize,fp);
  text[nchar] = '\0';
  fclose(fp);

  char *whitespace = " \t\n\f\r\0";
  char *word = strtok(text,whitespace);
  while (word) {
    kv->add(word,strlen(word)+1,NULL,0);
    word = strtok(NULL,whitespace);
  }

  delete [] text;
}

/* ----------------------------------------------------------------------
   count word occurrence
   emit key = word, value = # of multi-values
------------------------------------------------------------------------- */

void sum(char *key, int keybytes, char *multivalue,
	 int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
{
  int t = *(int *) key;
  double* s = (double *) multivalue;
  double sum = 0;
  for(int i = 0; i < nvalues; i++){
    sum += s[i];
  }
  kv->add(key,keybytes,(char *) &sum,sizeof(double));
}

/* ----------------------------------------------------------------------
   compare two counts
   order values by count, largest first
------------------------------------------------------------------------- */

int ncompare(char *p1, int len1, char *p2, int len2)
{
  int i1 = *(int *) p1;
  int i2 = *(int *) p2;
  if (i1 > i2) return -1;
  else if (i1 < i2) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   process a word and its count
   depending on flag, emit KV or print it, up to limit
------------------------------------------------------------------------- */
/*
void output(uint64_t itask, char *key, int keybytes, char *value,
	    int valuebytes, KeyValue *kv, void *ptr)
{
  Count *count = (Count *) ptr;
  count->n++;
  if (count->n > count->limit) return;

  int n = *(int *) value;
  if (count->flag) printf("%d %s\n",n,key);
  else kv->add(key,keybytes,(char *) &n,sizeof(int));
}*/

/*
Before compilation
 echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope

Compile
 make -f Makefile.pagerank
 mpic++ wordfreq.o ../src/libmrmpi_mpicc.a -o wordfreq

Run
mpirun -np 8 pagerank file1 file2

*/
