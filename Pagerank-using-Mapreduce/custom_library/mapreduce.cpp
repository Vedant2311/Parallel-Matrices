#include <bits/stdc++.h>
#include "mapreduce.h"
#include "mpi.h"
using namespace std;

// Initializing the class
MapReduce::MapReduce(int total_webpages)
{
	keyvalue.clear();
	keyvalue.resize(total_webpages, 0.0);
}

void MapReduce::map(vector<vector<int> > outgoing_links,vector<double> pageranks,int total_webpages){
	
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	for (int i=0; i<total_webpages/size; i++){
		int key = i + (rank * total_webpages/size);
		double page = pageranks[key];
		
		int children = outgoing_links[key].size();
		if (children!=0){
			for(int j=0; j<children; j++){
				int key_inner = outgoing_links[key][j];
				double value = (double)(page/children);
				
				// Partially performs the reduction during mapping to remove the complications related to the children of the nodes not being distributed in a contiguous manner
				keyvalue[key_inner]+=value;
			}
		}
		
	}

}

void MapReduce::reduce(int total_webpages){
	
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	if(rank==0){
		vector<double> temp(total_webpages);
		
		for (int j=1; j<size; j++){
			MPI_Recv(&temp[0], total_webpages, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			// Getting the sum value for all the different keyvalue vectors calculated by the different processes
			for (int i=0; i<keyvalue.size();i++){
				keyvalue[i] = keyvalue[i] + temp[i];
			}
		}
		

	}
	
	else{
		MPI_Send(&keyvalue[0], total_webpages, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&keyvalue[0],total_webpages,MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

	


