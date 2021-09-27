// 214101023_Kmeans.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"       
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

#define P 12
#define DELTA 0.0001

// global parameters
// assndRegIndices - to keep track of universe vector belonging to cluster
vector<vector<int>> assndRegIndices;
// CBSize - codebook size
int CBSize = 8;  

//tokhura weights
double tokhuraWt[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//function to read file named as fileName
bool read_csv(string fileName, vector<vector<long double>> &universe){
	fstream fin;
	fin.open(fileName);    

    //file does not exist
	if(!fin){
		fin.close();
		return false;
	}
	long double word;
	char delim;
    //until input is available
	while(fin >> word){
		vector<long double> line;
		line.push_back(word);

        //save whole row
		for(int i=1;i<P;i++){
			fin>>delim>>word;
			line.push_back(word);
		}
        //add row to universe vector
		universe.push_back(line);
	}

	fin.close();
	return true;
}

// calculating distance between the current region and current row of universe
long double calcTokhuraDistance(vector<long double> &currentRegion, vector<long double> &currentUniverse){
	long double ans = 0;
	//loop over each entry
	for(int i=0;i<P;i++){
		long double d = currentRegion[i]-currentUniverse[i];
		//using tokhura distance formula
		ans += (tokhuraWt[i] * d * d);  
	}
	return ans;
}

// function to assign clusters to vectors of universe, it uses tokhura distance
void assignClusters(vector<vector<long double>> &universe, vector<vector<long double>> &clusters){
	int M = universe.size();
	
	//loop over universe vectors
    for(int i=0;i<M;i++){
		int index = 0;
		vector<long double> current_universe = universe[i];   
		long double minDistance = LDBL_MAX;
		for(int j=0;j<CBSize;j++){
			vector<long double> current_region = clusters[j];
            //calculating the distance between current cluster and current row of universe
            long double distance = calcTokhuraDistance(current_region, current_universe);
			
            // save the nearest cluster's index
            if(distance < minDistance){       
				minDistance = distance;
				index = j;
			}
		}
        //add the current row index vector of universe to assigned cluster's vector 
		assndRegIndices[index].push_back(i);
	}
}

// calculating the avg total distortion present in given clusters and given universe
long double totalDistortion(vector<vector<long double>> &universe, vector<vector<long double>> &clusters){
	long double M = 0;
    long double total_disto = 0;
	//loop to all the clusters
	for(int i=0;i<CBSize;i++){
		int size = assndRegIndices[i].size();
		M += (long double)size;
		//loop to all the vectors of particular region
		for(int j=0;j<size;j++){
			int index = assndRegIndices[i][j];
			// total distortion = summation of tokhura distances between cluster and universe vectors
			total_disto += calcTokhuraDistance(clusters[i], universe[index]);
		}
	}
	// taking average of total distortion
	long double avg = total_disto/M;
	return avg;
}

// function which updates the clusters centroid
void updateClusters(vector<vector<long double>> &universe, vector<vector<long double>> &clusters){
	//loop over all region
	for(int i=0;i<CBSize;i++){
		//take temporary vector to sum up all the universe vector values
		vector<long double> t(P, 0.0);

		// take the size of the cluster
		int size = assndRegIndices[i].size();
		
		//loop over whole cluster
		for(int j=0;j<size;j++){
			//take index of the universe vector
			int index = assndRegIndices[i][j];

			//loop over P
			for(int k=0;k<P;k++){
				//add all the values of universe
				t[k] += universe[index][k];
			}
		}
		//take average for finding centroid
		for(int j=0;j<P;j++){
			clusters[i][j] = t[j]/((long double)size);
		}
	}
}

// lloyds algorithm to find distortion at each stage
void lloydsKmeansAlgo(vector<vector<long double>> &universe, vector<vector<long double>> &clusters){
	//assigning clusters to vectors of universe
	assignClusters(universe, clusters);
	long double current_dist = totalDistortion(universe, clusters);
	int round = 1;
	//storing old distortion
	long double old_dist = 0;  
    printf("Cycle %d has total distortion: %.8lf\n", round, current_dist);

	//loop until abs(current_dist - old_dist) becomes less than delta
	while(abs(current_dist - old_dist) > DELTA){
		//update clusters with new centroid values
		updateClusters(universe, clusters);
		
		//clear the indices saved into assigned clusters
		assndRegIndices.clear();
		assndRegIndices.resize(CBSize);

		//assign clusters for the next rounds
		assignClusters(universe, clusters);

		//save current distortion
		old_dist = current_dist;
		//find new distortion
		current_dist = totalDistortion(universe, clusters);
		
		round++;
        //printing on console
		printf("Cycle %d has total distortion: %.8lf\n", round, current_dist);
	}
}

//driver code
int _tmain(int argc, _TCHAR* argv[]) {
	system("MODE CON COLS=160");
	//setup global
	assndRegIndices.clear();
	assndRegIndices.resize(CBSize); 

	//file to be read
	string fileName = "Universe.csv";        

	//variables to hold the universe values and clusters
	vector<vector<long double>> universe;
	vector<vector<long double>> clusters;

	//read given universe file and save into universe vector
	if(!read_csv(fileName, universe)){
		printf("Failed to open file\n");
		system("pause");
		return 1;
	}

	//take 8 clusters
	int partition = universe.size()/CBSize;
	
	for(int i=0;i<CBSize;i++){
		int index = i*partition;
		clusters.push_back(universe[index]);
	}

	//printing initial codebook vectors
	printf("Initial Codebook - \n");
	if(universe.size() == 0) return 1;
	
	printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	for(int i=0;i<clusters.size();i++){
		for(int j=0;j<clusters[i].size();j++){
			printf("| %10.7lf ", clusters[i][j]);
		}
		printf("|\n");
	}
	printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");


	//applying lloyds
	printf("\nApplying Lloyds Algorithm - \n");
	lloydsKmeansAlgo(universe, clusters);
	printf("\n\n");
	
	//printing final codebook
	printf("\nFinal CodeBook (By Llyod's KMeans Algorithm):\n");

    //printing the codebook on console and saving in the file
	if(clusters.size() == 0) return 1;

	printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	for(int i=0;i<clusters.size();i++){
		for(int j=0;j<clusters[i].size();j++){
			printf("| %10.7lf ", clusters[i][j]);
		}
		printf("|\n");
	}
	printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	system("pause");
	return 0;
}