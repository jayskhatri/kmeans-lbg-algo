// 214101023_LBG.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"       
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

#define P 12
#define DELTA 0.0001 //delta value
#define EPS 0.03 //epsilon value

// global parameters
// assndRegIndices - to keep track of universe vector belonging to cluster
vector<vector<int>> assndRegIndices;
// CBSize - codebook size
int CBSize = 1;  

//tokhura weights
double tokhuraWt[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//function to read file named as fileName
bool read_csv(vector<vector<long double>> &universe, string fileName){
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
	for(int i=0;i<P;i++){
		long double d = currentRegion[i]-currentUniverse[i];
		ans += (tokhuraWt[i] * d * d);  
	}
	return ans;
}

// function to assign clusters to vectors of universe, it uses tokhura distance
void assignClusters(vector<vector<long double>> &universe, vector<vector<long double>> &clusters){
	int M = universe.size();
	//looping over all vectors of universe
    for(int i=0;i<M;i++){
		int index = 0;
		vector<long double> current_universe = universe[i];   
		long double minDistance = LDBL_MAX;
		//looping over all clusters
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
        //add the current row index vector of universe to assigned region vector 
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
	assndRegIndices.clear();
	assndRegIndices.resize(CBSize);
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

vector<long double> findCentroidOfUniverse(vector<vector<long double>> &universe){
    vector<long double> t(P, 0.0);

    int M = universe.size();
	//looping over all the vectors of universe
    for(int i=0; i<M; i++){
		//calculating the centroid of each entry of the universe
        for(int j = 0; j<P; j++){
            t[j] += universe[i][j];
        }
    }
	//taking average value by dividing the total number of vectors of universe
    for(int k = 0; k<P; k++){
        t[k] /= M;
    }
    return t;
}

//function to update codebook
void updateCodebook(vector<vector<long double>>&clusters, vector<vector<long double>> &copyCB){
	for(int i=0; i<copyCB.size(); i++){
        vector<long double> t1, t2;
		//calculating the new codebook entry
        for(int j=0; j<P; j++){
			long double v1 = ((1 + EPS) * copyCB[i][j]);
            t1.push_back(v1);
			long double v2 = ((1 - EPS) * copyCB[i][j]);
            t2.push_back(v2);
        }
        clusters.push_back(t1);
        clusters.push_back(t2);
    }
}

void LBGAlgo(vector<vector<long double>> &universe, vector<vector<long double>> &clusters){
    int round = 1;

    while(CBSize <= 8){
        printf("Codebook of size: %d\n", CBSize);

        lloydsKmeansAlgo(universe, clusters);
        printf("Codebook generated by LBG - \n");

        //printing codebook
        if(clusters.size() == 0){
            return;
        }
        printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        for(int i=0;i<clusters.size();i++){
            for(int j=0;j<clusters[i].size();j++){
                printf("| %10.7lf ", clusters[i][j]);
            }
            printf("|\n");
        }
        printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

        round++;
        if(CBSize >= 8) return;

        vector<vector<long double>> copyCodeBook(clusters);
        clusters.clear();
        CBSize *=2;
        updateCodebook(clusters, copyCodeBook);
    }
}

int _tmain(int argc, _TCHAR* argv[]){
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
	if(!read_csv(universe, fileName)){
		printf("Failed to open file\n");
		system("pause");
		return 1;
	}

    //start with one vector i.e. centroid of whole universe
    clusters.push_back(findCentroidOfUniverse(universe));

    //applying LBG algo to find vectors
    LBGAlgo(universe, clusters);
    system("pause");
	return 0;
}

