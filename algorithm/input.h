#ifndef INPUT_H
#define INPUT_H

#include <bits/stdc++.h>
using namespace std;

void input_basic(const string& fileName, int& taskN, int& workerN, int& Umax, int& sumC, int& seqN);
void input_weight(const string& path, vector<vector<double> >& weightArr);
void input_weight(const string& edgeFileName, const string& weightFileName, vector<vector<double> >& weightArr);

/**
0 w 57.437591 57.759672 20.000000 1 4 1.000000
0 t 52.688978 71.900817 4 1.000000
*/
void input_basic(const string& fileName, int& taskN, int& workerN, int& Umax, int& sumC, int& seqN) {
	ifstream fin(fileName.c_str(), ios::in);

	if (!fin.is_open()) {
		printf("Error openning FILE %s.\n", fileName.c_str());
		exit(1);
	}

	fin >> workerN >> taskN >> Umax >> sumC;
	seqN = workerN + taskN;

	fin.close();
}

void input_weight(const string& path, vector<vector<double> >& weightArr) {
	int taskN, workerN, sumC, seqN, Umax;

	// get basic parameter from any order file 
	string dataFileName = path + "/" + "order0.txt";
	input_basic(dataFileName, taskN, workerN, Umax, sumC, seqN);

	// get weight array
	string fileName = path + "/" + "weight.txt";
	ifstream fin(fileName.c_str(), ios::in);
	weightArr.clear();
	vector<double> weightRow(taskN, 0);

	if (!fin.is_open()) {
		printf("Error openning FILE %s.\n", fileName.c_str());
		exit(1);
	}

	weightArr.clear();
	for (int i=0; i<workerN; ++i) {
		for (int j=0; j<taskN; ++j) {
			fin >> weightRow[j];
		}
		weightArr.push_back(weightRow);
	}


	fin.close();
}

void input_weight(const string& edgeFileName, const string& weightFileName, vector<vector<double> >& weightArr) {
	int taskN, workerN, sumC, seqN, Umax;

	// get basic parameter from any order file 
	input_basic(edgeFileName, taskN, workerN, Umax, sumC, seqN);

	// get weight array
	string fileName = weightFileName;
	ifstream fin(fileName.c_str(), ios::in);
	weightArr.clear();
	vector<double> weightRow(taskN, 0);

	if (!fin.is_open()) {
		printf("Error openning FILE %s.\n", fileName.c_str());
		exit(1);
	}

	weightArr.clear();
	for (int i=0; i<workerN; ++i) {
		for (int j=0; j<taskN; ++j) {
			fin >> weightRow[j];
		}
		weightArr.push_back(weightRow);
	}


	fin.close();
}


#endif
