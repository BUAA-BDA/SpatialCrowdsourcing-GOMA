/**
	\author: Trasier
	\date: 2017.04.02
*/
#include <bits/stdc++.h>
using namespace std;
//#pragma comment(linker,"/STACK:102400000,1024000")
#include "monitor.h"
#include "input.h"


enum rule_t {
	worker, task
};

union W_un {
	double cost;
	double pay;
	double rate;
};
	
struct node_t {
	rule_t type;		// 0: task, 1: worker
	pair<double, double> loc;	// location
	int	cap;			// capacity
	int flow;			// flow
	double rad;			// radius	
	W_un cost;			// cost
	int begTime, endTime;	// time interval

	void print() {
		if (type == worker)
		 	printf("w: loc = (%.2lf, %.2lf), rad = %.2lf, cap = %d, time = (%d, %d), ratio = %.2lf\n",
		 			loc.first, loc.second, rad, cap, begTime, endTime, cost.rate);
		else
			printf("t: loc = (%.2lf, %.2lf), time = (%d, %d), pay = %.2lf\n",
		 			loc.first, loc.second, begTime, endTime, cost.pay);
	}
};

typedef long long LL;
int n, m, sumC;
double umax;
double utility;
int usedMemory = 0;

void init(int taskN, int workerN, double Umax, int sumc) {
	n = workerN;
	m = taskN;
	umax = Umax;
	utility = 0;
	usedMemory = 0;
	sumC = sumc;
}

void nextSeq(ifstream& fin, node_t& nd) {
	int timeId;
	string stype;

	fin >> nd.begTime >> stype;
	if (stype[0] == 'w') {
		nd.type = worker;
		fin >> nd.loc.first >> nd.loc.second >> nd.rad >> nd.cap >> nd.endTime >> nd.cost.pay;
		nd.endTime += nd.begTime;
	} else {
		nd.type = task;
		fin >> nd.loc.first >> nd.loc.second >> nd.endTime >> nd.cost.rate;
		nd.endTime += nd.begTime;
		nd.cap = 1;
	}

	nd.flow = 0;
}

inline double calcCost(const node_t& task, const node_t& worker) {
	return task.cost.rate * worker.cost.pay;
}

inline double Length(pair<double,double> pa, pair<double,double> pb) {
	return sqrt( (pa.first-pb.first)*(pa.first-pb.first) + (pa.second-pb.second)*(pa.second-pb.second) );
}

inline double Length2(pair<double,double> pa, pair<double,double> pb) {
	return (pa.first-pb.first)*(pa.first-pb.first) + (pa.second-pb.second)*(pa.second-pb.second);
}


inline bool satisfyLoc(const node_t& worker, const node_t& task) {
	// 4. condition of location
	if (Length2(worker.loc, task.loc) > worker.rad * worker.rad)
		return false;
	
	return true;
}

inline bool satisfyCap(const node_t& worker, const node_t& task) {
	// 2&3. capacity of worker & task
	if (worker.cap<=worker.flow || task.cap<=task.flow)
		return false;
	return true;
}

inline bool satisfyTime(const node_t& worker, const node_t& task) {
	// 1. condition of deadline
	if (!(worker.begTime<task.endTime && task.begTime<worker.endTime))
		return false;
	return true;
}

bool satisfy(const node_t& worker, const node_t& task) {
	return satisfyCap(worker, task) && satisfyTime(worker, task) && satisfyLoc(worker, task);
}

int chosenNextTask(const vector<node_t>& tasks, const node_t& worker, double costBound) {
	int taskN = tasks.size();
	double tmpCost;
	double mxCost = -1e8;
	int ret = -1;

	for (int i=0; i<taskN; ++i) {
		tmpCost = calcCost(tasks[i], worker);
		if (tmpCost>costBound && satisfy(worker, tasks[i])) {
			if (tmpCost > mxCost) {
				mxCost = tmpCost;
				ret = i;
			}
		}
	}
	return ret;
}

int chosenNextWorker(const vector<node_t>& workers, const node_t& task, double costBound) {
	int workerN = workers.size();
	double tmpCost;
	double mxCost = -1e8;
	int ret = -1;

	for (int i=0; i<workerN; ++i) {
		tmpCost = calcCost(task, workers[i]);
		if (tmpCost>costBound && satisfy(workers[i], task)) {
			if (tmpCost > mxCost) {
				mxCost = tmpCost;
				ret = i;
			}
		}
	}
	return ret;
}

void addOneMatch(node_t& task, node_t& worker) {
	// add cost to utility
	utility += calcCost(task, worker);
	// update the capacity of task & worker
	++task.flow;
	++worker.flow;
}

void Extend_Greedy_RT(ifstream& fin, int seqN, int k) {
	double costBound = (k==0) ? 0.0 : pow(exp(1.0), k);
	node_t node;
	vector<node_t> tasks, workers;
	int taskId, workerId;
	
	while (seqN--) {
		nextSeq(fin, node);
		int cap = node.cap;
		node.cap = 1;

		while (cap--) {

			if (node.type == task) { // node is task
				taskId = tasks.size();
				tasks.push_back(node);
				workerId = chosenNextWorker(workers, node, costBound);
			} else {
				workerId = workers.size();
				workers.push_back(node);
				taskId = chosenNextTask(tasks, node, costBound);
			}
			
			if (workerId>=0 && taskId>=0) {
				addOneMatch(tasks[taskId], workers[workerId]);
			}

		}
#ifdef WATCH_MEM
		watchSolutionOnce(getpid(), usedMemory);
#endif
	}
}

double calcUtility(const vector<double>& utilities) {
	double ret;
	const int sz = utilities.size();

	ret = 0.0;
	for (int i=0; i<sz; ++i)
		ret += utilities[i];
	ret /= sz;

	return ret;
}

void solve(string fileName) {
	int taskN, workerN, seqN, sumC;
	double Umax;
	{// get Umax to calculate theta
		ifstream fin(fileName.c_str(), ios::in);

		if (!fin.is_open()) {
			printf("Error openning FILE %s.\n", fileName.c_str());
			exit(1);
		}

		fin >> workerN >> taskN >> Umax >> sumC;
		fin.close();
	}
	
	int theta = ceil(log(Umax + 1.0));
	vector<double> utilities;
	
	for (int i=0; i<theta; ++i) {
		int k = (theta==0) ? 0 : i;
		ifstream fin(fileName.c_str(), ios::in);

		if (!fin.is_open()) {
			printf("Error openning FILE %s.\n", fileName.c_str());
			exit(1);
		}

		fin >> workerN >> taskN >> Umax >> sumC;
		seqN = taskN + workerN;
		init(taskN, workerN, Umax, sumC);
		Extend_Greedy_RT(fin, seqN, k);
		
		utilities.push_back(utility);
		
		printf("k = %d, utility = %.6lf\n", k, utility);
		
		fin.close();
	}
	
	utility = calcUtility(utilities);
}

int main(int argc, char* argv[]) {
	cin.tie(0);
	ios::sync_with_stdio(false);

	string edgeFileName;
	program_t begProg, endProg;

	if (argc > 1) {
		edgeFileName = string(argv[1]);
	} 

	save_time(begProg);
	solve(edgeFileName);
	save_time(endProg);

	double usedTime = calc_time(begProg, endProg);
	int theta = ceil(log(umax + 1.0));
#ifdef WATCH_MEM
	printf("Ext-GRT %.6lf %.6lf %d\n", utility, usedTime/theta, usedMemory/1024);
#else
	printf("Ext-GRT %.6lf %.6lf\n", utility, usedTime/theta);
#endif
	fflush(stdout);
	
	return 0;
}
