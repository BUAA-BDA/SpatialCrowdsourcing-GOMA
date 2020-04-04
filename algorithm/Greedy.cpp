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
int usedMemory;

void init(int taskN, int workerN, double Umax, int SumC) {
	n = workerN;
	m = taskN;
	umax = Umax;
	sumC = SumC;
	utility = 0;
	usedMemory = 0;
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
	//if (!(worker.begTime<=task.endTime && task.begTime<=worker.endTime))
	if (!(worker.begTime<task.endTime && task.begTime<worker.endTime))
		return false;
	return true;
}

bool satisfy(const node_t& worker, const node_t& task) {
	return satisfyCap(worker, task) && satisfyTime(worker, task) && satisfyLoc(worker, task);
	//return satisfyCap(worker, task) && satisfyLoc(worker, task);
}

inline double calcCost(const node_t& task, const node_t& worker) {
	return task.cost.rate * worker.cost.pay;
}

int chosenNextTask(const vector<node_t>& tasks, const node_t& worker) {
	int taskN = tasks.size();
	double tmpCost;
	double mxCost = 0.0;
	int ret = -1;

	for (int i=0; i<taskN; ++i) {
		tmpCost = calcCost(tasks[i], worker);
		if (satisfy(worker, tasks[i]) && tmpCost>mxCost) {
			mxCost = tmpCost;
			ret = i;
		}
	}
	return ret;
}

int chosenNextWorker(const vector<node_t>& workers, const node_t& task) {
	int workerN = workers.size();
	double tmpCost;
	double mxCost = 0.0;
	int ret = -1;

	for (int i=0; i<workerN; ++i) {
		tmpCost = calcCost(task, workers[i]);
		if (satisfy(workers[i], task) && tmpCost>mxCost) {
			mxCost = tmpCost;
			ret = i;
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

void Pure_Greedy(ifstream& fin, int seqN) {
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
				workerId = chosenNextWorker(workers, node);
			} else {
				workerId = workers.size();
				workers.push_back(node);
				taskId = chosenNextTask(tasks, node);
			}
			
			if (workerId>=0 && taskId>=0) {
				addOneMatch(tasks[taskId], workers[workerId]);
			}
		}
	}
	
#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif
}

void solve(string fileName) {
	int taskN, workerN, seqN, sumC;
	double Umax;
	ifstream fin(fileName.c_str(), ios::in);

	if (!fin.is_open()) {
		printf("Error openning FILE %s.\n", fileName.c_str());
		exit(1);
	}

	fin >> workerN >> taskN >> Umax >> sumC;
	seqN = taskN + workerN;
	init(taskN, workerN, Umax, sumC);
	Pure_Greedy(fin, seqN);
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
#ifdef WATCH_MEM
	printf("Greedy %.6lf %.6lf %d\n", utility, usedTime, usedMemory/1024);
#else
	printf("Greedy %.6lf %.6lf\n", utility, usedTime);
#endif
	fflush(stdout);
	
	return 0;
}
