#ifndef DATA_FORM_H
#define DATA_FORM_H

#include <iostream>
#include "Eigen/Dense"
#include <vector>

using namespace std;
using namespace Eigen;

typedef struct {
    double E;
    MatrixXd R;
    int num;
    MatrixXd t;
    bool mirror;
    MatrixXd P;
} Answer;

typedef struct {
    vector<int> assignments;
    double lower;
} Lower_candidate;

struct edge {int to,from;};

typedef struct {
    vector<edge> I;
    vector<edge> O;
    double second_largest;
    edge E;
    vector<int>  assignment;
    double cost;
} solution_space;

typedef struct {
    double cost;
    vector<int> froms;
} cost_from;

typedef struct {
    vector<int> assignment;
    vector<int> not_assign;
    vector<vector<int> > assign_candidate;
    double E;
} permutation_space;


bool operator< (const solution_space &sol1, const solution_space &sol2);

bool operator< (const cost_from &cost_from1, const cost_from &cost_from2);

bool operator< (const permutation_space &permutation_space1, const permutation_space &permutation_space2);



#endif