#ifndef CUBE_H
#define CUBE_H

#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include "../Eigen/Dense"
#include "data_form.h"

using namespace Eigen;

using namespace std;

typedef struct {
    vector<double> center;
    double diameter;
    double upper;
    double lower;
    int count; 
    bool mirror;
    vector<vector<int> > assignments;
} Cube;

typedef struct {
    vector<double> center_r;
    vector<double> center_t;
    double diameter_r;
    double diameter_t;
    double upper;
    double lower;
    int count; 
    bool mirror;
} Cube_Cube;


bool operator< (const Cube &cube1, const Cube &cube2);
bool operator< (const Cube_Cube &cube1, const Cube_Cube &cube2);

MatrixXd make_rotation(vector<double> r,bool mirror);
MatrixXd make_translation(vector<double> r);
double calculate_upper(Cube C,vector<MatrixXd> A,vector<MatrixXd> B,bool mirror);
vector<Cube> divide(Cube C);
bool check_in_sphere(Cube C);
bool calculate_lower_K_best(Cube C,vector<MatrixXd> As,vector<MatrixXd> Bs,bool mirror,Answer &pre,double &lower,int k,vector<vector<int> > a_assign,vector<vector<int> > b_assign);
bool calculate_lower_K_best_eps(Cube C,vector<MatrixXd> A,vector<MatrixXd> B,bool mirror,vector<Answer> &pre,double eps,int k,vector<vector<int> > a_assign,vector<vector<int> > b_assign,int len_a,int len_b);
Answer rotation_minimum_Best(vector<MatrixXd> As,vector<MatrixXd> Bs,Cube C_t,Answer &pre,bool permit_mirror,int K,int L,vector<vector<int> > a_assign,vector<vector<int> > b_assign);
#endif