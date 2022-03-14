#ifndef BASIC_FUNCTION_H
#define BASIC_FUNCTION_H

#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include "data_form.h"

using namespace std;
using namespace Eigen;

MatrixXd calculate_from_assignment(vector<MatrixXd> T,vector<MatrixXd> A,vector<MatrixXd> B,bool mirror);
MatrixXd calculate_from_assignment_not_vector(MatrixXd P,MatrixXd A,MatrixXd B,bool permit_mirror);
bool find_P(vector<Answer> ans_multiple,MatrixXd T);
double difference_second_largest_Dijkstra_edge_change_assignment(vector<int> &assignment,vector<vector<double> > M,edge &E,bool &cycle);
double calculate_ans(vector<MatrixXd> A,vector<MatrixXd> B,MatrixXd R,vector<vector<int> > assignments);
vector<vector<double> > calculate_lower_graph(MatrixXd A,MatrixXd B,MatrixXd R,double diameter);
int Matrix_vector(MatrixXd A,MatrixXd B,vector<MatrixXd> &As,vector<MatrixXd> &Bs,vector<int> label_A,vector<int> label_B,vector<vector<int> > &a_assign,vector<vector<int> > &b_assign, MatrixXd &Ag, MatrixXd &Bg);

#endif