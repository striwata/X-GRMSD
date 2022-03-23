#include <iostream>
#include "../Eigen/Dense"
#include "data_form.h"


#ifndef LOCALICP_H
#define LOCALICP_H

#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

Answer local_AO_same_size(vector<MatrixXd> As,vector<MatrixXd> Bs,MatrixXd initial_matrix,bool mirror,vector<vector<int> > a_assign,vector<vector<int> > b_assign);
Answer local_AO_diff_size(vector<MatrixXd> As, vector<MatrixXd> Bs,MatrixXd initial_rotation, bool mirror,MatrixXd initial_translation,vector<vector<int> > a_assign,vector<vector<int> > b_assign);
Answer local_TSR_same_size(vector<MatrixXd> As, vector<MatrixXd> Bs,MatrixXd initial_rotation, bool mirror,double alpha,vector<vector<int> > a_assign,vector<vector<int> > b_assign);
#endif