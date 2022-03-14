#include <iostream>
#include "Eigen/Dense"
#include "data_form.h"


#ifndef MAIN_FUNCTION_H
#define MAIN_FUNCTION_H

#include <iostream>
#include <vector>

Answer AO_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,vector<MatrixXd> Smentai,bool permit_mirror);
Answer TSR_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,vector<MatrixXd> Smentai,bool permit_mirror,double alpha);
Answer AO_diff_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,vector<MatrixXd> Smentai,bool permit_mirror);
Answer IsometryOpt_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror,double k ,double l);
Answer IsometryOpt_diff_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror,double k ,double l);
Answer MatchFastOpt_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror);
Answer MatchFastOpt_diff_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror);
vector<Answer> IsometrySearch(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror);
vector<Answer> MatchFPT(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror);
vector<Answer> Improved_MatchFPT(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror);

#endif