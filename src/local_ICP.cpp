#include "../Eigen/Dense"
#include <stdlib.h>
#include "local_ICP.h"
#include "Hungarian.h"
#include "basic_function.h"
#include "data_form.h"
// using namespace Eigen;

Answer local_AO_same_size(vector<MatrixXd> A,vector<MatrixXd> B,MatrixXd initial_matrix,bool mirror,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    HungarianAlgorithm HungAlgo;
    int N = A.size();
    double cost;
    MatrixXd R = initial_matrix;
    vector<MatrixXd> Ps;
    int len_a = 0;
    int len_b = 0;
    for(int iter = 0;iter<10000;iter++){
        vector<MatrixXd> temp_Ps;
        cost = 0;
        for(int i=0;i<N;i++){
            int a = A[i].cols();
            int b = B[i].cols();
            if(iter == 0){
                len_a += a;
                len_b += b;
            }
            if(a!=b){
                cout << "error" << endl;
            }
            MatrixXd P = MatrixXd::Zero(b,a);
            MatrixXd RB = R*B[i];
            vector<vector<double> >C(b, vector<double>(a));
            for(int i_b = 0;i_b<b;i_b++){
                for(int i_a=0;i_a<a;i_a++){
                    C[i_b][i_a] = 0;
                    for(int k=0;k<3;k++){
                        C[i_b][i_a]  += (A[i](k,i_a)-RB(k,i_b))*(A[i](k,i_a)-RB(k,i_b));
                    }
                }
            }
            vector<int> assignment;
            cost += HungAlgo.Solve(C, assignment);
            for(int i_t=0;i_t<a;i_t++){
                if(assignment[i_t]!=-1){
                    P(i_t,assignment[i_t]) = 1;
                }
            }
            temp_Ps.push_back(P);
        }
        R = calculate_from_assignment(temp_Ps,A,B,mirror);
        if(temp_Ps == Ps){
            break;
        }
        Ps = temp_Ps;
    }
    Answer ans;
    ans.E = cost;
    ans.R =R;
    MatrixXd P = MatrixXd::Zero(len_b,len_a);
    for(int i=0;i<N;i++){
        int a = A[i].cols();
        int b = B[i].cols();
        for(int i_a=0;i_a<a;i_a++){
            for(int i_b=0;i_b<b;i_b++){
                if(Ps[i](i_b,i_a) == 1){
                    P(b_assign[i][i_b],a_assign[i][i_a]) = 1;
                }
            }
        }
    }
    ans.P = P;
    return ans;
}

Answer local_TSR_same_size(vector<MatrixXd> As,vector<MatrixXd> Bs,MatrixXd initial_matrix,bool mirror,double alpha,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    HungarianAlgorithm HungAlgo;
    int N = As.size();
    double cost;
    MatrixXd R = initial_matrix;
    vector<MatrixXd> Ps;
    int len_a = 0;
    int len_b = 0;
    for(int iter = 0;iter<10000;iter++){
        vector<MatrixXd> temp_Ps;
        vector<MatrixXd> temp_Pss;
        cost = 0;
        for(int i=0;i<N;i++){
            int a = As[i].cols();
            int b = Bs[i].cols();
            if(iter == 0){
                len_a += a;
                len_b += b;
            }
            if(a!=b){
                cout << "error" << endl;
            }
            MatrixXd P = MatrixXd::Zero(b,a);
            MatrixXd RB = R*Bs[i];
            vector<vector<double> >C(b, vector<double>(a));
            for(int i_b = 0;i_b<b;i_b++){
                for(int i_a=0;i_a<a;i_a++){
                    C[i_b][i_a] = 0;
                    for(int k=0;k<3;k++){
                        C[i_b][i_a]  += (As[i](k,i_a)-RB(k,i_b))*(As[i](k,i_a)-RB(k,i_b));
                    }
                }
            }
            vector<int> assignment;
            cost += HungAlgo.Solve(C, assignment);
            for(int i_t=0;i_t<a;i_t++){
                if(assignment[i_t]!=-1){
                    P(i_t,assignment[i_t]) = 1;
                }
            }
            MatrixXd c = P*As[i].transpose()*R*Bs[i];
            MatrixXd S = MatrixXd::Zero(b,a);
            for(int i_a=0;i_a<a;i_a++){
                for(int i_b=0;i_b<b;i_b++){
                    S(i_a,i_b) = c(i_a,i_b)-c(i_b,i_a);
                }
            }
            double frobenius = 0;
            for(int i_a=0;i_a<a;i_a++){
                for(int i_b=0;i_b<b;i_b++){
                    frobenius += S(i_a,i_b)*S(i_a,i_b);
                }
            }
            frobenius = sqrt(frobenius);
            for(int i_a=0;i_a<a;i_a++){
                for(int i_b=0;i_b<b;i_b++){
                    S(i_a,i_b) = S(i_a,i_b)/frobenius;
                }
            }
            if(frobenius == 0){
                temp_Pss.push_back(P);
                temp_Ps.push_back(P);
            }
            if(frobenius != 0){
                temp_Pss.push_back(P);
                MatrixXd PP = P + alpha * S*P;
                temp_Ps.push_back(PP);
            }
            
        }
        R = calculate_from_assignment(temp_Ps,As,Bs,mirror);
        if(Ps == temp_Pss && temp_Pss.size() != 0){
            cost = 0;
            R = calculate_from_assignment(temp_Pss,As,Bs,mirror);
            for(int i=0;i<N;i++){
                int a = As[i].cols();
                int b = Bs[i].cols();
                if(a!=b){
                    cout << "error" << endl;
                }
                MatrixXd T = MatrixXd::Zero(b,a);
                MatrixXd RB = R*Bs[i];
                vector<vector<double> >C(b, vector<double>(a));
                for(int i_b = 0;i_b<b;i_b++){
                    for(int i_a=0;i_a<a;i_a++){
                        C[i_b][i_a] = 0;
                        for(int k=0;k<3;k++){
                            C[i_b][i_a]  += (As[i](k,i_a)-RB(k,i_b))*(As[i](k,i_a)-RB(k,i_b));
                        }
                    }
                }
                vector<int> assignment;
                cost += HungAlgo.Solve(C, assignment);
            }
            break;
        }
        Ps = temp_Pss;
    }
    Answer ans;
    ans.E = cost;
    ans.R = R;
    MatrixXd P = MatrixXd::Zero(len_b,len_a);
    for(int i=0;i<N;i++){
        int a = As[i].cols();
        int b = Bs[i].cols();
        for(int i_a=0;i_a<a;i_a++){
            for(int i_b=0;i_b<b;i_b++){
                if(Ps[i](i_b,i_a) == 1){
                    P(b_assign[i][i_b],a_assign[i][i_a]) = 1;
                }
            }
        }
    }
    ans.P = P;
    return ans;
}

Answer local_AO_diff_size(vector<MatrixXd> As, vector<MatrixXd> Bs,MatrixXd initial_rotation, bool mirror,MatrixXd initial_translation,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    HungarianAlgorithm HungAlgo;
    int N = As.size();
    double cost = 100;
    double cost_before = 100;
    MatrixXd R = initial_rotation;
    MatrixXd t = initial_translation;
    vector<MatrixXd> Ps;
    int len_a = 0;
    int len_b = 0;
    for(int iter = 0;iter<10000;iter++){
        int L = 0;
        vector<MatrixXd> temp_Ps;
        cost = 0;
        vector<MatrixXd> BBs = Bs;
        MatrixXd tt = MatrixXd::Zero(3,1);
        for(int i=0;i<N;i++){
            int a = As[i].cols();
            int b = Bs[i].cols();
            if(iter == 0){
                len_a += a;
                len_b += b;
            }
            MatrixXd P = MatrixXd::Zero(b,a);
            MatrixXd RB = R*Bs[i];
            vector<vector<double> >C(b, vector<double>(a));
            for(int i_b = 0;i_b<b;i_b++){
                for(int i_a=0;i_a<a;i_a++){
                    C[i_b][i_a] = 0;
                    for(int k=0;k<3;k++){
                        C[i_b][i_a]  += (As[i](k,i_a)-RB(k,i_b)-t(k,0))*(As[i](k,i_a)-RB(k,i_b)-t(k,0));
                    }
                }
            }
            vector<int> assignment;
            cost += HungAlgo.Solve(C, assignment);
            for(int i_t=0;i_t<b;i_t++){
                if(assignment[i_t]!=-1){
                    P(i_t,assignment[i_t]) = 1;
                }
            }

            MatrixXd BTs = Bs[i] * P;
            for(int k=0;k<3;k++){
                for(int j=0;j<a;j++){
                    tt(k,0) += BTs(k,j);
                }
            }
            
            temp_Ps.push_back(P);
        }
        cost_before = cost;
        for(int k=0;k<3;k++){
            tt(k,0) = -tt(k,0)/len_a;
            for(int i=0;i<N;i++){
                int b = BBs[i].cols();
                for(int j = 0;j<b;j++){
                    BBs[i](k,j) += tt(k,0);
                }
            }
        }
        R = calculate_from_assignment(temp_Ps,As,BBs,mirror);
        t = R*tt;
        if(temp_Ps == Ps){
            break;
        }
        Ps = temp_Ps;
    }
    Answer ans;
    ans.E = cost;
    ans.R =R;
    ans.t = t;
    MatrixXd P = MatrixXd::Zero(len_b,len_a);
    for(int i=0;i<N;i++){
        int a = As[i].cols();
        int b = Bs[i].cols();
        for(int i_a=0;i_a<a;i_a++){
            for(int i_b=0;i_b<b;i_b++){
                if(Ps[i](i_b,i_a) == 1){
                    P(b_assign[i][i_b],a_assign[i][i_a]) = 1;
                }
            }
        }
    }
    ans.P = P;
    return ans;

}
