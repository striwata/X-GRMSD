#include "../Eigen/Dense"
#include "../Eigen/SVD"
#include <stdlib.h>
#include "basic_function.h"
#include "Hungarian.h"
#include "data_form.h"
#include <set>
#define EPS 1e-10
#define PI 3.14159265359
using namespace Eigen;


// Judge whether a permutation matrix has been counted or not
bool find_P(vector<Answer> ans_multiple,MatrixXd P){
    int a = P.cols();
    int b = P.rows();
    for(int i=0;i<ans_multiple.size();i++){
        MatrixXd P0 = ans_multiple[i].P;
        bool find = true;
        for(int i=0;i<a;i++){
            bool found = false;
            for(int j=0;j<b;j++){
                if(P0(j,i) == 1){
                    for(int k=0;k<a;k++){
                        if(P(j,k) == 1){
                            found = true;
                            break;
                        }
                    }
                    break;
                }
            }
            if(!found){
                find = false;
                break;
            }
        }
        if(find){
            return true;
        }
        
    }
    return false;
}

// Caluculate the best rotation matrix with a fixed assignment
MatrixXd calculate_from_assignment(vector<MatrixXd> Ps,vector<MatrixXd> As,vector<MatrixXd> Bs,bool mirror){
    int N = Ps.size();
    MatrixXd BPA = MatrixXd::Zero(3,3);
    for(int i=0;i<N;i++){
        BPA += Bs[i]*Ps[i]*As[i].transpose();
    }
    MatrixXd umeyama = MatrixXd::Identity(3,3);
    umeyama(2,2) = -1;
    MatrixXd R;
    JacobiSVD< Matrix<double, 3, 3> > svd(BPA, Eigen::ComputeFullU |Eigen::ComputeFullV);
    if ((svd.matrixV().determinant()*svd.matrixU().determinant()>0) ^ (mirror)){
        R = svd.matrixV()*svd.matrixU().transpose();
    }else{
        R = svd.matrixV()*umeyama*svd.matrixU().transpose();
    }
    return R;
}

MatrixXd calculate_from_assignment_not_vector(MatrixXd P,MatrixXd A,MatrixXd B,bool permit_mirror){
    MatrixXd BPA = B*P*A.transpose();
    MatrixXd umeyama = MatrixXd::Identity(3,3);
    umeyama(2,2) = -1;
    MatrixXd R;
    JacobiSVD< Matrix<double, 3, 3> > svd(BPA, Eigen::ComputeFullU |Eigen::ComputeFullV);
    if ((svd.matrixV().determinant()*svd.matrixU().determinant()>0) || (permit_mirror)){
        R = svd.matrixV()*svd.matrixU().transpose();
    }else{
        R = svd.matrixV()*umeyama*svd.matrixU().transpose();
    }
    return R;
}

// Calculate the second largest matching
double difference_second_largest_Dijkstra_edge_change_assignment(vector<int> &assignment,vector<vector<double> > M,edge &E,bool &cycle){
    double inf = INFINITY;
    int N = assignment.size();
    int NN = M[0].size();
    vector<double> d(N+NN);
    vector<int> prev(N+NN);
    vector<double> ans(N);
    vector<int> assignment_inv(NN);
    for(int i=0;i<N;i++){
        if(assignment[i] != -1){
            assignment_inv[assignment[i]] = i;
        }
    }
    vector<bool> have_arrived(N+NN);
    double rep = inf;
    vector<int> path(N+NN);
    int change = 0;
    for(int i=0;i<N+NN;i++){
        path[i] = -1;
    }
    for(int i=0;i<N;i++){
        if(assignment[i] == -1){
            ans[i] = inf;
            int k;
            for(int j=0;j<NN;j++){
                if(ans[i] > M[i][j]){
                    ans[i] = M[i][j];
                    k = j;
                }
            }
            if(ans[i] < rep){
                cycle = false;
                rep = ans[i];
                E.from = assignment_inv[k];
                E.to = k;
            }
            continue;
        }
        for(int j=0;j<N+NN;j++){
            d[j] = inf;
            have_arrived[j] = false;
            prev[j] = -1;
        }
        d[i] = 0;
        while(!have_arrived[assignment[i]+N]){
            int v = -1;
            int k;
            for(int u=0;u<N+NN;u++){
                if(!have_arrived[u] && (v == -1 || d[u]<d[v])){
                    v = u;
                }
            }
            have_arrived[v] = true;
            if(v<N){
                for(int u=0;u<NN;u++){
                    if(u!=assignment[v]){
                        if(d[v] - d[u+N]+EPS<-M[v][u]){
                            d[u+N] = d[v] + M[v][u];
                            prev[u+N] = v;
                        }
                    }
                }
            }else{
                int u = assignment_inv[v-N];
                if(d[v] - M[u][v-N]-d[u]<-EPS){
                    d[u] = d[v] - M[u][v-N];
                    prev[u] = v;
                }
            }
        }
        ans[i] = d[assignment[i]+N]- M[i][assignment[i]];
        if(ans[i] < rep){
            cycle = true;
            rep = ans[i];
            E.from = i;
            E.to = assignment[i];
            for(int j=0;j<N+NN;j++){
                path[j] = prev[j];
            }
        }
    }
    if(rep == inf){
        return rep;
    }
    int t = E.to;
    if(!cycle){
        return rep;
    }
    while(true){
        int s = path[t+N];
        assignment_inv[t] = s;
        if(s==E.from){
            break;
        }
        t = path[s]-N;
    }
    for(int i=0;i<N;i++){
        bool found = false;
        for(int j=0;j<NN;j++){
            if(i == assignment_inv[j]){
                assignment[i] = j;
                found = true;
                break;
            }
        }
        if(!found){
            assignment[i] = -1;
        }
    }
    return rep;
}

// Calculate RMSD with a fixed rotation and assignment
double calculate_ans(vector<MatrixXd> As,vector<MatrixXd> Bs,MatrixXd R,vector<vector<int> > assignments){
    double ans = 0;
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        MatrixXd RB = R*Bs[iter];
        for(int i = 0;i<b;i++){
            if(assignments[iter][i]!=-1){
                for(int k=0;k<3;k++){
                    ans  += (As[iter](k,assignments[iter][i])-RB(k,i))*(As[iter](k,assignments[iter][i])-RB(k,i));
                }
            }
        }
    }
    return ans;
}

// Calculate lower bound of graph
vector<vector<double> > calculate_lower_graph(MatrixXd A,MatrixXd B,MatrixXd R,double diameter){
    int a = A.cols();
    int b = B.cols();
    vector<double> gamma(b);
    for(int i=0;i<b;i++){
        gamma[i] = 2*sin(min(sqrt(3)*diameter/2,PI/2))*B.col(i).norm();
    }
    MatrixXd RB = R*B;
    vector<vector<double> >M(b, vector<double>(a));
    for(int i = 0;i<b;i++){
        for(int j=0;j<a;j++){
            M[i][j] = 0;
            for(int k=0;k<3;k++){
                M[i][j]  += (A(k,j)-RB(k,i))*(A(k,j)-RB(k,i));
            }
            M[i][j] = sqrt(M[i][j]);
            M[i][j] = M[i][j]- gamma[i];
            M[i][j] = max(M[i][j],(double)0)*max(M[i][j],(double)0);
        }
    }
    return M;
}

// Turn two matrices into vectors of matrices
int Matrix_vector(MatrixXd A,MatrixXd B,vector<MatrixXd> &As,vector<MatrixXd> &Bs,vector<int> label_A,vector<int> label_B,vector<vector<int> > &a_assign,vector<vector<int> > &b_assign, MatrixXd &Ag, MatrixXd &Bg){
    int a = A.cols();
    int b = B.cols();	
	for(int i=0;i<3;i++){
		for(int j=0;j<a;j++){
			Ag(i,0) = Ag(i,0) + A(i,j);
		}
		for(int j=0;j<b;j++){
			Bg(i,0) = Bg(i,0) + B(i,j);
		}
		Ag(i,0) = Ag(i,0)/a;
		Bg(i,0) = Bg(i,0)/b;
		for(int j=0;j<a;j++){
			A(i,j) = A(i,j) - Ag(i,0);
		}
		for(int j=0;j<b;j++){
			B(i,j) = B(i,j) - Bg(i,0);
		}
	}
    set<int> labels;
    for(int i=0;i<label_B.size();i++){
        if(labels.find(label_B[i])==labels.end()){
            vector<int> a_assign_temp;
            vector<int> b_assign_temp;
            for(int j_b = i;j_b<label_B.size();j_b++){
                if(label_B[i] == label_B[j_b]){
                    b_assign_temp.push_back(j_b);
                }
            }
            for(int j_a = 0;j_a<label_A.size();j_a++){
                if(label_B[i] == label_A[j_a]){
                    a_assign_temp.push_back(j_a);
                }
            }
            labels.insert(label_B[i]);
            a_assign.push_back(a_assign_temp);
            b_assign.push_back(b_assign_temp);
        }
    }
    for(int i=0;i<b_assign.size();i++){
        int a = a_assign[i].size();
        int b = b_assign[i].size();
        if(a!=0){
            MatrixXd AA = MatrixXd::Zero(3,a);
            for(int i_a=0;i_a<a;i_a++){
                for(int k=0;k<3;k++){
                    AA(k,i_a) =  A(k,a_assign[i][i_a]);
                }
            } 
            As.push_back(AA);   
            MatrixXd BB = MatrixXd::Zero(3,b);
            for(int i_b=0;i_b<b;i_b++){
                for(int k=0;k<3;k++){
                    BB(k,i_b) =  B(k,b_assign[i][i_b]);
                }
            } 
            Bs.push_back(BB);
        }
    }
    return 0;
}
