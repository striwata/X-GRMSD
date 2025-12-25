#include <iostream>
#include "../Eigen/Dense"
#include "cube.h"
#include <queue>
#include <set>
#include "Hungarian.h"
#include "local_ICP.h"
#include "basic_function.h"
#include "time.h"
#define PI 3.141592653590

using namespace Eigen;

Answer AO_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,vector<MatrixXd> Smentai,bool permit_mirror){

    int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;
    
    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);

    Answer ans;
	ans.E = INFINITY;

    for(int i = 0;i<Smentai.size();i++){
		if(Smentai[i].determinant() > 0){
			Answer pre = local_AO_same_size(As,Bs,Smentai[i],0,a_assign,b_assign);
			if(pre.E < ans.E){
				ans = pre;
			}
		}else{
			if(!permit_mirror){
				continue;
			}else{
				Answer pre = local_AO_same_size(As,Bs,Smentai[i],1,a_assign,b_assign);
				if(pre.E < ans.E){
				    ans = pre;
			    }
			}
		}
	}
    ans.E = sqrt(ans.E/b);
	ans.t = Ag - ans.R*Bg;
	return ans;
}

Answer TSR_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,vector<MatrixXd> Smentai,bool permit_mirror,double alpha){

    int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;
    
    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);

    Answer ans;
	ans.E = INFINITY;

    for(int i = 0;i<Smentai.size();i++){
		if(Smentai[i].determinant() > 0){
			Answer pre = local_TSR_same_size(As,Bs,Smentai[i],0,alpha,a_assign,b_assign);
			if(pre.E < ans.E){
				ans = pre;
			}
		}else{
			if(!permit_mirror){
				continue;
			}else{
				Answer pre = local_TSR_same_size(As,Bs,Smentai[i],1,alpha,a_assign,b_assign);
				if(pre.E < ans.E){
				    ans = pre;
			    }
			}
		}
	}
    ans.E = sqrt(ans.E/b);
	ans.t = Ag - ans.R*Bg;
	return ans;
}


Answer AO_diff_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,vector<MatrixXd> Smentai,bool permit_mirror){
    int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;

    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);

    Answer ans;
	ans.E = INFINITY;
    for(int i_a = 0;i_a<a;i_a++){
		for(int i_b = 0;i_b<b;i_b++){
			MatrixXd t = MatrixXd::Zero(3,1);
			for(int k=0;k<3;k++){
				t(k,0) += Bs[0](k,i_b);
				t(k,0) -= As[0](k,i_a);
			}	
            for(int i = 0;i<120;i++){
                MatrixXd Rt = Smentai[i] * t;
                if(Smentai[i].determinant() > 0){
                    Answer pre = local_AO_diff_size(As,Bs,Smentai[i],0,Rt,a_assign,b_assign);
                    if(pre.E < ans.E){
                        ans = pre;
                    }
                }else{
                    if(!permit_mirror){
                        continue;
                    }else{
                        Answer pre = local_AO_diff_size(As,Bs,Smentai[i],1,Rt,a_assign,b_assign);
                        if(pre.E < ans.E){
                            ans = pre;
                        }
                    }
                }

            }
        }
    }
    ans.E = sqrt(ans.E/a);
	ans.t = ans.t + Ag - ans.R*Bg;
	return ans;

}

Answer IsometryOpt_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror,double k ,double l){
	int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;

    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);

    vector<double> r;
	for(int k=0;k<3;k++){
		r.push_back(0);
	}
	Cube Cr;
	Cr.center = r;
	Cr.diameter = PI;
	Cr.lower = 0;
	Cr.upper = INFINITY;
	Cr.mirror = 0;
	Cr.count = 0;
	priority_queue<Cube> que;
	que.push(Cr);
	if(permit_mirror){
		Cube Crr;
		Crr.count = 0;
		Crr.center = r;
		Crr.diameter = PI;
		Crr.lower = 0;
		Crr.upper = INFINITY;
		Crr.mirror = 1;
		que.push(Crr);
	}
	MatrixXd initial_matrix = make_rotation(r,0);
	Answer ans = local_AO_same_size(As,Bs,initial_matrix,0,a_assign,b_assign);
	for(int i=0;i<1000000;i++){
		if(que.size()==0){
			break;
		}
		Cube C = que.top();
		que.pop();
		vector<Cube> Cs = divide(C);
		for(int j=0;j<8;j++){
			Cube CC = Cs[j];
			if(!check_in_sphere(CC)){
				continue;
			}
			double E_upper = calculate_upper(CC,As,Bs,CC.mirror);
			CC.upper = E_upper;
			if(E_upper<ans.E){
				Answer temp_ans = local_AO_same_size(As,Bs,make_rotation(CC.center,CC.mirror),CC.mirror,a_assign,b_assign);
				if (temp_ans.E < ans.E){
					ans = temp_ans;
				}
			}
			int K = int(k*CC.count +l+1-k+0.01);
			bool fin = calculate_lower_K_best(CC,As,Bs,CC.mirror,ans,CC.lower,K,a_assign,b_assign);
			if(fin){
				continue;
			}
			que.push(CC);
		}
	}
    ans.E = sqrt(ans.E/a);
	ans.t = Ag - ans.R*Bg;
	return ans;

}

Answer IsometryOpt_same_size_eps(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror,double k ,double l, double eps){
	int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;

    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);

    vector<double> r;
	for(int k=0;k<3;k++){
		r.push_back(0);
	}
	Cube Cr;
	Cr.center = r;
	Cr.diameter = PI;
	Cr.lower = 0;
	Cr.upper = INFINITY;
	Cr.mirror = 0;
	Cr.count = 0;
	priority_queue<Cube> que;
	que.push(Cr);
	if(permit_mirror){
		Cube Crr;
		Crr.count = 0;
		Crr.center = r;
		Crr.diameter = PI;
		Crr.lower = 0;
		Crr.upper = INFINITY;
		Crr.mirror = 1;
		que.push(Crr);
	}
	MatrixXd initial_matrix = make_rotation(r,0);
	Answer ans = local_AO_same_size(As,Bs,initial_matrix,0,a_assign,b_assign);
	ans.E = eps*eps*a;
	for(int i=0;i<1000000;i++){
		if(que.size()==0){
			break;
		}
		Cube C = que.top();
		que.pop();
		vector<Cube> Cs = divide(C);
		for(int j=0;j<8;j++){
			Cube CC = Cs[j];
			if(!check_in_sphere(CC)){
				continue;
			}
			double E_upper = calculate_upper(CC,As,Bs,CC.mirror);
			CC.upper = E_upper;
			if(E_upper<ans.E){
				Answer temp_ans = local_AO_same_size(As,Bs,make_rotation(CC.center,CC.mirror),CC.mirror,a_assign,b_assign);
				if (temp_ans.E < ans.E){
					ans = temp_ans;
				}
			}
			int K = int(k*CC.count +l+1-k+0.01);
			bool fin = calculate_lower_K_best(CC,As,Bs,CC.mirror,ans,CC.lower,K,a_assign,b_assign);
			if(fin){
				continue;
			}
			que.push(CC);
		}
	}
    ans.E = sqrt(ans.E/a);
	ans.t = Ag - ans.R*Bg;
	cout << "finished" << endl;
	return ans;

}

Answer IsometryOpt_diff_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror,double k ,double l){
	int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;

    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);
	Cube Ct;
	Ct.diameter = 0;
	for(int iter=0;iter<As.size();iter++){
		for(int i=0;i<As[iter].cols();i++){
			double diameter_temp = 0;
            for(int k=0;k<3;k++){
				diameter_temp += As[iter](k,i)*As[iter](k,i);
			}
			diameter_temp = sqrt(diameter_temp);
			if(Ct.diameter < diameter_temp){
				Ct.diameter = diameter_temp;
			}
		}
		
	}
	priority_queue<Cube> que;
	vector<double> r;
	for(int k=0;k<3;k++){
		r.push_back(0);
	}
	Ct.center = r;
	Ct.lower = 0;
	Ct.upper = INFINITY;
	Ct.count = 1;
	que.push(Ct);
	MatrixXd initial_matrix = make_rotation(r,0);
	MatrixXd r_t = MatrixXd::Zero(3,1);
	Answer ans = local_AO_diff_size(As,Bs,initial_matrix,0,r_t,a_assign,b_assign);
	while(!que.empty()){
		Cube C = que.top();
		que.pop();
		vector<Cube> Cs = divide(C);
		for(int j=0;j<8;j++){
			Cube CC = Cs[j];
			Cube C_t;
			C_t.center = CC.center;
			C_t.lower = CC.lower;
			C_t.upper = CC.upper;
			C_t.mirror = CC.mirror;
			C_t.diameter = CC.diameter;
			Answer E_upper = local_AO_diff_size(As,Bs,initial_matrix,CC.mirror,make_translation(CC.center),a_assign,b_assign);
			CC.upper = E_upper.E;
			if(E_upper.E<ans.E){
				ans = E_upper;
			}
            int KK = int(k*CC.count +l+1-k+0.01);
			Answer E_lower = rotation_minimum_Best(As,Bs,CC,ans,permit_mirror,KK,10000000*KK,a_assign,b_assign);
			CC.lower = E_lower.E;
			if(E_lower.E>ans.E){
				continue;
			}
			que.push(CC);
		}
	}

	ans.E = sqrt(ans.E/a);
	ans.t = Ag - ans.R*(Bg + ans.t);
	return ans;
}

Answer MatchFastOpt_same_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror){
	
    int a = A.cols();
    int b = B.cols();

	MatrixXd Ag = MatrixXd::Zero(3,1);
    MatrixXd Bg = MatrixXd::Zero(3,1);

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

	Answer ans;
	ans.E = INFINITY;
	vector<permutation_space> Ps;
	permutation_space P0;
	P0.E = 0;
	for(int i=0;i<b;i++){
		P0.not_assign.push_back(i);
	}
	Ps.push_back(P0);
	while(!Ps.empty()){
		permutation_space P = Ps.back();
		Ps.pop_back();
		int N = P.assignment.size();
		for(int i=0;i<P.not_assign.size();i++){
			if(label_A[N] != label_B[P.not_assign[i]]){
				continue;
			}
			MatrixXd BBB = B;
            MatrixXd BB = MatrixXd::Zero(3,N+1);
			MatrixXd AA = MatrixXd::Zero(3,N+1);
			MatrixXd AAg = MatrixXd::Zero(3,1);
			MatrixXd BBg = MatrixXd::Zero(3,1);
			for(int j=0;j<N;j++){
                for(int k=0;k<3;k++){
					BB(k,j) = B(k,P.assignment[j]);
					AA(k,j) = A(k,j);
				}
			}
			for(int k=0;k<3;k++){
                BB(k,N) = B(k,P.not_assign[i]);
			    AA(k,N) = A(k,N);
			}
			MatrixXd R;
			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
			if(permit_mirror){
                 R = svd.matrixV()*svd.matrixU().transpose();
			}else{
				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
				    R = svd.matrixV()*svd.matrixU().transpose();
			    }else{
					MatrixXd umeyama = MatrixXd::Identity(3,3);
			        umeyama(2,2) = -1;
				    R = svd.matrixV()*umeyama*svd.matrixU().transpose();
				}
			}
			double pre_ans = 0;
			MatrixXd RBB = R*BB;
			for(int j=0;j<N+1;j++){
				for(int k=0;k<3;k++){
                    pre_ans += (AA(k,j)-RBB(k,j))*(AA(k,j)-RBB(k,j));
				}
			}
			if(pre_ans>ans.E){
				continue;
			}else{
				if(N+1==a){
					ans.E = pre_ans;
					ans.R = R;
					permutation_space new_per;
					MatrixXd T = MatrixXd::Zero(b,a); 
					new_per.assignment = P.assignment;
					new_per.assignment.push_back(P.not_assign[i]);
					for(int i=0;i<a;i++){
						T(new_per.assignment[i],i) = 1;
					}
					ans.P = T;
				}else{
					permutation_space new_per;
					new_per.E = pre_ans;
					new_per.assignment = P.assignment;
					new_per.assignment.push_back(P.not_assign[i]);
					new_per.not_assign = P.not_assign;
					new_per.not_assign.erase(new_per.not_assign.begin()+i);
					Ps.push_back(new_per);
				}
			}
		}
	}
	ans.E = sqrt(ans.E/a);
	ans.t = Ag - ans.R*Bg;
	return ans;
}

Answer MatchFastOpt_diff_size(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror){
    int a = A.cols();
    int b = B.cols();

	MatrixXd Ag = MatrixXd::Zero(3,1);
    MatrixXd Bg = MatrixXd::Zero(3,1);

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
	Answer ans;
	ans.E = INFINITY;
	vector<permutation_space> Ps;
	permutation_space P0;
	P0.E = 0;
	for(int i=0;i<b;i++){
		P0.not_assign.push_back(i);
	}
	Ps.push_back(P0);
	while(!Ps.empty()){
		permutation_space P = Ps.back();
		Ps.pop_back();
		int N = P.assignment.size();
		for(int i=0;i<P.not_assign.size();i++){
			if(label_A[N] != label_B[P.not_assign[i]]){
				continue;
			}
			MatrixXd BBB = B;
            MatrixXd BB = MatrixXd::Zero(3,N+1);
			MatrixXd AA = MatrixXd::Zero(3,N+1);
			MatrixXd AAg = MatrixXd::Zero(3,1);
			MatrixXd BBg = MatrixXd::Zero(3,1);
			for(int j=0;j<N;j++){
                for(int k=0;k<3;k++){
					BB(k,j) = B(k,P.assignment[j]);
					AA(k,j) = A(k,j);
				}
			}
			for(int k=0;k<3;k++){
                BB(k,N) = B(k,P.not_assign[i]);
			    AA(k,N) = A(k,N);
			}
			for(int k=0;k<3;k++){
				for(int j=0;j<N+1;j++){
					BBg(k,0) = BBg(k,0) + BB(k,j);
					AAg(k,0) = AAg(k,0) + AA(k,j);
				}
				BBg(k,0) = BBg(k,0)/(N+1);
				AAg(k,0) = AAg(k,0)/(N+1);
				for(int j=0;j<N+1;j++){
					AA(k,j) = AA(k,j) - AAg(k,0);
					BB(k,j) = BB(k,j) - BBg(k,0);
				}
				for(int j=0;j<b;j++){
					BBB(k,j) = B(k,j) - BBg(k,0);
				}
			}
			
			MatrixXd R;
			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
			if(permit_mirror){
                 R = svd.matrixV()*svd.matrixU().transpose();
			}else{
				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
				    R = svd.matrixV()*svd.matrixU().transpose();
			    }else{
					MatrixXd umeyama = MatrixXd::Identity(3,3);
			        umeyama(2,2) = -1;
				    R = svd.matrixV()*umeyama*svd.matrixU().transpose();
				}
			}
			double pre_ans = 0;
			MatrixXd RBB = R*BB;
			for(int j=0;j<N+1;j++){
				for(int k=0;k<3;k++){
                    pre_ans += (AA(k,j)-RBB(k,j))*(AA(k,j)-RBB(k,j));
				}
			}
			if(pre_ans>ans.E){
				continue;
			}else{
				if(N+1==a){
	                ans.E = pre_ans;
					ans.R = R;
					ans.t = AAg - ans.R*BBg;
					permutation_space new_per;
					MatrixXd T = MatrixXd::Zero(b,a); 
					new_per.assignment = P.assignment;
					new_per.assignment.push_back(P.not_assign[i]);
					for(int i=0;i<a;i++){
						T(new_per.assignment[i],i) = 1;
					}
					ans.P = T;		
				}else{
					permutation_space new_per;
					new_per.E = pre_ans;
					new_per.assignment = P.assignment;
					new_per.assignment.push_back(P.not_assign[i]);
					new_per.not_assign = P.not_assign;
					new_per.not_assign.erase(new_per.not_assign.begin()+i);
					Ps.push_back(new_per);
				}
			}
		}
	}
	ans.E = sqrt(ans.E/a);
	ans.t = ans.t + Ag - ans.R*Bg;
	return ans;
}

vector<Answer> IsometrySearch(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
	vector<Answer> ans_multiple;
	int a = A.cols();
    int b = B.cols();

    vector<MatrixXd> As;
    vector<MatrixXd> Bs;

	MatrixXd Ag = MatrixXd::Zero(3,1);
	MatrixXd Bg = MatrixXd::Zero(3,1);

	vector<vector<int> > a_assign;
	vector<vector<int> > b_assign;

    Matrix_vector(A,B,As,Bs,label_A,label_B,a_assign,b_assign,Ag,Bg);

	double max_diameter = 0;
	for(int j=0;j<As.size();j++){
		int aa = As[j].cols();
		for(int i=0;i<aa;i++){
	        double ell = 0;
			for(int k=0;k<3;k++){
				ell += As[j](k,i)*As[j](k,i);
			}
			if(max_diameter<ell){
				max_diameter = ell;
			}
		}
	}
	max_diameter = sqrt(max_diameter);
	for(int ja=0;ja<As.size();ja++){
        int aa = As[ja].cols();
		int bb = Bs[ja].cols();
		for(int ia=0;ia<aa;ia++){
			for(int ib=0;ib<bb;ib++){
				int NN = 0;
				vector<MatrixXd> Bs_temp;
				vector<MatrixXd> BBs;
				vector<MatrixXd> AAs;
				bool valid = true;
				for(int jb = 0;jb<Bs.size();jb++){
					MatrixXd B_temp = MatrixXd::Zero(3,Bs[jb].cols());
					for(int k=0;k<3;k++){
						for(int j=0;j<Bs[jb].cols();j++){
								B_temp(k,j) = Bs[jb](k,j) - Bs[ja](k,ib);
						}
					}
					Bs_temp.push_back(B_temp);
				}
				vector<vector<int> > bb_assign;
				int index_jb = 0;
				for(int jb = 0;jb<Bs.size();jb++){
					if(a_assign[index_jb].size() == 0){
						bb_assign.push_back(b_assign[index_jb]);
						index_jb += 1;
					}
					vector<int> assignment;
					for(int j=0;j<Bs_temp[jb].cols();j++){
						double b_norm = 0;
						for(int k=0;k<3;k++){
                            b_norm += (Bs_temp[jb](k,j))*(Bs_temp[jb](k,j));
						}
						b_norm = sqrt(b_norm);
						if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
							assignment.push_back(j);
							NN+=1;
						}
					}
					if(assignment.size() < a_assign[index_jb].size()){
						valid = false;
					}
			        MatrixXd BB = MatrixXd::Zero(3,assignment.size());
					vector<int> b_assign_temp;
					for(int j=0;j<assignment.size();j++){
						for(int k=0;k<3;k++){
							BB(k,j) = Bs_temp[jb](k,assignment[j]);
						}
						b_assign_temp.push_back(b_assign[index_jb][assignment[j]]);
					}
					bb_assign.push_back(b_assign_temp);
					BBs.push_back(BB);
					index_jb += 1;
				}
				if(!valid){
					continue;
				}
				for(int jb = 0;jb<As.size();jb++){
					MatrixXd AA = As[jb];
					for(int k=0;k<3;k++){
						for(int j=0;j<As[jb].cols();j++){
								AA(k,j) = As[jb](k,j) - As[ja](k,ia);
						}
					}
					AAs.push_back(AA);
				}
				if(NN<a){
					continue;
				}
				vector<double> r;
				for(int k=0;k<3;k++){
					r.push_back(0);
				}
				Cube Cr;
				Cr.center = r;
				Cr.diameter = PI;
				Cr.lower = 0;
				Cr.upper = INFINITY;
				Cr.count = 1;
				Cr.mirror = 0;
				priority_queue<Cube> que;
				que.push(Cr);
				if(permit_mirror){
					Cube Crr;
					Crr.center = r;
					Crr.diameter = PI;
					Crr.lower = 0;
					Crr.upper = INFINITY;
					Crr.mirror = 1;
					Crr.count = 1;
					que.push(Crr);
				}
				MatrixXd initial_matrix = make_rotation(r,0);
				for(int iteration=0;iteration<1000000;iteration++){
					if(que.size()==0){
						break;
					}
					Cube C = que.top();
					que.pop();
					vector<Cube> Cs = divide(C);
					for(int iter=0;iter<8;iter++){
						int K = ans_multiple.size();
						Cube CC = Cs[iter];
						if(!check_in_sphere(CC)){
							continue;
						}
						bool Flag = calculate_lower_K_best_eps(CC,AAs,BBs,CC.mirror,ans_multiple,eps,CC.count+K,a_assign,bb_assign,a,b);
						CC.lower = 0;
						if(!Flag){
							que.push(CC);
						}		
					}
				}
			}
		}
	}
	for(int i=0;i<ans_multiple.size();i++){
		MatrixXd translation = MatrixXd::Zero(3,1);
		for(int i_a = 0;i_a < a;i_a++){
			for(int i_b = 0;i_b < b;i_b++){
				if(ans_multiple[i].P(i_b,i_a) == 1){
                    translation = translation + B.col(i_b);
				}
			}
		}
		ans_multiple[i].t =  Ag - ans_multiple[i].R*(translation/a);
	}
	return ans_multiple;
}


vector<Answer> MatchFPT(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
	int a = A.cols();
    int b = B.cols();

	MatrixXd Ag = MatrixXd::Zero(3,1);
    MatrixXd Bg = MatrixXd::Zero(3,1);

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

	double max_diameter = 0;
	for(int j=0;j<a;j++){
		double ell = 0;
		for(int k=0;k<3;k++){
			ell += A(k,j)*A(k,j);
		}
		if(max_diameter<ell){
			max_diameter = ell;
		}
	}
	max_diameter = sqrt(max_diameter);
	vector<Answer> ans_multiple;
	for(int i=0;i<b;i++){
		vector<permutation_space> Ps;
		permutation_space P0;
        for(int j=0;j<b;j++){
			double b_norm = 0;
			for(int k=0;k<3;k++){
				b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
			}
			b_norm = sqrt(b_norm);
			if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
				P0.not_assign.push_back(j);
			}
		}
		Ps.push_back(P0);
		while(!Ps.empty()){
			permutation_space P = Ps.back();
			Ps.pop_back();
			int N = P.assignment.size();
			for(int i=0;i<P.not_assign.size();i++){
				if(label_A[N] != label_B[P.not_assign[i]]){
				    continue;
			    }
				MatrixXd BB = MatrixXd::Zero(3,N+1);
				MatrixXd AA = MatrixXd::Zero(3,N+1);
				for(int j=0;j<N;j++){
					for(int k=0;k<3;k++){
						BB(k,j) = B(k,P.assignment[j]);
						AA(k,j) = A(k,j);
					}
				}
				for(int k=0;k<3;k++){
					BB(k,N) = B(k,P.not_assign[i]);
					AA(k,N) = A(k,N);
				}
				MatrixXd BBg = MatrixXd::Zero(3,1);
				MatrixXd AAg = MatrixXd::Zero(3,1);
				for(int k=0;k<3;k++){
					for(int j=0;j<N+1;j++){
						BBg(k,0) = BBg(k,0) + BB(k,j);
						AAg(k,0) = AAg(k,0) + AA(k,j);
					}
					BBg(k,0) = BBg(k,0)/(N+1);
					AAg(k,0) = AAg(k,0)/(N+1);
					for(int j=0;j<N+1;j++){
						AA(k,j) = AA(k,j) - AAg(k,0);
						BB(k,j) = BB(k,j) - BBg(k,0);
					}
				}
				JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
				double R_sum = 0;
				if(permit_mirror){
					for(int k=0;k<3;k++){
						R_sum += svd.singularValues()(k,0);
					}
				}else{
					if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
						for(int k=0;k<3;k++){
							R_sum += svd.singularValues()(k,0);
						}
					}else{
						for(int k=0;k<2;k++){
							R_sum += svd.singularValues()(k,0);
						}
						R_sum -= svd.singularValues()(2,0);
					}
				}
				double pre_ans = 0;
				for(int j=0;j<N+1;j++){
					for(int k=0;k<3;k++){
						pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
					}
				}
				pre_ans -= R_sum*2;
				if(pre_ans < 0){
					pre_ans = 0;
				}
				pre_ans = sqrt(pre_ans/(a));
				if(pre_ans>eps){
					continue;
				}else{
					if(N+1==a){
						Answer ans;
						ans.E = pre_ans;
						MatrixXd T = MatrixXd::Zero(b,a);
						vector<int> assignment = P.assignment;
						assignment.push_back(P.not_assign[i]);
						for(int t=0;t<a;t++){
							T(assignment[t],t) = 1;
						}
						bool found = find_P(ans_multiple,T);
						if(!found){
							ans.P = T;
							ans.R = calculate_from_assignment_not_vector(T,A,B,permit_mirror);
							ans.t = (Ag + AAg) - ans.R*(Bg + BBg);
							ans_multiple.push_back(ans);
						}
					}else{
						permutation_space new_per;
						new_per.assignment = P.assignment;
						new_per.assignment.push_back(P.not_assign[i]);
						new_per.not_assign = P.not_assign;
						new_per.not_assign.erase(new_per.not_assign.begin()+i);
						Ps.push_back(new_per);
					}
				}
			}
	    }
	}
	return ans_multiple;
}

vector<Answer> MatchFPT_with_weight(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
	int a = A.cols();
    int b = B.cols();

	MatrixXd Ag = MatrixXd::Zero(3,1);
    MatrixXd Bg = MatrixXd::Zero(3,1);

	int total_label_A = 0;
	for(int i=0;i<3;i++){
		int sum_label_A = 0;
		int sum_label_B = 0;
		for(int j=0;j<a;j++){
			Ag(i,0) = Ag(i,0) + A(i,j)*label_A[j];
			sum_label_A += label_A[j];
		}
		for(int j=0;j<b;j++){
			Bg(i,0) = Bg(i,0) + B(i,j)*label_B[j];
			sum_label_B += label_B[j];
		}
		Ag(i,0) = Ag(i,0)/sum_label_A;
		Bg(i,0) = Bg(i,0)/sum_label_B;
		for(int j=0;j<a;j++){
			A(i,j) = A(i,j) - Ag(i,0);
		}
		for(int j=0;j<b;j++){
			B(i,j) = B(i,j) - Bg(i,0);
		}
		total_label_A = sum_label_A;
	}

	double max_diameter = 0;
	for(int j=0;j<a;j++){
		double ell = 0;
		for(int k=0;k<3;k++){
			ell += A(k,j)*A(k,j);
		}
		if(max_diameter<ell){
			max_diameter = ell;
		}
	}
	max_diameter = sqrt(max_diameter);
	vector<Answer> ans_multiple;
	for(int i=0;i<b;i++){
		vector<permutation_space> Ps;
		permutation_space P0;
        for(int j=0;j<b;j++){
			double b_norm = 0;
			for(int k=0;k<3;k++){
				b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
			}
			b_norm = sqrt(b_norm);
			b_norm = 0;
			if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
				P0.not_assign.push_back(j);
			}
		}
		Ps.push_back(P0);
		while(!Ps.empty()){
			permutation_space P = Ps.back();
			Ps.pop_back();
			int N = P.assignment.size();
			for(int i=0;i<P.not_assign.size();i++){
				if(label_A[N] != label_B[P.not_assign[i]]){
				    continue;
			    }
				MatrixXd BB = MatrixXd::Zero(3,N+1);
				MatrixXd AA = MatrixXd::Zero(3,N+1);
				for(int j=0;j<N;j++){
					for(int k=0;k<3;k++){
						BB(k,j) = B(k,P.assignment[j]);
						AA(k,j) = A(k,j);
					}
				}
				for(int k=0;k<3;k++){
					BB(k,N) = B(k,P.not_assign[i]);
					AA(k,N) = A(k,N);
				}
				MatrixXd BBg = MatrixXd::Zero(3,1);
				MatrixXd AAg = MatrixXd::Zero(3,1);
				for(int k=0;k<3;k++){
					double sum_label_A = 0;
					for(int j=0;j<N+1;j++){
						BBg(k,0) = BBg(k,0) + BB(k,j)*label_A[j];
						AAg(k,0) = AAg(k,0) + AA(k,j)*label_A[j];
						sum_label_A += label_A[j];
					}
					BBg(k,0) = BBg(k,0)/sum_label_A;
					AAg(k,0) = AAg(k,0)/sum_label_A;
					for(int j=0;j<N+1;j++){
						AA(k,j) = AA(k,j) - AAg(k,0);
						BB(k,j) = BB(k,j) - BBg(k,0);
					}
				}
				MatrixXd W = MatrixXd::Zero(N+1,N+1);
				for(int j=0;j<N+1;j++){
					W(j,j) = label_A[j];
				}
				JacobiSVD< Matrix<double, 3, 3> > svd(BB*W*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
				double R_sum = 0;
				if(permit_mirror){
					for(int k=0;k<3;k++){
						R_sum += svd.singularValues()(k,0);
					}
				}else{
					if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
						for(int k=0;k<3;k++){
							R_sum += svd.singularValues()(k,0);
						}
					}else{
						for(int k=0;k<2;k++){
							R_sum += svd.singularValues()(k,0);
						}
						R_sum -= svd.singularValues()(2,0);
					}
				}
				double pre_ans = 0;
				for(int j=0;j<N+1;j++){
					for(int k=0;k<3;k++){
						pre_ans += AA(k,j)*AA(k,j)*label_A[j] + BB(k,j)*BB(k,j)*label_A[j];
					}
				}
				pre_ans -= R_sum*2;
				if(pre_ans < 0){
					pre_ans = 0;
				}
				pre_ans = sqrt(pre_ans/(total_label_A));
				if(pre_ans>eps){
					continue;
				}else{
					if(N+1==a){
						Answer ans;
						ans.E = pre_ans;
						MatrixXd T = MatrixXd::Zero(b,a);
						vector<int> assignment = P.assignment;
						assignment.push_back(P.not_assign[i]);
						for(int i=0;i<a;i++){
							T(assignment[i],i) = 1;
						}
						bool found = find_P(ans_multiple,T);
						if(!found){
							ans.P = T;
							ans.R = calculate_from_assignment_not_vector_weight(T,A,B,permit_mirror,label_A);
							ans.t = (Ag + AAg) - ans.R*(Bg + BBg);
							ans_multiple.push_back(ans);
						}
					}else{
						permutation_space new_per;
						new_per.assignment = P.assignment;
						new_per.assignment.push_back(P.not_assign[i]);
						new_per.not_assign = P.not_assign;
						new_per.not_assign.erase(new_per.not_assign.begin()+i);
						Ps.push_back(new_per);
					}
				}
			}
	    }
	}
	return ans_multiple;
}

vector<Answer> Improved_MatchFPT(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
	int a = A.cols();
    int b = B.cols();

	MatrixXd Ag = MatrixXd::Zero(3,1);
    MatrixXd Bg = MatrixXd::Zero(3,1);

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
	vector<Answer> ans_multiple;
	vector<permutation_space> Ps;
	permutation_space P0;
	vector<int> assignment;
	for(int j=0;j<b;j++){
		assignment.push_back(j);
	}
	P0.assign_candidate.push_back(assignment);
	P0.not_assign = assignment;
	Ps.push_back(P0);
	while(!Ps.empty()){
		permutation_space P = Ps.back();
		Ps.pop_back();
		int N = P.assignment.size();
		if(P.assign_candidate[0].size() == 0){
			continue;
		}
		for(int i=0;i<P.assign_candidate[0].size();i++){
			if(label_A[N] != label_B[P.assign_candidate[0][i]]){
				continue;
			}
            MatrixXd BB = MatrixXd::Zero(3,N+1);
			MatrixXd AA = MatrixXd::Zero(3,N+1);
			for(int j=0;j<N;j++){
                for(int k=0;k<3;k++){
					BB(k,j) = B(k,P.assignment[j]);
					AA(k,j) = A(k,j);
				}
			}
			for(int k=0;k<3;k++){
                BB(k,N) = B(k,P.assign_candidate[0][i]);
			    AA(k,N) = A(k,N);
			}
			MatrixXd BBg = MatrixXd::Zero(3,1);
			MatrixXd AAg = MatrixXd::Zero(3,1);
			for(int k=0;k<3;k++){
				for(int j=0;j<N+1;j++){
					BBg(k,0) = BBg(k,0) + BB(k,j);
					AAg(k,0) = AAg(k,0) + AA(k,j);
				}
				BBg(k,0) = BBg(k,0)/(N+1);
				AAg(k,0) = AAg(k,0)/(N+1);
				for(int j=0;j<N+1;j++){
					AA(k,j) = AA(k,j) - AAg(k,0);
					BB(k,j) = BB(k,j) - BBg(k,0);
				}
			}
			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
			double R_sum = 0;
			if(permit_mirror){
				for(int k=0;k<3;k++){
					R_sum += svd.singularValues()(k,0);
				}
			}else{
				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
					for(int k=0;k<3;k++){
					    R_sum += svd.singularValues()(k,0);
				    }
			    }else{
					for(int k=0;k<2;k++){
					    R_sum += svd.singularValues()(k,0);
				    }
					R_sum -= svd.singularValues()(2,0);
				}
			}
			double pre_ans = 0;
			for(int j=0;j<N+1;j++){
				for(int k=0;k<3;k++){
                    pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
				}
			}
			pre_ans -= R_sum*2;
			if(pre_ans < 0){
				pre_ans = 0;
			}
			pre_ans = sqrt(pre_ans/(a));
			double a_b_norm = 0;
			for(int j=0;j<N+1;j++){
				double a_norm = 0;
				double b_norm = 0;
				for(int k=0;k<3;k++){
					a_norm += AA(k,j)*AA(k,j);
					b_norm += BB(k,j)*BB(k,j);
			    }
				a_norm = sqrt(a_norm);
				b_norm = sqrt(b_norm);
				a_b_norm += (b_norm - a_norm) * (b_norm - a_norm);
			}
			if(pre_ans>eps){
				continue;
			}else{
				if(N+1==a){
					Answer ans;
					ans.E = pre_ans;
					MatrixXd T = MatrixXd::Zero(b,a);
					vector<int> assignment = P.assignment;
					assignment.push_back(P.not_assign[i]);
					for(int i=0;i<a;i++){
						T(assignment[i],i) = 1;
					}
					bool found = find_P(ans_multiple,T);
					if(!found){
						ans.P = T;
						ans.R = calculate_from_assignment_not_vector(T,A,B,permit_mirror);
						ans.t = (Ag + AAg) - ans.R*(Bg + BBg);
                        ans_multiple.push_back(ans);
					}
				}else{
					permutation_space new_per;
					new_per.assignment = P.assignment;
					new_per.assignment.push_back(P.assign_candidate[0][i]);
					vector<vector<int> > assign_candidate(1);
					vector<int> not_assign;
					double max_diameter = 0;
						for(int j=0;j<a;j++){
							double ell = 0;
							for(int k=0;k<3;k++){
								ell += A(k,j)*A(k,j);
							}
							if(max_diameter<ell){
								max_diameter = ell;
							}
						}
						max_diameter = sqrt(max_diameter);
					for(int k1=0;k1<P.not_assign.size();k1++){
						if(P.assign_candidate[0][i]==P.not_assign[k1]){
							continue;
						}
						double distance_new_old = 0;
						for(int k=0;k<3;k++){
							distance_new_old += (B(k,P.not_assign[k1]) - B(k,P.assign_candidate[0][i]))*(B(k,P.not_assign[k1]) - B(k,P.assign_candidate[0][i]));
						}
						if(sqrt(distance_new_old)<2*(eps*sqrt(a) + max_diameter)){
							not_assign.push_back(P.not_assign[k1]);
						}
						if(label_A[N+1] != label_B[P.not_assign[k1]]){
							continue;
						}
						double value = 0;
						value += a_b_norm;
						double a_norm = 0;
						double b_norm = 0;
						for(int k=0;k<3;k++){
							a_norm += (A(k,N+1)-AAg(k,0))*(A(k,N+1)-AAg(k,0));
							b_norm += (B(k,P.not_assign[k1])-BBg(k,0))*(B(k,P.not_assign[k1])-BBg(k,0));
						}
						a_norm = sqrt(a_norm);
						b_norm = sqrt(b_norm);
						value += (b_norm - a_norm)*(b_norm - a_norm)*(N+1)/(N+2);
						if(value<a*eps*eps){
							assign_candidate[0].push_back(P.not_assign[k1]);
						}

					}
					new_per.assign_candidate = assign_candidate;
					new_per.not_assign = not_assign;
					Ps.push_back(new_per);
				}
			}
		}
	}
	return ans_multiple;
}

// vector<vector<int> > find_lattice(MatrixXd A, MatrixXd B, int index, double eps, permutation_space P, vector<double> lattice){
// 	vector<vector<int> > lattice_list;
// 	vector<int> assignment = P.assignment;
// 	int next_index = assignment.size();
// 	for(int i=0;i<assignment.size();i++){
// 		double length = 0;
// 		for(int j=0;j<3;j++){
// 			length += (A(j,i) - A(j,next_index))*(A(j,i) - A(j,next_index));
// 		}
// 		length = sqrt(length);
// 		cout << "length: " << length << endl;
// 		cout << next_index << " " << index << endl;
// 		if(i == 0){
// 			double r_1 = -eps*sqrt(A.cols()) + length;
// 			double r_2 = eps*sqrt(A.cols()) + length;
// 			cout << "r_1: " << r_1 << " r_2: " << r_2 << endl;
// 			// laticeの符号
// 			double k1_min = (B(0,assignment[i]) + P.lattice_count[i][0]*lattice[0] - B(0,index) - r_2)/abs(lattice[0]);
// 			double k1_max = (B(0,assignment[i]) + P.lattice_count[i][0]*lattice[0] - B(0,index) + r_2)/abs(lattice[0]);
// 			int k1_min_int = ceil(k1_min);
// 			int k1_max_int = floor(k1_max);
// 			if(lattice[0] < 0){
// 				int temp = -k1_min_int;
// 				k1_min_int = -k1_max_int;
// 				k1_max_int = temp;
// 			}
// 			for(int k1=k1_min_int;k1<=k1_max_int;k1++){
// 				double r_x = abs(B(0,assignment[i]) + (P.lattice_count[i][0] - k1)*lattice[0] - B(0,index));
// 				double r_2_x = r_2*r_2 - r_x*r_x;
// 				r_2_x = sqrt(r_2_x);
// 				double k2_min = (B(1,assignment[i]) + P.lattice_count[i][1]*lattice[1] - B(1,index) - r_2_x)/abs(lattice[1]);
// 				double k2_max = (B(1,assignment[i]) + P.lattice_count[i][1]*lattice[1] - B(1,index) + r_2_x)/abs(lattice[1]);
// 				int k2_min_int = ceil(k2_min);
// 				int k2_max_int = floor(k2_max);
// 				if(lattice[1] < 0){
// 					int temp = -k2_min_int;
// 					k2_min_int = -k2_max_int;
// 					k2_max_int = temp;
// 				}
// 				for(int k2=k2_min_int;k2<=k2_max_int;k2++){
// 					double r_y = abs(B(1,assignment[i]) + (P.lattice_count[i][1] - k2)*lattice[1] - B(1,index));
// 					double r_2_y = r_2_x*r_2_x - r_y*r_y;
// 					r_2_y = sqrt(r_2_y);
// 					double k3_min = (B(2,assignment[i]) + P.lattice_count[i][2]*lattice[2] - B(2,index) - r_2_y)/abs(lattice[2]);
// 					double k3_max = (B(2,assignment[i]) + P.lattice_count[i][2]*lattice[2] - B(2,index) + r_2_y)/abs(lattice[2]);
// 					int k3_min_int = ceil(k3_min);
// 					int k3_max_int = floor(k3_max);
// 					if(lattice[2] < 0){
// 						int temp = -k3_min_int;
// 						k3_min_int = -k3_max_int;
// 						k3_max_int = temp;
// 					}
// 					for(int k3=k3_min_int;k3<=k3_max_int;k3++){
// 						double r = 0;
// 						r += (B(0,assignment[i]) + (P.lattice_count[i][0] - k1)*lattice[0] - B(0,index))*(B(0,assignment[i]) + (P.lattice_count[i][0] - k1)*lattice[0] - B(0,index));
// 						r += (B(1,assignment[i]) + (P.lattice_count[i][1] - k2)*lattice[1] - B(1,index))*(B(1,assignment[i]) + (P.lattice_count[i][1] - k2)*lattice[1] - B(1,index));
// 						r += (B(2,assignment[i]) + (P.lattice_count[i][2] - k3)*lattice[2] - B(2,index))*(B(2,assignment[i]) + (P.lattice_count[i][2] - k3)*lattice[2] - B(2,index));
// 						r = sqrt(r);
// 						if (r > r_1){
// 							vector<int> lattice_temp;
// 							lattice_temp.push_back(k1);
// 							lattice_temp.push_back(k2);
// 							lattice_temp.push_back(k3);
// 							lattice_list.push_back(lattice_temp);
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return lattice_list;
// }



// vector<Answer> MatchFPT_lattice(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror, vector<double> lattice){
// 	int a = A.cols();
//     int b = B.cols();

// 	MatrixXd Ag = MatrixXd::Zero(3,1);
//     MatrixXd Bg = MatrixXd::Zero(3,1);

// 	for(int i=0;i<3;i++){
// 		for(int j=0;j<a;j++){
// 			Ag(i,0) = Ag(i,0) + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
// 			Bg(i,0) = Bg(i,0) + B(i,j);
// 		}
// 		Ag(i,0) = Ag(i,0)/a;
// 		Bg(i,0) = Bg(i,0)/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag(i,0);
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg(i,0);
// 		}
// 	}

// 	// double max_diameter = 0;
// 	// for(int j=0;j<a;j++){
// 	// 	double ell = 0;
// 	// 	for(int k=0;k<3;k++){
// 	// 		ell += A(k,j)*A(k,j);
// 	// 	}
// 	// 	if(max_diameter<ell){
// 	// 		max_diameter = ell;
// 	// 	}
// 	// }
// 	// max_diameter = sqrt(max_diameter);
// 	vector<Answer> ans_multiple;
// 	for(int i=0;i<b;i++){
// 		vector<permutation_space> Ps;
// 		permutation_space P0;
// 		P0.assignment.push_back(i);
// 		if(label_A[0] != label_B[i]){
// 			continue;
// 		}
// 		for(int j=0;j<b;j++){
// 			if(j == i){
// 				continue;
// 			}
// 			P0.not_assign.push_back(j);
// 		}
// 		P0.lattice_count.push_back(vector<int>(3,0));

//         // for(int j=0;j<b;j++){
// 		// 	double b_norm = 0;
// 		// 	for(int k=0;k<3;k++){
// 		// 		b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
// 		// 	}
// 		// 	b_norm = sqrt(b_norm);
// 		// 	if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
// 		// 		P0.not_assign.push_back(j);
// 		// 	}
// 		// }
// 		cout << "P0.not_assign.size() = " << P0.not_assign.size() << endl;
// 		Ps.push_back(P0);
// 		while(!Ps.empty()){
// 			permutation_space P = Ps.back();
// 			Ps.pop_back();
// 			int N = P.assignment.size();
// 			for(int i=0;i<P.not_assign.size();i++){
// 				if(label_A[N] != label_B[P.not_assign[i]]){
// 				    continue;
// 			    }
// 				cout << "P.not_assign[i] = " << P.not_assign[i] << endl;
// 				// P.assignmentを出力
// 				for (int i = 0; i < P.assignment.size(); i++){
// 					cout << P.assignment[i] << " ";
// 				}
// 				cout << endl;
// 				vector<vector<int> > lattice_list = find_lattice(A,B,P.not_assign[i],eps,P,lattice);
// 				// lattice_listを出力
// 				cout << "lattice_list.size() = " << lattice_list.size() << endl;
// 				for (int i = 0; i < lattice_list.size(); i++){
// 					for (int j = 0; j < lattice_list[i].size(); j++){
// 						cout << lattice_list[i][j] << " ";
// 					}
// 					cout << endl;
// 				}
// 				MatrixXd BB = MatrixXd::Zero(3,N+1);
// 				MatrixXd AA = MatrixXd::Zero(3,N+1);
// 				for(int j=0;j<N;j++){
// 					for(int k=0;k<3;k++){
// 						BB(k,j) = B(k,P.assignment[j]);
// 						AA(k,j) = A(k,j);
// 					}
// 				}
// 				for(int lattice_index = 0; lattice_index < lattice_list.size(); lattice_index++){
// 					for(int j=0;j<N;j++){
// 						for(int k=0;k<3;k++){
// 							BB(k,j) = B(k,P.assignment[j]) + P.lattice_count[j][k]*lattice[k];
// 							AA(k,j) = A(k,j);
// 						}
// 					}
// 					for(int k=0;k<3;k++){
// 						BB(k,N) = B(k,P.not_assign[i]) + lattice_list[lattice_index][k]*lattice[k];
// 						AA(k,N) = A(k,N);
// 					}
// 					MatrixXd BBg = MatrixXd::Zero(3,1);
// 					MatrixXd AAg = MatrixXd::Zero(3,1);
// 					for(int k=0;k<3;k++){
// 						for(int j=0;j<N+1;j++){
// 							BBg(k,0) = BBg(k,0) + BB(k,j);
// 							AAg(k,0) = AAg(k,0) + AA(k,j);
// 						}
// 						BBg(k,0) = BBg(k,0)/(N+1);
// 						AAg(k,0) = AAg(k,0)/(N+1);
// 						for(int j=0;j<N+1;j++){
// 							AA(k,j) = AA(k,j) - AAg(k,0);
// 							BB(k,j) = BB(k,j) - BBg(k,0);
// 						}
// 					}
// 					JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 					double R_sum = 0;
// 					if(permit_mirror){
// 						for(int k=0;k<3;k++){
// 							R_sum += svd.singularValues()(k,0);
// 						}
// 					}else{
// 						if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 							for(int k=0;k<3;k++){
// 								R_sum += svd.singularValues()(k,0);
// 							}
// 						}else{
// 							for(int k=0;k<2;k++){
// 								R_sum += svd.singularValues()(k,0);
// 							}
// 							R_sum -= svd.singularValues()(2,0);
// 						}
// 					}
// 					double pre_ans = 0;
// 					for(int j=0;j<N+1;j++){
// 						for(int k=0;k<3;k++){
// 							pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 						}
// 					}
// 					pre_ans -= R_sum*2;
// 					if(pre_ans < 0){
// 						pre_ans = 0;
// 					}
// 					pre_ans = sqrt(pre_ans/(a));
// 					if(pre_ans>eps){
// 						continue;
// 					}else{
// 						if(N+1==a){
// 							Answer ans;
// 							ans.E = pre_ans;
// 							MatrixXd T = MatrixXd::Zero(b,a);
// 							vector<int> assignment = P.assignment;
// 							assignment.push_back(P.not_assign[i]);
// 							for(int i=0;i<a;i++){
// 								T(assignment[i],i) = 1;
// 							}
// 							bool found = find_P(ans_multiple,T);
// 							if(!found){
// 								ans.P = T;
// 								vector<vector<int> > lattice_count = P.lattice_count;
// 								lattice_count.push_back(lattice_list[lattice_index]);
// 								ans.lattice_count = lattice_count;
// 								MatrixXd BBB = MatrixXd::Zero(3,b);
// 								for(int j=0;j<b;j++){
// 									for(int k=0;k<3;k++){
// 										BBB(k,j) = B(k,j) - BBg(k,0);
// 									}
// 								}
// 								cout << a << b << endl;
// 								for(int j=0;j<a;j++){
// 									for(int k=0;k<3;k++){
// 										BBB(k,assignment[j]) = B(k,assignment[j]) + lattice_count[j][k]*lattice[k] - BBg(k,0);
// 									}
// 								}
// 								ans.R = calculate_from_assignment_not_vector(T,A,BBB,permit_mirror);
// 								ans.t = (Ag + AAg) - ans.R*(Bg + BBg);
								
// 								ans_multiple.push_back(ans);
// 							}
// 						}else{
// 							permutation_space new_per;
// 							new_per.assignment = P.assignment;
// 							new_per.assignment.push_back(P.not_assign[i]);
// 							new_per.not_assign = P.not_assign;
// 							new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 							vector<vector<int> > lattice_count = P.lattice_count;
// 							lattice_count.push_back(lattice_list[lattice_index]);
// 							new_per.lattice_count = lattice_count;
// 							Ps.push_back(new_per);
// 						}
// 					}
// 				}
// 			}
// 	    }
// 	}
// 	return ans_multiple;
// }


vector<vector<int> > find_lattice(MatrixXd A, MatrixXd B, int index, double eps, permutation_space P, vector<vector<double> > lattice){
	vector<vector<int> > lattice_list;
	vector<int> assignment = P.assignment;
	int next_index = assignment.size();
	vector<double> lattice_norm;
	for(int i=0;i<3;i++){
		double norm = 0;
		for(int j=0;j<3;j++){
			norm += lattice[i][j]*lattice[i][j];
		}
		norm = sqrt(norm);
		lattice_norm.push_back(norm);
	}
	vector<double> lattice_dot;
	for(int i=0;i<2;i++){
		for(int j=i+1;j<3;j++){
			double dot = 0;
			for(int k=0;k<3;k++){
				dot += lattice[i][k]*lattice[j][k];
			}
			lattice_dot.push_back(dot);
		}
	}
	MatrixXd lattice_matrix = MatrixXd::Zero(3,3);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			lattice_matrix(i,j) = lattice[j][i];
		}
	}
	for(int i=0;i<assignment.size();i++){
		double length = 0;
		for(int j=0;j<3;j++){
			length += (A(j,i) - A(j,next_index))*(A(j,i) - A(j,next_index));
		}
		length = sqrt(length);
		MatrixXd B_temp = B.col(assignment[i]);
		for(int j=0;j<3;j++){
			B_temp = B_temp - P.lattice_count[i][j]*lattice_matrix.col(j);
		}
		MatrixXd C = lattice_matrix.inverse()*(B_temp - B.col(index));
		if(i == 0){
			double r_1 = -2*eps*sqrt(A.cols()) + length;
			double r_2 = 2*eps*sqrt(A.cols()) + length;
			vector<vector<int> > set_ans;
			vector<vector<int> > visited;
			for(int i_x = 0; i_x <=1; i_x++){
				for(int i_y = 0; i_y <=1; i_y++){
					for(int i_z = 0; i_z <=1; i_z++){
						vector<int> nint(3,0);
						nint[0] = -floor(C(0,0) + 1*i_x);
						nint[1] = -floor(C(1,0) + 1*i_y);
						nint[2] = -floor(C(2,0) + 1*i_z);
						set_ans.push_back(nint);
						visited.push_back(nint);
					}
				}
			}
			while(set_ans.size() != 0){
				vector<int> temp = set_ans.back();
				set_ans.pop_back();
				double r = 0;
				MatrixXd temp_mat = B_temp - B.col(index);
				for(int l=0;l<3;l++){
					temp_mat += lattice_matrix.col(l)*temp[l];
				}
				for(int l=0;l<3;l++){
					r += temp_mat(l,0)*temp_mat(l,0);
				}
				r = sqrt(r);
				// if(next_index == 1 && index == 3){
				// 	cout << "r_1: " << r_1 << " r_2: " << r_2 << endl;
				// 	cout << "temp: " << temp[0] << " " << temp[1] << " " << temp[2] << endl;
				// 	cout << "r: " << r << endl;
				// 	cout << "C: " << C << endl;
				// }
				if(r < r_2){
					for(int i1 = -1; i1 <= 1; i1++){
						for(int i2 = -1; i2 <= 1; i2++){
							for(int i3 = -1; i3 <= 1; i3++){
								vector<int> temp_temp = temp;
								temp_temp[0] += i1;
								temp_temp[1] += i2;
								temp_temp[2] += i3;
								if(find(visited.begin(), visited.end(), temp_temp) == visited.end()){
									set_ans.push_back(temp_temp);
									visited.push_back(temp_temp);
								}
							}
						}
					}
					// // 直線的なlatticeの探索
					// for(int i = 0; i <=2; i++){
					// 	if(-0.5 <= C(i,0) + temp[i] && C(i,0) + temp[i] <= 0.5){
					// 		for(int i1 = -1; i1 <= 1; i1 += 2){
					// 			for(int i2 = -1; i2 <= 1; i2 += 2){
					// 				vector<int> temp_temp = temp;
					// 				// i以外の要素をそれぞれi1,i2だけ動かす
					// 				bool first = true;
					// 				for(int j = 0; j < 3; j++){
					// 					if(j != i && first){
					// 						temp_temp[j] += i1;
					// 						first = false;
					// 					}else if(j != i && !first){
					// 						temp_temp[j] += i2;
					// 					}
					// 				}
					// 				if(find(visited.begin(), visited.end(), temp_temp) == visited.end()){
					// 					set_ans.push_back(temp_temp);
					// 					visited.push_back(temp_temp);
					// 				}
					// 			}
					// 		}
					// 		for(int j = 0; j < i; j++){
					// 			if(-0.5 <= C(j,0) + temp[j] && C(j,0) + temp[j] <= 0.5){
					// 				for(int i1 = -1; i1 <= 1; i1 += 2){
					// 					vector<int> temp_temp = temp;
					// 					for(int k = 0; k < 3; k++){
					// 						if(k != i && k != j){
					// 							temp_temp[k] += i1;
					// 						}
					// 					}
					// 					if(find(visited.begin(), visited.end(), temp_temp) == visited.end()){
					// 						set_ans.push_back(temp_temp);
					// 						visited.push_back(temp_temp);
					// 					}
					// 				}
					// 			}
					// 		}
					// 	}
					// }
					// vector<int> temp_temp = temp;
					// for(int i1 = -1; i1 <= 1; i1 += 2){
					// 	for(int i2 = -1; i2 <= 1; i2 += 2){
					// 		for(int i3 = -1; i3 <= 1; i3 += 2){
					// 			vector<int> temp_temp = temp;
					// 			temp_temp[0] += i1;
					// 			temp_temp[1] += i2;
					// 			temp_temp[2] += i3;
					// 			if(find(visited.begin(), visited.end(), temp_temp) == visited.end()){
					// 				set_ans.push_back(temp_temp);
					// 				visited.push_back(temp_temp);
					// 			}
					// 		}
					// 	}
					// }
				}
				if (r > r_1 && r < r_2){
					lattice_list.push_back(temp);
				}

			}
		}

		// if(P.assignment[0] == 0 && P.assignment[1] == 1 && index == 2){
		// 	cout << P.lattice_count[i][0] << " " << P.lattice_count[i][1] << " " << P.lattice_count[i][2] << endl;
		// 	cout << "B_temp: " << B_temp << endl;
		// 	cout << "B_temp - B.col(index): " << B_temp - B.col(index) << endl;
		// }
		// MatrixXd C = lattice_matrix.inverse()*(B_temp - B.col(index));
		// if(P.assignment[0] == 0 && P.assignment[1] == 1 && index == 2){
		// 	cout << "C: " << C << endl;
		// }

		// if(i == 0){
			
		// 	// laticeの符号
		// 	double k1_min = - r_2/lattice_norm[0] - C(0,0);
		// 	double k1_max = r_2/lattice_norm[0] - C(0,0);
		// 	int k1_min_int = ceil(k1_min);
		// 	int k1_max_int = floor(k1_max);
		// 	if(P.assignment[0] == 0 && P.assignment[1] == 1 && index == 2){
		// 		cout << r_2 << endl;
		// 		cout << lattice_norm[0] << endl;
		// 		cout << k1_min << " " << k1_max << endl;
		// 	}
		// 	for(int k1=k1_min_int;k1<=k1_max_int;k1++){
		// 		double d_1 = C(0,0) + k1;
		// 		double r_2_x = r_2*r_2 + lattice_dot[0]*lattice_dot[0]*d_1*d_1/(lattice_norm[1]*lattice_norm[1]) - lattice_norm[0]*lattice_norm[0]*d_1*d_1;
		// 		r_2_x = sqrt(r_2_x);
		// 		double k2_min = - r_2_x/lattice_norm[1] - C(1,0) - lattice_dot[0]*d_1/(lattice_norm[1]*lattice_norm[1]);
		// 		double k2_max = r_2_x/lattice_norm[1] - C(1,0) - lattice_dot[0]*d_1/(lattice_norm[1]*lattice_norm[1]);
		// 		int k2_min_int = ceil(k2_min);
		// 		int k2_max_int = floor(k2_max);
				
		// 		for(int k2=k2_min_int;k2<=k2_max_int;k2++){
		// 			double d_2 = C(1,0) + k2;
		// 			double temp = d_1*lattice_dot[1] + d_2*lattice_dot[2];
		// 			temp = temp/(lattice_norm[2]*lattice_norm[2]);
		// 			double r_2_y = r_2*r_2 - d_1*d_1*lattice_norm[0]*lattice_norm[0] - d_2*d_2*lattice_norm[1]*lattice_norm[1] + temp*temp*lattice_norm[2]*lattice_norm[2] - 2*d_1*d_2*lattice_dot[0];
		// 			r_2_y = sqrt(r_2_y);
		// 			double k3_min = - r_2_y/lattice_norm[2] - C(2,0) - temp/lattice_norm[2];
		// 			double k3_max = r_2_y/lattice_norm[2] - C(2,0) - temp/lattice_norm[2];
		// 			int k3_min_int = ceil(k3_min);
		// 			int k3_max_int = floor(k3_max);
		// 			for(int k3=k3_min_int;k3<=k3_max_int;k3++){
		// 				double r = 0;
		// 				MatrixXd temp_mat = B_temp - B.col(index) + lattice_matrix.col(0)*k1 + lattice_matrix.col(1)*k2 + lattice_matrix.col(2)*k3;
		// 				for(int l=0;l<3;l++){
		// 					r += temp_mat(l,0)*temp_mat(l,0);
		// 				}
		// 				r = sqrt(r);
		// 				if (r > r_1 && r < r_2){
		// 					vector<int> lattice_temp;
		// 					lattice_temp.push_back(k1);
		// 					lattice_temp.push_back(k2);
		// 					lattice_temp.push_back(k3);
		// 					lattice_list.push_back(lattice_temp);
		// 				}
		// 			}
		// 		}
		// 	}
		// }
	}
	return lattice_list;
}

bool check_lattice_normalize(vector<vector<double> > &lattice){
	bool flag = true;
	vector<double> lattice_norm(3);
	for(int i=0;i<3;i++){
		double norm = 0;
		for(int j=0;j<3;j++){
			norm += lattice[i][j]*lattice[i][j];
		}
		norm = sqrt(norm);
		lattice_norm[i] = norm;
	}
	for(int i=0;i<3;i++){
		for(int j=i+1;j<3;j++){
			double dot = 0;
			for(int k=0;k<3;k++){
				dot += lattice[i][k]*lattice[j][k];
			}
			// cout << i << " " << j << endl;
			// cout << "dot " << dot << endl;
			// cout << "lattice_norm[i] " << lattice_norm[i] << endl;
			// cout << "lattice_norm[j] " << lattice_norm[j] << endl; 
			if(abs(dot) > lattice_norm[i]*lattice_norm[i]/2 || abs(dot) > lattice_norm[j]*lattice_norm[j]/2){
				flag = false;
				return flag;
			}
		}
	}
	return flag;
}

MatrixXd P_select(MatrixXd lattice_vec, int index){
	if(0 <= index && index <= 2){
		return lattice_vec.row(index);
	}else if(index == 3){
		return -(lattice_vec.row(0) + lattice_vec.row(1) + lattice_vec.row(2));
	}else if(index == 4){
		return lattice_vec.row(0) + lattice_vec.row(1);
	}else if(index == 5){
		return lattice_vec.row(0) + lattice_vec.row(2);
	}else if(index == 6){
		return lattice_vec.row(1) + lattice_vec.row(2);
	}
	cout << "index is out of range" << endl;
	return lattice_vec.row(0);
}


vector<Answer> MatchFPT_lattice(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror, vector<vector<double> > lattice){
	int a = A.cols();
    int b = B.cols();

	MatrixXd Ag = MatrixXd::Zero(3,1);
    MatrixXd Bg = MatrixXd::Zero(3,1);

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

	MatrixXd lattice_mat = MatrixXd::Zero(3,3);
	MatrixXd P_ans = MatrixXd::Identity(3,3);


	if(!check_lattice_normalize(lattice)){
		// cout << "lattice is not normalized" << endl;
		bool flag = false;
		// latticeを行列に変換
		
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				lattice_mat(i,j) = lattice[i][j];
			}
		}

		MatrixXd lattice_before = lattice_mat;

		// cout << "lattice_mat" << endl;
		// cout << lattice_mat << endl;
		MatrixXd lattice_vec = MatrixXd::Identity(3,3);
		VectorXd sum_vec_neg = -(lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
		while(flag == false){
			flag = true;
			if(lattice_mat.row(0).dot(lattice_mat.row(1)) > 0){
				MatrixXd P_mat = MatrixXd::Identity(3,3);
				P_mat(0,0) = -1;
				P_mat(2,0) = 1;
				lattice_vec = P_mat*lattice_vec;
				lattice_mat = P_mat*lattice_mat;
				sum_vec_neg = - (lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
				flag = false;
			}
			if(lattice_mat.row(0).dot(lattice_mat.row(2)) > 0){
				MatrixXd P_mat = MatrixXd::Identity(3,3);
				P_mat(2,2) = -1;
				P_mat(1,2) = 1;
				lattice_vec = P_mat*lattice_vec;
				lattice_mat = P_mat*lattice_mat;
				sum_vec_neg = - (lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
				flag = false;
			}
			if(lattice_mat.row(1).dot(lattice_mat.row(2)) > 0){
				MatrixXd P_mat = MatrixXd::Identity(3,3);
				P_mat(0,1) = 1;
				P_mat(1,1) = -1;
				lattice_vec = P_mat*lattice_vec;
				lattice_mat = P_mat*lattice_mat;
				sum_vec_neg = - (lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
				flag = false;
			}
			if(lattice_mat.row(0).dot(sum_vec_neg) > 0){
				MatrixXd P_mat = MatrixXd::Identity(3,3);
				P_mat(0,0) = -1;
				P_mat(1,0) = 1;
				P_mat(2,0) = 1;
				lattice_vec = P_mat*lattice_vec;
				lattice_mat = P_mat*lattice_mat;
				sum_vec_neg = - (lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
				flag = false;
			}
			if(lattice_mat.row(1).dot(sum_vec_neg) > 0){
				MatrixXd P_mat = MatrixXd::Identity(3,3);
				P_mat(0,1) = 1;
				P_mat(1,1) = -1;
				P_mat(2,1) = 1;
				lattice_vec = P_mat*lattice_vec;
				lattice_mat = P_mat*lattice_mat;
				sum_vec_neg = - (lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
				flag = false;
			}
			if(lattice_mat.row(2).dot(sum_vec_neg) > 0){
				MatrixXd P_mat = MatrixXd::Identity(3,3);
				P_mat(0,2) = 1;
				P_mat(1,2) = 1;
				P_mat(2,2) = -1;
				lattice_vec = P_mat*lattice_vec;
				lattice_mat = P_mat*lattice_mat;
				sum_vec_neg = - (lattice_mat.row(0) + lattice_mat.row(1) + lattice_mat.row(2));
				flag = false;
			}
		}

		// cout << "lattice_mat" << endl;
		// cout << lattice_mat << endl;

		// cout << "lattice_vec" << endl;
		// cout << lattice_vec << endl;

		// cout << "lattice_vec @ lattice_before" << endl;
		// cout << lattice_vec * lattice_before << endl;

		MatrixXd new_vec_array = MatrixXd::Zero(7,3);
		new_vec_array.row(0) = lattice_mat.row(0);
		new_vec_array.row(1) = lattice_mat.row(1);
		new_vec_array.row(2) = lattice_mat.row(2);
		new_vec_array.row(3) = sum_vec_neg;
		new_vec_array.row(4) = lattice_mat.row(0) + lattice_mat.row(1);
		new_vec_array.row(5) = lattice_mat.row(0) + lattice_mat.row(2);
		new_vec_array.row(6) = lattice_mat.row(1) + lattice_mat.row(2);


		vector<vector<vector<double> > > lattice_candidate;
		vector<MatrixXd> P_candidate;
		for(int i=0;i<7;i++){
			for(int j=i+1;j<7;j++){
				for(int k=j+1;k<7;k++){
					MatrixXd P_new = MatrixXd::Zero(3,3);
					bool flag = true;
					vector<vector<double> > lattice_new(3, vector<double>(3,0));
					lattice_new[0] = {new_vec_array(i,0),new_vec_array(i,1),new_vec_array(i,2)};
					lattice_new[1] = {new_vec_array(j,0),new_vec_array(j,1),new_vec_array(j,2)};
					lattice_new[2] = {new_vec_array(k,0),new_vec_array(k,1),new_vec_array(k,2)};
					if(check_lattice_normalize(lattice_new)){
						P_new.row(0) = P_select(lattice_vec,i);
						P_new.row(1) = P_select(lattice_vec,j);
						P_new.row(2) = P_select(lattice_vec,k);
						P_candidate.push_back(P_new);
						lattice_candidate.push_back(lattice_new);
					}
				}
			}
		}
		double length_min = 1000000;
		for(int i = 0; i < P_candidate.size(); i++){
			vector<vector<double> > lattice_temp = lattice_candidate[i];
			double length = 0;
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){
					length += lattice_temp[j][k]*lattice_temp[j][k];
				}
			}
			length = sqrt(length);
			if(length < length_min){
				length_min = length;
				P_ans = P_candidate[i];
				lattice = lattice_candidate[i];
			}
		}

		if(!check_lattice_normalize(lattice)){
			cout << "lattice is not normalized after normalization" << endl;
		}
		// cout << "lattice is normalized" << endl;
	}

	// for(int i=0;i<3;i++){
	// 	for(int j=0;j<3;j++){
	// 		cout << lattice[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

	// double max_diameter = 0;
	// for(int j=0;j<a;j++){
	// 	double ell = 0;
	// 	for(int k=0;k<3;k++){
	// 		ell += A(k,j)*A(k,j);
	// 	}
	// 	if(max_diameter<ell){
	// 		max_diameter = ell;
	// 	}
	// }
	// max_diameter = sqrt(max_diameter);
	vector<Answer> ans_multiple;
	for(int iter=0;iter<b;iter++){
		vector<permutation_space> Ps;
		permutation_space P0;
		P0.assignment.push_back(iter);
		if(label_A[0] != label_B[iter]){
			continue;
		}
		for(int j=0;j<b;j++){
			if(j == iter){
				continue;
			}
			P0.not_assign.push_back(j);
		}
		P0.lattice_count.push_back(vector<int>(3,0));

        // for(int j=0;j<b;j++){
		// 	double b_norm = 0;
		// 	for(int k=0;k<3;k++){
		// 		b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
		// 	}
		// 	b_norm = sqrt(b_norm);
		// 	if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
		// 		P0.not_assign.push_back(j);
		// 	}
		// }
		Ps.push_back(P0);
		while(!Ps.empty()){
			permutation_space P = Ps.back();
			Ps.pop_back();
			int N = P.assignment.size();
			for(int i=0;i<P.not_assign.size();i++){
				if(label_A[N] != label_B[P.not_assign[i]]){
				    continue;
			    }
				
				// cout << endl;
				vector<vector<int> > lattice_list = find_lattice(A,B,P.not_assign[i],eps,P,lattice);
				// lattice_listを出力
				// if(iter == 19){
				// 	cout << i << " " << P.not_assign[i] << endl;
				// 	cout << "lattice_list.size() = " << lattice_list.size() << endl;
				// 	for (int i = 0; i < lattice_list.size(); i++){
				// 		for (int j = 0; j < lattice_list[i].size(); j++){
				// 			cout << lattice_list[i][j] << " ";
				// 		}
				// 		cout << endl;
				// 	}
				// }
				
				MatrixXd BB = MatrixXd::Zero(3,N+1);
				MatrixXd AA = MatrixXd::Zero(3,N+1);
				for(int j=0;j<N;j++){
					for(int k=0;k<3;k++){
						BB(k,j) = B(k,P.assignment[j]);
						AA(k,j) = A(k,j);
					}
				}
				for(int lattice_index = 0; lattice_index < lattice_list.size(); lattice_index++){
					for(int j=0;j<N;j++){
						for(int k=0;k<3;k++){
							BB(k,j) = B(k,P.assignment[j]);
							AA(k,j) = A(k,j);
							for(int l=0;l<3;l++){
								BB(k,j) = BB(k,j) - P.lattice_count[j][l]*lattice[l][k];
							}
						}
					}
					for(int k=0;k<3;k++){
						BB(k,N) = B(k,P.not_assign[i]);
						AA(k,N) = A(k,N);
						for(int l=0;l<3;l++){
							BB(k,N) = BB(k,N) - lattice_list[lattice_index][l]*lattice[l][k];
						}
					}
					MatrixXd BBg = MatrixXd::Zero(3,1);
					MatrixXd AAg = MatrixXd::Zero(3,1);
					for(int k=0;k<3;k++){
						for(int j=0;j<N+1;j++){
							BBg(k,0) = BBg(k,0) + BB(k,j);
							AAg(k,0) = AAg(k,0) + AA(k,j);
						}
						BBg(k,0) = BBg(k,0)/(N+1);
						AAg(k,0) = AAg(k,0)/(N+1);
						for(int j=0;j<N+1;j++){
							AA(k,j) = AA(k,j) - AAg(k,0);
							BB(k,j) = BB(k,j) - BBg(k,0);
						}
					}
					JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
					double R_sum = 0;
					if(permit_mirror){
						for(int k=0;k<3;k++){
							R_sum += svd.singularValues()(k,0);
						}
					}else{
						if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
							for(int k=0;k<3;k++){
								R_sum += svd.singularValues()(k,0);
							}
						}else{
							for(int k=0;k<2;k++){
								R_sum += svd.singularValues()(k,0);
							}
							R_sum -= svd.singularValues()(2,0);
						}
					}
					double pre_ans = 0;
					for(int j=0;j<N+1;j++){
						for(int k=0;k<3;k++){
							pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
						}
					}
					pre_ans -= R_sum*2;
					if(pre_ans < 0){
						pre_ans = 0;
					}
					pre_ans = sqrt(pre_ans/(a));
					// if(P.assignment[0] == 0 && P.assignment[1] == 1){
					// 	cout << "lattice_list: " << endl;
					// 	for(int j=0;j<P.lattice_count.size();j++){
					// 		for(int k=0;k<3;k++){
					// 			cout << P.lattice_count[j][k] << " ";
					// 		}
					// 		cout << endl;
					// 	}
					// 	cout << endl;
					// }
					if(pre_ans>eps){
						continue;
					}else{
						if(N+1==a){
							Answer ans;
							ans.E = pre_ans;
							MatrixXd T = MatrixXd::Zero(b,a);
							vector<int> assignment = P.assignment;
							assignment.push_back(P.not_assign[i]);
							for(int i=0;i<a;i++){
								T(assignment[i],i) = 1;
							}
							bool found = find_P(ans_multiple,T);
							if(!found){
								ans.P = T;
								vector<vector<int> > lattice_count = P.lattice_count;


								lattice_count.push_back(lattice_list[lattice_index]);

								vector<vector<int> > lattice_count_temp(a,vector<int>(3,0));
								// cout << "lattice_count: " << endl;
								for(int j=0;j<a;j++){
									for(int k=0;k<3;k++){
										// cout << lattice_count[j][k] << " ";
										for(int l=0;l<3;l++){
											lattice_count_temp[j][k] += P_ans(l,k)*lattice_count[j][l];
										}
									}
									// cout << endl;
								}
								// cout << "P_ans: " << endl;
								// cout << P_ans << endl;
								ans.lattice_count = lattice_count_temp;
								
								MatrixXd BBB = MatrixXd::Zero(3,b);
								for(int j=0;j<b;j++){
									for(int k=0;k<3;k++){
										BBB(k,j) = B(k,j) - BBg(k,0);
									}
								}
								// cout << a << b << endl;
								for(int j=0;j<a;j++){
									for(int k=0;k<3;k++){
										for(int l=0;l<3;l++){
											BBB(l,assignment[j]) -= lattice_count[j][k]*lattice[k][l];
										}
									}
								}
								ans.R = calculate_from_assignment_not_vector(T,A,BBB,permit_mirror);
								ans.t = (Ag + AAg) - ans.R*(Bg + BBg);
								
								ans_multiple.push_back(ans);
							}
						}else{
							permutation_space new_per;
							new_per.assignment = P.assignment;
							new_per.assignment.push_back(P.not_assign[i]);
							new_per.not_assign = P.not_assign;
							new_per.not_assign.erase(new_per.not_assign.begin()+i);
							vector<vector<int> > lattice_count = P.lattice_count;
							lattice_count.push_back(lattice_list[lattice_index]);
							// if(new_per.assignment[0] == 0 && new_per.assignment[1] == 1){
							// 	cout << "lattice_list[lattice_index]: ";
							// 	for(int j=0;j<3;j++){
							// 		cout << lattice_list[lattice_index][j] << " ";
							// 	}
							// 	cout << endl;
							// }
							new_per.lattice_count = lattice_count;
							Ps.push_back(new_per);
						}
					}
				}
			}
	    }
	}
	return ans_multiple;
}
