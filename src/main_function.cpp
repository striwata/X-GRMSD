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
				t(k,0) += B(k,i_b);
				t(k,0) -= A(k,i_a);
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
				// cout << BB << endl;
				// cout << AA << endl;
				// cout << pre_ans << endl;
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
				for(int jb = 0;jb<Bs.size();jb++){
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
					if (assignment.size() < aa){
						valid = false;
					}
			        MatrixXd BB = MatrixXd::Zero(3,assignment.size());
					vector<int> b_assign_temp;
					for(int j=0;j<assignment.size();j++){
						for(int k=0;k<3;k++){
							BB(k,j) = Bs_temp[jb](k,assignment[j]);
						}
						b_assign_temp.push_back(b_assign[jb][assignment[j]]);
					}
					bb_assign.push_back(b_assign_temp);
					BBs.push_back(BB);
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
						// cout << i << endl;
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
						bool Flag = calculate_lower_K_best_eps(CC,AAs,BBs,CC.mirror,ans_multiple,eps,CC.count+K,a_assign,bb_assign);
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
						P.assignment.push_back(P.not_assign[i]);
						for(int i=0;i<a;i++){
							T(P.assignment[i],i) = 1;
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
					P.assignment.push_back(P.not_assign[i]);
					for(int i=0;i<a;i++){
						T(P.assignment[i],i) = 1;
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


// Answer Go_ICP_best_different_atoms_64(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror,double k ,double l){
// 	int a = A.cols();
//     int b = B.cols();
//     vector<MatrixXd> As;
//     vector<MatrixXd> Bs;
//     Matrix_vector(A,B,As,Bs,label_A,label_B);
// 	Cube_Cube C;
// 	C.diameter_t = 0;
// 	for(int iter=0;iter<As.size();iter++){
// 		for(int i=0;i<As[iter].cols();i++){
// 			double diameter_temp = 0;
//             for(int k=0;k<3;k++){
// 				diameter_temp += As[iter](k,i)*As[iter](k,i);
// 			}
// 			diameter_temp = sqrt(diameter_temp);
// 			if(C.diameter_t < diameter_temp){
// 				C.diameter_t = diameter_temp;
// 			}
// 		}
		
// 	}
// 	priority_queue<Cube_Cube> que;
// 	vector<double> r = {0,0,0};
// 	C.center_r = r;
// 	C.lower = 0;
// 	C.upper = INFINITY;
// 	C.count = 1;
// 	C.center_t = r;
// 	C.diameter_r = PI;
// 	C.mirror = 0;
// 	que.push(C);
// 	if(permit_mirror){
// 		Cube_Cube CC;
// 		CC.center_r = r;
// 		CC.lower = 0;
// 		CC.upper = INFINITY;
// 		CC.count = 1;
// 		CC.center_t = r;
// 		CC.diameter_r = PI;
// 		CC.diameter_t = C.diameter_t;
// 		CC.mirror = 1;
// 		que.push(CC);
// 	}
// 	MatrixXd initial_matrix = make_rotation(r,0);
// 	Answer pre = local_ICP_diff_atom(As,Bs,initial_matrix,1,r);
// 	while(!que.empty()){
// 		Cube_Cube C = que.top();
// 		que.pop();
// 		vector<Cube_Cube> Cs = divide_divide(C);
// 		cout << C.lower << endl;
// 		cout << C.diameter_r << endl;
// 		cout << C.diameter_t << endl;
// 		cout << pre.E << endl;
// 		for(int j=0;j<64;j++){
// 			Cube_Cube CC = Cs[j];
// 			if(!check_in_sphere_64(CC)){
// 				continue;
// 			}
// 			Answer E_upper = local_ICP_diff_atom(As,Bs,make_rotation(CC.center_r,CC.mirror),CC.mirror,CC.center_t);
// 			CC.upper = E_upper.E;
// 			if(E_upper.E<pre.E){
// 				pre = E_upper;
// 			}
//             int KK = int(k*CC.count +l+1-k+0.01);
// 			bool fin = calculate_lower_K_best_64(CC,As,Bs,permit_mirror,pre,CC.lower,KK);
// 			if(fin){
// 				continue;
// 			}
// 			que.push(CC);
// 		}
// 	}
// 	pre.E = sqrt(pre.E/a);
// 	return pre;
	
// }



// Answer permutation_best_same_atom(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
//             Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	vector<MatrixXd> As;
// 	vector<MatrixXd> Bs;
// 	Matrix_vector(A,B,As,Bs,label_A,label_B);
// 	Answer pre;
// 	pre.E = INFINITY;
// 	priority_queue<permutation_space> Ps;
// 	permutation_space P0;
// 	P0.assignment = {};
// 	P0.not_assign = {};
// 	P0.E = 0;
// 	for(int i=0;i<b;i++){
// 		P0.not_assign.push_back(i);
// 	}
// 	Ps.push(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.top();
// 		Ps.pop();
// 		int N = P.assignment.size();
// 		for(int i=0;i<P.not_assign.size();i++){
// 			if(label_A[N] != label_B[P.not_assign[i]]){
// 				continue;
// 			}
// 			count += 1;
// 			MatrixXd BBB = B;
//             MatrixXd BB = MatrixXd::Zero(3,N+1);
// 			MatrixXd AA = MatrixXd::Zero(3,N+1);
// 			MatrixXd AAg = MatrixXd::Zero(3,1);
// 			MatrixXd BBg = MatrixXd::Zero(3,1);
// 			for(int j=0;j<N;j++){
//                 for(int k=0;k<3;k++){
// 					BB(k,j) = B(k,P.assignment[j]);
// 					AA(k,j) = A(k,j);
// 				}
// 			}
// 			for(int k=0;k<3;k++){
//                 BB(k,N) = B(k,P.not_assign[i]);
// 			    AA(k,N) = A(k,N);
// 			}
// 			MatrixXd R;
// 			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			if(permit_mirror){
//                  R = svd.matrixV()*svd.matrixU().transpose();
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 				    R = svd.matrixV()*svd.matrixU().transpose();
// 			    }else{
// 					MatrixXd umeyama = MatrixXd::Identity(3,3);
// 			        umeyama(2,2) = -1;
// 				    R = svd.matrixV()*umeyama*svd.matrixU().transpose();
// 				}
// 			}
// 			double pre_ans = 0;
// 			MatrixXd RBB = R*BB;
// 			for(int j=0;j<N+1;j++){
// 				for(int k=0;k<3;k++){
//                     pre_ans += (AA(k,j)-RBB(k,j))*(AA(k,j)-RBB(k,j));
// 				}
// 			}
// 			if(pre_ans>pre.E){
// 				// cout << BB << endl;
// 				// cout << AA << endl;
// 				// cout << pre_ans << endl;
// 				continue;
// 			}else{
// 				Answer ans = local_ICP(As,Bs,R,permit_mirror);
// 				if(ans.E<pre.E){
// 					pre = ans;
// 				}
// 				if(N+1==a){
// 					Answer ans;
// 					ans.E = pre_ans;
// 					if(ans.E<pre.E){
// 					    pre = ans;
// 				    }
// 				}else{
// 					permutation_space new_per;
// 					new_per.E = pre_ans;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.not_assign[i]);
// 					new_per.not_assign = P.not_assign;
// 					new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 					Ps.push(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	pre.E = sqrt(pre.E/a);
// 	return pre;
// }

// Answer permutation_best_same_atom_opt_fast(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
//             Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	vector<MatrixXd> As;
// 	vector<MatrixXd> Bs;
// 	Matrix_vector(A,B,As,Bs,label_A,label_B);
// 	Answer pre;
// 	pre.E = INFINITY;
// 	vector<permutation_space> Ps;
// 	permutation_space P0;
// 	P0.assignment = {};
// 	P0.not_assign = {};
// 	P0.E = 0;
// 	for(int i=0;i<b;i++){
// 		P0.not_assign.push_back(i);
// 	}
// 	Ps.push_back(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.back();
// 		Ps.pop_back();
// 		int N = P.assignment.size();
// 		for(int i=0;i<P.not_assign.size();i++){
// 			if(label_A[N] != label_B[P.not_assign[i]]){
// 				continue;
// 			}
// 			count += 1;
// 			MatrixXd BBB = B;
//             MatrixXd BB = MatrixXd::Zero(3,N+1);
// 			MatrixXd AA = MatrixXd::Zero(3,N+1);
// 			MatrixXd AAg = MatrixXd::Zero(3,1);
// 			MatrixXd BBg = MatrixXd::Zero(3,1);
// 			for(int j=0;j<N;j++){
//                 for(int k=0;k<3;k++){
// 					BB(k,j) = B(k,P.assignment[j]);
// 					AA(k,j) = A(k,j);
// 				}
// 			}
// 			for(int k=0;k<3;k++){
//                 BB(k,N) = B(k,P.not_assign[i]);
// 			    AA(k,N) = A(k,N);
// 			}
// 			MatrixXd R;
// 			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			if(permit_mirror){
//                  R = svd.matrixV()*svd.matrixU().transpose();
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 				    R = svd.matrixV()*svd.matrixU().transpose();
// 			    }else{
// 					MatrixXd umeyama = MatrixXd::Identity(3,3);
// 			        umeyama(2,2) = -1;
// 				    R = svd.matrixV()*umeyama*svd.matrixU().transpose();
// 				}
// 			}
// 			double pre_ans = 0;
// 			MatrixXd RBB = R*BB;
// 			for(int j=0;j<N+1;j++){
// 				for(int k=0;k<3;k++){
//                     pre_ans += (AA(k,j)-RBB(k,j))*(AA(k,j)-RBB(k,j));
// 				}
// 			}
// 			if(pre_ans>pre.E){
// 				// cout << BB << endl;
// 				// cout << AA << endl;
// 				// cout << pre_ans << endl;
// 				continue;
// 			}else{
// 				// Answer ans = local_ICP(As,Bs,R,permit_mirror);
// 				// if(ans.E<pre.E){
// 				// 	pre = ans;
// 				// }
// 				if(N+1==a){
// 					Answer ans;
// 					ans.E = pre_ans;
// 					if(ans.E<pre.E){
// 					    pre = ans;
// 				    }
// 				}else{
// 					permutation_space new_per;
// 					new_per.E = pre_ans;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.not_assign[i]);
// 					new_per.not_assign = P.not_assign;
// 					new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 					Ps.push_back(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	pre.E = sqrt(pre.E/a);
// 	return pre;
// }

// vector<Answer> permutation_best_multiple_like_thesis_apply(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
//             Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	vector<Answer> ans_multiple;
// 	vector<permutation_space> Ps;
// 	permutation_space P0;
// 	P0.assignment = {};
// 	P0.assign_candidate = {};
// 	for(int i=0;i<1;i++){
// 		vector<int> assignment;
// 		for(int j=0;j<b;j++){
//             assignment.push_back(j);
// 		}
// 		P0.assign_candidate.push_back(assignment);
// 	}
// 	Ps.push_back(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.back();
// 		Ps.pop_back();
// 		int N = P.assignment.size();
// 		if(P.assign_candidate.size() == 0){
// 			continue;
// 		}
// 		for(int i=0;i<P.assign_candidate[0].size();i++){
// 			if(label_A[N] != label_B[P.assign_candidate[0][i]]){
// 				continue;
// 			}
// 			count += 1;
//             MatrixXd BB = MatrixXd::Zero(3,N+1);
// 			MatrixXd AA = MatrixXd::Zero(3,N+1);
// 			for(int j=0;j<N;j++){
//                 for(int k=0;k<3;k++){
// 					BB(k,j) = B(k,P.assignment[j]);
// 					AA(k,j) = A(k,j);
// 				}
// 			}
// 			for(int k=0;k<3;k++){
//                 BB(k,N) = B(k,P.assign_candidate[0][i]);
// 			    AA(k,N) = A(k,N);
// 			}
// 			vector<double> BBg(3);
// 			vector<double> AAg(3);
// 			for(int k=0;k<3;k++){
// 				AAg[k] = 0;
// 				BBg[k] = 0;
				
// 				for(int j=0;j<N+1;j++){
// 					BBg[k] = BBg[k] + BB(k,j);
// 					AAg[k] = AAg[k] + AA(k,j);
// 				}
// 				BBg[k] = BBg[k]/(N+1);
// 				AAg[k] = AAg[k]/(N+1);
// 				for(int j=0;j<N+1;j++){
// 					AA(k,j) = AA(k,j) - AAg[k];
// 					BB(k,j) = BB(k,j) - BBg[k];
// 				}
// 			}
			
// 			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			double R_sum = 0;
// 			if(permit_mirror){
// 				for(int k=0;k<3;k++){
// 					R_sum += svd.singularValues()(k,0);
// 				}
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 					for(int k=0;k<3;k++){
// 					    R_sum += svd.singularValues()(k,0);
// 				    }
// 			    }else{
// 					for(int k=0;k<2;k++){
// 					    R_sum += svd.singularValues()(k,0);
// 				    }
// 					R_sum -= svd.singularValues()(2,0);
// 				}
// 			}
// 			double pre_ans = 0;
// 			for(int j=0;j<N+1;j++){
// 				for(int k=0;k<3;k++){
//                     pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 				}
// 			}
// 			pre_ans -= R_sum*2;
// 			pre_ans = sqrt(pre_ans/(a));
			
// 			if(pre_ans>eps){
// 				continue;
// 			}else{
// 				if(N+1==a){
// 					Answer ans;
// 					ans.E = pre_ans;
// 					ans_multiple.push_back(ans);
// 				}else{
// 					double max_diameter = 0;
// 					for(int j=0;j<a;j++){
// 						double ell = 0;
// 						for(int k=0;k<3;k++){
// 							ell += A(k,j)*A(k,j);
// 						}
// 						if(max_diameter<ell){
// 							max_diameter = ell;
// 						}
// 					}
// 					max_diameter = sqrt(max_diameter);
// 					permutation_space new_per;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.assign_candidate[0][i]);
// 					vector<vector<int> > assign_candidate(1);
// 					for(int k2=0;k2<P.assign_candidate[0].size();k2++){
// 						if(i==k2){
// 							continue;
// 						}
// 						double distance_new_old = 0;
// 						for(int k=0;k<3;k++){
// 							distance_new_old += (B(k,P.assign_candidate[0][k2]) - B(k,P.assign_candidate[0][i]))*(B(k,P.assign_candidate[0][k2]) - B(k,P.assign_candidate[0][i]));
// 						}
// 						if(sqrt(distance_new_old)<2*(eps*sqrt(a) + max_diameter)){
// 							assign_candidate[0].push_back(P.assign_candidate[0][k2]);
// 						}
// 					}
// 					new_per.assign_candidate = assign_candidate;
// 					Ps.push_back(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	return ans_multiple;
// }


// vector<Answer> permutation_best_multiple_like_thesis_apply_fastest_a(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
//             Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	double max_diameter = 0;
// 	for(int j=0;j<a;j++){
// 		double ell = 0;
// 		for(int k=0;k<3;k++){
// 			ell += A(k,j)*A(k,j);
// 		}
// 		if(max_diameter<ell){
// 			max_diameter = ell;
// 		}
// 	}
// 	max_diameter = sqrt(max_diameter);
// 	vector<Answer> ans_multiple;
// 	vector<permutation_space> Ps;
// 	permutation_space P0;
// 	P0.assignment = {};
// 	P0.assign_candidate = {};
// 	vector<int> assignment;
// 	for(int j=0;j<b;j++){
// 		assignment.push_back(j);
// 	}
// 	P0.assign_candidate.push_back(assignment);
// 	P0.not_assign = assignment;
// 	P0.ABt = MatrixXd::Zero(3,3);
// 	P0.Ag = MatrixXd::Zero(3,1);
// 	P0.Bg = MatrixXd::Zero(3,1);
// 	P0.AtA = 0;
// 	P0.BtB = 0;
// 	Ps.push_back(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.back();
// 		Ps.pop_back();
// 		int N = P.assignment.size();
// 		if(P.not_assign.size() == 0){
// 			continue;
// 		}
// 		for(int i=0;i<P.not_assign.size();i++){
// 			if(label_A[N] != label_B[P.not_assign[i]]){
// 				continue;
// 			}
// 			count += 1;
// 			MatrixXd Ag = MatrixXd::Zero(3,1);
// 			MatrixXd Bg = MatrixXd::Zero(3,1);
// 			MatrixXd Ag_diff = MatrixXd::Zero(3,1);
// 			MatrixXd Bg_diff = MatrixXd::Zero(3,1);
// 			for(int k=0;k<3;k++){
// 				Ag(k,0) = (P.Ag(k,0)*N + A(k,N))/(N+1);
// 				Ag_diff(k,0) = Ag(k,0)-P.Ag(k,0);
// 				Bg(k,0) = (P.Bg(k,0)*N + B(k,P.not_assign[i]))/(N+1);
// 				Bg_diff(k,0) = Bg(k,0)-P.Bg(k,0);
// 			}

// 			double AtA = P.AtA + (N*Ag_diff.transpose()*Ag_diff)(0,0);
// 			double BtB = P.BtB + (N*Bg_diff.transpose()*Bg_diff)(0,0);
// 			MatrixXd ABt = P.ABt + N*Ag_diff*Bg_diff.transpose();
// 			for(int k=0;k<3;k++){
// 				AtA += (A(k,N)-Ag(k,0))*(A(k,N)-Ag(k,0));
// 				BtB += (B(k,P.not_assign[i])-Bg(k,0))*(B(k,P.not_assign[i])-Bg(k,0));
// 			}
// 			for(int k1=0;k1<3;k1++){
// 				for(int k2=0;k2<3;k2++){
// 					ABt(k1,k2) += (A(k1,N)-Ag(k1,0))*(B(k2,P.not_assign[i])-Bg(k2,0));
// 				}
// 			}
// 			JacobiSVD< Matrix<double, 3, 3> > svd(ABt.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			double R_sum = 0;
// 			if(permit_mirror){
// 				for(int k=0;k<3;k++){
// 					R_sum += svd.singularValues()(k,0);
// 				}
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 					for(int k=0;k<3;k++){
// 					    R_sum += svd.singularValues()(k,0);
// 				    }
// 			    }else{
// 					for(int k=0;k<2;k++){
// 					    R_sum += svd.singularValues()(k,0);
// 				    }
// 					R_sum -= svd.singularValues()(2,0);
// 				}
// 			}
// 			double pre_ans = AtA + BtB;
// 			// for(int j=0;j<N+1;j++){
// 			// 	for(int k=0;k<3;k++){
//             //         pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 			// 	}
// 			// }
// 			pre_ans -= R_sum*2;
// 			pre_ans = sqrt(pre_ans/(a));
// 			double a_b_norm = 0;
// 			// for(int j=0;j<N+1;j++){
// 			// 	double a_norm = 0;
// 			// 	double b_norm = 0;
// 			// 	for(int k=0;k<3;k++){
// 			// 		a_norm += AA(k,j)*AA(k,j);
// 			// 		b_norm += BB(k,j)*BB(k,j);
// 			//     }
// 			// 	a_norm = sqrt(a_norm);
// 			// 	b_norm = sqrt(b_norm);
// 			// 	a_b_norm += (b_norm - a_norm) * (b_norm - a_norm);
// 			// }
// 			if(pre_ans>eps){
// 				continue;
// 			}else{
// 				if(N+1==a){
// 					Answer ans;
// 					ans.E = pre_ans;
// 					ans_multiple.push_back(ans);
// 				}else{
// 					permutation_space new_per;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.not_assign[i]);
// 					vector<vector<int> > assign_candidate(1);
// 					vector<int> not_assign;
// 					for(int k1=0;k1<P.not_assign.size();k1++){
// 						if(i==k1){
// 							continue;
// 						}
// 						double distance_new_old = 0;
// 						for(int k=0;k<3;k++){
// 							distance_new_old += (B(k,P.not_assign[k1]) - B(k,P.not_assign[i]))*(B(k,P.not_assign[k1]) - B(k,P.not_assign[i]));
// 						}
// 						if(sqrt(distance_new_old)<2*(eps*sqrt(a) + max_diameter)){
// 							not_assign.push_back(P.not_assign[k1]);
// 						}
// 						if(label_A[N+1] != label_B[P.not_assign[k1]]){
// 							continue;
// 						}
// 						// double value = 0;
// 						// value += a_b_norm;
// 						// double a_norm = 0;
// 						// double b_norm = 0;
// 						// for(int k=0;k<3;k++){
// 						// 	a_norm += (A(k,N+1)-AAg[k])*(A(k,N+1)-AAg[k]);
// 						// 	b_norm += (B(k,P.not_assign[k1])-BBg[k])*(B(k,P.not_assign[k1])-BBg[k]);
// 						// }
// 						// a_norm = sqrt(a_norm);
// 						// b_norm = sqrt(b_norm);
// 						// value += (b_norm - a_norm)*(b_norm - a_norm)*(N+1)/(N+2);
// 						// if(value<a*eps*eps){
// 						// 	assign_candidate[0].push_back(P.not_assign[k1]);
// 						// }

// 					}
// 					new_per.ABt = ABt;
// 					new_per.Ag = Ag;
// 					new_per.Bg = Bg;
// 					new_per.AtA = AtA;
// 					new_per.BtB = BtB;
// 					new_per.assign_candidate = assign_candidate;
// 					new_per.not_assign = not_assign;
// 					Ps.push_back(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	return ans_multiple;
// }

// vector<Answer> permutation_best_multiple_like_thesis_apply_fastest_b(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
// 			Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	double max_diameter = 0;
// 	for(int j=0;j<a;j++){
// 		double ell = 0;
// 		for(int k=0;k<3;k++){
// 			ell += A(k,j)*A(k,j);
// 		}
// 		if(max_diameter<ell){
// 			max_diameter = ell;
// 		}
// 	}
// 	max_diameter = sqrt(max_diameter);
// 	vector<Answer> ans_multiple;
// 	for(int i=0;i<b;i++){
// 		vector<permutation_space> Ps ={};
// 		permutation_space P0;
//         for(int j=0;j<b;j++){
// 			double b_norm = 0;
// 			for(int k=0;k<3;k++){
// 				b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
// 			}
// 			b_norm = sqrt(b_norm);
// 			if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
// 				P0.not_assign.push_back(j);
// 			}
// 		}
// 		P0.assignment = {};
// 		P0.assignment.push_back(i);
// 		Ps.push_back(P0);
// 		while(!Ps.empty()){
// 			permutation_space P = Ps.back();
// 			Ps.pop_back();
// 			int N = P.assignment.size();
// 			double a_b_norm = 0;
// 			for(int j=0;j<N;j++){
// 				double a_norm = 0;
// 				double b_norm = 0;
// 				for(int k=0;k<3;k++){
// 					a_norm += A(k,j)*A(k,j);
// 					b_norm += B(k,P.assignment[j])*B(k,P.assignment[j]);
// 				}
// 				a_norm = sqrt(a_norm);
// 				b_norm = sqrt(b_norm);
// 				a_b_norm += (b_norm - a_norm) * (b_norm - a_norm);
				
// 			}
// 			for(int i=0;i<P.not_assign.size();i++){
// 				if(label_A[N] != label_B[P.not_assign[i]]){
// 				    continue;
// 			    }
// 				MatrixXd BB = MatrixXd::Zero(3,N+1);
// 				MatrixXd AA = MatrixXd::Zero(3,N+1);
// 				for(int j=0;j<N;j++){
// 					for(int k=0;k<3;k++){
// 						BB(k,j) = B(k,P.assignment[j]);
// 						AA(k,j) = A(k,j);
// 					}
// 				}
// 				for(int k=0;k<3;k++){
// 					BB(k,N) = B(k,P.not_assign[i]);
// 					AA(k,N) = A(k,N);
// 				}
// 				vector<double> BBg(3);
// 				vector<double> AAg(3);
// 				for(int k=0;k<3;k++){
					
// 					for(int j=0;j<N+1;j++){
// 						BBg[k] = BBg[k] + BB(k,j);
// 						AAg[k] = AAg[k] + AA(k,j);
// 					}
// 					BBg[k] = BBg[k]/(N+1);
// 					AAg[k] = AAg[k]/(N+1);
// 					for(int j=0;j<N+1;j++){
// 						AA(k,j) = AA(k,j) - AAg[k];
// 						BB(k,j) = BB(k,j) - BBg[k];
// 					}
// 				}
// 				double value = 0;
// 				value += a_b_norm;
// 				double a_norm = 0;
// 				double b_norm = 0;
// 				for(int k=0;k<3;k++){
// 					a_norm += AAg[k]*AAg[k];
// 					b_norm += BBg[k]*BBg[k];
// 				}
// 				a_norm = sqrt(a_norm);
// 				b_norm = sqrt(b_norm);
// 				value -= (b_norm - a_norm)*(b_norm - a_norm)*(N+1);
// 				double a_norm_2 = 0;
// 				double b_norm_2 = 0;
// 				for(int k=0;k<3;k++){
// 					a_norm_2 += A(k,N)*A(k,N);
// 					b_norm_2 += B(k,P.not_assign[i])*B(k,P.not_assign[i]);
// 				}
// 				a_norm_2 = sqrt(a_norm_2);
// 				b_norm_2 = sqrt(b_norm_2);
// 				value += (b_norm_2-a_norm_2)*(b_norm_2-a_norm_2);
// 				if(value>a*eps*eps){
// 					continue;
// 				}
// 				count += 1;
// 				JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 				double R_sum = 0;
// 				if(permit_mirror){
// 					for(int k=0;k<3;k++){
// 						R_sum += svd.singularValues()(k,0);
// 					}
// 				}else{
// 					if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 						for(int k=0;k<3;k++){
// 							R_sum += svd.singularValues()(k,0);
// 						}
// 					}else{
// 						for(int k=0;k<2;k++){
// 							R_sum += svd.singularValues()(k,0);
// 						}
// 						R_sum -= svd.singularValues()(2,0);
// 					}
// 				}
// 				double pre_ans = 0;
// 				for(int j=0;j<N+1;j++){
// 					for(int k=0;k<3;k++){
// 						pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 					}
// 				}
// 				pre_ans -= R_sum*2;
// 				if(pre_ans < 0){
// 					pre_ans = 0;
// 				}
// 				pre_ans = sqrt(pre_ans/(a));
// 				if(pre_ans>eps){
// 					continue;
// 				}else{
// 					if(N+1==a){
// 						bool found = false;
// 						for(int i=0;i<ans_multiple.size();i++){
// 							if(abs(ans_multiple[i].E  - pre_ans)<1e-10){
// 								found = true;
// 								continue;
// 							}                
// 						}
// 						if(!found){
// 							permutation_space new_per;
// 							vector<MatrixXd> Ts;
// 							MatrixXd T = MatrixXd::Zero(b,a); 
// 							new_per.assignment = P.assignment;
// 						    new_per.assignment.push_back(P.not_assign[i]);
// 							for(int i=0;i<a;i++){
// 								T(new_per.assignment[i],i) = 1;
// 							}
// 							Ts.push_back(T);
// 							Answer ans;
// 							ans.E = pre_ans;
// 							ans.Ts = Ts;
// 							ans_multiple.push_back(ans);
// 						}
// 					}else{
// 						permutation_space new_per;
// 						new_per.assignment = P.assignment;
// 						new_per.assignment.push_back(P.not_assign[i]);
// 						new_per.not_assign = P.not_assign;
// 						new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 						Ps.push_back(new_per);
// 					}
// 				}
// 			}
// 	    }
// 	}
// 	cout << count << endl;
// 	return ans_multiple;
// }

// vector<Answer> permutation_best_multiple_original(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
//             Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	vector<Answer> ans_multiple;
// 	vector<permutation_space> Ps;
// 	permutation_space P0;
// 	P0.assignment = {};
// 	P0.assign_candidate = {};
// 	for(int i=0;i<a;i++){
// 		vector<int> assignment;
// 		for(int j=0;j<b;j++){
//             assignment.push_back(j);
// 		}
// 		P0.assign_candidate.push_back(assignment);
// 	}
// 	Ps.push_back(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.back();
// 		Ps.pop_back();
// 		int N = P.assignment.size();
// 		if(P.assign_candidate[0].size() == 0){
// 			continue;
// 		}
// 		for(int i=0;i<P.assign_candidate[0].size();i++){
// 			if(label_A[N] != label_B[P.assign_candidate[0][i]]){
// 				continue;
// 			}
// 			count += 1;
//             MatrixXd BB = MatrixXd::Zero(3,N+1);
// 			MatrixXd AA = MatrixXd::Zero(3,N+1);
// 			for(int j=0;j<N;j++){
//                 for(int k=0;k<3;k++){
// 					BB(k,j) = B(k,P.assignment[j]);
// 					AA(k,j) = A(k,j);
// 				}
// 			}
// 			for(int k=0;k<3;k++){
//                 BB(k,N) = B(k,P.assign_candidate[0][i]);
// 			    AA(k,N) = A(k,N);
// 			}
// 			vector<double> BBg(3);
// 			vector<double> AAg(3);
// 			for(int k=0;k<3;k++){
// 				AAg[k] = 0;
// 				BBg[k] = 0;
				
// 				for(int j=0;j<N+1;j++){
// 					BBg[k] = BBg[k] + BB(k,j);
// 					AAg[k] = AAg[k] + AA(k,j);
// 				}
// 				BBg[k] = BBg[k]/(N+1);
// 				AAg[k] = AAg[k]/(N+1);
// 				for(int j=0;j<N+1;j++){
// 					AA(k,j) = AA(k,j) - AAg[k];
// 					BB(k,j) = BB(k,j) - BBg[k];
// 				}
// 			}
			
// 			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			double R_sum = 0;
// 			if(permit_mirror){
// 				for(int k=0;k<3;k++){
// 					R_sum += svd.singularValues()(k,0);
// 				}
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 					for(int k=0;k<3;k++){
// 					    R_sum += svd.singularValues()(k,0);
// 				    }
// 			    }else{
// 					for(int k=0;k<2;k++){
// 					    R_sum += svd.singularValues()(k,0);
// 				    }
// 					R_sum -= svd.singularValues()(2,0);
// 				}
// 			}
// 			double pre_ans = 0;
// 			for(int j=0;j<N+1;j++){
// 				for(int k=0;k<3;k++){
//                     pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 				}
// 			}
// 			pre_ans -= R_sum*2;
// 			pre_ans = sqrt(pre_ans/(a));
// 			double a_b_norm = 0;
// 			for(int j=0;j<N+1;j++){
// 				double a_norm = 0;
// 				double b_norm = 0;
// 				for(int k=0;k<3;k++){
// 					a_norm += AA(k,j)*AA(k,j);
// 					b_norm += BB(k,j)*BB(k,j);
// 			    }
// 				a_norm = sqrt(a_norm);
// 				b_norm = sqrt(b_norm);
// 				a_b_norm += (b_norm - a_norm) * (b_norm - a_norm);
// 			}
// 			if(pre_ans>eps){
// 				continue;
// 			}else{
// 				if(N+1==a){
// 					Answer ans;
// 					ans.E = pre_ans;
// 					ans_multiple.push_back(ans);
// 				}else{
// 					permutation_space new_per;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.assign_candidate[0][i]);
// 					vector<vector<int> > assign_candidate(a-N-1);
// 					for(int k1=1;k1<P.assign_candidate.size();k1++){
// 						for(int k2=0;k2<P.assign_candidate[k1].size();k2++){
// 							if(label_A[k1] != label_B[P.assign_candidate[k1][k2]]){
// 								continue;
// 							}
// 							if(P.assign_candidate[0][i]==P.assign_candidate[k1][k2]){
// 								continue;
// 							}
// 							double value = 0;
// 							value += a_b_norm;
// 							double a_norm = 0;
// 							double b_norm = 0;
// 							for(int k=0;k<3;k++){
// 								a_norm += (A(k,k1+N)-AAg[k])*(A(k,k1+N)-AAg[k]);
// 								b_norm += (B(k,P.assign_candidate[k1][k2])-BBg[k])*(B(k,P.assign_candidate[k1][k2])-BBg[k]);
// 							}
// 							a_norm = sqrt(a_norm);
// 							b_norm = sqrt(b_norm);
// 							value += (b_norm - a_norm)*(b_norm - a_norm)*(N+1)/(N+2);
// 							if(value<eps*a*eps){
// 								assign_candidate[k1-1].push_back(P.assign_candidate[k1][k2]);
// 							}
// 						}
// 					}
// 					new_per.assign_candidate = assign_candidate;
// 					Ps.push_back(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	return ans_multiple;
// }


// vector<Answer> permutation_best_multiple_like_thesis_a(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
// 			Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	double max_diameter = 0;
// 	for(int j=0;j<a;j++){
// 		double ell = 0;
// 		for(int k=0;k<3;k++){
// 			ell += A(k,j)*A(k,j);
// 		}
// 		if(max_diameter<ell){
// 			max_diameter = ell;
// 		}
// 	}
// 	max_diameter = sqrt(max_diameter);
// 	vector<Answer> ans_multiple;
// 	for(int i=0;i<b;i++){
// 		vector<permutation_space> Ps ={};
// 		permutation_space P0;
//         for(int j=0;j<b;j++){
// 			double b_norm = 0;
// 			for(int k=0;k<3;k++){
// 				b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
// 			}
// 			b_norm = sqrt(b_norm);
// 			if(b_norm < 2*(eps*sqrt(a) + max_diameter) && i != j){
// 				P0.not_assign.push_back(j);
// 			}
// 		}
// 		P0.assignment = {};
// 		if(label_A[0] != label_B[i]){
// 			continue;
// 		}
// 		P0.assignment.push_back(i);
// 		Ps.push_back(P0);
// 		while(!Ps.empty()){
// 			permutation_space P = Ps.back();
// 			Ps.pop_back();
// 			int N = P.assignment.size();
// 			for(int i=0;i<P.not_assign.size();i++){
// 				count += 1;
// 				if(label_A[N] != label_B[P.not_assign[i]]){
// 				    continue;
// 			    }
// 				MatrixXd BB = MatrixXd::Zero(3,N+1);
// 				MatrixXd AA = MatrixXd::Zero(3,N+1);
// 				for(int j=0;j<N;j++){
// 					for(int k=0;k<3;k++){
// 						BB(k,j) = B(k,P.assignment[j]);
// 						AA(k,j) = A(k,j);
// 					}
// 				}
// 				for(int k=0;k<3;k++){
// 					BB(k,N) = B(k,P.not_assign[i]);
// 					AA(k,N) = A(k,N);
// 				}
// 				for(int k=0;k<3;k++){
// 					double BBg = 0;
// 					double AAg = 0;
// 					for(int j=0;j<N+1;j++){
// 						BBg = BBg + BB(k,j);
// 						AAg = AAg + AA(k,j);
// 					}
// 					BBg = BBg/(N+1);
// 					AAg = AAg/(N+1);
// 					for(int j=0;j<N+1;j++){
// 						AA(k,j) = AA(k,j) - AAg;
// 						BB(k,j) = BB(k,j) - BBg;
// 					}
// 				}
// 				JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 				double R_sum = 0;
// 				if(permit_mirror){
// 					for(int k=0;k<3;k++){
// 						R_sum += svd.singularValues()(k,0);
// 					}
// 				}else{
// 					if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 						for(int k=0;k<3;k++){
// 							R_sum += svd.singularValues()(k,0);
// 						}
// 					}else{
// 						for(int k=0;k<2;k++){
// 							R_sum += svd.singularValues()(k,0);
// 						}
// 						R_sum -= svd.singularValues()(2,0);
// 					}
// 				}
// 				double pre_ans = 0;
// 				for(int j=0;j<N+1;j++){
// 					for(int k=0;k<3;k++){
// 						pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 					}
// 				}
// 				pre_ans -= R_sum*2;
// 				if(pre_ans < 0){
// 					pre_ans = 0;
// 				}
// 				pre_ans = sqrt(pre_ans/(a));
// 				if(pre_ans>eps){
// 					continue;
// 				}else{
// 					if(N+1==a){
// 						bool found = false;
// 						for(int i=0;i<ans_multiple.size();i++){
// 							if(abs(ans_multiple[i].E  - pre_ans)<1e-10){
// 								found = true;
// 								continue;
// 							}                
// 						}
// 						if(!found){
// 							permutation_space new_per;
// 							vector<MatrixXd> Ts;
// 							MatrixXd T = MatrixXd::Zero(b,a); 
// 							new_per.assignment = P.assignment;
// 						    new_per.assignment.push_back(P.not_assign[i]);
// 							for(int i=0;i<a;i++){
// 								T(new_per.assignment[i],i) = 1;
// 							}
// 							Ts.push_back(T);
// 							Answer ans;
// 							ans.E = pre_ans;
// 							ans.Ts = Ts;
// 							ans_multiple.push_back(ans);
// 						}
// 					}else{
// 						permutation_space new_per;
// 						new_per.assignment = P.assignment;
// 						new_per.assignment.push_back(P.not_assign[i]);
// 						new_per.not_assign = P.not_assign;
// 						new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 						Ps.push_back(new_per);
// 					}
// 				}
// 			}
// 	    }
// 	}
// 	cout << count << endl;
// 	return ans_multiple;
// }

// // vector<Answer> permutation_best_multiple_like_thesis_a(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// // 	int a = A.cols();
// //     int b = B.cols();
// // 	int count = 0;
// // 	for(int i=0;i<3;i++){
// // 		double Ag = 0;
// // 		double Bg = 0;
// // 		for(int j=0;j<a;j++){
// // 			Ag = Ag + A(i,j);
// // 		}
// // 		for(int j=0;j<b;j++){
// // 			Bg = Bg + B(i,j);
// // 		}
// // 		Ag = Ag/a;
// // 		Bg = Bg/b;
// // 		for(int j=0;j<a;j++){
// // 			A(i,j) = A(i,j) - Ag;
// // 		}
// // 		for(int j=0;j<b;j++){
// // 			B(i,j) = B(i,j) - Bg;
// // 		}
// // 	}
// // 	double max_diameter = 0;
// // 	for(int j=0;j<a;j++){
// // 		double ell = 0;
// // 		for(int k=0;k<3;k++){
// // 			ell += A(k,j)*A(k,j);
// // 		}
// // 		if(max_diameter<ell){
// // 			max_diameter = ell;
// // 		}
// // 	}
// // 	max_diameter = sqrt(max_diameter);
// // 	vector<Answer> ans_multiple;
// // 	for(int i=0;i<b;i++){
// // 		vector<permutation_space> Ps ={};
// // 		permutation_space P0;
// //         for(int j=0;j<b;j++){
// // 			double b_norm = 0;
// // 			if(i==b){
// // 				continue;
// // 			}
// // 			for(int k=0;k<3;k++){
// // 				b_norm += (B(k,j) - B(k,i))*(B(k,j) - B(k,i));
// // 			}
// // 			b_norm = sqrt(b_norm);
// // 			if(b_norm < 2*(eps*sqrt(a) + max_diameter)){
// // 				P0.not_assign.push_back(j);
// // 			}
// // 		}
// // 		P0.assignment = {};
// // 		P0.assignment.push_back(i);
// // 		P0.ABt = MatrixXd::Zero(3,3);
// // 	P0.Ag = MatrixXd::Zero(3,1);
// // 	P0.Bg = MatrixXd::Zero(3,1);
// // 	P0.AtA = 0;
// // 	P0.BtB = 0;
// // 	if(P0.not_assign.size()<A.cols()){
// // 		continue;
// // 	}
// // 	Ps.push_back(P0);
// // 	while(!Ps.empty()){
// // 		permutation_space P = Ps.back();
// // 		Ps.pop_back();
// // 		int N = P.assignment.size();
// // 		if(P.not_assign.size() == 0){
// // 			continue;
// // 		}
// // 		for(int i=0;i<P.not_assign.size();i++){
// // 			if(label_A[N] != label_B[P.not_assign[i]]){
// // 				continue;
// // 			}
// // 			count += 1;
// // 			MatrixXd Ag = MatrixXd::Zero(3,1);
// // 			MatrixXd Bg = MatrixXd::Zero(3,1);
// // 			MatrixXd Ag_diff = MatrixXd::Zero(3,1);
// // 			MatrixXd Bg_diff = MatrixXd::Zero(3,1);
// // 			for(int k=0;k<3;k++){
// // 				Ag(k,0) = (P.Ag(k,0)*N + A(k,N))/(N+1);
// // 				Ag_diff(k,0) = Ag(k,0)-P.Ag(k,0);
// // 				Bg(k,0) = (P.Bg(k,0)*N + B(k,P.not_assign[i]))/(N+1);
// // 				Bg_diff(k,0) = Bg(k,0)-P.Bg(k,0);
// // 			}

// // 			double AtA = P.AtA + (N*Ag_diff.transpose()*Ag_diff)(0,0);
// // 			double BtB = P.BtB + (N*Bg_diff.transpose()*Bg_diff)(0,0);
// // 			MatrixXd ABt = P.ABt + N*Ag_diff*Bg_diff.transpose();
// // 			for(int k=0;k<3;k++){
// // 				AtA += (A(k,N)-Ag(k,0))*(A(k,N)-Ag(k,0));
// // 				BtB += (B(k,P.not_assign[i])-Bg(k,0))*(B(k,P.not_assign[i])-Bg(k,0));
// // 			}
// // 			for(int k1=0;k1<3;k1++){
// // 				for(int k2=0;k2<3;k2++){
// // 					ABt(k1,k2) += (A(k1,N)-Ag(k1,0))*(B(k2,P.not_assign[i])-Bg(k2,0));
// // 				}
// // 			}
// // 			JacobiSVD< Matrix<double, 3, 3> > svd(ABt.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// // 			double R_sum = 0;
// // 			if(permit_mirror){
// // 				for(int k=0;k<3;k++){
// // 					R_sum += svd.singularValues()(k,0);
// // 				}
// // 			}else{
// // 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// // 					for(int k=0;k<3;k++){
// // 					    R_sum += svd.singularValues()(k,0);
// // 				    }
// // 			    }else{
// // 					for(int k=0;k<2;k++){
// // 					    R_sum += svd.singularValues()(k,0);
// // 				    }
// // 					R_sum -= svd.singularValues()(2,0);
// // 				}
// // 			}
// // 			double pre_ans = AtA + BtB;
// // 			// for(int j=0;j<N+1;j++){
// // 			// 	for(int k=0;k<3;k++){
// //             //         pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// // 			// 	}
// // 			// }
// // 			pre_ans -= R_sum*2;
// // 			pre_ans = sqrt(pre_ans/(a));
// // 			double a_b_norm = 0;
// // 			// for(int j=0;j<N+1;j++){
// // 			// 	double a_norm = 0;
// // 			// 	double b_norm = 0;
// // 			// 	for(int k=0;k<3;k++){
// // 			// 		a_norm += AA(k,j)*AA(k,j);
// // 			// 		b_norm += BB(k,j)*BB(k,j);
// // 			//     }
// // 			// 	a_norm = sqrt(a_norm);
// // 			// 	b_norm = sqrt(b_norm);
// // 			// 	a_b_norm += (b_norm - a_norm) * (b_norm - a_norm);
// // 			// }
// // 			if(pre_ans>eps){
// // 				continue;
// // 				}else{
// // 					if(N+1==a){
// // 						bool found = false;
// // 						for(int i=0;i<ans_multiple.size();i++){
// // 							if(abs(ans_multiple[i].E  - pre_ans)<1e-10){
// // 								found = true;
// // 								continue;
// // 							}                
// // 						}
// // 						if(!found){
// // 							permutation_space new_per;
// // 							vector<MatrixXd> Ts;
// // 							MatrixXd T = MatrixXd::Zero(b,a); 
// // 							new_per.assignment = P.assignment;
// // 						    new_per.assignment.push_back(P.not_assign[i]);
// // 							for(int i=0;i<a;i++){
// // 								T(new_per.assignment[i],i) = 1;
// // 							}
// // 							Ts.push_back(T);
// // 							Answer ans;
// // 							ans.E = pre_ans;
// // 							ans.Ts = Ts;
// // 							ans_multiple.push_back(ans);
// // 						}
// // 					}else{
// // 						permutation_space new_per;
// // 						new_per.assignment = P.assignment;
// // 						new_per.assignment.push_back(P.not_assign[i]);
// // 						new_per.not_assign = P.not_assign;
// // 						new_per.not_assign.erase(new_per.not_assign.begin()+i);
// // 						Ps.push_back(new_per);
// // 					}
// // 				}
// // 			}
// // 	    }
// // 	}
// // 	cout << count << endl;
// // 	return ans_multiple;
// // }

// Answer permutation_best_multiple_like_thesis_opt(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
// 			Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	double max_diameter = 0;
// 	for(int j=0;j<a;j++){
// 		double ell = 0;
// 		for(int k=0;k<3;k++){
// 			ell += A(k,j)*A(k,j);
// 		}
// 		if(max_diameter<ell){
// 			max_diameter = ell;
// 		}
// 	}
// 	max_diameter = sqrt(max_diameter);
// 	Answer ans;
// 	ans.E = INFINITY;
// 	vector<permutation_space> Ps ={};
// 	permutation_space P0;
// 	for(int j=0;j<b;j++){
// 		P0.not_assign.push_back(j);
// 	}
// 	P0.assignment = {};
// 	Ps.push_back(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.back();
// 		Ps.pop_back();
// 		int N = P.assignment.size();
// 		for(int i=0;i<P.not_assign.size();i++){
// 			count += 1;
// 			if(label_A[N] != label_B[P.not_assign[i]]){
// 				continue;
// 			}
// 			MatrixXd BB = MatrixXd::Zero(3,N+1);
// 			MatrixXd AA = MatrixXd::Zero(3,N+1);
// 			for(int j=0;j<N;j++){
// 				for(int k=0;k<3;k++){
// 					BB(k,j) = B(k,P.assignment[j]);
// 					AA(k,j) = A(k,j);
// 				}
// 			}
// 			for(int k=0;k<3;k++){
// 				BB(k,N) = B(k,P.not_assign[i]);
// 				AA(k,N) = A(k,N);
// 			}
// 			for(int k=0;k<3;k++){
// 				double BBg = 0;
// 				double AAg = 0;
// 				for(int j=0;j<N+1;j++){
// 					BBg = BBg + BB(k,j);
// 					AAg = AAg + AA(k,j);
// 				}
// 				BBg = BBg/(N+1);
// 				AAg = AAg/(N+1);
// 				for(int j=0;j<N+1;j++){
// 					AA(k,j) = AA(k,j) - AAg;
// 					BB(k,j) = BB(k,j) - BBg;
// 				}
// 			}
// 			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			double R_sum = 0;
// 			if(permit_mirror){
// 				for(int k=0;k<3;k++){
// 					R_sum += svd.singularValues()(k,0);
// 				}
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 					for(int k=0;k<3;k++){
// 						R_sum += svd.singularValues()(k,0);
// 					}
// 				}else{
// 					for(int k=0;k<2;k++){
// 						R_sum += svd.singularValues()(k,0);
// 					}
// 					R_sum -= svd.singularValues()(2,0);
// 				}
// 			}
// 			double pre_ans = 0;
// 			for(int j=0;j<N+1;j++){
// 				for(int k=0;k<3;k++){
// 					pre_ans += AA(k,j)*AA(k,j) + BB(k,j)*BB(k,j);
// 				}
// 			}
// 			pre_ans -= R_sum*2;
// 			if(pre_ans < 0){
// 				pre_ans = 0;
// 			}
// 			pre_ans = sqrt(pre_ans/(a));
// 			if(pre_ans>ans.E){
// 				continue;
// 			}else{
// 				if(N+1==a){
// 					ans.E = pre_ans;
// 					permutation_space new_per;
// 					vector<MatrixXd> Ts;
// 					MatrixXd T = MatrixXd::Zero(b,a); 
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.not_assign[i]);
// 					for(int i=0;i<a;i++){
// 						T(new_per.assignment[i],i) = 1;
// 					}
// 					Ts.push_back(T);
// 					Answer ans;
// 					ans.E = pre_ans;
// 					ans.Ts = Ts;
// 				}else{
// 					permutation_space new_per;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.not_assign[i]);
// 					new_per.not_assign = P.not_assign;
// 					new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 					Ps.push_back(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	return ans;
// }

// vector<Answer> permutation_best_multiple(MatrixXd A,MatrixXd B,vector<int> label_A,vector<int> label_B,double eps,bool permit_mirror){
// 	int a = A.cols();
//     int b = B.cols();
// 	int count = 0;
// 	for(int i=0;i<3;i++){
// 		double Ag = 0;
// 		double Bg = 0;
// 		for(int j=0;j<a;j++){
// 			Ag = Ag + A(i,j);
// 		}
// 		for(int j=0;j<b;j++){
//             Bg = Bg + B(i,j);
// 		}
// 		Ag = Ag/a;
// 		Bg = Bg/b;
// 		for(int j=0;j<a;j++){
// 			A(i,j) = A(i,j) - Ag;
// 		}
// 		for(int j=0;j<b;j++){
// 			B(i,j) = B(i,j) - Bg;
// 		}
// 	}
// 	vector<Answer> ans_multiple;
// 	vector<permutation_space> Ps;
// 	permutation_space P0;
// 	P0.assignment = {};
// 	P0.not_assign = {};
// 	for(int i=0;i<b;i++){
// 		P0.not_assign.push_back(i);
// 	}
// 	Ps.push_back(P0);
// 	while(!Ps.empty()){
// 		permutation_space P = Ps.back();
// 		Ps.pop_back();
// 		int N = P.assignment.size();
// 		for(int i=0;i<P.not_assign.size();i++){
// 			count += 1;
//             MatrixXd BB = MatrixXd::Zero(3,N+1);
// 			MatrixXd AA = MatrixXd::Zero(3,N+1);
// 			for(int j=0;j<N;j++){
//                 for(int k=0;k<3;k++){
// 					BB(k,j) = B(k,P.assignment[j]);
// 					AA(k,j) = A(k,j);
// 				}
// 			}
// 			for(int k=0;k<3;k++){
//                 BB(k,N) = B(k,P.not_assign[i]);
// 			    AA(k,N) = A(k,N);
// 			}
// 			for(int k=0;k<3;k++){
// 				double BBg = 0;
// 				double AAg = 0;
// 				for(int j=0;j<N+1;j++){
// 					BBg = BBg + BB(k,j);
// 					AAg = AAg + AA(k,j);
// 				}
// 				BBg = BBg/(N+1);
// 				AAg = AAg/(N+1);
// 				for(int j=0;j<N+1;j++){
// 					AA(k,j) = AA(k,j) - AAg;
// 					BB(k,j) = BB(k,j) - BBg;
// 				}
// 			}
			
// 			MatrixXd R;
// 			JacobiSVD< Matrix<double, 3, 3> > svd(BB*AA.transpose(), Eigen::ComputeFullU |Eigen::ComputeFullV);
// 			if(permit_mirror){
//                  R = svd.matrixV()*svd.matrixU().transpose();
// 			}else{
// 				if((svd.matrixV().determinant()*svd.matrixU().determinant()>0)){
// 				    R = svd.matrixV()*svd.matrixU().transpose();
// 			    }else{
// 					MatrixXd umeyama = MatrixXd::Identity(3,3);
// 			        umeyama(2,2) = -1;
// 				    R = svd.matrixV()*umeyama*svd.matrixU().transpose();
// 				}
// 			}
// 			double pre_ans = 0;
// 			MatrixXd RBB = R*BB;
// 			for(int j=0;j<N+1;j++){
// 				for(int k=0;k<3;k++){
//                     pre_ans += (AA(k,j)-RBB(k,j))*(AA(k,j)-RBB(k,j));
// 				}
// 			}
// 			pre_ans = sqrt(pre_ans/(a));
// 			if(pre_ans>eps){
// 				// cout << BB << endl;
// 				// cout << AA << endl;
// 				// cout << pre_ans << endl;
// 				continue;
// 			}else{
// 				if(N+1==a){
// 					Answer ans;
// 					ans.E = pre_ans;
// 					ans_multiple.push_back(ans);
// 				}else{
// 					permutation_space new_per;
// 					new_per.assignment = P.assignment;
// 					new_per.assignment.push_back(P.not_assign[i]);
// 					new_per.not_assign = P.not_assign;
// 					new_per.not_assign.erase(new_per.not_assign.begin()+i);
// 					Ps.push_back(new_per);
// 				}
// 			}
// 		}
// 	}
// 	cout << count << endl;
// 	return ans_multiple;
// }





