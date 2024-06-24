#include "../Eigen/Dense"
#include <stdlib.h>
#include "cube.h"
#include "Hungarian.h"
#include <math.h>
#include "data_form.h"
#include "basic_function.h"
#include <cmath> 

#define PI 3.141592653590
using namespace Eigen;
using namespace std;

bool operator< (const Cube &cube1, const Cube &cube2){
    return cube1.lower > cube2.lower;
};

bool operator< (const Cube_Cube &cube1, const Cube_Cube &cube2){
    return cube1.lower > cube2.lower;
};

Cube make_new_cube(vector<double> center,double diameter, double upper, double lower){
    Cube c = {center,diameter,upper,lower};
    return c;
};

bool check_in_sphere(Cube C){
    double r = 0;
    for(int i=0;i<3;i++){
        r += C.center[i]*C.center[i];
    }
    if(sqrt(r) < C.diameter + PI){
        return true;
    }else{
        return false;
    }
}

MatrixXd make_rotation(vector<double> r,bool mirror){
    MatrixXd r_t = MatrixXd::Zero(3,3);
    r_t(0,1) = -r[2];
    r_t(0,2) = r[1];
    r_t(1,0) = r[2];
    r_t(1,2) = -r[0];
    r_t(2,0) = -r[1];
    r_t(2,1) = r[0];
    double r_n = 0;
    for(int i=0;i<3;i++){
        r_n += r[i]*r[i];
    }
    r_n = sqrt(r_n);
    if(r_n == 0){
        return MatrixXd::Identity(3,3);
    }else{
        MatrixXd R = MatrixXd::Identity(3,3) + sin(r_n)*r_t/r_n + (1-cos(r_n))*r_t*r_t/(r_n*r_n);
        if(!mirror){
            return R;
        }else{
            R(2,2) = -R(2,2);
            R(2,1) = -R(2,1);
            R(2,0) = -R(2,0); 
            return R;
        }
        
    }
}

MatrixXd make_translation(vector<double> r){
    MatrixXd r_t = MatrixXd::Zero(3,1);
    for(int k=0;k<3;k++){
        r_t(k,0) = r[k];
    }
    return r_t;
}

vector<Cube> divide(Cube C){
    vector<Cube> Cs(8);
    vector< vector<double> > direction = {{1.0,1.0,1.0},{1.0,1.0,-1.0},{1.0,-1.0,1.0},{1.0,-1.0,-1.0},{-1.0,1.0,1.0},{-1.0,1.0,-1.0},{-1.0,-1.0,1.0},{-1.0,-1.0,-1.0}};
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){
            direction[i][j] = direction[i][j]*C.diameter/2 + C.center[j];
        }
        Cube c;
        c.diameter = C.diameter/2;
        c.count = C.count + 1; // 後で消す
        c.center = direction[i];
        c.assignments = C.assignments;
        c.mirror = C.mirror;
        Cs[i] = c;
    }
    return Cs;
}

double calculate_upper(Cube C,vector<MatrixXd> As,vector<MatrixXd> Bs,bool mirror){
    HungarianAlgorithm HungAlgo;
    double cost = 0;
    MatrixXd R = make_rotation(C.center,mirror);
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        MatrixXd RB = R*Bs[iter];
        vector<vector<double> >M(b, vector<double>(a));
        for(int i = 0;i<b;i++){
            for(int j=0;j<a;j++){
                M[i][j] = 0;
                for(int k=0;k<3;k++){
                    M[i][j]  += (As[iter](k,j)-RB(k,i))*(As[iter](k,j)-RB(k,i));
                }
            }
        }
        vector<int> assignment;
        cost += HungAlgo.Solve(M, assignment);
    }
    return cost;
}

bool calculate_lower_K_best(Cube C,vector<MatrixXd> As,vector<MatrixXd> Bs,bool mirror,Answer &pre,double &lower,int k,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    HungarianAlgorithm HungAlgo;
    vector<vector<vector<double> > > Ms;
    MatrixXd R = make_rotation(C.center,mirror);
    int len_a = 0;
    int len_b = 0;
    for(int iter=0;iter<As.size();iter++){
        vector<vector<double> > M = calculate_lower_graph(As[iter],Bs[iter],R,C.diameter);
        Ms.push_back(M);
    }
    vector<priority_queue<solution_space> > ques;
    vector<vector<Lower_candidate> > Lower_candidates(As.size());
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        len_a += a;
        len_b += b;
        vector<int> assignment(b);
        vector<vector<double> > MM = Ms[iter];
        double cost = HungAlgo.Solve(MM, assignment);
        Lower_candidate Lower;
        Lower.lower = cost;
        Lower.assignments = assignment;
        Lower_candidates[iter].push_back(Lower);

        priority_queue<solution_space> que;
        vector<int> assignment2 = assignment;
        solution_space Sol;
        edge E;
        bool cycle = true;
        Sol.second_largest = cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment2,MM,E,cycle);
        Sol.assignment = assignment2;
        Sol.E = E;
        Sol.I = {};
        Sol.O = {};

        Lower_candidate Lower2;
        Lower2.lower = Sol.second_largest;
        Lower2.assignments =  assignment2;

        Lower_candidates[iter].push_back(Lower2);

        que.push(Sol);
        ques.push_back(que);
    }

    priority_queue<cost_from> Q;
    vector<int> froms;
    for(int iter=0;iter<As.size();iter++){
        froms.push_back(0);
    }

    double cost_sum = 0;
    for(int iter=0;iter<As.size();iter++){
        cost_sum += Lower_candidates[iter][froms[iter]].lower;
    }

    cost_from Cost;
    Cost.cost = cost_sum;
    Cost.froms = froms;
    Q.push(Cost);
    for(int s=1;s<k;s++){
        cost_from Cost = Q.top();
        Q.pop();
        if(Cost.cost>pre.E){
            return true;
        }
        lower = Cost.cost;
        vector<MatrixXd> Ps;
        vector<vector<int> > assignments;
        for(int iter=0;iter<As.size();iter++){
            int a = As[iter].cols();
            int b = Bs[iter].cols();
            assignments.push_back(Lower_candidates[iter][froms[iter]].assignments);
            MatrixXd P = MatrixXd::Zero(b,a);
            for(int i=0;i<b;i++){
                if(Lower_candidates[iter][froms[iter]].assignments[i] != -1){
                    P(i,Lower_candidates[iter][froms[iter]].assignments[i]) = 1;
                }
            }
            Ps.push_back(P);
        }
        MatrixXd R = calculate_from_assignment(Ps,As,Bs,mirror);
        double ans_candidate = calculate_ans(As,Bs,R,assignments);
        if(ans_candidate < pre.E){
            pre.E = ans_candidate;
            pre.R = R;
            MatrixXd P = MatrixXd::Zero(len_b,len_a);
            for(int i=0;i<As.size();i++){
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
            pre.P = P;
        }
        for(int iter=0;iter<As.size();iter++){
            cost_from Cost_new;
            Cost_new.froms = Cost.froms;
            Cost_new.froms[iter] += 1;

            if(Cost_new.froms[iter]<Lower_candidates[iter].size()){
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }else{
                int a = As[iter].cols();
                int b = Bs[iter].cols();
                solution_space Sol = ques[iter].top();
                ques[iter].pop();
                solution_space Sol1;
                Sol1.I = Sol.I;
                Sol1.I.push_back(Sol.E);
                Sol1.O = Sol.O;
                vector<vector<double> > M1 = Ms[iter];
    
                for(int i=0;i<Sol1.I.size();i++){
                    for(int j=0;j<a;j++){
                        if(Sol1.I[i].to != j){
                            M1[Sol1.I[i].from][j] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol1.O.size();j++){
                    M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
                }

                vector<int> assignment_2(b);
                double cost = HungAlgo.Solve(M1, assignment_2);
                edge E;
                bool cycle = true;
                Sol1.second_largest = cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment_2,M1,E,cycle); 
                Sol1.E = E;
                Sol1.assignment = assignment_2;
                ques[iter].push(Sol1);
                solution_space Sol2;
                Sol2.I = Sol.I;
                Sol2.O = Sol.O;
                Sol2.O.push_back(Sol.E);
                vector<vector<double> > M2 = Ms[iter];


                for(int i=0;i<Sol2.I.size();i++){
                    for(int j=0;j<a;j++){
                        if(Sol2.I[i].to != j){
                            M2[Sol2.I[i].from][j] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol2.O.size();j++){
                    M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
                }

                vector<int> assignment_3(a);
                cost = HungAlgo.Solve(M2, assignment_3);
                Sol2.second_largest = cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment_3,M2,E,cycle);
                Sol2.E = E;
                Sol2.assignment = assignment_3;
                ques[iter].push(Sol2);
                solution_space Sol_f = ques[iter].top();
                Lower_candidate Lower;
                Lower.assignments = Sol_f.assignment;
                Lower.lower = Sol_f.second_largest;
                Lower_candidates[iter].push_back(Lower);
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }
        }
    }
    return false;
    
}

vector<Cube_Cube> divide_divide(Cube_Cube C){
    vector<Cube_Cube> Cs(64);
    vector< vector<double> > direction = {{1.0,1.0,1.0},{1.0,1.0,-1.0},{1.0,-1.0,1.0},{1.0,-1.0,-1.0},{-1.0,1.0,1.0},{-1.0,1.0,-1.0},{-1.0,-1.0,1.0},{-1.0,-1.0,-1.0}};
    for(int i=0;i<8;i++){
        for(int ii = 0;ii<8;ii++){
            vector<double> direct_t(3);
            for(int j=0;j<3;j++){
                direct_t[j] = direction[i][j]*C.diameter_t/2 + C.center_t[j];
            }
            vector<double> direct_r(3);
            for(int j=0;j<3;j++){
                direct_r[j] = direction[ii][j]*C.diameter_r/2 + C.center_r[j];
            }
            Cube_Cube c;
            c.diameter_r = C.diameter_r/2;
            c.diameter_t = C.diameter_t/2;
            c.count = C.count + 1; // 後で消す
            c.center_r = direct_r;
            c.center_t = direct_t;
            c.mirror = C.mirror;
            Cs[i*8 + ii] = c;
        }
    }
    return Cs;
}

double K_Best_matching_rotation_upper(vector<MatrixXd> As,vector<MatrixXd> Bs,Answer pre,int K,bool mirror,Cube C_r,Cube C_t,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    vector<vector<vector<double> > >Ms;
    MatrixXd R = make_rotation(C_r.center,mirror);
    HungarianAlgorithm HungAlgo;
    int len_a = 0;
    int len_b = 0;
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        len_a += a;
        len_b += b;
        vector<vector<double> > M(b, vector<double>(a));
        MatrixXd RB = R*Bs[iter];
        for(int i = 0;i<b;i++){
            for(int j=0;j<a;j++){
                M[i][j] = 0;
                for(int k=0;k<3;k++){
                    M[i][j]  += (As[iter](k,j)-RB(k,i)-C_t.center[k])*(As[iter](k,j)-RB(k,i)-C_t.center[k]);
                }
                M[i][j] = sqrt(M[i][j]);
                M[i][j] = M[i][j] - C_t.diameter;
                M[i][j] = max(M[i][j],(double)0)*max(M[i][j],(double)0);
            }
        }
        Ms.push_back(M);
    }
    vector<priority_queue<solution_space> > ques;
    vector<vector<Lower_candidate> > Lower_candidates(As.size());
    for(int iter=0;iter<As.size();iter++){

        int a = As[iter].cols();
        int b = Bs[iter].cols();
        vector<int> assignment(b);
        vector<vector<double> > MM = Ms[iter];
        double cost = HungAlgo.Solve(MM, assignment);
        Lower_candidate Lower;
        Lower.lower = cost;
        Lower.assignments = assignment;
        Lower_candidates[iter].push_back(Lower);

        priority_queue<solution_space> que;
        vector<int> assignment2 = assignment;
        solution_space Sol;
        edge E;
        bool cycle1 = true;
        Sol.second_largest = cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment2,MM,E,cycle1);
		if(!cycle1){
			vector<vector<double> > MM = Ms[iter];
			MM[E.from][E.to] = INFINITY;
			HungAlgo.Solve(MM, assignment2);
        }
        Sol.assignment = assignment2;
        Sol.E = E;
        Sol.I = {};
        Sol.O = {};

        Lower_candidate Lower2;
        Lower2.lower = Sol.second_largest;
        Lower2.assignments =  assignment2;

        Lower_candidates[iter].push_back(Lower2);

        que.push(Sol);
        ques.push_back(que);
    }

    priority_queue<cost_from> Q;
    vector<int> froms;
    for(int iter=0;iter<As.size();iter++){
        froms.push_back(0);
    }

    double cost_sum = 0;
    for(int iter=0;iter<As.size();iter++){
        cost_sum += Lower_candidates[iter][froms[iter]].lower;
    }

    cost_from Cost;
    Cost.cost = cost_sum;
    Cost.froms = froms;
    Q.push(Cost);
    double lower;
    
    
    for(int s=1;s<K;s++){
        cost_from Cost = Q.top();
        Q.pop();
        lower = Cost.cost;
        if(Cost.cost == INFINITY){
            return INFINITY;
        }
        vector<MatrixXd> Ps;
        vector<vector<int> > assignments;
        vector<MatrixXd> BBs = Bs;
        vector<double> Bg(3);
        Bg[0] = 0;
        Bg[1] = 0;
        Bg[2] = 0;
        int N = 0;
        for(int iter=0;iter<As.size();iter++){
            int a = As[iter].cols();
            N += a;
            int b = Bs[iter].cols();
            assignments.push_back(Lower_candidates[iter][Cost.froms[iter]].assignments);
            MatrixXd P = MatrixXd::Zero(b,a); 
            for(int i=0;i<b;i++){
                if(Lower_candidates[iter][Cost.froms[iter]].assignments[i] != -1){
                    P(i,Lower_candidates[iter][Cost.froms[iter]].assignments[i]) = 1;
                    for(int k=0;k<3;k++){
                        Bg[k] += Bs[iter](k,i);
                    }
                }
            }
            Ps.push_back(P);
        }
        for(int iter=0;iter<As.size();iter++){
            int b = Bs[iter].cols();
            for(int i=0;i<b;i++){
                for(int k=0;k<3;k++){
                    BBs[iter](k,i) -= Bg[k]/N;
                }
            }
        }


        MatrixXd R = calculate_from_assignment(Ps,As,BBs,mirror);
        double ans_candidate = calculate_ans(As,BBs,R,assignments);
        if(ans_candidate < pre.E){
            pre.E = ans_candidate;
            pre.R = R;
            MatrixXd P = MatrixXd::Zero(len_b,len_a);
            pre.t = make_translation(Bg);
            pre.t = pre.t / N;
            for(int i=0;i<As.size();i++){
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
            pre.P = P;
        }

        for(int iter=0;iter<As.size();iter++){
            cost_from Cost_new;
            Cost_new.froms = Cost.froms;
            Cost_new.froms[iter] += 1;

            if(Cost_new.froms[iter]<Lower_candidates[iter].size()){
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }else{
                int a = As[iter].cols();
                int b = Bs[iter].cols();
                solution_space Sol = ques[iter].top();
                ques[iter].pop();
                solution_space Sol1;
                Sol1.I = Sol.I;
                Sol1.I.push_back(Sol.E);
                Sol1.O = Sol.O;
                vector<vector<double> > M1 = Ms[iter];
    
                for(int i=0;i<Sol1.I.size();i++){
                    for(int j=0;j<b;j++){
                        if(Sol1.I[i].from != j){
                            M1[j][Sol1.I[i].to] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol1.O.size();j++){
                    M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
                }


                vector<int> assignment_2(b);
                double cost = HungAlgo.Solve(M1, assignment_2);
                edge E;
                bool cycle2 = true;
                Sol1.second_largest = cost +  difference_second_largest_Dijkstra_edge_change_assignment(assignment_2,M1,E,cycle2); 
				if( !cycle2){
					vector<vector<double> > M1 = Ms[iter];
					for(int i=0;i<Sol1.I.size();i++){
						for(int j=0;j<b;j++){
							if(Sol1.I[i].from != j){
								M1[j][Sol1.I[i].to] = INFINITY;
							}
						}
					}
					for(int j=0;j<Sol1.O.size();j++){
						M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
					}
					M1[E.from][E.to] = INFINITY;
					HungAlgo.Solve(M1, assignment_2);
				}
                Sol1.E = E;
                Sol1.assignment = assignment_2;
                ques[iter].push(Sol1);
                solution_space Sol2;
                Sol2.I = Sol.I;
                Sol2.O = Sol.O;
                Sol2.O.push_back(Sol.E);
                vector<vector<double> > M2 = Ms[iter];
        
                for(int i=0;i<Sol2.I.size();i++){
                    for(int j=0;j<b;j++){
                        if(Sol2.I[i].from != j){
                            M2[j][Sol2.I[i].to] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol2.O.size();j++){
                    M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
                }

                vector<int> assignment_3(b);
                cost = HungAlgo.Solve(M2, assignment_3);
                bool cycle3 = true;
                Sol2.second_largest =cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment_3,M2,E,cycle3);
				if(!cycle3){
					vector<vector<double> > M2 = Ms[iter];
					for(int i=0;i<Sol2.I.size();i++){
						for(int j=0;j<b;j++){
							if(Sol2.I[i].from != j){
								M2[j][Sol2.I[i].to] = INFINITY;
							}
						}
					}
					for(int j=0;j<Sol2.O.size();j++){
						M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
					}
					M2[E.from][E.to] = INFINITY;
					HungAlgo.Solve(M2, assignment_3);
				}
                Sol2.E = E;
                Sol2.assignment = assignment_3;
                ques[iter].push(Sol2);
                solution_space Sol_f = ques[iter].top();
                Lower_candidate Lower;
                Lower.assignments = Sol_f.assignment;
                Lower.lower = Sol_f.second_largest;
                Lower_candidates[iter].push_back(Lower);
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }
        }
    }
    Cost = Q.top();
    return Cost.cost;
}

Answer K_Best_matching_rotation(vector<MatrixXd> As, vector<MatrixXd> Bs,double pro_solution,double pro_lower,int K,double &lower,bool &Flag,bool mirror,Cube C_r,Cube C_t,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    Answer Ans;
    Ans.E = 0;
    int len_a = 0;
    int len_b = 0;
    vector<vector<vector<double> > >Ms;
    MatrixXd R = make_rotation(C_r.center,mirror);
    HungarianAlgorithm HungAlgo;
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        len_a += a;
        len_b += b;
        vector<double> gamma(b);
        for(int i=0;i<b;i++){
            gamma[i] = 2*sin(min(sqrt(3)*C_r.diameter/2,PI/2))*Bs[iter].col(i).norm();
        }
        vector<vector<double> > M(b, vector<double>(a));
        MatrixXd RB = R*Bs[iter];
        for(int i = 0;i<b;i++){
            for(int j=0;j<a;j++){
                M[i][j] = 0;
                for(int k=0;k<3;k++){
                    M[i][j]  += (As[iter](k,j)-RB(k,i)-C_t.center[k])*(As[iter](k,j)-RB(k,i)-C_t.center[k]);
                }
                M[i][j] = sqrt(M[i][j]);
                M[i][j] = M[i][j] - C_t.diameter - gamma[i];
                M[i][j] = max(M[i][j],(double)0)*max(M[i][j],(double)0);
            }
        }
        Ms.push_back(M);
    }
    vector<priority_queue<solution_space> > ques;
    vector<vector<Lower_candidate> > Lower_candidates(As.size());
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        vector<int> assignment(b);
        vector<vector<double> > MM = Ms[iter];
        double cost = HungAlgo.Solve(MM, assignment);
        Lower_candidate Lower;
        Lower.lower = cost;
        Lower.assignments = assignment;
        Lower_candidates[iter].push_back(Lower);

        priority_queue<solution_space> que;
        vector<int> assignment2 = assignment;
        solution_space Sol;
        edge E;
        bool cycle1 = true;
        Sol.second_largest = cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment2,MM,E,cycle1);
		if(!cycle1){
			vector<vector<double> > MM = Ms[iter];
			MM[E.from][E.to] = INFINITY;
			HungAlgo.Solve(MM, assignment2);
        }
        Sol.assignment = assignment2;
        Sol.E = E;
        Sol.I = {};
        Sol.O = {};

        Lower_candidate Lower2;
        Lower2.lower = Sol.second_largest;
        Lower2.assignments =  assignment2;

        Lower_candidates[iter].push_back(Lower2);

        que.push(Sol);
        ques.push_back(que);
    }

    priority_queue<cost_from> Q;
    vector<int> froms;
    for(int iter=0;iter<As.size();iter++){
        froms.push_back(0);
    }

    double cost_sum = 0;
    for(int iter=0;iter<As.size();iter++){
        cost_sum += Lower_candidates[iter][froms[iter]].lower;
    }

    cost_from Cost;
    Cost.cost = cost_sum;
    Cost.froms = froms;
    Q.push(Cost);

    
    for(int s=0;s<K;s++){
        cost_from Cost = Q.top();
        Q.pop();
        
        if(Cost.cost>=pro_lower){
            return Ans;
        }
        lower = Cost.cost;
        vector<MatrixXd> Ps;
        vector<vector<int> > assignments;
        vector<MatrixXd> BBs = Bs;
        vector<double> Bg(3);
        Bg[0] = 0;
        Bg[1] = 0;
        Bg[2] = 0;
        int N = 0;
        for(int iter=0;iter<As.size();iter++){
            int a = As[iter].cols();
            N += a;
            int b = Bs[iter].cols();
            assignments.push_back(Lower_candidates[iter][Cost.froms[iter]].assignments);
            MatrixXd P = MatrixXd::Zero(b,a); 
            for(int i=0;i<b;i++){
                if(Lower_candidates[iter][Cost.froms[iter]].assignments[i] != -1){
                    P(i,Lower_candidates[iter][Cost.froms[iter]].assignments[i]) = 1;
                    for(int k=0;k<3;k++){
                        Bg[k] += Bs[iter](k,i);
                    }
                }
            }
            Ps.push_back(P);
        }
        for(int iter=0;iter<As.size();iter++){
            int b = Bs[iter].cols();
            for(int i=0;i<b;i++){
                for(int k=0;k<3;k++){
                    BBs[iter](k,i) -= Bg[k]/N;
                }
            }
        }


        MatrixXd R = calculate_from_assignment(Ps,As,BBs,mirror);
        double ans_candidate = calculate_ans(As,BBs,R,assignments);
        if(ans_candidate < pro_solution){
            pro_solution = ans_candidate;
            Ans.E = ans_candidate;
            Ans.R = R;
            Ans.t = make_translation(Bg);
            Ans.t = Ans.t / N;
            MatrixXd P = MatrixXd::Zero(len_b,len_a);
            for(int i=0;i<As.size();i++){
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
            Ans.P = P;
        }

        for(int iter=0;iter<As.size();iter++){
            cost_from Cost_new;
            Cost_new.froms = Cost.froms;
            Cost_new.froms[iter] += 1;

            if(Cost_new.froms[iter]<Lower_candidates[iter].size()){
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }else{
                int a = As[iter].cols();
                int b = Bs[iter].cols();
                solution_space Sol = ques[iter].top();
                ques[iter].pop();
                solution_space Sol1;
                Sol1.I = Sol.I;
                Sol1.I.push_back(Sol.E);
                Sol1.O = Sol.O;
                vector<vector<double> > M1 = Ms[iter];
                for(int i=0;i<Sol1.I.size();i++){
                    for(int j=0;j<b;j++){
                        if(Sol1.I[i].from != j){
                            M1[j][Sol1.I[i].to] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol1.O.size();j++){
                    M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
                }

                vector<int> assignment_2(b);
                double cost = HungAlgo.Solve(M1, assignment_2);
                edge E;
                bool cycle2 = true;
                Sol1.second_largest = cost +  difference_second_largest_Dijkstra_edge_change_assignment(assignment_2,M1,E,cycle2); 
				if( !cycle2){
					vector<vector<double> > M1 =Ms[iter];
					for(int i=0;i<Sol1.I.size();i++){
						for(int j=0;j<b;j++){
							if(Sol1.I[i].from != j){
								M1[j][Sol1.I[i].to] = INFINITY;
							}
						}
					}
					for(int j=0;j<Sol1.O.size();j++){
						M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
					}
					M1[E.from][E.to] = INFINITY;
					HungAlgo.Solve(M1, assignment_2);
				}
                Sol1.E = E;
                Sol1.assignment = assignment_2;
                ques[iter].push(Sol1);
                solution_space Sol2;
                Sol2.I = Sol.I;
                Sol2.O = Sol.O;
                Sol2.O.push_back(Sol.E);
                vector<vector<double> > M2 = Ms[iter];


        
                for(int i=0;i<Sol2.I.size();i++){
                    for(int j=0;j<b;j++){
                        if(Sol2.I[i].from != j){
                            M2[j][Sol2.I[i].to] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol2.O.size();j++){
                    M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
                }


                vector<int> assignment_3(b);
                cost = HungAlgo.Solve(M2, assignment_3);
                bool cycle3 = true;
                Sol2.second_largest =cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment_3,M2,E,cycle3);
				if(!cycle3){
					vector<vector<double> > M2 =Ms[iter];
					for(int i=0;i<Sol2.I.size();i++){
						for(int j=0;j<b;j++){
							if(Sol2.I[i].from != j){
								M2[j][Sol2.I[i].to] = INFINITY;
							}
						}
					}
					for(int j=0;j<Sol2.O.size();j++){
						M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
					}
					M2[E.from][E.to] = INFINITY;
					HungAlgo.Solve(M2, assignment_3);
				}
                Sol2.E = E;
                Sol2.assignment = assignment_3;
                ques[iter].push(Sol2);


                solution_space Sol_f = ques[iter].top();
                Lower_candidate Lower;
                Lower.assignments = Sol_f.assignment;
                Lower.lower = Sol_f.second_largest;
                Lower_candidates[iter].push_back(Lower);
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }
        }
    }
    Flag = false;
    return Ans;
}

Answer rotation_minimum_Best(vector<MatrixXd> As,vector<MatrixXd> Bs,Cube C_t,Answer &pre,bool permit_mirror,int K,int L,vector<vector<int> > a_assign,vector<vector<int> > b_assign){
    priority_queue<Cube> Q;
    Cube Cr;
    vector<double> r = {0,0,0};
	Cr.center = r;
	Cr.diameter = PI;
	Cr.lower = 0;
	Cr.upper = INFINITY;
	Cr.mirror = 0;
    Cr.count = 1;
    Q.push(Cr);
	if(permit_mirror){
		Cube Crr;
		Crr.center = r;
		Crr.diameter = PI;
		Crr.lower = 0;
		Crr.upper = INFINITY;
		Crr.mirror = 1;
		Crr.count = 1; // 後で消す
		Q.push(Crr);
	}
    Answer ans;
    ans.E = INFINITY;
    bool fin = false;
    double E_upper = INFINITY;
    for(int iter=0;iter<L;iter++){
        if(Q.empty()){
            fin = true;
            break;
        }
        Cube c_r = Q.top();
        Q.pop();
        vector<Cube> Cs = divide(c_r);
		for(int j=0;j<8;j++){
			Cube CC = Cs[j];
			MatrixXd R = make_rotation(CC.center,CC.mirror);
            double E_upper= K_Best_matching_rotation_upper(As,Bs,pre,K,CC.mirror,CC,C_t,a_assign,b_assign);
			CC.upper = E_upper;
			if(E_upper<ans.E){
				ans.E = E_upper;
                ans.t = make_translation(CC.center);
                ans.mirror = CC.mirror;
			}
			bool Flag = true;
			Answer Ans= K_Best_matching_rotation(As,Bs,pre.E,ans.E,K+CC.count,CC.lower,Flag,CC.mirror,CC,C_t,a_assign,b_assign);
			if(Ans.E<pre.E  && Ans.E != 0){
				pre = Ans;
			}
			if(!Flag){
                Q.push(CC);
			}
        }

    }
    if(!fin){
        ans.E = Q.top().lower;
    }
    
    return ans;
}

bool calculate_lower_K_best_eps(Cube C,vector<MatrixXd> As,vector<MatrixXd> Bs,bool mirror,vector<Answer> &ans_multiple,double eps,int k,vector<vector<int> > a_assign,vector<vector<int> > b_assign,int len_a,int len_b){
    HungarianAlgorithm HungAlgo;
    vector<vector<vector<double> > > Ms;
    int N = 0;
    MatrixXd R = make_rotation(C.center,mirror);
    for(int iter=0;iter<As.size();iter++){
        N+=As[iter].cols();
        vector<vector<double> > M = calculate_lower_graph(As[iter],Bs[iter],R,C.diameter);
        Ms.push_back(M);
    }
    vector<priority_queue<solution_space> > ques;
    vector<vector<Lower_candidate> > Lower_candidates(As.size());
    for(int iter=0;iter<As.size();iter++){
        int a = As[iter].cols();
        int b = Bs[iter].cols();
        vector<int> assignment(b);
        vector<vector<double> > MM = Ms[iter];
        double cost = HungAlgo.Solve(MM, assignment);
        Lower_candidate Lower;
        Lower.lower = cost;
        Lower.assignments = assignment;
        Lower_candidates[iter].push_back(Lower);

        priority_queue<solution_space> que;
        vector<int> assignment2 = assignment;
        solution_space Sol;
        edge E;
        bool cycle1 = true;
        Sol.second_largest = cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment2,MM,E,cycle1);
        if(!cycle1){
			vector<vector<double> > MM = Ms[iter];
			MM[E.from][E.to] = INFINITY;
			HungAlgo.Solve(MM, assignment2);
        }
        Sol.assignment = assignment2;
        Sol.E = E;
        Sol.I = {};
        Sol.O = {};

        Lower_candidate Lower2;
        Lower2.lower = Sol.second_largest;
        Lower2.assignments =  assignment2;

        Lower_candidates[iter].push_back(Lower2);

        que.push(Sol);
        ques.push_back(que);
    }

    priority_queue<cost_from> Q;
    vector<int> froms;
    for(int iter=0;iter<As.size();iter++){
        froms.push_back(0);
    }

    double cost_sum = 0;
    for(int iter=0;iter<As.size();iter++){
        cost_sum += Lower_candidates[iter][froms[iter]].lower;
    }

    cost_from Cost;
    Cost.cost = cost_sum;
    Cost.froms = froms;
    Q.push(Cost);
    
    for(int s=1;s<k;s++){
        cost_from Cost = Q.top();
        Q.pop();
        if(sqrt(Cost.cost/N)>sqrt(2)*eps){
            return true;
        }
        vector<MatrixXd> Ps;
        vector<vector<int> > assignments;
        vector<MatrixXd> BBs = Bs;
        vector<double> Bg(3);
        Bg[0] = 0;
        Bg[1] = 0;
        Bg[2] = 0;
        vector<MatrixXd> AAs = As;
        vector<double> Ag(3);
        Ag[0] = 0;
        Ag[1] = 0;
        Ag[2] = 0;
        for(int iter=0;iter<As.size();iter++){
            int a = As[iter].cols();
            int b = Bs[iter].cols();
            assignments.push_back(Lower_candidates[iter][Cost.froms[iter]].assignments);
            MatrixXd P = MatrixXd::Zero(b,a); 
            for(int i=0;i<b;i++){
                if(Lower_candidates[iter][Cost.froms[iter]].assignments[i] != -1){
                    P(i,Lower_candidates[iter][Cost.froms[iter]].assignments[i]) = 1;
                    for(int k=0;k<3;k++){
                        Bg[k] += Bs[iter](k,i);
                    }
                    for(int k=0;k<3;k++){
                        Ag[k] += As[iter](k,Lower_candidates[iter][Cost.froms[iter]].assignments[i]);
                    }
                }
            }
            Ps.push_back(P);
        }
        for(int iter=0;iter<Bs.size();iter++){
            int b = Bs[iter].cols();
            for(int i=0;i<b;i++){
                for(int k=0;k<3;k++){
                    BBs[iter](k,i) = Bs[iter](k,i) - Bg[k]/N;
                }
            }
        }
        for(int iter=0;iter<As.size();iter++){
            int a = As[iter].cols();
            for(int i=0;i<a;i++){
                for(int k=0;k<3;k++){
                    AAs[iter](k,i) = As[iter](k,i) - Ag[k]/N;
                }
            }
        }
        MatrixXd R = calculate_from_assignment(Ps,AAs,BBs,mirror);
        double ans_candidate = calculate_ans(AAs,BBs,R,assignments);
        MatrixXd R_2 = calculate_from_assignment(Ps,AAs,BBs,!mirror);
        double ans_candidate_2 = calculate_ans(AAs,BBs,R_2,assignments);
        if(ans_candidate > ans_candidate_2){
            ans_candidate = ans_candidate_2;
        }
        if(sqrt(ans_candidate/N) < eps){
            MatrixXd P = MatrixXd::Zero(len_b,len_a);
            int index_i = 0;
            for(int i=0;i<As.size();i++){
                int a = As[i].cols();
                int b = Bs[i].cols();
                if(a_assign[index_i].size() == 0){
                    index_i++;
                }
                for(int i_a=0;i_a<a;i_a++){
                    for(int i_b=0;i_b<b;i_b++){
                        if(Ps[i](i_b,i_a) == 1){
                            P(b_assign[index_i][i_b],a_assign[index_i][i_a]) = 1;
                        }
                    }
                }
                index_i++;
            }
            bool found = find_P(ans_multiple,P);
            if(!found){
                Answer ans;
                ans.E = sqrt(ans_candidate/N);
                ans.R = R;
                ans.P = P;
                ans_multiple.push_back(ans);
            }
            
        }
        for(int iter=0;iter<As.size();iter++){
            cost_from Cost_new;
            Cost_new.froms = Cost.froms;
            Cost_new.froms[iter] += 1;

            if(Cost_new.froms[iter]<Lower_candidates[iter].size()){
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                // Cost_new.cost = Lower_candidates[iter][Cost_new.froms[iter]].lower
                Q.push(Cost_new);
            }else{
                int a = As[iter].cols();
                int b = Bs[iter].cols();
                solution_space Sol = ques[iter].top();
                ques[iter].pop();
                solution_space Sol1;
                Sol1.I = Sol.I;
                Sol1.I.push_back(Sol.E);
                Sol1.O = Sol.O;
                vector<vector<double> > M1 = Ms[iter];
    
                for(int i=0;i<Sol1.I.size();i++){
                    for(int j=0;j<b;j++){
                        if(Sol1.I[i].from != j){
                            M1[j][Sol1.I[i].to] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol1.O.size();j++){
                    M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
                }

                vector<int> assignment_2(b);
                double cost = HungAlgo.Solve(M1, assignment_2);
                edge E;
                bool cycle2 = true;
                Sol1.second_largest = cost +  difference_second_largest_Dijkstra_edge_change_assignment(assignment_2,M1,E,cycle2); 
				if( !cycle2){
					vector<vector<double> > M1 = Ms[iter];
					for(int i=0;i<Sol1.I.size();i++){
						for(int j=0;j<b;j++){
							if(Sol1.I[i].from != j){
								M1[j][Sol1.I[i].to] = INFINITY;
							}
						}
					}
					for(int j=0;j<Sol1.O.size();j++){
						M1[Sol1.O[j].from][Sol1.O[j].to] = INFINITY;
					}
					M1[E.from][E.to] = INFINITY;
					HungAlgo.Solve(M1, assignment_2);
				}
                Sol1.E = E;
                Sol1.assignment = assignment_2;
                ques[iter].push(Sol1);
                solution_space Sol2;
                Sol2.I = Sol.I;
                Sol2.O = Sol.O;
                Sol2.O.push_back(Sol.E);
                vector<vector<double> > M2 = Ms[iter];

        
                for(int i=0;i<Sol2.I.size();i++){
                    for(int j=0;j<b;j++){
                        if(Sol2.I[i].from != j){
                            M2[j][Sol2.I[i].to] = INFINITY;
                        }
                    }
                }
                for(int j=0;j<Sol2.O.size();j++){
                    M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
                }

                vector<int> assignment_3(b);
                cost = HungAlgo.Solve(M2, assignment_3);
                bool cycle3 = true;
                Sol2.second_largest =cost + difference_second_largest_Dijkstra_edge_change_assignment(assignment_3,M2,E,cycle3);
				if(!cycle3){
					vector<vector<double> > M2 = Ms[iter];
					for(int i=0;i<Sol2.I.size();i++){
						for(int j=0;j<b;j++){
							if(Sol2.I[i].from != j){
								M2[j][Sol2.I[i].to] = INFINITY;
							}
						}
					}
					for(int j=0;j<Sol2.O.size();j++){
						M2[Sol2.O[j].from][Sol2.O[j].to] = INFINITY;
					}
					M2[E.from][E.to] = INFINITY;
					HungAlgo.Solve(M2, assignment_3);
				}
                Sol2.E = E;
                Sol2.assignment = assignment_3;
                ques[iter].push(Sol2);
                solution_space Sol_f = ques[iter].top();
                Lower_candidate Lower;
                Lower.assignments = Sol_f.assignment;
                Lower.lower = Sol_f.second_largest;
                Lower_candidates[iter].push_back(Lower);
                Cost_new.cost = Cost.cost - Lower_candidates[iter][Cost_new.froms[iter]-1].lower + Lower_candidates[iter][Cost_new.froms[iter]].lower;
                Q.push(Cost_new);
            }
        }
    }
    return false;
}