#include <iostream>
#include "../Eigen/Dense"
#include "main_function.h"
#include "basic_function.h"
#include "Hungarian.h"
#include "cube.h"
#include <time.h>
#include <fstream>

vector<double> split(string& input, char delimiter){
    istringstream stream(input);
    string field;
    vector<double> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(stod(field));
    }
    return result;
}

int main(int argc, char *argv[]){
	string function_name = argv[1]; 
	ofstream outputfile("../output/" + function_name + ".txt");
    ifstream ifs_csv_file_model(argv[2]); 
    ifstream ifs_csv_file_data(argv[3]);
	bool permit_mirror = atoi(argv[4]);
    double r = stod(argv[5]);
	string line;
    getline(ifs_csv_file_data, line);
	vector<double> strvec = split(line, ',');
    strvec = split(line, ',');
	int N_data = (int)strvec[0];
	vector<MatrixXd> AA;
	vector<vector<int> > pen;
	int M=0;
	MatrixXd A = MatrixXd::Zero(3,N_data);
	vector<int> p(N_data);
	while(getline(ifs_csv_file_data, line)) {
		strvec = split(line, ',');
		A(0,M%N_data) = strvec[0];
		A(1,M%N_data) = strvec[1];
		A(2,M%N_data) = strvec[2];
		p[M%N_data] = (int)strvec[3];
		M = M + 1;
		if(M%N_data==0){
			AA.push_back(A);
			pen.push_back(p);
			MatrixXd A = MatrixXd::Zero(3,N_data);
			vector<int> p(N_data);
		}
	}

    getline(ifs_csv_file_model, line);
    strvec = split(line, ',');
    int N_model = (int)strvec[0];
    MatrixXd B = MatrixXd::Zero(3,N_model);
    while(getline(ifs_csv_file_model, line)){
        strvec = split(line, ',');
		B(0,M%N_model) = strvec[0];
		B(1,M%N_model) = strvec[1];
		B(2,M%N_model) = strvec[2];
		p[M%N_model] = (int)strvec[3];
		M = M + 1;
	}	

    if(function_name == "IsometrySearch"){
        for(int i=0;i<AA.size();i++){
            vector<Answer> ans_multiple = IsometrySearch(B,AA[i],p,pen[i],r,permit_mirror);
            for(int i=0;i<ans_multiple.size();i++){
                Answer ans = ans_multiple[i];
                outputfile << ans.E << endl;
                outputfile << ans.t << endl;
                outputfile << ans.R << endl;
                outputfile << ans.P << endl;
            }
        }
    }

    if(function_name == "MatchFPT"){
        for(int i=0;i<AA.size();i++){
            vector<Answer> ans_multiple = IsometrySearch(B,AA[i],p,pen[i],r,permit_mirror);
            for(int i=0;i<ans_multiple.size();i++){
                Answer ans = ans_multiple[i];
                outputfile << ans.E << endl;
                outputfile << ans.t << endl;
                outputfile << ans.R << endl;
                outputfile << ans.P << endl;
            }
        }
    }

    if(function_name == "Improved-MatchFPT"){
        for(int i=0;i<AA.size();i++){
            cout << i << endl;
            vector<Answer> ans_multiple = IsometrySearch(B,AA[i],p,pen[i],r,permit_mirror);
            for(int i=0;i<ans_multiple.size();i++){
                Answer ans = ans_multiple[i];
                outputfile << ans.E << endl;
                outputfile << ans.t << endl;
                outputfile << ans.R << endl;
                outputfile << ans.P << endl;
            }
        }
    }

	
	
	outputfile.close();
    return 0;
}