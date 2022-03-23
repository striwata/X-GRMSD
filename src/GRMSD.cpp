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
	ifstream ifs_csv_file("../input/s20mentai.csv");
	string function_name = argv[1]; 
	ofstream outputfile("../output/" + function_name + ".txt");
    ifstream ifs_csv_file_model(argv[2]); 
    ifstream ifs_csv_file_data(argv[3]);
	bool permit_mirror = atoi(argv[4]);
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
	
    if(function_name == "AO" || function_name == "TSR"){
        vector<MatrixXd> Smentai;
		Smentai.reserve(120);
		for(int i=0;i<120;i++){
			Smentai.push_back(MatrixXd::Zero(3,3));
		}
		int j = 0;
		while(getline(ifs_csv_file, line)) {    
			istringstream i_stream(line);
			vector<double> strvec = split(line, ',');
			for (int i=0; i<strvec.size()/3;i++){
				Smentai[i](j,0) = strvec[3*i];
				Smentai[i](j,1) = strvec[3*i+1];
				Smentai[i](j,2) = strvec[3*i+2];
			}
			j = j + 1;
		}

		if(N_data == N_model){
            if(function_name == "AO"){
				for(int i=0;i<AA.size();i++){
					Answer ans = AO_same_size(B,AA[i],p,pen[i],Smentai,permit_mirror);
					outputfile << ans.E << endl;
					outputfile << ans.t << endl;
					outputfile << ans.R << endl;
					outputfile << ans.P << endl;
				}
			}

			if(function_name == "TSR"){
				for(int i=0;i<AA.size();i++){
					Answer ans = TSR_same_size(B,AA[i],p,pen[i],Smentai,permit_mirror,0.3);
					outputfile << ans.E << endl;
					outputfile << ans.t << endl;
					outputfile << ans.R << endl;
					outputfile << ans.P << endl;
				}
			}
		}
        
		if(N_data != N_model){
            if(function_name == "AO"){
				for(int i=0;i<AA.size();i++){
					Answer ans = AO_diff_size(B,AA[i],p,pen[i],Smentai,permit_mirror);
					outputfile << ans.E << endl;
					outputfile << ans.t << endl;
					outputfile << ans.R << endl;
					outputfile << ans.P << endl;
				}
			}

			if(function_name == "TSR"){
				cout << "TSR is not defined with molecules of different sizes" << endl;
			}
		}
		

		
	}
    if(N_data == N_model){
		if(function_name == "IsometryOpt"){
			for(int i=0;i<AA.size();i++){
				Answer ans = IsometryOpt_same_size(B,AA[i],p,pen[i],permit_mirror,1,1);
				outputfile << ans.E << endl;
				outputfile << ans.t << endl;
				outputfile << ans.R << endl;
				outputfile << ans.P << endl;
			}
		}

		if(function_name == "MatchFastOpt"){
			for(int i=0;i<AA.size();i++){
				Answer ans = MatchFastOpt_same_size(B,AA[i],p,pen[i],permit_mirror);
				outputfile << ans.E << endl;
				outputfile << ans.t << endl;
				outputfile << ans.R << endl;
				outputfile << ans.P << endl;
			}
		}
	}

	if(N_data != N_model){
		if(function_name == "IsometryOpt"){
			for(int i=0;i<AA.size();i++){
				Answer ans = IsometryOpt_diff_size(B,AA[i],p,pen[i],permit_mirror,1,1);
				outputfile << ans.E << endl;
				outputfile << ans.t << endl;
				outputfile << ans.R << endl;
				outputfile << ans.P << endl;
			}
		}

		if(function_name == "MatchFastOpt"){
			for(int i=0;i<AA.size();i++){
				Answer ans = MatchFastOpt_diff_size(B,AA[i],p,pen[i],permit_mirror);
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