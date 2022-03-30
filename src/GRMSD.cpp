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

void output(Answer ans, ofstream &outputfile){
	outputfile << "Value of Minimized RMSD: " << ans.E << endl;
	outputfile << endl;
	outputfile << "Translation vector" << endl;
	outputfile << ans.t << endl;
	outputfile << endl;
	outputfile << "Orthogonal matrix" << endl;
	outputfile << ans.R << endl;
	outputfile << endl;
	outputfile << "Permutation matrix" << endl;
	outputfile << ans.P << endl;
	outputfile << endl;
}

int main(int argc, char *argv[]){
	ifstream ifs_csv_file("../input/csv/s20mentai.csv");
	string function_name = argv[1]; 
	string model_name = argv[2];
	string data_name = argv[3];
	bool permit_mirror = atoi(argv[4]);
	ofstream outputfile(argv[5]);
    ifstream ifs_csv_file_model(model_name); 
    ifstream ifs_csv_file_data(data_name);
	string line;
	vector<MatrixXd> AA;
	vector<vector<int> > pen_data;
	while(getline(ifs_csv_file_data, line)) {
		vector<double> strvec = split(line, ',');
		int N_data = (int)strvec[0];
		MatrixXd A = MatrixXd::Zero(3,N_data);
	    vector<int> p(N_data);
		for(int i=0;i<N_data;i++){
			getline(ifs_csv_file_data, line);
			vector<double> strvec = split(line, ',');
			A(0,i) = strvec[0];
		    A(1,i) = strvec[1];
		    A(2,i) = strvec[2];
			p[i] = (int)strvec[3];
		}
		AA.push_back(A);
		pen_data.push_back(p);
	}

    vector<MatrixXd> BB;
	vector<vector<int> > pen_model;
	while(getline(ifs_csv_file_model, line)) {
		vector<double> strvec = split(line, ',');
		int N_model = (int)strvec[0];
		MatrixXd B = MatrixXd::Zero(3,N_model);
	    vector<int> p(N_model);
		for(int i=0;i<N_model;i++){
			getline(ifs_csv_file_model, line);
			vector<double> strvec = split(line, ',');
			B(0,i) = strvec[0];
		    B(1,i) = strvec[1];
		    B(2,i) = strvec[2];
			p[i] = (int)strvec[3];
		}
		BB.push_back(B);
		pen_model.push_back(p);
	}
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

	for(int i=0;i<AA.size();i++){
		if(AA[i].size() == BB[i].size()){
			if(function_name == "AO"){
				Answer ans = AO_same_size(BB[i],AA[i],pen_model[i],pen_data[i],Smentai,permit_mirror);
				output(ans,outputfile);
			}
			if(function_name == "TSR"){
				Answer ans = TSR_same_size(BB[i],AA[i],pen_model[i],pen_data[i],Smentai,permit_mirror,0.3);
				output(ans,outputfile);
			}
			if(function_name == "IsometryOpt"){
				Answer ans = IsometryOpt_same_size(BB[i],AA[i],pen_model[i],pen_data[i],permit_mirror,1,1);
				output(ans,outputfile);
			}
			if(function_name == "MatchFastOpt"){
				Answer ans = MatchFastOpt_same_size(BB[i],AA[i],pen_model[i],pen_data[i],permit_mirror);
				output(ans,outputfile);
			}
		}

		if(AA[i].size() != BB[i].size()){
			if(function_name == "AO"){
				Answer ans = AO_diff_size(BB[i],AA[i],pen_model[i],pen_data[i],Smentai,permit_mirror);
				output(ans,outputfile);
			}
			if(function_name == "TSR"){
				cout << "TSR is not defined with molecules of different sizes" << endl;
			}
			if(function_name == "IsometryOpt"){
				Answer ans = IsometryOpt_diff_size(BB[i],AA[i],pen_model[i],pen_data[i],permit_mirror,1,1);
				output(ans,outputfile);
			}
			if(function_name == "MatchFastOpt"){
				Answer ans = MatchFastOpt_diff_size(BB[i],AA[i],pen_model[i],pen_data[i],permit_mirror);
				output(ans,outputfile);
			}
		}
	}
	outputfile.close();
    return 0;
}