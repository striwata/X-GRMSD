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

void output(vector<Answer> ans_multiple, ofstream &outputfile){
    outputfile << "Number of substructure: " << ans_multiple.size() << endl;
    for(int i=0;i<ans_multiple.size();i++){
        outputfile << "Value of Minimized RMSD: " << ans_multiple[i].E << endl;
        outputfile << endl;
        outputfile << "Translation vector" << endl;
        outputfile << ans_multiple[i].t << endl;
        outputfile << endl;
        outputfile << "Orthogonal matrix" << endl;
        outputfile << ans_multiple[i].R << endl;
        outputfile << endl;
        outputfile << "Permutation matrix" << endl;
        outputfile << ans_multiple[i].P << endl;
        outputfile << endl;
    }
	
}

int main(int argc, char *argv[]){
	string function_name = argv[1]; 
	string model_name = argv[2];
	string data_name = argv[3];
    ifstream ifs_csv_file_model(model_name); 
    ifstream ifs_csv_file_data(data_name);
	bool permit_mirror = atoi(argv[4]);
    double r = stod(argv[5]);
    ofstream outputfile(argv[6]);
    bool multiple_query = atoi(argv[7]);
    if(multiple_query){
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

        if(function_name == "IsometrySearch"){
            for(int i=0;i<AA.size();i++){
                vector<Answer> ans_multiple = IsometrySearch(BB[i],AA[i],pen_model[i],pen_data[i],r,permit_mirror);
                output(ans_multiple,outputfile);
            }
        }

        if(function_name == "MatchFPT"){
            for(int i=0;i<AA.size();i++){
                vector<Answer> ans_multiple = MatchFPT(BB[i],AA[i],pen_model[i],pen_data[i],r,permit_mirror);
                output(ans_multiple,outputfile);
            }
        }

        if(function_name == "Improved-MatchFPT"){
            for(int i=0;i<AA.size();i++){
                vector<Answer> ans_multiple = Improved_MatchFPT(BB[i],AA[i],pen_model[i],pen_data[i],r,permit_mirror);
                output(ans_multiple,outputfile);
            }
        }
        outputfile.close();
    }else{
        string line;
        while(getline(ifs_csv_file_model, line)) {
            vector<double> strvec = split(line, ',');
            int N_model = (int)strvec[0];
            MatrixXd B = MatrixXd::Zero(3,N_model);
            vector<int> pen_model(N_model);
            for(int i=0;i<N_model;i++){
                getline(ifs_csv_file_model, line);
                vector<double> strvec = split(line, ',');
                B(0,i) = strvec[0];
                B(1,i) = strvec[1];
                B(2,i) = strvec[2];
                pen_model[i] = (int)strvec[3];
            }
            cout << "Model loaded" << endl;
            while(getline(ifs_csv_file_data, line)) {
                vector<double> strvec = split(line, ',');
                int N_data = (int)strvec[0];
                MatrixXd A = MatrixXd::Zero(3,N_data);
                vector<int> pen_data(N_data);
                for(int i=0;i<N_data;i++){
                    getline(ifs_csv_file_data, line);
                    vector<double> strvec = split(line, ',');
                    A(0,i) = strvec[0];
                    A(1,i) = strvec[1];
                    A(2,i) = strvec[2];
                    pen_data[i] = (int)strvec[3];
                }
                if(function_name == "IsometrySearch"){
                    vector<Answer> ans = IsometrySearch(B,A,pen_model,pen_data,r,permit_mirror);
                    output(ans,outputfile);
                }else if(function_name == "MatchFPT"){
                    vector<Answer> ans = MatchFPT(B,A,pen_model,pen_data,r,permit_mirror);
                    output(ans,outputfile);
                }else if(function_name == "Improved-MatchFPT"){
                    vector<Answer> ans = Improved_MatchFPT(B,A,pen_model,pen_data,r,permit_mirror);
                    output(ans,outputfile);
                }
            }
        }
        outputfile.close();
    }

	
    return 0;
}