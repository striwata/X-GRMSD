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
    ofstream outputfile("test.txt");
	ifstream ifs_csv_file("s20mentai.csv"); 
	vector<MatrixXd> Smentai;
	Smentai.reserve(120);
	for(int i=0;i<120;i++){
		Smentai.push_back(MatrixXd::Zero(3,3));
	}
    string line;
	int j = 0;
	while(getline(ifs_csv_file, line)) {    
                // 「,」区切りごとにデータを読み込むためにistringstream型にする
		istringstream i_stream(line);
		vector<double> strvec = split(line, ',');
		for (int i=0; i<strvec.size()/3;i++){
			Smentai[i](j,0) = strvec[3*i];
			Smentai[i](j,1) = strvec[3*i+1];
			Smentai[i](j,2) = strvec[3*i+2];
		}
		j = j + 1;
	}
    ofstream outputfile_3("test_3.txt");

    int NN = 9;
    int M = 5;
    int KK =1;
    int L = 2;
    srand(17); 
    for(int k=2;k<=36;k++){
		int N = k;
		outputfile << N << endl;
		for(int i=0;i<M;i++){
			vector<Answer> ans;
			clock_t start;
			clock_t end;
			MatrixXd A = MatrixXd::Random(3,N);
			MatrixXd B = MatrixXd::Random(3,N);
			MatrixXd initial_matrix = MatrixXd::Identity(3,3);
			A = A/2;
			B = A;
            vector<int> label_A;
            vector<int> label_B;
            for(int j=N/2;j<N;j++){
				label_A.push_back(0);
			}
            for(int j=0;j<N/2;j++){
				label_A.push_back(1);
			}
            for(int j=N/2;j<N;j++){
				label_B.push_back(0);
			}
			for(int j=0;j<N/2;j++){
				label_B.push_back(1);
			}
            
            cout << i << endl;
			start = clock();
            outputfile_3 << "start_Iso" << endl;
			ans = IsometrySearch(A,B,label_A,label_B,0.1,0);
			end = clock();
			// cout << (double)(end - start)/ CLOCKS_PER_SEC << endl;	
            for(int iter=0;iter<ans.size();iter++){
            outputfile_3 << ans[iter].E << endl;
			outputfile_3 << ans[iter].P << endl;
            outputfile_3 << ans[iter].R << endl;
            outputfile_3 << ans[iter].t << endl;
			// outputfile_3 << ans[iter].R << endl;
			outputfile_3  << (double)(end - start)/ CLOCKS_PER_SEC << endl;
            MatrixXd RBP = ans[iter].R*B*ans[iter].P;
            double value = 0;
            for(int i=0;i<N;i++){
                for(int k=0;k<3;k++){
                    value += (A(k,i) - RBP(k,i) - ans[iter].t(k,0))*(A(k,i) - RBP(k,i) - ans[iter].t(k,0));
                }
            }
            outputfile_3 << sqrt(value/N) << endl;
            }	
			

            outputfile_3 << "start_Match" << endl;
            start = clock();
			ans = Improved_MatchFPT(A,B,label_A,label_B,0.1,0);
			end = clock();
            
			// cout << (double)(end - start)/ CLOCKS_PER_SEC << endl;
            for(int iter=0;iter<ans.size();iter++){
            outputfile_3 << ans[iter].E << endl;
			outputfile_3 << ans[iter].P << endl;
            outputfile_3 << ans[iter].R << endl;
            outputfile_3 << ans[iter].t << endl;
			// outputfile_3 << ans[iter].R << endl;
			outputfile_3  << (double)(end - start)/ CLOCKS_PER_SEC << endl;
            MatrixXd RBP = ans[iter].R*B*ans[iter].P;
            double value = 0;
            for(int i=0;i<N;i++){
                for(int k=0;k<3;k++){
                    value += (A(k,i) - RBP(k,i) - ans[iter].t(k,0))*(A(k,i) - RBP(k,i) - ans[iter].t(k,0));
                }
            }
            outputfile_3 << sqrt(value/N) << endl;
            }		
			
        }
    }
	 outputfile_3.close();
     return 0;
}
