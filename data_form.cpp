#include "Eigen/Dense"
#include <stdlib.h>
#include "data_form.h"

using namespace Eigen;
using namespace std;

bool operator< (const solution_space &sol1, const solution_space &sol2){
    return sol1.second_largest > sol2.second_largest;
};

bool operator< (const cost_from &cost_from1, const cost_from &cost_from2){
    return cost_from1.cost > cost_from2.cost;
};

bool operator< (const permutation_space &permutation_space1, const permutation_space &permutation_space2){
    return permutation_space1.E > permutation_space2.E;
}