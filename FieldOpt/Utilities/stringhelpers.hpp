/// This file contains helper methods for working with time.
#ifndef STRINGHELPERS_FUNCTIONS_H
#define STRINGHELPERS_FUNCTIONS_H

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <Eigen/Core>
#include <iomanip>

using namespace std;

template <typename T>
inline string vec_to_str(vector<T> vec) {
    stringstream str;
    str << std::fixed << scientific << std::setprecision(3);
    if (vec.size() == 0)
        return "";
    str << vec[0];
    if (vec.size() > 1) {
        for (int i = 1; i < vec.size(); ++i) {
            str << ", " << vec[i];
        }
    }
    return str.str();
}

inline string eigenvec_to_str(Eigen::VectorXd vec) {
    string str = boost::lexical_cast<string>(vec(0));
    for (int i = 1; i < vec.size(); ++i) {
        str = str + ", " + boost::lexical_cast<string>(vec(i));
    }
    return str;
}

#endif // STRINGHELPERS_FUNCTIONS_H
