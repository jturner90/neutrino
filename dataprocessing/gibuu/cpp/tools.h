#pragma once

#include <unordered_map>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "Eigen/Core"

std::map<std::string, std::vector<int> >  indexMap(std::vector<std::string> const & input) {
    auto work = input[0];
    std::map<std::string, std::vector<int> > IDX;

    std::vector<int> temp;
    temp.clear();
    for (unsigned int i=0;i<input.size()-1;++i) {
        temp.push_back(i);
        if (input[i+1]!=work) {
            IDX[work] = temp;
            work = input[i+1];
            temp.clear();
        }
    }
    IDX[work] = temp;
    return IDX;
}

std::vector< std::vector<std::string> > readCSV(std::string fname) {
    std::ifstream file(fname);
    std::vector<std::vector<std::string> > data;
    std::string line = "";
    getline(file, line);
    while (getline(file, line)) {
        std::vector<std::string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        data.push_back(vec);
    }
    file.close();
    return data;
}

std::tuple<std::vector<std::string>, Eigen::ArrayXXd> convert(std::vector<std::vector<std::string> > const & data, int ncols=5) {
  std::vector<std::string> eventids;
    Eigen::ArrayXXd ret(data.size(), ncols); // rows, cols
    for (unsigned int r=0;r<ret.rows();++r) {
        for (unsigned int c=1;c<ret.cols()+1;++c) {
            ret(r,c-1) = std::stod(data[r][c]);
        }

        eventids.push_back(data[r][0]);
    }
    return std::make_tuple(eventids, ret);
}
