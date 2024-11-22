#ifndef PLOT_SORT_H
#define PLOT_SORT_H


#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std; 

namespace plt = matplotlibcpp;

void plotSortingPerformance(
    const vector<int>& sizes,  
    const ector<string>& colNames, // name of the sort alg 
    const vector<ector<double>>& sortedArrayPerf, //te vector of vector which has the perf 
    const string& savePath = "sorting_performance.png") //name of the images it needs to be 
{
    // vonv sizes to log10(sizes) for the x-axis
    vector<double> logSizes;
    for (int size : sizes) {
        logSizes.push_back(log10(size));
    }

    // Plot each algorithm's performance
    for (size_t i = 0; i < colNames.size(); ++i) {
        vector<double> perf;
        for (size_t j = 0; j < sortedArrayPerf.size(); ++j) {
            perf.push_back(sortedArrayPerf[j][i]);
        }
        plt::named_plot(colNames[i], logSizes, perf);
    }

    // Configure the plot
    plt::title("Sorting Algorithm Performance");
    plt::xlabel("Array Size (log scale)");
    plt::ylabel("Time Performance (s)");
    plt::legend();
    plt::grid(true);

    // Save the plot with the specified file name
    plt::save(savePath);
}
#endif // PLOT_SORT_H
