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
    const vector<string>& colNames, // name of the sort alg 
    const vector<vector<double>>& sortedArrayPerf, //the vector of vector which has the perf 
    const string& savePath = "sorting_performance.png",
    string plot_title = "Sorting Algorithm Performance") //name of the images it needs to save under
{
    if (sizes.size() != sortedArrayPerf.size()) {
        cerr << "Error: sizesV and sortedArrayPerf do not have the same number of elements." << endl;
        cerr << "sizesV size: " << sizes.size() << ", sortedArrayPerf size: " << sortedArrayPerf.size() << endl;
        return;
    }

    
    // converst sizes to log10(size) for the x-axis

    vector<double> logSizes;
    for (int size : sizes) {
        logSizes.push_back(log10(size));
    }


    plt::clf(); // Clear the current plot
    //plt::subplot(1, 2, 1);

    // plot each algorithm's perf
    for (size_t i = 0; i < colNames.size(); ++i) {
        vector<double> perf;
        for (size_t j = 0; j < sortedArrayPerf.size(); ++j) {
            perf.push_back(sortedArrayPerf[j][i]);
        }
        plt::named_plot(colNames[i], logSizes, perf);
    }
    // configure the plot
    plt::title(plot_title);
    plt::xlabel("Array Size (log scale)");
    plt::ylabel("Time Performance (ms)");
    plt::legend();
    plt::grid(true);
    plt::save("full" + savePath);

    //zoomed in
    plt::clf();// 1 row, 2 columns, 1st plot
    for (size_t i = 0; i < colNames.size(); ++i) {
        vector<double> perf;
        for (size_t j = 0; j < sortedArrayPerf.size(); ++j) {
            perf.push_back(sortedArrayPerf[j][i]);
        }
        plt::named_plot(colNames[i], logSizes, perf);
    }
    plt::title(plot_title + " Zoomed In");
    plt::xlabel("Array Size (log scale)");
    plt::ylabel("Time Performance (ms)");
    plt::legend();
    plt::grid(true);
    plt::ylim(0, 40);
    // save the plot with the specified file name
    plt::save(savePath);
}


void printVec(vector<double>& vec)
{
    cout << "{";
    for(double ele:vec)
    {
        cout << ele << " ";
    }
    cout << "}";
}


//plot for indiv algorithem
void plotComparisonPerformance(
    const vector<int>& sizesV, //
    const vector<string>& colNames, 
    const vector<vector<double>>& sortedArrayPerf,
    const vector<vector<double>>& revSortedPerf,
    const vector<vector<double>>& randomArrayPerf)
{
    // dimentions check
    if (sortedArrayPerf.size() != revSortedPerf.size() ||
        sortedArrayPerf.size() != randomArrayPerf.size() ||
        sortedArrayPerf[0].size() != revSortedPerf[0].size() ||
        sortedArrayPerf[0].size() != randomArrayPerf[0].size()) {
        cerr << "Error: Performance datasets have inconsistent dimensions." << endl;
        return;
    }

    vector<double> logSizes;
    for (int size : sizesV) {
        logSizes.push_back(log10(size)); // log conversion
    }

    plt::clf(); // Clear the current plot
    // plot per algorithm 
    for (size_t algIndex = 0; algIndex < colNames.size(); ++algIndex) { // go through each alg/ column
        plt::clf(); // Clear the current plot

        // extract performance data for the current algorithm
        vector<double> bestPerf, worstPerf, avgPerf;
        for (size_t i = 0; i < sortedArrayPerf.size(); ++i) { //going through each row/size and get the corresponding alg row

            bestPerf.push_back(sortedArrayPerf[i][algIndex]);
            worstPerf.push_back(revSortedPerf[i][algIndex]);
            avgPerf.push_back(randomArrayPerf[i][algIndex]);
        } //here for a given alg, the data fro all the sizes are pushed into the vectors
        
        plt::named_plot("Best Case", logSizes, bestPerf); //on the same plot
        plt::named_plot("Worst Case", logSizes, worstPerf);
        plt::named_plot("Average Case", logSizes, avgPerf);
        
        // Configure the plot
        plt::title("Performance Comparison: " + colNames[algIndex]);
        plt::xlabel("Array Size (log scale)");
        plt::ylabel("Time Performance (ms)");
        plt::legend();
        plt::grid(true);
        //plt::show();
        // save the plot to a file
        string savePath = colNames[algIndex] + "_comparison.png";
        plt::save(savePath);
    }
}

#endif // PLOT_SORT_H
