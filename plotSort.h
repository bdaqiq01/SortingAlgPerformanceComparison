#ifndef PLOT_SORT_H
#define PLOT_SORT_H

#include <matplotlibcpp.h>
#include <vector>
#include <string>
#include <cmath>

namespace plt = matplotlibcpp;

void plotSortingPerformance(
    const std::vector<int>& sizes,
    const std::vector<std::string>& colNames,
    const std::vector<std::vector<double>>& sortedArrayPerf,
    const std::string& savePath = "sorting_performance.png")
{
    // Convert sizes to log10(sizes) for the x-axis
    std::vector<double> logSizes;
    for (int size : sizes) {
        logSizes.push_back(std::log10(size));
    }

    // Plot each algorithm's performance
    for (size_t i = 0; i < colNames.size(); ++i) {
        std::vector<double> perf;
        for (size_t j = 0; j < sortedArrayPerf.size(); ++j) {
            perf.push_back(sortedArrayPerf[j][i]);
        }
        plt::plot(logSizes, perf, {{"label", colNames[i]}});
    }

    // Configure the plot
    plt::title("Sorting Algorithm Performance");
    plt::xlabel("Array Size (log scale)");
    plt::ylabel("Time Performance (s)");
    plt::legend();
    plt::grid(true);

    // Save the plot
    plt::save(savePath);
}

#endif // PLOT_SORT_H
