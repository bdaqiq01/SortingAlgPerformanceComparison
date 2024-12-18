# SortingAlgPerformanceComparison

# Sorting Algorithm Performance Comparison

This project compares the performance of various sorting algorithms (Bubble Sort, Insertion Sort, Selection Sort, Merge Sort, Quick Sort, and Quick Sort with Median of Medians) under different input conditions (sorted, reversed, and random arrays). It measures the execution time of each algorithm and plots the results. When the main program is executed, the user can input the name of a CSV file containing arrays. The expected CSV format includes the following columns: the array size in the first column, the array type ("sorted," "reversed," or "random") in the second column, and the array data (separated by spaces) in the third column. The program calculates the sorting time for each array using all six algorithms and plots data size versus execution time for each algorithm. It also saves the timing information in a tabular format.

If the user does not provide a CSV file, the program will automatically generate sorted, reverse-sorted, and randomly sorted arrays of sizes ranging from 10,000 to 150,000, in increments of 10,000. It will then measure the performance of all six algorithms. At the end of the program execution, the following outputs are generated and saved to the project directory:

Two plots showing the performance of all six algorithms on sorted arrays versus the array size.
Two plots showing the performance of all six algorithms on reverse-sorted arrays versus the array size.
   Note: In this case, for Quick Sort, the pivot is set as the minimum element in each partition.
Two plots showing the performance of all six algorithms on randomly sorted arrays.
For each of the above scenarios, one plot is a zoomed-in version of the other.
Additionally, six plots are generated to show the performance of each sorting algorithm across all three array types.

---

## Prerequisites

1. **Python**: Ensure Python 3.10 or higher is installed.
   - Install the development libraries if not already installed:
     ```bash
     sudo apt-get install python3.x-dev
     ```
   replace python3.x-dev with the corresponding Python library version installed on your system
2. **C++ Compiler**: A C++17-compatible compiler is required (e.g., `g++`).

---

## Setup and Installation

1. Clone the repository or download the files.

2. Make sure the required libraries are installed:
   - Install Matplotlib for Python:
     ```bash
     pip install matplotlib
     ```

3. Ensure the `matplotlibcpp` header file is available in the same directory as `main.cpp`.

---

## Compiling the Code

1. Navigate to the project directory:
   ```bash
   cd SortingAlgPerformanceComparison

option 1. run the following command: g++ -g main.cpp -o main $(python3-config --cflags --ldflags) -lpython3.X

- replace python3-config with the appropriate version of the python-config command for your Python installation
- replace -lpython3.X with the corresponding Python library version installed on your system

option 2: using the provided makefile
depending on your system you can just run make and it might work, otherwise, modify "lpython3.10" on the following line on the makefile:
LDFLAGS = $(shell python3-config --ldflags) -lpython3.10 
to match your python version. 

This directory also contents the outputs of a run as a sample output files for refrences. 





