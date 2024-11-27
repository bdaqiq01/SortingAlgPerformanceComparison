#include <iostream>
#include <algorithm> 
#include <cstdlib>   
#include <vector>
#include <random>
#include <numeric>
#include <chrono>
#include <time.h>
#include <limits.h>
#include<math.h>
#include<fstream>
#include "plotSort.h"
#include <random>  
#include <fstream>
#include <iomanip> // 
#include <sstream>
#include <filesystem>
using namespace std;




bool isSorted(const int arr[], int n)
{
    for(int i = 1; i <n; i++)
    {
        if(arr[i-1] > arr[i]) { return false; }
    }
    return true; 
}

void printArray(int arr[], int size)
{
    for(int i = 0; i <size; i++)
    {
        cout << arr[i] << " ";
    }
    cout << endl; 
}

void bubbleSort(int arr[], int n) //in c++ an array is automatically passed by refeence, therefore the function directly operates on the original array
{
    int k = n; //k is the black region 
    while(k >0) // k stop expanding when it has reached index 0 
    {
        bool swapped = false;
        //printArray(arr, n);
        for( int i = 0; i < k-1; i++) // no cpmparison made with the actual element in k position I, 
        {                           //assuming k is the size so at the very first round the last element index is k-1
            if(arr[i] > arr[i+1])
            { 
                swap(arr[i], arr[i+1]);
                swapped = true; 
            }
        }
        //if no swap has happend break 
        if(!swapped)
        {
            break;
        }
        k--; 
    }

}


void insertionSort2(int arr[], int size) //
{
    for( int i = 1; i < size; i++ ) //i should stop at the last element that needs sorintg
    {
        int temp = arr[i]; // the element in the region we are cheching out, sicne I have the hold of this element I can replace it in the array
        int j = i-1; // the last element in the sorted region
        while (j >= 0 && temp < arr[j])
        {
            arr[j+1] = arr[j]; //move the value of j at the emptied spot 
            j--; //decrement
        } //this loop breaks when temp is no longer smaller and the 
        arr[j+1] = temp; 
        
    }
}


void insertionSort(int arr[], int size) //
{
    int i = 1; 
    while (i < size)
    {
        int j = i-1; //the resting part of the j should be 0
        while (j >= 0 )
        {
            if(arr[j+1] < arr[j]) //the new element is in j+1 
            {
                swap(arr[j+1], arr[j]); // this brings the elemtn to j position 
                j--;
            }
            else // if no sswap is needed break from the inner loop
            {
                break;
            }
        }    
        i++; 
    } 
}


/*Selection sort: USER API
look at the entire the array and find the minum put it in index 0, then we assue the second element is out 2nd minum iterate through the 2 to n-1 element 
if we fina another min we replace that

*/
void SelectionSort(int arr[], int size) //n^2
{
    int i =0;
    while(i < size) //
    {
        int MinInde = i;
        int j = i+1;
        while( j< size) // iterate through the rest and change the need 
        {
            if(arr[j] < arr[MinInde])
            {
                MinInde = j; 
            }
            j++;
        }
        if(MinInde != i ) { swap(arr[i], arr[MinInde]); }
        i++; //only swap onces u have found the min 
    }   
}



void merge(int arr[], int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temp arrays
    int L[n1], R[n2];

    int i = 0;
    while (i < n1) {
        L[i] = arr[l + i];
        i++;
    }
    int j = 0;
    while (j < n2) {
        R[j] = arr[m + 1 + j];
        j++;
    }

    // Merge the temp arrays back into arr[l..r]
    i = 0;
    j = 0;
    int k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of L[], if any
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of R[], if any
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
//USER API
void mergeSort(int arr[], int l, int r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}

void moveZeros(int num[], int size)
{
    int j = 0;
    int i = size; 
    while(j < i)
    {
        if(num[j] != 0)
        {
            j++;
        }
        else 
        {
            int k = j+1; //have to move the question mark to the left and then empty one spot at the i-1 
            while (k<i) //we have exhausted all the question mark 
            {
                num[k-1] = num[k]; // all the question mark have been moved to the left and there is a hole in i-1 
            }
            num[i-1] = 0; 
            i--; 
                     
        }
    }
}

void moveZeros2(int num[], int size)
{
    int i = 0;
    int j = -1; 
    while( i < size)
    {
        if(num[j] == 0)
        {
            i++;
        }
        else 
        {
            swap(num[j+1], num[i]);
            j++;
            i++;            
        }
    }
}

int partitionNew(int arr[], int low, int high) {
    // Randomly choose a pivot
    int randomIndex = rand() % ((high - low)+1);
    int pivot=arr[randomIndex+low];
    //cout <<pivot<<endl;
    int j=low;
    int k=high+1;
    int i=low;
    while(j<k)
    {
        if (arr[j]==pivot)
        {
            j++;
        }
        else if (arr[j]<pivot)
        {
            swap(arr[j],arr[i]);
            i++;
            j++;
        }
        else
        {
            swap(arr[k-1],arr[j]);
            k--;
        }
    }
    return i;
}

/* user API:
Worst case T(n) = T(n-1) + n 
escribes the worst-case scenario for QuickSort when a deterministic strategy is used, 
like always selecting the first element or the last element as the pivot, and the input list 
is already sorted (or reverse sorted). I
*/
void quickSortRandomPivot(int arr[], int low, int high) {
    if (low < high) {
        int pi = partitionNew(arr, low, high);

        quickSortRandomPivot(arr, low, pi - 1);
        quickSortRandomPivot(arr, pi + 1, high);
    }
}



int partitionWorst(int arr[], int low, int high) {
    // find the min element and set as pivot
    int minIndex = low;
    for (int i = low + 1; i <= high; i++) {
        if (arr[i] < arr[minIndex]) {
            minIndex = i;
        }
    }
    // Place the minimum element at the start of the range to use it as the pivot
    int pivot = arr[minIndex];
    
    int j = low;
    int k = high + 1;
    int i = low;
    while (j < k) {
        if (arr[j] == pivot) {
            j++;
        } else if (arr[j] < pivot) {
            std::swap(arr[j], arr[i]);
            i++;
            j++;
        } else {
            std::swap(arr[k - 1], arr[j]);
            k--;
        }
    }
    return i;
}

void quickSortWorstCase(int arr[], int low, int high) {
    if (low < high) {
        int pi = partitionWorst(arr, low, high);

        // Recursively sort elements before and after partition
        quickSortWorstCase(arr, low, pi - 1);
        quickSortWorstCase(arr, pi + 1, high);
    }
}


//The below functions are for quicksort MOM 
// Helper function to find the median of five numbers
int medianOfFive(int arr[], int start) {
    sort(arr + start, arr + start + 5);
    return arr[start + 2];
}

/*
T(n) = T(n/5) + n 
*/
int medianOfMedians(int arr[], int low, int high) {
    int n = high - low + 1;
    int medians[(n + 4) / 5];
    int i;
    for (i = 0; i < n / 5; i++) {
        medians[i] = medianOfFive(arr, low + i * 5);
    }
    if (i * 5 < n) {
        sort(arr + low + i * 5, arr + high + 1);  
        //medians[i] = arr[low + i * 5 + ((high - (low + i * 5)) / 2)]; 
        medians[i] = arr[low + i * 5]; 
        i++;
    }
    if (i == 1)
        return medians[0];
    else
        return medianOfMedians(medians, 0, i - 1);
}

/*
T(n) = n 
*/

int partitionWithGivenPivotNew(int arr[], int low, int high, int pivot) {
    int j=low;
    int k=high+1;
    int i=low;
    while(j<k)
    {
        if (arr[j]==pivot)
        {
            j++;
        }
        else if (arr[j]<pivot)
        {
            swap(arr[j],arr[i]);
            i++;
            j++;
        }
        else
        {
            swap(arr[k-1],arr[j]);
            k--;
        }
    }
    return i;

}

/* QuicksortMOM user API
T(n) = n + T(3n/10)+T(7n/10)

*/
void quickSortMedianOfMedians(int arr[], int low, int high) {
    if (low < high) {
        int pivotValue = medianOfMedians(arr, low, high);
        int pi = partitionWithGivenPivotNew(arr, low, high, pivotValue);

        quickSortMedianOfMedians(arr, low, pi - 1);
        quickSortMedianOfMedians(arr, pi + 1, high);
    }
}




vector<double> measureSortTime(vector<int>& vector_array, int size, bool worstBubble = false)
{
    vector<double> time(6);
    //int* arr = new int[size]; //
    vector<int> arr(size);
    // 1. Bubble Sort
    copy(vector_array.begin(), vector_array.end(), arr.begin());
    auto bubbleStart = chrono::high_resolution_clock::now();
    bubbleSort(arr.data(), size);
    auto bubbleEnd = chrono::high_resolution_clock::now();
    //time[0] = chrono::duration_cast<chrono::duration<double>>(bubbleEnd - bubbleStart).count();
    time[0] = chrono::duration_cast<chrono::milliseconds>(bubbleEnd - bubbleStart).count();

    // 2. Insertion Sort
    copy(vector_array.begin(), vector_array.end(), arr.begin());
    auto insertionStart = chrono::high_resolution_clock::now();
    insertionSort(arr.data(), size);
    auto insertionEnd = chrono::high_resolution_clock::now();
    //time[1] = chrono::duration_cast<chrono::duration<double>>(insertionEnd - insertionStart).count();
    time[1] = chrono::duration_cast<chrono::milliseconds>(insertionEnd - insertionStart).count();

    // 3. Selection Sort
    copy(vector_array.begin(), vector_array.end(), arr.begin());
    auto selectionStart = chrono::high_resolution_clock::now();
    SelectionSort(arr.data(), size);
    auto selectionEnd = chrono::high_resolution_clock::now();
    //time[2] = chrono::duration_cast<chrono::duration<double>>(selectionEnd - selectionStart).count();
    time[2] = chrono::duration_cast<chrono::milliseconds>(selectionEnd - selectionStart).count();

    // 4. Merge Sort
    copy(vector_array.begin(), vector_array.end(), arr.begin());
    auto mergeStart = chrono::high_resolution_clock::now();
    mergeSort(arr.data(), 0, size - 1);
    auto mergeEnd = chrono::high_resolution_clock::now();
    //time[3] = chrono::duration_cast<chrono::duration<double>>(mergeEnd - mergeStart).count();
    time[3] = chrono::duration_cast<chrono::milliseconds>(mergeEnd - mergeStart).count();

    // 5. quick with random or worst case 
    copy(vector_array.begin(), vector_array.end(), arr.begin());
    auto quickStart = chrono::high_resolution_clock::now();
    if (worstBubble) {
        quickSortWorstCase(arr.data(), 0, size);
    } else {
        quickSortRandomPivot(arr.data(), 0, size);
    }
    auto quickEnd = chrono::high_resolution_clock::now();
    //time[4] = chrono::duration_cast<chrono::duration<double>>(quickEnd - quickStart).count();
    time[4] = chrono::duration_cast<chrono::milliseconds>(quickEnd - quickStart).count();

    // 6. quick sort with MOM
    copy(vector_array.begin(), vector_array.end(), arr.begin());
    auto medianStart = chrono::high_resolution_clock::now();
    quickSortMedianOfMedians(arr.data(), 0, size);
    auto medianEnd = chrono::high_resolution_clock::now();
    //time[5] = chrono::duration_cast<std::chrono::duration<double>>(medianEnd - medianStart).count();
    time[5] = chrono::duration_cast<chrono::milliseconds>(medianEnd - medianStart).count();
    return time;
}

struct ArrayData {
    int size;
    string type;
    vector<int> elements;
};

vector<ArrayData> readInputCSV(const string& filename) {
    vector<ArrayData> arrays; //  read from the file
    ifstream file(filename);

    if (!file.is_open()) {
        throw runtime_error("Error: Unable to open file " + filename);
    }

    string line;
    getline(file, line); // skip the colnames 

    while (getline(file, line)) {
        stringstream ss(line);
        string sizeStr, typeStr, elementsStr;

        // get the array info
        getline(ss, sizeStr, ',');
        getline(ss, typeStr, ',');
        getline(ss, elementsStr, ',');

        int size = stoi(sizeStr); // convert to int
        vector<int> elements;
        stringstream elementsStream(elementsStr);
        int num;
        while (elementsStream >> num) { // convert to vect
            elements.push_back(num);
        }

        arrays.push_back({size, typeStr, elements}); // add to the vector of arrays
    }

    file.close();
    return arrays; // Return all parsed arrays
}


void saveToCSVfile(
    const string& filename,
    const vector<int>& sizesV,
    const vector<string>& colNames,
    const vector<vector<double>>& perfData)
{
    filesystem::path filePath = filesystem::current_path() / filename;

    // open in write into the file
    ofstream file(filePath);

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filePath << " for writing." << endl;
        return;
    }

    // setting the col names 
    file << "Size";
    for (const string& colName : colNames) {
        file << "," << colName;
    }
    file << "\n";

    // Write the data rows
    for (size_t i = 0; i < sizesV.size(); ++i) {
        file << sizesV[i]; // size in the first col 
        for (double perf : perfData[i]) {
            file << "," << std::fixed << std::setprecision(6) << perf; // Write performance values
        }
        file << "\n"; // End the row
    }

    // Close the file
    file.close();

    cout << "file saved to: " << filePath << endl;
}





int main()
{
    // assessing the perf in already sorted arrays aka best case
    int max_size = 150000;
    int step = 10000;
    int rows = max_size / step;
    int cols = 6;
    vector<int> sizesV;
    vector<string> colNames = {"Bubble", "Insertion", "Selection", "Merge", "QuickRandom", "QuickMOM"};
    
    // performance vectors
    vector<vector<double>> sortedArrayPerf(rows, vector<double>(cols, 0.0));
    vector<vector<double>> revSortedPerf(rows, vector<double>(cols, 0.0));
    vector<vector<double>> randomArrayPerf(rows, vector<double>(cols, 0.0));

    vector<int> sortedArray(max_size);
    random_device rd;  // Random seed
    mt19937 g(rd());   // Random number generator

    // if the user does not put a valid file name, it will create array sizes and test the sorting runtime
    string inputFilename;
    cout << "Enter an input filename (CSV format with arrays) or press Enter to test auto-generated arrays of diffrent sizes to test: ";
    getline(cin, inputFilename);

    if (!inputFilename.empty() && filesystem::exists(inputFilename)) {
        // ff a file is provided, read arrays and measure their performance
        try {
            vector<ArrayData> arrays = readInputCSV(inputFilename);

            vector<int> fileSizes;
            vector<vector<double>> filePerfData(arrays.size(), vector<double>(cols, 0.0));

            for (size_t i = 0; i < arrays.size(); ++i) {
                fileSizes.push_back(arrays[i].size);
                filePerfData[i] = measureSortTime(arrays[i].elements, arrays[i].size);
            }

            // Save results to a CSV
            saveToCSVfile("InputFilePerformance.csv", fileSizes, colNames, filePerfData);
            cout << "Performance results for input file saved to InputFilePerformance.csv" << endl;

            // Plot the performance results
            plotSortingPerformance(fileSizes, colNames, filePerfData, "input_file_performance.png", "Performance from Input File");
            cout << "Performance plot for input file saved to input_file_performance.png" << endl;
        } catch (const std::exception& e) {
            cerr << "Error processing input file: " << e.what() << endl;
            return 1;
        }
    } else if (!inputFilename.empty()) {
        cerr << "Error: File " << inputFilename << " does not exist." << endl;
        return 1;
    } else {
        for (int i = 0, size = step; size <= max_size; size += step, i++) {
            sizesV.push_back(size);

            sortedArray.resize(size);
            iota(sortedArray.begin(), sortedArray.end(), 0); // Fill with increasing integers

            // Measure sorting performance for the current array size
            vector<double> times = measureSortTime(sortedArray, size);
            // Store the resulting times in sortedArrayPerf
            sortedArrayPerf[i] = times; // times is the a vector of 6 time for a specific size of array here

            // Create a reversed array
            vector<int> reversedArray = sortedArray;
            reverse(reversedArray.begin(), reversedArray.end()); // reverse the order
            vector<double> reversedTimes = measureSortTime(reversedArray, size, true); // worst case of the bubble sort choosing the min as pivot
            // Store the resulting times in the rev sorted perf
            revSortedPerf[i] = reversedTimes; // Store reversed array times

            // avg case perf of three runs
            vector<double> avgTimes(cols, 0.0); // Initialize averages to 0
            for (int run = 0; run < 3; ++run) {
                vector<int> randomArray = reversedArray;
                shuffle(randomArray.begin(), randomArray.end(), g); // create a copy of shuffled
                vector<double> randomTimes = measureSortTime(randomArray, size); // random time of the six algorithms
                for (int j = 0; j < cols; ++j) {
                    avgTimes[j] += randomTimes[j]; // random time added together
                }
            }
            for (int j = 0; j < cols; ++j) {
                avgTimes[j] /= 3.0; // the accumulated average
            }
            randomArrayPerf[i] = avgTimes; // store the average times for a given size
        } // for loop

            // Plot the sorting performance
        string savePath = "sorting_performance_best_case.png";
        plotSortingPerformance(sizesV, colNames, sortedArrayPerf, savePath, "Sorting Performance on Sorted Array");

        // plot the rev sorted perf
        string revSavePath = "sorting_performance_rev_case.png";
        plotSortingPerformance(sizesV, colNames, revSortedPerf, revSavePath, "Sorting Performance on Reversed Array");

        // plot the Avg case
        string AvgSavePath = "sorting_performance_Avg_case.png";
        plotSortingPerformance(sizesV, colNames, randomArrayPerf, AvgSavePath, "Sorting Performance on Random Array");


        saveToCSVfile("SortedArrayPerformance.csv", sizesV, colNames, sortedArrayPerf);
        saveToCSVfile("RevArrayPerformance.csv", sizesV, colNames, revSortedPerf);
        saveToCSVfile("AvgArrayPerformance.csv", sizesV, colNames, randomArrayPerf);
        
        // each algorithm separately
        plotComparisonPerformance(sizesV, colNames, sortedArrayPerf, revSortedPerf, randomArrayPerf);

    } // else

    return 0;
}
