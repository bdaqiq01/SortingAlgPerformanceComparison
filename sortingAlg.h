#ifndef sortingALg_H
#define sortingALg_H


#include <iostream>
#include <algorithm> 
#include <cstdlib>   
#include <vector>
#include <random>
#include <numeric>
#include <chrono>
//#include "plotSort.h"
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



#endif //sortingALg_H