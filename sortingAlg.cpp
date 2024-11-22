#include <iostream>
#include <algorithm> 
#include <cstdlib>   
#include <vector>
#include <random>
#include <numeric>
# include <chrono>
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

/* bublesort (user api) start with the first item compare to the next item if one is bigger than the other then swap, this puts the biggest number in its place, 
after each pass the largest item of the unsorted portion bubbles up to its correct position at the end 
so at the worst case T(n) t(n-1) + n, proof that this is O(n)
expanding T(n)
T(n)  = T(n-2) + (n-1) + n 
    = T(n-3) + (n-2)+ (n-1) + n 
    = T(1) +2+3 +... +(n-1)+ n // this expands to the sum of n intergers aka 1+2+3+ ....+ n 
    = T(1) + [ from k = 2 to n of k, T(1) is constant we can simplify 
    [ from k = 2 to n of k = [ from k = 1 to n of k - 1, cause we are adding one and susbstracting 1 
    Based on the sum of arithmatic series( addition of series of number) the sum of the first n numbers is equal to 
    k = 1 to k = n, n(n+1)/2 here as N gets larger the dominent factor is the n^2 not the other ones the variable with the highest power 
    therefore we say it is O(n)
    proof of arithmatic ser
    Sn = 1 + 2 + .... + (n-1) + n 
    sn = n +(n-1) + ....+ 2 + 1 
    add the two equation togather 
    2Sn = (n+1) + (2+(n-1)) + (3+(n-2)) + ...+ (n-1+2) + (n+1) , each pair sums to n+1 and there are n terms 
    2Sn = n(n+1)  -> Sn = n(n+1)/2

    for the bubble sort we initally everything is unsorted lets say indicated by color black
    at 0+1 step is the is all is unsorted except the last element lets say indicated by color green. 
    lets k is keeping track of the sorted rgion since the green region is expanding. 
    in each pass the last swap that can possibly happen is between k -1 and k -2 
    So I want iterate to swap between every item up until it is smaller than k, then if no swap is needed
    substract k by 1 and start from index 0 one more time

    Best case: if the array is already sorted, the alg still will run the same amount of time 
    so lets pick up if the array is already sorted to improve its performence. if realize in the first pass 
    mp swaps happend then we can determine that the array is already sorted 

*/
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

/* insertion sort:
start with assumption that the first element is already sorted, and the remaining elements arenot
pick the next element in the unsorted part of the list and insert it in the correct location in the sorted part
by comparing it to every lement in the sorted list, moving from right to left
shift all the element that are greater to make room for the new element
at the very begining our region looks like the first element is green the rest are black 
the middle time looks like the half way are green and the last time it looks like all are green 
so lets make i keep track of the region which is sorted and let k be the last element that needs sortin 
then j be the rgion 

*/

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

/* insertion sort: USER API
we have two regions of the array the left side starte with 1 item and it is sorted but the right is not
the left side continuous to be sorted at all times 
we look at he first element and since one element is sorted it is already sorted. 
and we extend from the left region. index i keeps track of unsorted region
the finale value of sho
i starte from 1 and goes up tp n-1, n-1 -1 +1 = n-1 checks that there are n-1 items in the unsorted region 
what keeps track of the sorted region is 0 to i-1 
end index - start index +1 = the number of element inbetween 
*/
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


/* ,erge sort

*/

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

/*lets say we are solving the problems where in a array we want all the non-zero number in the left side mainting the order 
and all the zeros in the right side 
at the beining there is all question mark 
at some point it t 
we have three regions the left size has some non-zero numbers the middle region are the element we havent looked at yet
and the right region are some zeros, lets say j is keeping track of the start of the regions with question marks
and I keeps track of the region with zeros 
lets look at the item at j, it can be zero or non zero, the moment we look at the new element regardless the j should increment by 1
if j is non zero the j should just increment by 1 and no moving is needed, if it zero it should move the i-1 location to maintain the 
order. and j goes up until the value of i = j is equal, at which point we have looked at all the question marks
so the inner loop shoudl run until j is smaller than i 

more efficinet application: if we again have three region where the first region is non-zero, the secon region is zereos and the third
region is question mark. We keep track of the begining of the question mark region using i and the end of the non-zero region with j 
then we look at the ith item if it is zero then we just increment i, and it it is non-zero then we swap j+1 item with the i item, and then
increment both j and i
*/

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


/* quicksort: look at the array select pivot element them move the object smaller to the left and the element bigger, 
then choose another pivot do the same and repeat until the array is sorted 

*/

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
    int* arr = new int[size];

    // 1. Bubble Sort
    copy(vector_array.begin(), vector_array.end(), arr);
    auto bubbleStart = chrono::high_resolution_clock::now();
    bubbleSort(arr, size);
    auto bubbleEnd = chrono::high_resolution_clock::now();
    time[0] = chrono::duration_cast<chrono::duration<double>>(bubbleEnd - bubbleStart).count();

    // 2. Insertion Sort
    copy(vector_array.begin(), vector_array.end(), arr);
    auto insertionStart = chrono::high_resolution_clock::now();
    insertionSort(arr, size);
    auto insertionEnd = chrono::high_resolution_clock::now();
    time[1] = chrono::duration_cast<chrono::duration<double>>(insertionEnd - insertionStart).count();

    // 3. Selection Sort
    copy(vector_array.begin(), vector_array.end(), arr);
    auto selectionStart = chrono::high_resolution_clock::now();
    SelectionSort(arr, size);
    auto selectionEnd = chrono::high_resolution_clock::now();
    time[2] = chrono::duration_cast<chrono::duration<double>>(selectionEnd - selectionStart).count();

    // 4. Merge Sort
    copy(vector_array.begin(), vector_array.end(), arr);
    auto mergeStart = chrono::high_resolution_clock::now();
    mergeSort(arr, 0, size - 1);
    auto mergeEnd = chrono::high_resolution_clock::now();
    time[3] = chrono::duration_cast<chrono::duration<double>>(mergeEnd - mergeStart).count();

    // 5. Quick Sort with Random Pivot
    copy(vector_array.begin(), vector_array.end(), arr);
    auto quickStart = chrono::high_resolution_clock::now();
    if (worstBubble) {
        quickSortWorstCase(arr, 0, size);
    } else {
        quickSortRandomPivot(arr, 0, size);
    }
    auto quickEnd = chrono::high_resolution_clock::now();
    time[4] = chrono::duration_cast<chrono::duration<double>>(quickEnd - quickStart).count();

    // 6. Quick Sort with Median of Medians Pivot
    copy(vector_array.begin(), vector_array.end(), arr);
    auto medianStart = chrono::high_resolution_clock::now();
    quickSortMedianOfMedians(arr, 0, size);
    auto medianEnd = chrono::high_resolution_clock::now();
    time[5] = chrono::duration_cast<std::chrono::duration<double>>(medianEnd - medianStart).count();
    cout << "line 512 here" <<endl;

    if (arr != nullptr) {
    delete[] arr;
    arr = nullptr;
    }
    return time;
}


int main()
{
    //assesing the perf in alredy sorted arrays aka best case
    int rows = 40;
    int cols = 6; 
    int max_size = 1000;
    int step = 500; 
    vector<string> colNames = {"Bubble", "Insertion", "Selection", "Merge", "QuickRandom", "QuickMOM"};
    vector<vector<double>> sortedArrayPerf(rows, vector<double>(cols, 0.0));
    //perf for the best time
    for (int i = 0, size = step; size <= max_size; size += step, i++) {
        vector<int> sortedArray(size);
        iota(sortedArray.begin(), sortedArray.end(), 0); // Fill with increasing integers

        // Measure sorting performance for the current array size
        vector<double> times = measureSortTime(sortedArray, size);
        // Store the resulting times in sortedArrayPerf
        cout << "test " << endl; 
        sortedArrayPerf[i] = times;
    }
    
    cout <<  " end of the loop" << endl;

    return 0;
}