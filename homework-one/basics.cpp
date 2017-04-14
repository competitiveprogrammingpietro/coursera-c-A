/* 
   Coursera "C++ as a better C" 
   Peer Graded Assignment, Homework 1: Convert a C program to C++

   September 26, 2016
   Pietro Paolini
 */
#include <iostream>
#include <vector>
using namespace std;
const int LENGTH = 40;

// These functions mirror the ones defined in the C implementation but they take
// advantages of some of the C++ features discussed during the first week lectures.
template <class T>
inline void sum (T& accumulator, int length, vector<T>& data)
{
  accumulator = 0;
  for (vector<int>::size_type i = 0; i < data.size(); ++i)
    accumulator += data[i];
}


int main()
{
  int accumulator;
  vector<int> data(LENGTH);

  for (vector<int>::size_type i = 0; i < data.size(); ++i)
    data[i] = i;
  sum(accumulator, LENGTH, data);
  cout << "sum is " << accumulator << endl;
}

