/* 
   Coursera "C++ as a better C" 
   TODO
   Pietro Paolini
 */
#include <iostream>
#include <vector>
#include <iostream>
#include <stdlib.h>     /* srand, rand */

using namespace std;

template <class T> // int, float, double etc..
class Graph {
public:

  // Documentation, it has to be between [0..100]
  Graph(int size = 10, int density = 10, T range = 0):
    m_size(size),
    m_density(density),
    m_range(range) {
    generate();
  }

  ~Graph() {
     for (int i = 0; i < m_size; i++)
  	delete[] m_graph[i];
      delete[] m_graph;
  }
     
  void generate() {

    // // Sanitize parameters

    // Allocate all the memory first
    m_graph = new T*[m_size];
    for (int i = 0; i < m_size; i++)
      m_graph[i] = new T[m_size];

    
    // Initilize all the edges
    for (int i = 0; i < m_size; i++) {
      srand(time(NULL));
      for (int j = 0; j < m_size; j++) {
	
	// No loop
	if (i == j) {
	  m_graph[i][j] = m_graph[i][j] = 0;
	  continue;
	}
	bool edge = (rand() % 100 + 1) < m_density;
	
	// No connection between i <-> j
	if (!edge)
	  continue;

	// Assign a cost to the edge
	T cost = rand() % m_range + 1;
	m_graph[i][j] = m_graph[j][i] = cost; 
      }
    }
  }
  
  friend ostream& operator << (std::ostream& out, const Graph<T>& object) {
    out << "Graph size, density, range = ["
	<< object.m_size << "]["
	<< object.m_density << "]["
	<< object.m_range << "]" << endl;
    out << "Matrix: " << endl;
    for (int i = 0; i < object.m_size; i++) {
      for (int j = 0; j < object.m_size; j++)
	out << "[" << i << "][" << j << "] = " << object.m_graph[i][j] << '\t';
      out << endl;
    }
    return out;
  }
  
private:
  T** m_graph;
  const float m_density;
  const int m_range;
  const int m_size;
};

int main() {
  {
    Graph<int> test(10, 40, 10);
    cout << test;
  }
}
