/* 
   Coursera "C++ as a better C" 
   Homework 2: Implement Dijkstra's Algorithm
   Pietro Paolini
 */
#include <iostream>
#include <vector>
#include <stack>
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

  T getRange() {
    return m_range;
  }

  int getSize() {
    return m_size;
  }

  T* getAdjances(int node) {
    return m_graph[node];
  }
  
  inline bool adjacent(int x, int y) {
    return m_graph[x][y] != 0;
  }

  // Not needed
  vector<int>* neighbors(int x) {
    vector<T>* neighbors = new vector<T>();
    T* x_star = m_graph[x];
    for (int i = 0; i < m_size; i++)
      if (x_star[i] != 0)
	neighbors->push_back(i);
    return neighbors;
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
  const T m_range;
  const int m_size;

  void generate() {
    srand(time(NULL));
    // Allocate all the memory first
    m_graph = new T*[m_size];
    for (int i = 0; i < m_size; i++)
      m_graph[i] = new T[m_size];
    
    // Initilize all the edges
    for (int i = 0; i < m_size; i++) {
    
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
};

template <class T>
class Dijkstra {
public:
  Dijkstra(Graph<T> *graph)
    :m_graph(graph){ }

 T compute_path(int root, int dest) {
    stack<int> Q;
    T *cost = new int[m_graph->getSize()];
    int *pred = new int[m_graph->getSize()];
    
    for (int i = 0; i < m_graph->getSize(); i++)
      cost[i] = m_graph->getRange() + 1;
    for (int i = 0; i < m_graph->getSize(); i++) {
      if (i == root)
	pred[i] = root;
      else 
	pred[i] = 0;
    }
    Q.push(root);
    while (Q.size() > 0) {
      int node = Q.top();
      Q.pop();
      T* edges = m_graph->getAdjances(node);
     
      // Outgoing start
      for (int i = 0; i < m_graph->getSize(); i++) {
	T icost = cost[node] + edges[i];

	cout << "star(" << node << ")" <<  icost << "<" << cost[i] << endl;
	//Bellman-Ford condition
	if (icost < cost[i]) {
	  pred[i] = node;
	  cost[i] = icost;
	  cout << "pushing" << i;
	  Q.push(i);
	}
      }
    }

    int tmp = dest;
    while (tmp != root) {
      cout << tmp << ",";
      tmp = pred[tmp];
    }
    cout << root << endl;
      
    T value = cost[dest];
    delete[] cost;
    delete[] pred;
    return cost[dest]; 
  }
  
private:
  Graph<T>* m_graph; 
  
};


int main() {
  {
    Graph<int> test(40, 80, 5);
    cout << test;
    vector<int>* neighbors = test.neighbors(1);
    for (vector<int>::iterator it = neighbors->begin();
	 it != neighbors->end(); ++it) {
      cout << *it;
    }
    cout << endl;

    Dijkstra<int> d(&test);
    d.compute_path(0, 3);
  }
}

