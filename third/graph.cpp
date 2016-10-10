/* 
   Coursera "C++ as a better C" 
   Homework 2: Implement Dijkstra's Algorithm
   Pietro Paolini
 */
#include <iostream>
#include <vector>
#include <stack>
#include <cstring>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <fstream>
#include <string>
#include <limits>

using namespace std;

template <class T> // int, float, double etc..
class Graph {
public:

  // 
  Graph(int size = 10, int density = 10, T range = 0, bool verbose = false):
    m_size(size),
    m_verbose(verbose) {
    generate(density, range);
  }

  Graph(string filename, bool verbose = false):
    m_verbose(verbose) {
    ifstream graph_data(filename.c_str());

    // The first line gives the graph's size
    graph_data >> m_size;
    m_graph = new T*[m_size];
    for (int i = 0; i < m_size; i++) {
      m_graph[i] = new T[m_size];

      // No loop
      m_graph[i][0] = numeric_limits<T>::max();
    }

    int i, j;
    T cost, max = numeric_limits<T>::max();  
    while (graph_data >> i >> j >> cost)
      m_graph[i][j] = cost;
  }

  ~Graph() {
     for (int i = 0; i < m_size; i++)
  	delete[] m_graph[i];
      delete[] m_graph;
  }

  int getSize() {
    return m_size;
  }

  T* getAdjances(int node) {
    return &m_graph[node][0];
  }

  inline bool adjacent(int x, int y) {
    return m_graph[x][y] != 0;
  }

  vector<int>* neighbors(int x) {
    vector<T>* neighbors = new vector<T>();
    T* x_star = m_graph[x];
    for (int i = 0; i < m_size; i++)
      if (x_star[i] != 0)
	neighbors->push_back(i);
    return neighbors;
  }

  
 T SPT(int root, int dest = -1) {
    stack<int> Q;
    T *cost = new int[m_size];
    int *pred = new int[m_size];
    bool *visited = new bool[m_size]();

    // The cost vector is initializated with range plus one value in order
    // to accept any path at the first iterations.
    for (int i = 0; i < m_size; i++)
      cost[i] = numeric_limits<T>::max();
    cost[root] = 0;
    for (int i = 0; i < m_size; i++)
      pred[i] = root;
    Q.push(root);
    while (Q.size() > 0) {
      int node = Q.top();
      Q.pop();
      visited[node] = true;
      T* edges = getAdjances(node);
     
      // Outwards star
      for (int i = 0; i < m_size; i++) {

	// No loop
	if (node == i)
	  continue;

	// Skip not existing edges
	if (edges[i] == numeric_limits<T>::max())
	  continue;
	
	T icost = cost[node] + edges[i];
	if (m_verbose)
	  cout << "star(" << node << "," << i << ") "
	       <<  cost[node] << " + " << edges[i] << "(" << icost << ")"
	       << " < " << cost[i] << endl;
	
	//Bellman-Ford condition
	if (icost < cost[i]) {
	  pred[i] = node;
	  cost[i] = icost;
	  if (!visited[i])
	    Q.push(i);
	}
      }
    }

    // Print vector cost
    if (m_verbose) {
      cout << "[";
      for (int i = 0; i < m_size; i++)
	cout << cost[i] << ",";
      cout << '\b' << "]" << endl;
    }		     

    cout << "[";
    for(int i = 0; i < m_size; i++)
      cout << cost[i] << "," ;
    cout << '\b' << "]" << endl;
    
    // Compute average as requested and return it
    if (dest == -1) {
      T accumulator = 0;
      for (int i = 0; i < m_size; i++)
	accumulator += cost[i];
      delete[] cost;
      delete[] pred;
      delete[] visited;
      return accumulator / m_size;
    }

    // Print SPT for destination node and return its cost, this can be useful
    // when debugging the corrctness of the algorithm.
    int node = pred[dest];
    cout << dest << ",";
    while (node != root) {
      cout << node << ",";
      node = pred[node];
    }
    cout << root << endl;
    T value = cost[dest];
    delete[] cost;
    delete[] pred;
    delete[] visited;
    return value;
  }  

  static Graph<T>* custom(T* cgraph, int size, T range, bool verbose = false) {
    Graph<T> *graph = new Graph<T>(size, 10, range, verbose);
    for (int i = 0; i < size; i++)
      delete[] graph->m_graph[i];
    delete[] graph->m_graph;

    // Allocate and copy it over
    T** copy = new T*[size];
    for (int i = 0; i < size; i++) {
      copy[i] = new T[size];
      memcpy(&copy[i][0], &cgraph[i * size], sizeof(T) * size);
    }
    graph->m_graph = copy;
  }
    
  friend ostream& operator << (std::ostream& out, const Graph<T>& object) {
    out << "Graph size  = [" << object.m_size << "]" << endl;
    out << "Matrix: " << endl;
    for (int i = 0; i < object.m_size; i++) {
      for (int j = 0; j < object.m_size; j++)
	out << "[" << i << "][" << j << "] = " << object.m_graph[i][j] << '\t';
      out << endl;
    }
    return out;
  }
  
private:
  T**         m_graph;
  int         m_size;
  const bool  m_verbose;

  void generate(int density, T range) {
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
	  m_graph[i][j] = m_graph[i][j] = numeric_limits<T>::max();
	  continue;
	}
	bool edge = (rand() % 100 + 1) < density;
	
	// No connection between i <-> j
	if (!edge) {
	  m_graph[i][j] = m_graph[i][j] = numeric_limits<T>::max();
	  continue;
	}

	// Assign a cost to the edge
	T cost = rand() % range + 1;
	m_graph[i][j] = m_graph[j][i] = cost; 
      }
    }
  }
};



template <class T>
class Dijkstra {
public:
  Dijkstra(Graph<T> *graph, bool verbose = false)
    :m_graph(graph), m_verbose(verbose){ }
  
private:
  Graph<T>*  m_graph;
  const bool m_verbose;
};


int nodes[4][4] = {
  { 6, 1, 5, 5 },
  { 5, 6, 3, 3 },
  { 1, 3, 6, 1 },
  { 5, 3, 1, 6 }
};
  
int main() {
    // {
    //   // First run uses the custom graph, useful for validation
    //   Graph<int> *custom = Graph<int>::custom((int *) &nodes[0], 4, 5, true);
    //   cout << *custom;
    //   cout << custom->SPT(0) << endl;
    //   delete custom;
    // }
    // {
    //   // Second attempt with much bigger graph, broken
    //   Graph<int> generated(60, 80, 50, true);
    //   cout << generated;
    //   cout << generated.SPT(0) << endl;
    // }
    {
      Graph<int> graph_from_file = Graph<int>("graph_data.txt");
      cout << graph_from_file << endl;
      cout << graph_from_file.SPT(0) << endl;
    }
}

