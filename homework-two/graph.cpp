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
#include <time.h>

using namespace std;

template <class T> // int, float, double etc..
class Graph {
public:

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

  static Graph<T>* custom(T* cgraph, int size, T range) {
    Graph<T> *graph = new Graph<T>(size, 10, range);
    for (int i = 0; i < size; i++)
      delete[] graph->m_graph[i];
    delete[] graph->m_graph;
    
    T** copy = new T*[size];
    for (int i = 0; i < size; i++) {
      copy[i] = new T[size]; 
      for (int j = 0; j < size; j++)
	copy[i][j] = cgraph[i*size + j];
    }
    graph->m_graph = copy;
    return graph;
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
  const int m_size;
  const float m_density;
  const T m_range;

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
	  m_graph[i][j] = m_graph[j][i] = m_range + 1;
	  continue;
	}
	bool edge = (rand() % 100 + 1) < m_density;
	
	// No connection between i <-> j
	if (!edge) {
	  m_graph[i][j] = m_graph[j][i] = m_range + 1;
	  continue; 
	}

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
  Dijkstra(Graph<T>* graph, bool verbose = false)
    :m_graph(graph), m_verbose(verbose){ }  
   
  ~Dijkstra() { }

 // The Shortest Path Tree (SPT) method computes the shortest paths between
 // the root node given as argument and all the other nodes.
 // If a destination node is give the method computes the shortest path
 // between the root and the destination node only.
 T SPT(int root, int dest = -1) {

    // Stack used for Breadth First Search 
    stack<int> Q;
    T *cost = new T[m_graph->getSize()];
    int *pred = new int[m_graph->getSize()];
    bool *visited = new bool[m_graph->getSize()]();

    // The cost vector is initializated with range plus one value in order
    // to accept any path at the first iterations.
    for (int i = 0; i < m_graph->getSize(); i++)
      cost[i] = m_graph->getRange() + 1;
    cost[root] = 0;
    
    // Default predecessor is root for all nodes
    for (int i = 0; i < m_graph->getSize(); i++)
      pred[i] = root;
    Q.push(root);
    while (!Q.empty()) {
      int node = Q.top();
      Q.pop();
      visited[node] = true;
      T* edges = m_graph->getAdjances(node);
     
      // Examining the outwards star
      for (int i = 0; i < m_graph->getSize(); i++) {

	// Do not consider loops
	if (node == i)
	  continue;

	// There is no edge in this case
	if (edges[i] == m_graph->getRange() + 1)
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
      cout << "Costs: [";
      for (int i = 0; i < m_graph->getSize(); i++)
	cout << cost[i] << ",";
      cout << '\b' << "]" << endl;
    }		     
    
    // Compute average as requested and return it
    if (dest == -1) {
      T accumulator = 0;
      for (int i = 0; i < m_graph->getSize(); i++)
	accumulator += cost[i];
      delete[] cost;
      delete[] pred;
      delete[] visited;
      return accumulator / m_graph->getSize();
    }

    // Print SPT for destination node and return its cost, this can be useful
    // when debugging the correctness of the algorithm.
    int node = pred[dest];
    cout << "Path: " << dest << ",";
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
  
private:
  Graph<T>*  m_graph;
  const bool m_verbose;
};

// Custom graph used to validate the results
int nodes[4][4] = {
  { 6, 1, 5, 5 },
  { 1, 6, 3, 1 },
  { 5, 3, 6, 3 },
  { 5, 1, 3, 6 }
};

int main() {
    {
      
      // First run uses the custom graph, useful for validation
      Graph<int> *custom = Graph<int>::custom((int *) &nodes[0], 4, 5);
      Dijkstra<int> spt(custom, false);
      cout << "Average SPT from root node(0) : " << spt.SPT(0) << endl;
      delete custom;
    }
    {
      
      // Second attempt with much bigger graph
      Graph<int> generated(50, 10, 10);
      Dijkstra<int> spt(&generated, false);
      cout << "Average SPT from root node(0) : " << spt.SPT(0) << endl;
    }
}

