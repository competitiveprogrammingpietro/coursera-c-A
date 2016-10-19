/* 
   Coursera "C++ as a better C" 
   Homework 2: Implement Dijkstra's Algorithm
   Pietro Paolini
 */
#include <vector>
#include <stack>
#include <cstring>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <fstream>
#include <string>
#include <limits>
#include <queue>

using namespace std;


// The custom comparator class compare two node ids using the cost vector given
// as a parameter in the constructor, the less the cost the higher the priority
template <class T> // int, float, double etc..
class NodeComparison  {
public:
  NodeComparison(T*& cost = 0):m_cost(cost) { }
  bool operator() (const int i, const int j) {
    return  (m_cost[i] > m_cost[j]);
  }
private:
  const T* m_cost;
};

// Class holding the MST result
template <class T>
class MST {
public:
  int*  m_pred;
  const T   m_cost;
  const int m_size;
  
  MST(T*& pred, int size, T cost):
    m_cost(cost), m_size(size) {
    m_pred = new int[size];
    memcpy(m_pred, pred, sizeof(int) * size);
  }
  
  MST(const MST<T>& obj):
    m_cost(obj.m_cost), m_size(obj.m_size) {
    m_pred = new int[obj.m_size];
    memcpy(m_pred, obj.m_pred, sizeof(int) * obj.m_size);
  }

  ~MST() {
    delete [] m_pred;
  }
  
  friend ostream& operator << (std::ostream& out, const MST<T>& object) {
     out << "MST cost  = [" << object.m_cost << "]" << endl;
     out << "MST Predecessor (edges): " << '\n' << "[";
     for (int i = 0; i < object.m_size; i++)
       out << object.m_pred[i] << ",";
     out << '\b' << "]" << endl;
     return out;
   }

};

// General Graph class contains the MST, SPT method plus a bunch of
// accessories methods.
template <class T>
class Graph {
public:
  Graph(int size = 10, int density = 10, T range = 0, bool verbose = false):
    m_size(size),
    m_verbose(verbose) {
    generate(density, range);
  }

  Graph(string filename, bool verbose = false):
    m_verbose(verbose) {
    ifstream graph_data(filename.c_str());
    T cost, max = numeric_limits<T>::max();  

    // The first line gives the graph's size
    graph_data >> m_size;
    m_graph = new T*[m_size];
    for (int i = 0; i < m_size; i++) {
      m_graph[i] = new T[m_size];
      
      // Initalize the array, this eliminates loops as well
      // as the max value means no edge between the two nodes
      for (int j = 0; j < m_size; j++)
	m_graph[i][j] = max;
    }

    // Read the tuples (src, dst, cost) and store them
    int i, j;
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

  /* Dijkstra's algorithm, the algorithm is structured following the page 105 of:
   * 
   * http://www.di.unipi.it/optimize/Courses/ROM/1314/Appunti/Appunti1314.pdf
   * 
   * The document is in Italian but the pseudocode should be well understanble 
   * without knowing the language though.
   * 
   * The algorithm returns :
   * 1) the average value for the SPT computed if no destination is given
   * 2) the cost to get to the given destination otherwise
   */ 
  T Dijkstra_SPT(int root, int dest = -1) {
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
      if (m_verbose)
	cout << accumulator << "/" << m_size << endl;
      return accumulator / m_size;
    }

    // Print SPT for destination node and return its cost, this can be useful
    // when debugging the correctness of the algorithm.
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

    // Return the cost to the destination for our SPT
    return value;
  }  

  /* Prim's algorithm, the algorithm is structured following the page 254 of:
   * 
   * http://www.di.unipi.it/optimize/Courses/ROM/1314/Appunti/Appunti1314.pdf
   * 
   * The document is in Italian, the pseudocode should be understanble though.
   * 
   * The algorithm returns the MST under the form of a vector of predecessors
   * and a root node where to start is required as a parameter.
   */ 
  MST<T> Prim_MST(int root) {
    int *pred = new int[m_size];
    T *cost = new int[m_size];
    bool *visited = new bool[m_size] ();
    priority_queue<int, vector<int>, NodeComparison<T> > Q(cost);
    
    for (int i = 0; i < m_size; i++) {
      pred[i] = root;
      cost[i] = numeric_limits<T>::max();
    }
    cost[root] = 0;
    Q.push(root);
    do {
      int node = Q.top();
      Q.pop();
      visited[node] = true;
      T* edges = getAdjances(node);
      for (int i = 0; i < m_size; i++) {
	if (edges[i] >= cost[i]) 
	  continue;
	cost[i] = edges[i];
	pred[i] = node;
	if (visited[i])
	  continue;
	Q.push(i);
      }
    } while (Q.size() > 0);

    T total = 0;
    for (int i = 0; i < m_size; i++)
      total += cost[i];
    MST<T> result = MST<T>(pred, m_size, total);
    
    delete [] visited;
    delete [] cost;
    delete [] pred;
    return result;
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
      for (int j = 0; j < size; j++)
	copy[i][j] = cgraph[i*size + j];
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

int nodes[4][4] = {
  { 2147483647, 1, 5, 5 },
  { 1, 2147483647, 3, 1 },
  { 5, 3, 2147483647, 3 },
  { 5, 1, 3, 2147483647}
};
  
int main(int argn, char **argv) {
  bool verbose = false;
  
  if (argv[1]) {
    verbose = true;
  }
    {
      // First run uses the custom graph, useful for validation
      Graph<int> *custom = Graph<int>::custom((int *) &nodes[0], 4, 5, true);
      cout << *custom;
      if (verbose)
	cout << "SPT cost average " << custom->Dijkstra_SPT(0) << endl;
      MST<int> result = custom->Prim_MST(0);
      cout << result << endl;
      delete custom;
    }
    {
      // Second attempt with much bigger graph
      Graph<int> generated(60, 80, 50, false);
      if (verbose)
	cout << generated;
      cout << "SPT cost average " << generated.Dijkstra_SPT(0) << endl;
      MST<int> result = generated.Prim_MST(0);
      cout << result << endl;
    }
    {
      Graph<int> graph_from_file = Graph<int>("graph_data.txt", false);
      if (verbose)
	cout << graph_from_file << endl;
      cout << "SPT cost average " << graph_from_file.Dijkstra_SPT(0) << endl;
      cout << graph_from_file.Prim_MST(0) << endl;
    }
}

