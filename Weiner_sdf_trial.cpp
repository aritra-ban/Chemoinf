#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <map>

using namespace std;

class Graph {
private:
    int number_vertices;
    vector<vector<int> > adjacency_matrix;
	public:
    void printGraph() const {
        cout << "Graph adjacency matrix: " << endl;
        for (int i = 0; i < number_vertices; ++i) {
            for (int j = 0; j < number_vertices; ++j) {
                if (adjacency_matrix[i][j] == numeric_limits<int>::max()) {
                    cout << "INF ";
                } else {
                    cout << adjacency_matrix[i][j] << " ";
                }
            }
            cout << endl;
        }
    }
public:
    Graph(int vertices) : number_vertices(vertices), adjacency_matrix(vertices, vector<int>(vertices, numeric_limits<int>::max())) {
        for (int i = 0; i < number_vertices; i++) {
            adjacency_matrix[i][i] = 0;
        }
    }

    void addEdge(int u, int v) {
        adjacency_matrix[u][v] = 1;
        adjacency_matrix[v][u] = 1;
    }

    void floydWarshall() {
        for (int k = 0; k < number_vertices; k++) {
            for (int i = 0; i < number_vertices; i++) {
                for (int j = 0; j < number_vertices; j++) {
                    if (adjacency_matrix[i][k] < numeric_limits<int>::max() && adjacency_matrix[k][j] < numeric_limits<int>::max()) {
                        adjacency_matrix[i][j] = min(adjacency_matrix[i][j], adjacency_matrix[i][k] + adjacency_matrix[k][j]);
                    }
                }
            }
        }
    }

    long long calculateWienerIndex() const {
        long long wienerIndex = 0;
        for (int i = 0; i < number_vertices; i++) {
            for (int j = i + 1; j < number_vertices; j++) {
                if (adjacency_matrix[i][j] != numeric_limits<int>::max()) {
                    wienerIndex += adjacency_matrix[i][j];
                }
            }
        }
        return wienerIndex;
    }
};

Graph readSdfFile(const string& filename) {
    ifstream file(filename.c_str());  
    if (!file.is_open()) {
        cout << "Could not open file: " << filename << endl;  
        exit(1); 
    }
    string line;
    getline(file, line);  
    getline(file, line);
    getline(file, line);

    getline(file, line);
    istringstream countsLine(line);
    int numAtoms, numBonds;
    countsLine >> numAtoms >> numBonds;

    vector<bool> isHeavyAtom(numAtoms, false);
    int heavyAtomCount = 0;
    map<int, int> atomIndexToVertex;

    for (int i = 0; i < numAtoms; i++) {
        getline(file, line);
        if (line[31] != 'H') {  
            isHeavyAtom[i] = true;
            atomIndexToVertex[i] = heavyAtomCount++;
        }
    }

    Graph graph(heavyAtomCount);

    for (int i = 0; i < numBonds; i++) {
        getline(file, line);
        istringstream bondLine(line);
        int atom1, atom2;
        bondLine >> atom1 >> atom2;
        atom1--; atom2--;  

        if (isHeavyAtom[atom1] && isHeavyAtom[atom2]) {
            graph.addEdge(atomIndexToVertex[atom1], atomIndexToVertex[atom2]);
        }
    }

    return graph;
}

int main() {
    string filename;
    cout << "Enter the path to the SDF file: ";
    cin >> filename;

    try {
        Graph moleculeGraph = readSdfFile(filename);
        moleculeGraph.floydWarshall();

        cout << "Constructed Graph:" << endl;
        moleculeGraph.printGraph();

        cout << "Calculating Wiener Index..." << endl;
        long long wienerIndex = moleculeGraph.calculateWienerIndex();
        cout << "Wiener Index: " << wienerIndex << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1; 
    }

    return 0; 
}

