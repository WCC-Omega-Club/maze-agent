#include <cstdlib>
#include <iostream>
#include <cmath>
#include <deque>
#include <string>
#include <fstream>
#include <algorithm>
#include "Matrix.hpp"
using namespace std;

 


/// <summary>
/// Used to define different blocks on the map
/// </summary>
enum class tileType {
/// <summary>
/// Exit in top right corner
/// </summary>
	CORNER = 0,
/// <summary>
/// Empty space
/// </summary>
	FREE = 1,
/// <summary>
/// Starting position
/// </summary>
	START = 4,
/// <summary>
/// Represents a maze wall
/// </summary>
	BLOCK = 5,
/// <summary>
///  Node added to open queue but not expanded yet
/// </summary>
	OPEN = 6,
/// <summary>
/// After node is expanded and explores other nodes around it
/// </summary>
	CLOSED = 7,
/// <summary>
/// Exit marked E1
/// </summary>
	E1 = 8,
/// <summary>
/// Exit marked E2
/// </summary>
	E2 = 9,

};

/// Flags for keeping track of which search technique we are using
/// <summary>
/// Breadth First Search -> Start using this search
/// </summary>
static bool BFS = true;
/// <summary>
/// Depth First Search
/// </summary>
static bool DFS = false; 
/// <summary>
///  A * Search
/// </summary>
static bool AS = false;  	

							
/// <summary>
/// Represents a coordinate in the maze. As well as search attributes of the coordinate.
/// </summary>
struct Coord 								
{
private:	
/// <summary>
/// The position of the agent
/// </summary>
	int row, col;      	
/// <summary>
/// The agent's current depth
/// </summary>
	int depth;		
/// <summary>
/// // f(n) = g(n) + h(n) , this is the estimated total path cost
/// </summary>
	int functionF;	
/// <summary>
/// Keeps track of the parents address
/// </summary>
	Coord * parent;					

public:

/// <summary>
/// Initializes a new instance of the <see cref="Coord"/> class.
/// for use with BFS and DFS
/// </summary>
/// <param name="row">The row.</param>
/// <param name="col">The col.</param>
/// <param name="depth">The depth.</param>
/// <param name="parent">The parent.</param>
	Coord(int row, int col, int depth, Coord * parent) {
		this->row = row;
		this->col = col;
		this->depth = depth;
		this->parent = parent;
	}

	
/// <summary>
/// Initializes a new instance of the <see cref="Coord"/> class.
///	Constructor used for A* search
/// </summary>
/// <param name="row">The row.</param>
/// <param name="col">The col.</param>
/// <param name="depth">The depth.</param>
/// <param name="functionF">The function f.</param>
/// <param name="parent">The parent.</param>
	Coord(int row, int col, int depth, int functionF, Coord * parent) {
		this->row = row;
		this->col = col;
		this->depth = depth;
		this->functionF = functionF;
		this->parent = parent;
	}

	 	
/// <summary>
/// Finalizes an instance of the <see cref="Coord"/> class.
///	Delete the pointer to the parent node
/// </summary>
	~Coord() { delete parent; }

		
#pragma region Accessor Methods



	/// <summary>
	/// Gets the row.
	/// </summary>
	/// <returns></returns>
	int getRow() { return row; };	
	/// <summary>
	/// Gets the col.
	/// </summary>
	/// <returns></returns>
	int getCol() { return col; };	
	/// <summary>
	/// Gets the depth.
	/// </summary>
	/// <returns></returns>
	int getDepth() { return depth; };	
	/// <summary>
	/// Gets the current search algorithm.
	/// </summary>
	/// <returns></returns>
	int getFunctionF() { return functionF; };	
	/// <summary>
	/// Gets the parent.
	/// </summary>
	/// <returns></returns>
	Coord * getParent() { return parent; };
#pragma endregion
};


/// <summary>
/// Represents an agent in an enviroment represented by a 2D maze.
/// </summary>
class Agent
{
private:	
/// <summary>
/// Two dimensional matrix for mapping the original maze
/// </summary>
	numeric_lib::matrix<int>* Maze;
/// <summary>
/// Two dimensional matrix for keeping track of open and closed nodes
/// </summary>
	numeric_lib::matrix<int>* Route; 					
/// <summary>
/// Store the starting position of the agent
/// </summary>
	int startRow, startCol;     		
/// <summary>
/// Store the starting position for Agenting for E1 or E2
/// </summary>
	int startE1E2Row, startE1E2Col;	 	
/// <summary>
/// Store exit coordinates, either E1 or E2
/// </summary>
	int exitRow, exitCol;				
/// <summary>
/// Store coordinates of E1
/// </summary>
	int exitE1Row, exitE1Col;			
/// <summary>
/// Store coordinates of E2
/// </summary>
	int exitE2Row, exitE2Col;		 	
/// <summary>
/// Total number of moves made
/// </summary>
	int cost;						 	
/// <summary>
///  Keep track of maximum open queue size to report on memory performance
/// </summary>
	int maxOpenQSize; 				

									
/// <summary>
/// Queue of opened nodes	
/// </summary>
	deque<Coord *> openDeque;	
/// <summary>
/// Queue of closed nodes.
/// </summary>
	deque<Coord *> closedDeque;

public:	
/// <summary>
/// Initializes a new instance of the <see cref="Agent"/> class.
/// </summary>
	Agent();

	~Agent() {

		delete Maze;
		delete Route;
	}

/// <summary>
/// Executes a search through this agents enviroment for the specified tile type
/// </summary>
/// <param name="exit">The tile type.</param>
	void StartSearch(tileType tile);
	
/// <summary>
/// Executes a search of the coordinates neighbooring the specified row and column.
/// </summary>
/// <param name="row">The row.</param>
/// <param name="col">The col.</param>
	void SearchNeighboors(int, int);
/// <summary>
/// Executes an an insertion sort for use with the A* search algorithm. Similar to Djstrika's algorithm.
/// For increased search speed a 4-ary heap should be implemented instead of using an insertion sort.
/// </summary>
/// <param name="">The .</param>
/// <param name="">The .</param>
/// <param name="">The .</param>
/// <param name="">The .</param>
	void AStarSort(int, int, int, Coord *);	

/// <summary>
/// Prints the maze.
/// </summary>
	void Print();
/// <summary>
/// Cleans up allocated heap memory.
/// </summary>
	void CleanUp();	
/// <summary>
/// Switches the search type.
/// </summary>
	void SwitchSearch();
};



Agent::Agent() :cost(0), maxOpenQSize(0) {
	ifstream mapFile;
	mapFile.open("map.txt", ios::in);
	if (mapFile.is_open())
	{
		string line;
		getline(mapFile, line);
		std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
		line.erase(end_pos, line.end());
		mapFile.seekg(0, ios::beg);
																// Initialize both maps
		int size = line.size();
		Maze = new numeric_lib::matrix<int>(size, size);
		Route = new numeric_lib::matrix<int>(size, size);
		for (auto i = 0; i < size; ++i) {
			for (auto j = 0; j <size; ++j) {
				mapFile >> (*Maze)(i,j);
				(*Route)(i,j) =(*Maze)(i,j);

																// Mark the beginning and exit coordinates
				if ((*Maze)(i,j) == (int)tileType::START) {
					startE1E2Row = i;
					startE1E2Col = j;
				}
				else if ((*Maze)(i,j) == (int)tileType::E1) {
					exitE1Row = i;
					exitE1Col = j;
				}
				else if ((*Maze)(i,j) == (int)tileType::E2) {
					exitE2Row = i;
					exitE2Col = j;
				}
			}
		}
		mapFile.close();
	}
	else {
		cout << "Unable to open file";
	}
}



inline void Agent::StartSearch(tileType exit) {
																   // Define the search goal and set the starting position to the coordinates of the goal.
	if (exit == tileType::E1) {
		startRow = startE1E2Row;
		startCol = startE1E2Col;
		exitRow = exitE1Row;
		exitCol = exitE1Col;
	}
	else if (exit == tileType::E2) {
		startRow = startE1E2Row;
		startCol = startE1E2Col;
		exitRow = exitE2Row;
		exitCol = exitE2Col;
	}
	else if (exit == tileType::OPEN) {
		startRow = Maze->size() -1;
		startCol = 0;
		exitRow = 0;
		exitCol = Maze->size() - 1;
	}

																	// Load the starting coordinates into the open queue
	if (BFS) {
		openDeque.push_back(new Coord(startRow, startCol, 0, 0));
	}
	else if (DFS) {
		openDeque.push_front(new Coord(startRow, startCol, 0, 0));
	}
	else if (AS) {
		int functionF = abs(startRow - exitRow) + abs(startCol - exitCol);
		openDeque.push_back(new Coord(startRow, startCol, 0, functionF, 0));
	}
	maxOpenQSize = 1;

	(*Route)(startRow, startCol) = (int)tileType::OPEN;

																	// Keep track of where we are in the maze
	int row = 0;
	int col = 0;
	bool win = false;												// Did we solve the maze?

																	// Keep searching an exit is determined or no solution is found
	while (openDeque.size() != 0) {

																	// Get the row and column of the current position
		row = openDeque.front()->getRow();
		col = openDeque.front()->getCol();

		(*Route)(row, col) = (int)tileType::CLOSED; 				// Current position has now been opened and explored
		cost++; 													// Increase cost for each node explored

		if ((*Maze)(row, col) == (int)exit) { 						// Check if we have found the goal

			closedDeque.push_back(openDeque.front());
			openDeque.pop_front();
			win = true;
			break; 													// Goal has been found

		}
		else { 														// Goal has not been found so we must check surrounding nodes

			closedDeque.push_back(openDeque.front());
			openDeque.pop_front(); 									// Get rid of the first element in the deque
			SearchNeighboors(row, col);
		}
	}

	if (!win) 
	{ 
		cout << "This maze has no solution!" << endl; 
	}

	Print();
}


inline void Agent::SearchNeighboors(int row, int col) {
	// Use this method to push_back nodes for BFS and call insertion sort for A*
	int depth = (closedDeque.back()->getDepth()) + 1; 	// Depth of the next node will be +1 of the current node recently added to the closed deque
	Coord * parent = closedDeque.back();					// Parent address is the previous node

	// Try below the current position
	if ((col - 1 > -1) && (*Route)(row, col - 1) != (int)tileType::BLOCK && (*Route)(row, col - 1) != (int)tileType::CLOSED && (*Route)(row, col - 1) != (int)tileType::OPEN) {

		(*Route)(row, col - 1) = (int)tileType::OPEN;

		if (BFS) {
			openDeque.push_back(new Coord(row, col - 1, depth, parent));
		}
		else if (DFS) {
			openDeque.push_front(new Coord(row, col - 1, depth, parent));
		}
		else if (AS) {
			AStarSort(row, col - 1, depth, parent);
		}
	}

	// Try to the left of the current position
	if ((row - 1 > -1) && (*Route)(row - 1, col) != (int)tileType::BLOCK && (*Route)(row - 1, col) != (int)tileType::CLOSED && (*Route)(row - 1,col) != (int)tileType::OPEN) {

		(*Route)(row - 1, col) = (int)tileType::OPEN;

		if (BFS) {
			openDeque.push_back(new Coord(row - 1, col, depth, parent));
		}
		else if (DFS) {
			openDeque.push_front(new Coord(row - 1, col, depth, parent));
		}
		else if (AS) {
			AStarSort(row - 1, col, depth, parent);
		}
	}
	// Try above the current position
	if ((col + 1 <Maze->cols) && (*Route)(row , col +1) != (int)tileType::BLOCK && (*Route)(row, col + 1) != (int)tileType::CLOSED && (*Route)(row, col + 1) != (int)tileType::OPEN) {

		(*Route)(row, col + 1 )= (int)tileType::OPEN;

		if (BFS) {
			openDeque.push_back(new Coord(row, col + 1, depth, parent));
		}
		else if (DFS) {
			openDeque.push_front(new Coord(row, col + 1, depth, parent));
		}
		else if (AS) {
			AStarSort(row, col + 1, depth, parent);
		}
	}
	// Try to the right of the current position
	if ((row + 1 <Maze->rows) && (*Route)(row + 1, col) != (int)tileType::BLOCK && (*Route)(row + 1, col) != (int)tileType::CLOSED && (*Route)(row + 1, col) != (int)tileType::OPEN) {

		(*Route)(row + 1, col) = (int)tileType::OPEN;

		if (BFS) {
			openDeque.push_back(new Coord(row + 1, col, depth, parent));
		}
		else if (DFS) {
			openDeque.push_front(new Coord(row + 1, col, depth, parent));
		}
		else if (AS) {
			AStarSort(row + 1, col, depth, parent);
		}
	}

	maxOpenQSize = (openDeque.size() > maxOpenQSize) ? openDeque.size() : maxOpenQSize;
}



void Agent::AStarSort(int row, int col, int depth, Coord * parent) {
	// A * Search uses a priority queue
	// Insert elements in the deque according to priority
	int gN = depth;
	int hN = abs(exitRow - row) + abs(exitCol - col);
	int functionF = gN + hN;
	bool insertSuccess = false;
	// Perform insertion sort to insert element into deque according to highest priority ( ie. lowest functionF value )
	if (openDeque.size() == 0) 
		openDeque.push_back(new Coord(row, col, depth, functionF, parent));
	else {
		deque<Coord *>::iterator it;
		int i = 0;
		for (it = openDeque.begin(); it != openDeque.end(); it++) {
			if (functionF < openDeque[i]->getFunctionF()) {
				openDeque.insert(it, new Coord(row, col, depth, functionF, parent));
				insertSuccess = true;
				break;
			}
			i++;
		}
		if (!insertSuccess)
			openDeque.push_back(new Coord(row, col, depth, functionF, parent));

	}
}


void Agent::Print() {
	// Building the path taken by the agent by climbing the tree from the exit state
	Coord * treeIterator;
	deque<Coord *> path;
	for (treeIterator = closedDeque.back(); treeIterator->getParent() != 0; treeIterator = treeIterator->getParent()) {
		path.push_front(treeIterator);
	}

	path.push_front(treeIterator);

	// Display the the path found in the Maze
	// Used tokens that are easy to read
	cout << "\nPath Taken" << endl;
	for (int i = 0; i < Maze->rows; i++) {
		for (int j = 0; j < Maze->cols; j++) {

			if ((*Maze)(i,j) == (int)tileType::BLOCK) 
				cout << "¦";
			else if ((*Maze)(i,j) == (int)tileType::START) 
				cout << "S";
			else if ((*Maze)(i,j) == (int)tileType::E1) 
				cout << "E";
			else if ((*Maze)(i,j) == (int)tileType::E2) 
				cout << "F";
			else {
				bool on_the_path = false;
				for (int k = 0; k < path.size(); k++) {
					if (path[k]->getRow() == i && path[k]->getCol() == j) 
						on_the_path = true;
				}
				if (on_the_path) 
					cout << ".";
				else 
					cout << " ";
			}
		}
		cout << "\n";
	}

	cout << "\nComplete path: " << endl;
	for (int k = 0; k < path.size(); k++) {
		cout << "(" << Maze->rows - (path[k]->getRow()) - 1 << "," << (path[k]->getCol()) << ")";
		if (k < path.size() - 1) 
			cout << " -> ";
	}
	cout << endl << endl;

	cout << "Path Cost: " << path.size() << endl;
	cout << "Total Cost: " << cost - 1 << endl; // Initial state is not counted
	cout << "Maximum Size of Open Queue (fringe): " << maxOpenQSize << endl;
	cout << "Final Size of Open Queue : " << openDeque.size() << endl;
	cout << "Final Size of Closed Queue (expanded states): " << closedDeque.size() << endl;
	cout << "Net States Explored: " << openDeque.size() + closedDeque.size() << endl;
}

void Agent::CleanUp() {
	// Reinitialize the maze
	for (int i = 0; i < Maze->rows; i++) {
		for (int j = 0; j < Maze->cols; j++) {
			(*Route)(i,j) = (*Maze)(i,j);
		}
	}

	// Clear out the deque and reset variables
	openDeque.clear();
	closedDeque.clear();
	cost = 0;
	maxOpenQSize = 0;
}

void Agent::SwitchSearch() {
	// Specify the next search to be completed
	if (BFS) {
		BFS = false;
		DFS = true;
	}
	else if (DFS) {
		DFS = false;
		AS = true;
	}
}


