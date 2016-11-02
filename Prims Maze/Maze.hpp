#ifndef __MazeGenerator__MazeGenerator__
#define __MazeGenerator__MazeGenerator__
#include <iostream>
#include <vector>

enum class CellState
{
	in,
	out,
	frontier,
};


enum class Direction
{
	North,
	South,
	West,
	East
};

enum class Border
{
	Space,

	HLine,
	UpperHLine,
	LowerHLine,

	VLine,
	UpperVLine,
	LowerVLine,

	TLCorner,
	TRCorner,
	BLCorner,
	BRCorner,

	HalfBVH,
	HalfUVH,
	HalfRHV,
	HalfLHV,

	HV
};

class Cell
{

public:
	bool N_wall;
	bool S_wall;
	bool W_wall;
	bool E_wall;
	CellState inMaze;

	//cell coordinates
	int x;
	int y;

	Cell()
	{
		inMaze = CellState::out;
		N_wall = true;
		S_wall = true;
		W_wall = true;
		E_wall = true;
	}

	Cell& operator = (const Cell& that)
	{
		if (this != &that)
		{
			(*this).N_wall = that.N_wall;
			(*this).S_wall = that.S_wall;
			(*this).W_wall = that.W_wall;
			(*this).E_wall = that.E_wall;
			(*this).inMaze = that.inMaze;
			(*this).x = that.x;
			(*this).y = that.y;
		}

		return *this;
	}
};


class Maze
{
private:
	int width;
	int height;
	Cell** grid;
	std::vector<Cell> frontier;

public:

	Maze(int rows, int cols);

	~Maze();

	//This method adds the frontier cells to the frontier vector
	void addFrontier(int x, int y);

	//This method marks a cell in the maze as "in" and add the cells around to the frontier set
	void expandMaze(int x, int y);

	//This method creates a vector that holds all the "in" state neighbor cells and 
	//returns the vector 
	std::vector<Cell> addNeighbour(int x, int y);

	Direction getDirection(int x, int y, int new_x, int new_y);

	//This method connects cell (x,y) with cell(new_x,new_y)
	void connectCells(int x, int y, int new_x, int new_y, Direction dir);

	void generateMaze();

	void printMaze();
};

#endif 
