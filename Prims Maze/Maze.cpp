#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <time.h>
#include "Maze.hpp"


//default constructor
Maze::Maze(int rows, int cols)
{
	width = cols;
	height = rows;
	grid = new Cell*[rows];
	for (int i = 0; i < rows; ++i)
	{
		grid[i] = new Cell[cols];
		for (int j = 0; j < cols; ++j)//initialize all cells as "out"
		{
			grid[i][j].inMaze = CellState::out;
			grid[i][j].x = i;
			grid[i][j].y = j;
		}
	}
	std::vector<Cell> frontier(cols*rows);
}

//default destructor
Maze::~Maze()
{
	for (int i = 0; i < height; ++i)
	{
		delete grid[i];
	}
	delete[] grid;//or delete grid 
}

//This method adds the frontier cells to the frontier vector
void Maze::addFrontier(int x, int y)
{
	if (x >= 0 &&
		y >= 0 &&
		x <= height - 1 &&
		y <= width - 1 &&
		grid[x][y].inMaze == CellState::out)
	{
		grid[x][y].inMaze = CellState::frontier;
		frontier.push_back(grid[x][y]);
	}
}

//This method marks a cell in the maze as "in" and add the cells around to the frointer set
void Maze::expandMaze(int x, int y)
{
	grid[x][y].inMaze = CellState::in;
	addFrontier(x - 1, y);
	addFrontier(x + 1, y);
	addFrontier(x, y - 1);
	addFrontier(x, y + 1);
}

//This method creates a vector that holds all the "in" state neighbor cells and 
//returns the vector 
std::vector<Cell> Maze::addNeighbour(int x, int y)
{
	std::vector<Cell> neighbour;
	if (x > 0 && grid[x - 1][y].inMaze == CellState::in)
	{
		neighbour.push_back(grid[x - 1][y]);
	}
	if (x < height - 1 && grid[x + 1][y].inMaze == CellState::in)
	{
		neighbour.push_back(grid[x + 1][y]);
	}
	if (y > 0 && grid[x][y - 1].inMaze == CellState::in)
	{
		neighbour.push_back(grid[x][y - 1]);
	}
	if (y < width - 1 && grid[x][y + 1].inMaze == CellState::in)
	{
		neighbour.push_back(grid[x][y + 1]);
	}

	return neighbour;
}

Direction Maze::getDirection(int x, int y, int new_x, int new_y)
{
	if (x > new_x)
	{
		return Direction::South;
	}
	else if (x < new_x)
	{
		return Direction::North;
	}
	else if (y > new_y)
	{
		return Direction::East;
	}
	else //if( y < new_y)
	{
		return Direction::West;
	}
}


//This method "connects" grid[x][y] to grid[new_x][new_y]
//Note that dir is the direction going from (new_x, new_y) to (x,y)
//
//e.g dir == Direction::North, that means (x,y) is north to (new_x, new_y)

void Maze::connectCells(int x, int y, int new_x, int new_y, Direction dir)
{
	//grid[x][y] to grid[new_x][new_y] is opposite to dir 
	if (dir == Direction::North)
	{
		grid[new_x][new_y].N_wall = false;
		grid[x][y].S_wall = false;
	}
	else if (dir == Direction::South)
	{
		grid[new_x][new_y].S_wall = false;
		grid[x][y].N_wall = false;
	}
	else if (dir == Direction::West)
	{
		grid[new_x][new_y].W_wall = false;
		grid[x][y].E_wall = false;
	}
	else if (dir == Direction::East)
	{
		grid[new_x][new_y].E_wall = false;
		grid[x][y].W_wall = false;
	}
}


//1.remove an item at random "f_index" from frontier vector 
//2.addNeighbour(new_x, new_y);
//3.randomly choose a neighbour cell(the cell that is "in" the maze)
//4.connect the neighbour cell to the frontier cell 
//5.remove the frontier cell from the frontier vector
//6.repeat the process until no more frontier cell is in the vector
void Maze::generateMaze()
{

	srand(static_cast<int>(time(nullptr)));

	//randomly create the starting coordinate
	int x = rand() % (height);
	int y = rand() % (width);

	expandMaze(x, y);

	unsigned long f_size = frontier.size();
	while (f_size > 0)
	{
		//randomly choose a cell from the frontier set
		int f_index = rand() % (f_size);
		Cell f_cell = frontier[f_index];
		int f_new_x = f_cell.x;
		int f_new_y = f_cell.y;

		std::vector<Cell> neighbour = addNeighbour(f_new_x, f_new_y);
		unsigned long n_size = neighbour.size();

		int  n_index = rand() % (n_size);
		Cell n_cell = neighbour[n_index];
		int  n_new_x = n_cell.x;
		int  n_new_y = n_cell.y;

		//dir is the direction going from n_cell to f_cell;
		Direction dir = getDirection(f_new_x, f_new_y, n_new_x, n_new_y);
		connectCells(f_new_x, f_new_y, n_new_x, n_new_y, dir);

		//Remove the chosen frontier cell from the vector
		frontier.erase(frontier.begin() + f_index);

		expandMaze(f_new_x, f_new_y);

		f_size = frontier.size();//update f_size and check if the frontier vector is empty
	}


}


//This method prints the maze into the console
void Maze::printMaze()
{
	std::cout << "Initialization" << std::endl;
	Border border[2 * height + 1][2 * width + 1];

	for (int i = 0; i < 2 * height + 1; ++i)
	{
		for (int j = 0; j< 2 * width + 1; ++j)
		{
			//Doubles the horizontal size of a printed cell
			const unsigned int number_of_prints = (j % 2 == 1) ? 2 : 1;
			for (int k = 0; k < number_of_prints; ++k)
			{

				//4 corners 
				if (i == 0 && j == 0)
				{
					border[i][j] = Border::TLCorner;
					std::cout << u8"┌";
				}

				else if (i == 0 && j == 2 * width)
				{
					border[i][j] = Border::TRCorner;
					std::cout << u8"┐" << std::endl;
				}

				else if (i == 2 * height && j == 0)
				{
					border[i][j] = Border::BLCorner;
					std::cout << u8"└";
				}

				else if (i == 2 * height && j == 2 * width)
				{
					border[i][j] = Border::BRCorner;
					std::cout << u8"┘" << std::endl;
				}

				//borders 
				else if (i % 2 == 0 && j % 2 != 0)
				{
					border[i][j] = Border::HLine;  //u8"─"
					std::cout << u8"─";
				}
				else if (i % 2 != 0 && j % 2 == 0)
				{
					border[i][j] = Border::VLine;  //u8"│"
					std::cout << u8"│";
					if (j == 2 * width)
						std::cout << std::endl;
				}
				else if (i == 0 && j % 2 == 0)
				{
					border[i][j] = Border::HalfBVH;//u8"┬"
					std::cout << u8"┬";
				}
				else if (i == 2 * height && j % 2 == 0)
				{
					border[i][j] = Border::HalfUVH;//u8"┴"
					std::cout << u8"┴";
				}
				else if (i % 2 == 0 && j == 0)
				{
					border[i][j] = Border::HalfRHV;//u8"├"
					std::cout << u8"├";
				}
				else if (i % 2 == 0 && j == 2 * width)
				{
					border[i][j] = Border::HalfLHV;//u8"┤"
					std::cout << u8"┤" << std::endl;
				}
				else if (i % 2 == 0 && j % 2 == 0)
				{
					border[i][j] = Border::HV; //u8"┼
					std::cout << u8"┼";
				}
				//cells
				else
				{
					border[i][j] = Border::Space;
					std::cout << " ";
				}
			}
		}
	}
	//--------------------------------------------------------------
	//This modifies the border 2d array based on the cell grid
	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{

			//North Border
			if (grid[y][x].N_wall == false)
			{
				border[2 * y][2 * x + 1] = Border::Space;//break the north border

				switch (border[2 * y][2 * x])//top left corner of the current cell
				{
				case Border::HalfRHV://"├"
					border[2 * y][2 * x] = Border::VLine;
					break;
				case Border::HV://"┼
					border[2 * y][2 * x] = Border::HalfLHV;
					break;
				case Border::HalfBVH:
					border[2 * y][2 * x] = Border::TRCorner;
					break;
				case Border::HalfUVH:
					border[2 * y][2 * x] = Border::BRCorner;
					break;
				case Border::BLCorner:
					border[2 * y][2 * x] = Border::UpperVLine;
					break;
				case Border::TLCorner:
					border[2 * y][2 * x] = Border::LowerVLine;
					break;
				case Border::HLine:
					border[2 * y][2 * x] = Border::LowerHLine;
					break;
				case Border::UpperHLine:
					border[2 * y][2 * x] = Border::Space;
					break;
				default: break;
				}

				switch (border[2 * y][2 * x + 2])//top right corner of the current cell
				{
				case Border::HalfLHV://"┤"
					border[2 * y][2 * x + 2] = Border::VLine;
					break;
				case Border::HV://"┼
					border[2 * y][2 * x + 2] = Border::HalfRHV;
					break;
				case Border::HalfBVH://"┬"
					border[2 * y][2 * x + 2] = Border::TLCorner;
					break;
				case Border::HalfUVH: //"┴"
					border[2 * y][2 * x + 2] = Border::BLCorner;
					break;
				case Border::BRCorner:
					border[2 * y][2 * x + 2] = Border::UpperVLine;
					break;
				case Border::TRCorner:
					border[2 * y][2 * x + 2] = Border::LowerVLine;
					break;
				case Border::HLine:
					border[2 * y][2 * x + 2] = Border::UpperHLine;
					break;
				case Border::LowerHLine:
					border[2 * y][2 * x + 2] = Border::Space;
					break;
				default: break;

				}

			}

			//East Border 
			if (grid[y][x].E_wall == false)
			{
				border[2 * y + 1][2 * x + 2] = Border::Space;//break the east wall

				switch (border[2 * y][2 * x + 2])//top right corner
				{
				case Border::HalfBVH:
					border[2 * y][2 * x + 2] = Border::HLine;
					break;
				case Border::HV:
					border[2 * y][2 * x + 2] = Border::HalfUVH;
					break;
				case Border::HalfLHV:
					border[2 * y][2 * x + 2] = Border::BRCorner;
					break;
				case Border::HalfRHV:
					border[2 * y][2 * x + 2] = Border::BLCorner;
					break;
				case Border::TLCorner:
					border[2 * y][2 * x + 2] = Border::UpperHLine;
					break;
				case Border::TRCorner:
					border[2 * y][2 * x + 2] = Border::LowerHLine;
					break;
				default: break;
				}

				switch (border[2 * y + 2][2 * x + 2])//bottom right corner
				{
				case Border::HalfUVH:
					border[2 * y + 2][2 * x + 2] = Border::HLine;
					break;
				case Border::HV:
					border[2 * y + 2][2 * x + 2] = Border::HalfBVH;
					break;
				case Border::HalfLHV:
					border[2 * y + 2][2 * x + 2] = Border::TRCorner;
					break;
				case Border::HalfRHV:
					border[2 * y + 2][2 * x + 2] = Border::TLCorner;
					break;
				case Border::BLCorner:
					border[2 * y + 2][2 * x + 2] = Border::UpperHLine;
					break;
				case Border::BRCorner:
					border[2 * y + 2][2 * x + 2] = Border::LowerHLine;
					break;
				default: break;
				}
			}
		}
	}
	//-------------------------------------------------------

	std::cout << "Maze" << std::endl;
	for (int y = 0; y < 2 * height + 1; ++y)
	{
		for (int x = 0; x < 2 * width + 1; ++x)
		{
			//Doubles the horizontal size of a printed cell
			const unsigned int number_of_prints = (x % 2 == 1) ? 2 : 1;
			for (int z = 0; z < number_of_prints; ++z)
			{

				switch (border[y][x])
				{
				case Border::Space:
					std::cout << " ";
					break;
				case Border::HLine:
					std::cout << u8"─";
					break;
				case Border::UpperHLine:
					std::cout << u8"╶";
					break;
				case Border::LowerHLine:
					std::cout << u8"╴";
					break;
				case Border::VLine:
					std::cout << u8"│";
					break;
				case Border::UpperVLine:
					std::cout << u8"╵";
					break;
				case Border::LowerVLine:
					std::cout << u8"╷";
					break;
				case Border::TLCorner:
					std::cout << u8"┌";
					break;
				case Border::TRCorner:
					std::cout << u8"┐";
					break;
				case Border::BLCorner:
					std::cout << u8"└";
					break;
				case Border::BRCorner:
					std::cout << u8"┘";
					break;
				case Border::HalfBVH:
					std::cout << u8"┬";
					break;
				case Border::HalfUVH:
					std::cout << u8"┴";
					break;
				case Border::HalfRHV:
					std::cout << u8"├";
					break;
				case Border::HalfLHV:
					std::cout << u8"┤";
					break;
				case Border::HV:
					std::cout << u8"┼";
					break;
				}
			}
		}
		std::cout << std::endl;
	}
}