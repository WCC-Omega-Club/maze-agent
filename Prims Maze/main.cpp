#include "MazeAgent.hpp"
#include <chrono>
#include <iostream>
#include <iterator>


void PrintAnalytics(std::chrono::time_point<std::chrono::steady_clock> t1, std::chrono::time_point<std::chrono::steady_clock> t2);
typedef std::chrono::high_resolution_clock Clock;
int main() {
	Agent * agent = new Agent();
	
	
	cout << "			    Begin Breadth First Search				      " << endl;
	
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for E1 (first exit)" << endl << endl;
	std::chrono::time_point<std::chrono::steady_clock> sbf1 = Clock::now();
	agent->StartSearch(tileType::E1);	
	std::chrono::time_point<std::chrono::steady_clock> ebf1 = Clock::now();
	PrintAnalytics(sbf1, ebf1);
	agent->CleanUp();
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for E2 (second exit)" << endl << endl;
	std::chrono::time_point<std::chrono::steady_clock> sbf2 = Clock::now();
	agent->StartSearch(tileType::E2);		
	std::chrono::time_point<std::chrono::steady_clock> ebf2 = Clock::now();
	PrintAnalytics(sbf2, ebf2);
	agent->CleanUp();
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for Corner Exit" << endl << endl;
	std::chrono::time_point<std::chrono::steady_clock> sbf3 = Clock::now();
	agent->StartSearch(tileType::CORNER);
	std::chrono::time_point<std::chrono::steady_clock> ebf3 = Clock::now();
	PrintAnalytics(sbf3, ebf3);
	agent->CleanUp();
	agent->SwitchSearch();		

	cout << "\n=========================================================\n" << endl;
	cout << "				Begin Depth First Search					" << endl;

	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for E1 (first exit)" << endl << endl;
	std::chrono::time_point<std::chrono::steady_clock> sdf1 = Clock::now();
	agent->StartSearch(tileType::E1);
	std::chrono::time_point<std::chrono::steady_clock> edf1 = Clock::now();
	PrintAnalytics(sdf1, edf1);
	agent->CleanUp();
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for E2 (second exit)" << endl << endl;
	std::chrono::time_point<std::chrono::steady_clock> sdf2 = Clock::now();
	agent->StartSearch(tileType::E2);
	std::chrono::time_point<std::chrono::steady_clock> edf2 = Clock::now();
	PrintAnalytics(sdf2, edf2);
	agent->CleanUp();
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for Corner Exit" << endl << endl;
	std::chrono::time_point<std::chrono::steady_clock> sdf3 = Clock::now();
	agent->StartSearch(tileType::CORNER);
	std::chrono::time_point<std::chrono::steady_clock> edf3 = Clock::now();
	PrintAnalytics(sdf3, edf3);
	agent->CleanUp();
	agent->SwitchSearch();	
	

	cout << "\n=========================================================\n" << endl;
	cout << "				Begin A* Search								  " << endl;

	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for E1 (first exit)" << endl;
	std::chrono::time_point<std::chrono::steady_clock> sa1 = Clock::now();
	agent->StartSearch(tileType::E1);
	std::chrono::time_point<std::chrono::steady_clock> ea1 = Clock::now();
	PrintAnalytics(sa1, ea1);
	agent->CleanUp();
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for second exit (E2)" << endl;
	std::chrono::time_point<std::chrono::steady_clock> sa2 = Clock::now();
	agent->StartSearch(tileType::E2);
	std::chrono::time_point<std::chrono::steady_clock> ea2 = Clock::now();
	PrintAnalytics(sa2, ea2);
	agent->CleanUp();
	cout << "\n---------------------------------------------------------\n" << endl;
	cout << "Search for Corner Exit" << endl;
	std::chrono::time_point<std::chrono::steady_clock>	sa3 = Clock::now();
	agent->StartSearch(tileType::CORNER);
	std::chrono::time_point<std::chrono::steady_clock> ea3 = Clock::now();
	PrintAnalytics(sa3, ea3);
	agent->CleanUp();

	cout << "\n\n Challenge: Write a faster search than the current A* search implementation for a 1000x1000 maze" << endl;

	
	delete agent;

	return 3;
}

void PrintAnalytics(std::chrono::time_point<std::chrono::steady_clock> t1, std::chrono::time_point<std::chrono::steady_clock> t2) {
	
	cout << "\n"
	<< std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
	<< " milliseconds" << std::endl;
}