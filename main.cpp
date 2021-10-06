// main.cpp : Defines the entry point for the console application.
//

#include "iostream"
#include "pgl_functs.hpp"
#include "RI.hpp"
#include "tinyxml2.hpp"

using namespace std;
using namespace PGL;

int main(int argc, char* argv[])
{
	std::cerr <<"WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;
	
	std::cerr<<Functs::DetectExisting("E:\\ece")<<std::endl;
	std::cerr << Functs::DetectExisting("E:\\ece\\") << std::endl;

	system("pause");
	return 0;
}








