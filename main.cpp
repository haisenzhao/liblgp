// main.cpp : Defines the entry point for the console application.
//

#include "iostream"
#include "pgl_functs.hpp"

using namespace std;
using namespace PGL;

int main(int argc, char* argv[])
{
	std::cerr << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << Functs::WinGetUserName() << std::endl;
	Functs::MAssert("Test pgl library...");
	system("pause");
	return 0;
}

