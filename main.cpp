// main.cpp : Defines the entry point for the console application.
//

#include "iostream"
#include "pgl_functs.hpp"

using namespace std;
using namespace PGL;

int main(int argc, char* argv[])
{
	std::cerr <<"WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;
	
	std::string path("E:/Dropbox/Mold/microstructures");

	auto files = Functs::GetFilesInDirectory("E:/Dropbox/Mold/microstructures");
	std::cerr << path << std::endl;
	std::cerr << files.size() << std::endl;
	for (auto file : files)
	{
		std::cerr << file << std::endl;
	}

	Functs::MAssert("Test PGL library...");
	system("pause");
	return 0;
}








