// main.cpp : Defines the entry point for the console application.
//

#include "iostream"
#include "pgl_functs.hpp"
#include "RI.hpp"
#include "tinyxml2.hpp"

//#include <sys/stat.h>
//#include <string>
//#include <fstream>

using namespace std;
using namespace PGL;

int main(int argc, char* argv[])
{
	std::cerr <<"WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;

   // CC().LoadConfig(argc, argv);
    
	auto a = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible");
	auto b = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible2");
	Functs::CerrLine(std::to_string(a));
	Functs::CerrLine(std::to_string(b));

	system("pause");
	return 0;
}








