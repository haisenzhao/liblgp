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
	
	auto a = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible");
	auto b = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible2");
	Functs::CerrLine(a);
	Functs::CerrLine(a);

	system("pause");
	return 0;
}








