// main.cpp : Defines the entry point for the console application.
//

#include "iostream"
#include "pgl_functs.hpp"
#include "RI.hpp"
#include "tinyxml2.hpp"

using namespace std;
using namespace PGL;

//#include <sys/stat.h>
//#include <string>
//#include <fstream>

inline bool exists_test3(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}


int main(int argc, char* argv[])
{
	std::cerr <<"WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;
	
	auto a = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible");
	auto b = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible2");
	Functs::CerrLine(std::to_string(a));
	Functs::CerrLine(std::to_string(b));

	auto c = exists_test3("C:\\Working");
	auto d = exists_test3("C:\\expert.txt");
	Functs::CerrLine(std::to_string(c));
	Functs::CerrLine(std::to_string(d));

	system("pause");
	return 0;
}








