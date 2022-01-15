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

void CClearFolder(const std::string& path)
{
	if (!FF::DetectExisting(path))
	{
		auto folders = FF::SplitStr(FF::StringReplace(path, "\\", "/"), "/");
		if (folders.back() == "")
			folders.erase(folders.begin()+folders.size()-1);
		std::string str;
		for (int i = 0; i < folders.size(); i++)
		{
			str += folders[i] + "/";
			if (!FF::DetectExisting(str))
			{
				if (_mkdir(str.c_str())) {};
			}
		}
	}
	else
	{
		std::string del_cmd = "del /f/s/q " + path + " > nul";
		system(del_cmd.c_str());
		std::string rmdir_cmd = "rmdir /s/q " + path;
		system(rmdir_cmd.c_str());

		if (_mkdir(path.c_str())) {};
	}
}


int main(int argc, char* argv[])
{
	std::cerr <<"WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;
	
	CClearFolder("C:\\123\\456");

   // CC().LoadConfig(argc, argv);
    
	auto a = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible");
	auto b = Functs::DetectExisting("E:\\Dropbox\\Mold\\microstructures\\demo\\pipeline\\3_accessible2");
	Functs::CerrLine(std::to_string(a));
	Functs::CerrLine(std::to_string(b));

	system("pause");
	return 0;
}








