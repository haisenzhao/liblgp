#ifndef _RI_H
#define _RI_H

#include <iostream>
#include <map>
#include "pgl_functs.hpp"
#include <random>
#include <unordered_map>
#include <unordered_set>

using namespace PGL;


//Resume IO Info
struct RI
{
public:

	int numb = 0;

	std::ofstream ofile;
	std::ifstream ifile;
	bool ol;

	RI(const bool ol_, const string path);
	void Close();
	int Check();
	void Read(const string t, const string l, const string& name);

	int NB(const string& name, const int t);
	void Out(const string& name, int& t);
	void Out(const string& name, double& t);
	void Out(const string& name, bool& t);
	void Out(const string& name, string& t);
	void Out(const string& name, const string& pt, string& t);
	void Out(const string& name, VectorTI3& t);
	void Out(const string& name, Vector1i1& t);
	void Out(const string& name, Vector1d1& t);
	void Out(const string& name, Vector1i2& t);
	void Out(const string& name, Vector1d2& t);
	void Out(const string& name, VectorPI1& t);
	void Out(const string& name, std::vector<string>& t);
	void Out(const string& name, const string& pt, std::vector<string>& t);
	void Out(const string& name, std::vector <std::vector<string>>& t);
	void Out(const string& name, VectorPI2& t);
	void Out(const string& name, std::map<int, int>& t);
	void Out(const string& name, std::map<string, string>& t);
	void Out(const string& name, std::map<string, int>& t);
	void Out(const string& name, std::vector<std::pair<string, int>>& t);
	void Out(const string& name, std::map<string, bool>& t);
	void Out(const string& name, std::map<int, std::tuple<int, int, int, int, int>>& t);
	void Out(const string& name, std::map <string, std::map<string, bool>>& t);
	void Out(const string& name, std::map <string, Vector1i1>& t);
	void Out(const string& name, std::map<int, TI3>& t);
	void Out(const string& name, std::vector<std::map <string, int>>& t);
	void Out(const string& name, std::vector<std::map<int, TI3>>& t);
	void Out(const string& name, std::vector<std::tuple<int, int, Vector1i1>>& t);
	void Out(const string& name, std::map<std::pair<int, int>, double>& t);
	void Out(const string& name, Vector2d& t);
	void Out(const string& name, Vector3d& t);
	void Out(const string& name, glm::dmat4& t);
	void Out(const string& name, Vector3d2& t);
	void Out(const string& name, Vector2d1& t);
	void Out(const string& name, Vector3d1& t);
	void Out(const string& name, std::vector<std::pair<Vector3d, Vector3d>>& t);
	void Out(const string& name, std::vector<std::tuple<Vector3d, Vector3d, Vector3d, Vector3d>>& t);
	void Out(const string& name, std::vector<std::pair<Vector3d1, double>>& t);
};

#endif // _PCE_COMBINE_H
