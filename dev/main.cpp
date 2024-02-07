#define GLM_ENABLE_EXPERIMENTAL
#include "iostream"
#include "liblgp.hpp"
#include "RI.hpp"
#include "tinyxml2.hpp"

using namespace std;
using namespace liblgp;

int main(int argc, char* argv[])
{
	std::cerr <<"WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;


	Vector3d bb = FF::RotationAxis(Vector3d(-1,0,1), FF::Angle2Radian(135.0), Vector3d(0,1,0));

	auto m = FF::RotationMatrix(Vector3d(0, 1, 0), FF::Angle2Radian(135.0));
	Vector3d bb1 = FF::PosApplyMatrix(Vector3d(-1, 0, 1), m);

	auto mm = glm::rotate(FF::Angle2Radian(135.0), Vector3d(0, 1, 0));
	auto rtv = mm* glm::dvec4(Vector3d(-1, 0, 1), 1.0);
	auto bb2 = Vector3d(rtv[0], rtv[1], rtv[2]);

	auto a = Functs::DetectExisting("E:\\file.txt");
	Functs::CerrLine(std::to_string(a));

	system("pause");
	return 0;
}








