This code is header-only library developed by Haisen Zhao.

# Dependency

Depend on [glm](https://github.com/g-truc/glm.git) and [eigen](https://github.com/libigl/eigen.git) but you don't need to install them explicitly.


# Call liblgp

## Call liblgp in your project

- First, download this  repository "[liblgp](https://github.com/haisenzhao/liblgp)" to your computer.
- Open the "liblgp" folder and copy all files and sub-folders.<br> <img src="dev/images/1.png" width = "50%" />
- Open the folder where your own project is located, then paste the files.<br> <img src="dev/images/2.png" width = "50%" />
- To call liblgp in your code, you can run the following code to check if it can be called.


```cpp
#include "iostream"
#include "liblgp.hpp"
#include "RI.hpp"
#include "tinyxml2.hpp"

using namespace std;
using namespace liblgp;

int main(int argc, char* argv[])
{
	std::cerr << "WinGetCurDirectory: " << Functs::WinGetCurDirectory() << std::endl;
	std::cerr << "WinGetUserName: " << Functs::WinGetUserName() << std::endl;
	system("pause");
	return 0;
}
```


# Usage in Cmake


## Cmake option
- When you use Cmake, you can choose third-party library or self-content library.<br>
 <img src="dev/images/cmake_option.png" width = "80%" />

# License
All rights about the program are reserved by the authors of this project. The programs can only be used for research purposes. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program.
