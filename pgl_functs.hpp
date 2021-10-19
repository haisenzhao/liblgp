#ifndef pgl_hpp
#define pgl_hpp
#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <map>
#include <sstream>
#include <iostream>
#include <math.h>
#include <random>
#include <fstream>
#ifdef __APPLE__
#include <sys/uio.h>
#include <unistd.h>
#else
#include <io.h>
#include <direct.h>
#include <windows.h>
#include <tchar.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <locale>
#include <codecvt>
#include <chrono>
#include <thread>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/random.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>


//glm modules
//http://glm.g-truc.net/0.9.8/api/modules.html

template <typename datum>
using Vector1 = std::vector<datum>;

template <typename datum>
using Vector2 = std::vector<std::vector<datum>>;

template <typename datum>
using Vector3 = std::vector<std::vector<std::vector<datum>>>;

typedef glm::highp_dvec2 Vector2d;
typedef glm::highp_dvec3 Vector3d;
typedef glm::highp_ivec2 Vector2i;
typedef glm::highp_ivec3 Vector3i;

typedef Vector1<Vector2d> Vector2d1;
typedef Vector2<Vector2d> Vector2d2;
typedef Vector3<Vector2d> Vector2d3;

typedef Vector1<Vector3d> Vector3d1;
typedef Vector2<Vector3d> Vector3d2;
typedef Vector3<Vector3d> Vector3d3;

typedef Vector1<bool> Vector1b1;
typedef Vector2<bool> Vector1b2;
typedef Vector3<bool> Vector1b3;

typedef Vector1<int> Vector1i1;
typedef Vector2<int> Vector1i2;
typedef Vector3<int> Vector1i3;

typedef Vector1<double> Vector1d1;
typedef Vector2<double> Vector1d2;
typedef Vector3<double> Vector1d3;

typedef Vector1<std::string> VectorStr1;
typedef Vector2<std::string> VectorStr2;
typedef Vector3<std::string> VectorStr3;


typedef Vector1<Vector2i> Vector2i1;
typedef Vector2<Vector2i> Vector2i2;
typedef Vector3<Vector2i> Vector2i3;

typedef Vector1<Vector3i> Vector3i1;
typedef Vector2<Vector3i> Vector3i2;
typedef Vector3<Vector3i> Vector3i3;

typedef Vector1<std::pair<int, int>> VectorPI1;
typedef Vector2<std::pair<int, int>> VectorPI2;
typedef Vector3<std::pair<int, int>> VectorPI3;

typedef Vector1<std::pair<bool, bool>> VectorPB1;
typedef Vector2<std::pair<bool, bool>> VectorPB2;
typedef Vector3<std::pair<bool, bool>> VectorPB3;

typedef std::tuple<int, int, int> TI3;
typedef Vector1<std::tuple<int, int, int>> VectorTI3;


using namespace std;

namespace PGL {

	static double DOUBLE_EPSILON = 1.0E-05;
	static double Math_PI = 3.14159265358979323846;
	static double SINGLE_EPSILON = 1.0E-05f;
	static double MAXDOUBLE = 100000000000.0;
	static std::random_device MATHRD;  //Will be used to obtain a seed for the random number engine
	static std::mt19937 MATHGEN(0); //Standard mersenne_twister_engine seeded with rd()
	static std::string CERR_ITER = "  ";

	struct TimeClock
	{
	public:
		TimeClock():start(0), end(0), duration(0){ start = clock(); };
		void StartClock(){ start = clock();};
		double EndClock()
		{
			end = clock();
			duration = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
			return duration;
		};
		double start, end, duration;
	};


	class Functs
	{
	public:
	
#pragma region StatisticsCombinationSet

		template <class Type>
		static Type Variance(const std::vector<Type>& resultSet)
		{
			Type sum = std::accumulate(std::begin(resultSet), std::end(resultSet), 0.0);
			Type mean = sum / resultSet.size(); //均值
			Type accum = 0.0;
			std::for_each(std::begin(resultSet), std::end(resultSet), [&](const Type d) {accum += (d - mean) * (d - mean); });
			Type stdev = sqrt(accum / (resultSet.size() - 1)); //方差

			return stdev;
		}

		//nb>=1
		static int Factorial(const int& n)
		{
			if (n > 1)
				return n * Factorial(n - 1);
			else
				return 1;
		};

		//arrangement and combination
		static void Combination(const Vector1i3& combs, const int& nb, const Vector1i2& sequence, const int& max_nb, Vector1i3& output)
		{
			if (output.size() > max_nb)return;

			if (nb == combs.size())
			{
				output.emplace_back(sequence);
				return;
			}

			for (int i = 0; i < combs[nb].size(); i++)
			{
				Vector1i2 seq = sequence;
				seq.emplace_back(combs[nb][i]);
				Combination(combs, nb + 1, seq, max_nb, output);
			}
		}

		//arrangement and combination
		//goups={3,2,3}
		//0,0,0 //0,0,1 //0,0,2 //0,1,0 //0,1,1 //0,1,2
		//1,0,0 //1,0,1 //1,0,2 //1,1,0 //1,1,1 //1,1,2
		//2,0,0 //2,0,1 //2,0,2 //2,1,0 //2,1,1 //2,1,2
		static Vector1i2 Selection(const Vector1i1& groups)
		{
			int nb = 1;
			for (auto& group : groups) nb = nb * group;

			if (nb == 0)
			{
				Functs::MAssert("nb == 0 in Selection(const Vector1i1& groups)");
			}

			Vector1i1 nbs(1, nb);
			for (auto& group : groups)
			{
				nbs.emplace_back(nbs.back() / group);
			}
			nbs.erase(nbs.begin());

			Vector1i2 combs;
			for (int i = 0; i < nb; i++)
			{
				int id = i;
				combs.emplace_back(Vector1i1());
				for (auto nbs_ : nbs)
				{
					int a = (id - id % nbs_) / nbs_;
					combs.back().emplace_back(a);
					id = id - a * nbs_;
				}
			}
			return combs;
		}

		//arrangement and combination
		//with repeat selection
		//n=3,m=2
		//1,2
		//1,3
		//2,3
		static Vector1i2 CombNonRepeat(const int& n, const int& m)
		{
			std::vector<std::vector<int>> combs;
			if (n > m)
			{
				vector<int> p, set;
				p.insert(p.end(), m, 1);
				p.insert(p.end(), static_cast<int64_t>(n) - m, 0);
				for (int i = 0; i != p.size(); ++i)
					set.push_back(i + 1);
				vector<int> vec;
				size_t cnt = 0;
				do {
					for (int i = 0; i != p.size(); ++i)
						if (p[i])
							vec.push_back(set[i]);
					combs.emplace_back(vec);
					cnt++;
					vec.clear();
				} while (prev_permutation(p.begin(), p.end()));
			}
			else
			{
				combs.emplace_back(std::vector<int>());
				for (int i = 1; i <= n; i++)
					combs.back().emplace_back(i);
			}

			return combs;
		}

		//arrangement and combination
		//N=3,K=2
		//0,0
		//0,1
		//0,2
		//1,1
		//1,2
		//2,2
		static std::vector<std::vector<int>> CombRepeat(const int& N, const int& K)
		{
			auto Combination = [](const int& N, const int& K) {
				std::vector<std::vector<int>> temps;
				std::string bitmask(K, 1); // K leading 1's
				bitmask.resize(N, 0); // N-K trailing 0's
				// print integers and permute bitmask
				do {
					std::vector<int> temp;
					for (int i = 0; i < N; ++i) // [0..N-1] integers
						if (bitmask[i]) temp.emplace_back(i);
					if (!temp.empty()) temps.emplace_back(temp);
				} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
				return temps;
			};

			std::vector<std::vector<int>> combs = Combination(N + K - 1, K);

			std::vector<std::vector<int>> repeatCombs;
			for (auto& comb : combs) {
				std::vector<int> nbs((int)(static_cast<int64_t>(N) + static_cast<int64_t>(K) - 1), 1);
				for (auto& i : comb) nbs[i] = 0;
				std::vector<int> temp;
				int sum = 0;
				for (auto& nb : nbs) {
					if (nb == 0)temp.emplace_back(sum);
					sum += nb;
				}
				if (!temp.empty()) repeatCombs.emplace_back(temp);
			}
			return repeatCombs;
		}

		//remove duplicated elements
		//vec={3,2,1,2,3,4}
		//output:{1,2,3,4}
		static Vector1i1 UniqueSet(const Vector1i1& vec)
		{
			Vector1i1 s;
			for (auto& v : vec)
			{
				if (std::find(s.begin(), s.end(), v) == s.end()) s.emplace_back(v);
			}
			std::sort(s.begin(), s.end());
			return s;
		};

		static Vector1i1 SetUnion(const Vector1i1& first, const Vector1i1& second)
		{
			Vector1i1 v = first;
			for (auto& s : second)
				if (std::find(first.begin(), first.end(), s) == first.end())
					v.emplace_back(s);
			return v;
		}

		static Vector1i1 SetUnion(Vector1i2& sets)
		{
			Vector1i1 start = sets[0];

			for (int i = 1; i < sets.size(); i++)
			{
				start = SetUnion(start, sets[i]);
			}

			return start;
		}

		static Vector1i1 SetIntersection(const Vector1i1& first, const Vector1i1& second)
		{
			Vector1i1 v;
			for (auto& s : second)
				if (std::find(first.begin(), first.end(), s) != first.end())
					v.emplace_back(s);
			return v;
		}
		static Vector1i1 SetSubtraction(const Vector1i1& first, const Vector1i1& second)
		{
			Vector1i1 v;
			for (auto& s : first)
				if (std::find(second.begin(), second.end(), s) == second.end())
					v.emplace_back(s);
			return v;
		}


		static Vector1i1 FindSetCombination(Vector1i2& input_)
		{
			auto FSC = [](vector<set<int>>& input, set<int>& target, vector<int>& output)
			{
				set<int> full;
				for (auto it : input) {
					full.insert(it.begin(), it.end());
				}

				if (!includes(full.begin(), full.end(), target.begin(), target.end())) {
					return;
				}

				for (int i = static_cast<int>(input.size()) - 1; i > 0; --i) {
					vector<bool> vec(input.size(), false);
					fill(vec.begin() + i, vec.end(), true);
					set<int> comb;

					do {
						for (int j = 0; j < vec.size(); ++j) {
							if (vec[j]) {
								comb.insert(input[j].begin(), input[j].end());
							}
						}

						if (includes(comb.begin(), comb.end(), target.begin(), target.end())) {
							for (int j = 0; j < vec.size(); ++j) {
								if (vec[j]) {
									output.push_back(j);
								}
							}
							return;
						}
						comb.clear();

					} while (next_permutation(vec.begin(), vec.end()));
				}
			};


			vector<set<int>> input;
			set<int> target;
			for (auto& a : input_)
			{
				for (int x : a) target.insert(x);
				input.emplace_back(ConvertToSet(a));
			}
			Vector1i1 output;
			FSC(input, target, output);
			return output;
		}



#pragma endregion


#pragma region StringDataStructure

		static int INe(const int& i)
		{
			return i + 1;
		}

		static int IPr(const int& i)
		{
			return i - 1;
		}

		static bool StringContain(const string& str, const string& sub)
		{
			return str.find(sub) != std::string::npos;
		}

		static std::string StringReplace(const string& source, const string& toReplace,const string& replaceWith)
		{
			size_t pos = 0;
			size_t cursor = 0;
			int repLen = (int)toReplace.length();
			stringstream builder;

			do
			{
				pos = source.find(toReplace, cursor);

				if (string::npos != pos)
				{
					//copy up to the match, then append the replacement
					builder << source.substr(cursor, pos - cursor);
					builder << replaceWith;

					// skip past the match 
					cursor = pos + repLen;
				}
			} while (string::npos != pos);

			//copy the remainder
			builder << source.substr(cursor);

			return (builder.str());
		}

		
		template <class Type>
		static std::string IntString(Type& i)
		{
			std::stringstream ss;
			std::string str;
			ss << i;
			ss >> str;
			return str;
		}

		template <class Type>
		static std::string IntString(const std::vector<Type>& vecs, bool order = false, std::string insert_str = "")
		{
			std::vector<Type> a = vecs;
			if (order)sort(a.begin(), a.end());

			std::string str;
			for (int i = 0; i < a.size(); i++)
			{
				str += IntString(a[i]);
				if (i != a.size() - 1) str += insert_str;
			}
			return str;

		}

		template <class Type>
		static std::string IntString(const std::vector<std::vector<Type>>& vecs, bool order = false, std::string insert_str_0 = "", std::string insert_str_1 = "")
		{
			std::string str;

			for (int i = 0; i < vecs.size(); i++)
			{
				str += IntString(vecs[i], order, insert_str_0);
				if (i != vecs.size() - 1) str += insert_str_1;
			}

			return str;
		}
		template <class Type>
		static std::string IntString(const std::vector<std::vector<std::vector<Type>>>& vecs, bool order = false, std::string insert_str_0 = "", std::string insert_str_1 = "", std::string insert_str_2 = "")
		{
			std::string str;

			for (int i = 0; i < vecs.size(); i++)
			{
				str += IntString(vecs[i], order, insert_str_0, insert_str_1);
				if (i != vecs.size() - 1) str += insert_str_2;
			}

			return str;
		}

		template <class Type>
		static std::string VectorString(const Type& v, const string insert_str="", int p = 3)
		{
			std::string str;
			for (int i = 0; i < v.length(); i++)
			{
				if(i== v.length()-1)
					str += DoubleString(v[i], p);
				else
					str += DoubleString(v[i], p) + insert_str;
			}
			return str;
		};

		static std::string VectorString(const Vector3d1& vecs)
		{
			auto Comp = [](Vector3d& v_0, Vector3d& v_1)
			{
				if (IsAlmostZero(abs(v_1[0] - v_0[0])))
				{
					if (IsAlmostZero(abs(v_1[1] - v_0[1])))
					{
						if (IsAlmostZero(abs(v_1[2] - v_0[2])))
							return true;
						return v_0[2] < v_1[2];
					}
					return v_0[1] < v_1[1];
				}
				return v_0[0] < v_1[0];
			};

			Vector3d1 vecs_1 = vecs;
			std::sort(vecs_1.begin(), vecs_1.end(), Comp);

			std::string str;
			for (auto& p : vecs_1)
			{
				double x = floor(p[0] * 10.0f + 0.5) / 10.0f;
				double y = floor(p[1] * 10.0f + 0.5) / 10.0f;
				double z = floor(p[2] * 10.0f + 0.5) / 10.0f;
				str += DoubleString(x);
				str += DoubleString(y);
				str += DoubleString(z);
			}

			return str;
		}

		static std::string VectorString(const Vector3d3& vecs_3)
		{
			Vector3d1 vecs_1;
			for (int i = 0; i < vecs_3.size(); i++)
				for (int j = 0; j < vecs_3[i].size(); j++)
					for (int k = 0; k < vecs_3[i][j].size(); k++)
						vecs_1.emplace_back(vecs_3[i][j][k]);

			return VectorString(vecs_1);
		};

		static std::string IntString(const VectorPI2& vecs, bool order = false,
			const std::string insert_str_0 = "", const std::string insert_str_1 = "", const std::string insert_str_2 = "")
		{
			std::string str;
			for (auto& vec : vecs)
			{
				str += IntString(vec, order, insert_str_0, insert_str_1) + insert_str_2;
			}
			return str;
		}

		static std::string IntString(const VectorPI1& vecs, bool order = false, const std::string insert_str_0 = "", const std::string insert_str_1 = "")
		{
			if (order)
			{
				std::vector<std::pair<int, int>> a = vecs;
				sort(a.begin(), a.end());
				std::string str;
				for (int i = 0; i < a.size(); i++)
				{
					str += IntString(a[i].first) + insert_str_0 + IntString(a[i].second);
					if (i != a.size() - 1)
						str += insert_str_1;
				}
				return str;
			}
			else
			{
				std::string str;
				for (int i = 0; i < vecs.size(); i++)
				{
					str += IntString(vecs[i].first) + insert_str_0 + IntString(vecs[i].second);
					if (i != vecs.size() - 1) str += insert_str_1;
				}
				return str;
			}
		}

		static std::string DoubleString(const double& d, int p = 8)
		{
			//double d_ = abs(d);
			//double d0 = (int)d_;
			//double d1 = d_ - d0;

			//std::string str;
			//if (Functs::IsAlmostZero(d1))
			//	str = std::to_string((int)d_);
			//else
			//{
			//	str = std::to_string((int)d_) + ".";
			//	for (int i = 1; i <= p; i++)
			//	{
			//		int a = (int)(d1 * pow(10, i));
			//		str += std::to_string(a);
			//		d1 = d1 - (double)a / (double)pow(10, i);
			//	}
			//	return str;
			//}

			//if (d >= 0)
			//	return str;
			//else
			//	return "-" + str;



			double d_ = abs(d);
			double d0 = (int)d_;
			double d1 = d_ - d0;

			{
				double d1_ = d1;
				double d2_ = 0.0;
				if (!Functs::IsAlmostZero(d1_))
				{
					for (int i = 1; i <= p; i++)
					{
						int a = (int)(d1_ * pow(10, i));
						d1_ = d1_ - (double)a / (double)pow(10, i);
						d2_ += (double)a / (double)pow(10, i);
					}

					auto temp_0 = std::pow(10, -p);
					auto temp_1 = abs(temp_0 - d1_);

					if (Functs::IsAlmostZero_Double(temp_1, pow(10, -4)))
					{
						d1 = d2_ + temp_0;
						d_ += (int)d1;
						d1 = d1 - (int)d1;
					}
				}
			}

			std::string str;
			if (Functs::IsAlmostZero(d1))
				str = std::to_string((int)d_);
			else
			{
				str = std::to_string((int)d_) + ".";
				for (int i = 1; i <= p; i++)
				{
					int a = (int)(d1 * pow(10, i));
					str += std::to_string(a);
					d1 = d1 - (double)a / (double)pow(10, i);
				}
				return str;
			}

			if (d >= 0)
				return str;
			else
				return "-" + str;
		}

		static std::string DoubleString(const Vector1d1& ds_, int p = 8, bool order = false, const std::string insert_str_0 = "")
		{
			Vector1d1 ds = ds_;
			if (order)std::sort(ds.begin(), ds.end());
			std::string str;
			for (int i = 0; i < ds.size(); i++)
				str += DoubleString(ds[i], p) + (i == ds.size() - 1 ? "" : insert_str_0);
			return str;
		}

		static double StringToDouble(const string& str)
		{
			istringstream iss(str);
			double num;
			iss >> num;
			return num;
		}
		
		template <class Type>
		static Type StringToNum(const string& str)
		{
			istringstream iss(str);
			Type num;
			iss >> num;
			return num;
		}
	

		static vector<string> SplitStr(const string& str, const char& delimiter)
		{
			vector<string> internal_strs;
			stringstream ss(str); // Turn the string into a stream.
			string tok;
			while (getline(ss, tok, delimiter))
				internal_strs.emplace_back(tok.c_str());
			return internal_strs;
		}

		static vector<string> SplitStr(const string& str, const string& delimiter)
		{
			vector<string> internal_strs;
			std::string s = str;
			size_t pos = 0;
			std::string token;
			while ((pos = s.find(delimiter)) != std::string::npos) {
				token = s.substr(0, pos);
				internal_strs.push_back(token);
				s.erase(0, pos + delimiter.length());
			}
			internal_strs.push_back(s);
			return internal_strs;
		}
		static vector<double> SplitD(const string& str, const char& delimiter)
		{
			vector<double> internal_d;
			stringstream ss(str); // Turn the string into a stream.
			string tok;
			while (getline(ss, tok, delimiter))
				internal_d.emplace_back(atof(tok.c_str()));
			return internal_d;
		};

		static vector<int> SplitI(const string& str, const char& delimiter) {
			vector<int> internal_d;
			stringstream ss(str); // Turn the string into a stream.
			string tok;
			while (getline(ss, tok, delimiter))
				internal_d.emplace_back(atoi(tok.c_str()));
			return internal_d;
		};

		static Vector1i1 ShuffleVector(const int& size)
		{
			Vector1i1 shuffle_vec;
			for (int i = 0; i < size; i++)
				shuffle_vec.emplace_back(i);
			std::random_shuffle(shuffle_vec.begin(), shuffle_vec.end());
			return shuffle_vec;
		};

		static set<int> ConvertToSet(const vector<int>& v)
		{
			set<int> s;
			for (int x : v) s.insert(x);
			return s;
		};

		static vector<int> ConvertToVector(const set<int>& s)
		{
			vector<int> v;
			for (auto x : s)v.emplace_back(x);
			return v;
		}

		template <class T1, class T2>
		static T2 MapFind(const std::map<T1, T2>& mt, const T1& t1)
		{
			if (mt.find(t1) == mt.end())
				MAssert("if (mt.find(t1) == mt.end())");
			return mt[t1];
		}

		template <class T1, class T2>
		static bool MapContain(const std::map<T1, T2>& mt, const T1& t1)
		{
			return mt.find(t1) != mt.end();
		}

		static Vector1i1 RemoveDuplicate(const Vector1i1& vec_)
		{
			Vector1i1 vec = vec_;
			sort(vec.begin(), vec.end());
			vec.erase(unique(vec.begin(), vec.end()), vec.end());
			return vec;
		}

		static Vector3d EigenVector(const Eigen::Vector3d& vec)
		{
			return Vector3d(vec[0],vec[1],vec[2]);
		};

		static Vector3d1 EigenVector(const std::vector<Eigen::Vector3d>& vecs)
		{
			Vector3d1 vs;
			for (auto& vec : vecs)
				vs.push_back(EigenVector(vec));
			return vs;
		};

#pragma endregion

#pragma region BasicGeomFunctions

		template <class Type>
		static double GetLength(const Type& v) {
			return glm::length(v);
		}

		template <class Type>
		static double GetAngleBetween(const Type& v1, const Type& v2) {
			double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));
			if (IsAlmostZero(d - 1.0))
				return 0.0;
			if (IsAlmostZero(d + 1.0))
				return Math_PI;
			return glm::acos(d);
		}
		
		static double RadiantoAngle(const double& r)
		{
			return r / Math_PI*180.0;
		}
		
		static double Radian2Angle(const double& radian)
		{
			return radian / Math_PI * 180.0;
		}

		static double Angle2Radian(const double& angle)
		{
			return angle / 180.0 * Math_PI;
		}

		static Vector3d SetVectorLength(Vector3d& v, const double& length)
		{
			double l = GetLength(v);

			v[0] = v[0] / l * length;
			v[1] = v[1] / l * length;
			v[2] = v[2] / l * length;

			return v;
		}

		static Vector2d SetVectorLength(Vector2d& v, const double& length)
		{
			double l = GetLength(v);
			v[0] = v[0] / l * length;
			v[1] = v[1] / l * length;
			return v;
		}

		//existing bugs in this function
		static Vector3d Vector3dBase(const Vector3d& v)
		{
			Vector3d n(1.0, 1.0, 1.0);
			if (!IsAlmostZero(v[0])) {
				n[0] = -(v[1] + v[2]) / v[0];
				return n;
			}
			if (!IsAlmostZero(v[1])) {
				n[1] = -(v[0] + v[2]) / v[1];
				return n;
			}
			if (!IsAlmostZero(v[2])) {
				n[2] = -(v[0] + v[1]) / v[2];
				return n;
			}
			return n;
		}

		static Vector3d GetCenter(const Vector3d1& points)
		{
			Vector3d center(0.0, 0.0, 0.0);
			for (int i = 0; i < points.size(); i++)
				center += points[i];
			center = center / (double)points.size();
			return center;
		}

		static Vector3d GetCenter(const Vector3d2& points)
		{
			Vector3d center(0.0, 0.0, 0.0);
			int nb = 0;
			for (int i = 0; i < points.size(); i++)
			{
				for (int j = 0; j < points[i].size(); j++)
					center += points[i][j];
				nb += static_cast<int>(points[i].size());
			}
			center = center / (double)nb;
			return center;
		}

		static Vector3d GetCenter(const Vector3d3& points)
		{
			Vector3d center(0.0, 0.0, 0.0);
			int nb = 0;
			for (int i = 0; i < points.size(); i++)
			{
				for (int j = 0; j < points[i].size(); j++)
				{
					for (int k = 0; k < points[i][j].size(); k++)
					{
						center += points[i][j][k];
						nb++;
					}
				}
			}
			center = center / (double)nb;
			return center;
		}

		static Vector2d GetCenter(const std::vector<Vector2d>& points)
		{
			Vector2d center(0.0, 0.0);
			for (int i = 0; i < points.size(); i++)
				center += points[i];
			center = center / (double)points.size();
			return center;
		}

		static Vector2d GetCenter(const std::vector<std::vector<Vector2d>>& points)
		{
			Vector2d center(0.0, 0.0);
			int nb = 0;
			for (int i = 0; i < points.size(); i++)
			{
				for (int j = 0; j < points[i].size(); j++)
					center += points[i][j];
				nb += static_cast<int>(points[i].size());
			}
			center = center / (double)nb;
			return center;
		}

		static void GetBoundingBox(const Vector3d2& points, Vector3d& minimal_corner, Vector3d& maximal_corner)
		{
			minimal_corner = Vector3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
			maximal_corner = Vector3d(-MAXDOUBLE, -MAXDOUBLE, -MAXDOUBLE);
			for (int i = 0; i < points.size(); i++)
			{
				for (int j = 0; j < points[i].size(); j++)
				{
					minimal_corner[0] = min(minimal_corner[0], points[i][j][0]);
					minimal_corner[1] = min(minimal_corner[1], points[i][j][1]);
					minimal_corner[2] = min(minimal_corner[2], points[i][j][2]);
					maximal_corner[0] = max(maximal_corner[0], points[i][j][0]);
					maximal_corner[1] = max(maximal_corner[1], points[i][j][1]);
					maximal_corner[2] = max(maximal_corner[2], points[i][j][2]);
				}
			}
		}

		static void GetBoundingBox(const Vector3d1& points, Vector3d& minimal_corner, Vector3d& maximal_corner)
		{
			minimal_corner = Vector3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
			maximal_corner = Vector3d(-MAXDOUBLE, -MAXDOUBLE, -MAXDOUBLE);

			for (int i = 0; i < points.size(); i++)
			{
				minimal_corner[0] = min(minimal_corner[0], points[i][0]);
				minimal_corner[1] = min(minimal_corner[1], points[i][1]);
				minimal_corner[2] = min(minimal_corner[2], points[i][2]);
				maximal_corner[0] = max(maximal_corner[0], points[i][0]);
				maximal_corner[1] = max(maximal_corner[1], points[i][1]);
				maximal_corner[2] = max(maximal_corner[2], points[i][2]);
			}
		}

		static void GetBoundingBox(const std::vector<Vector2d>& points, Vector2d& minimal_corner, Vector2d& maximal_corner)
		{
			minimal_corner = Vector2d(MAXDOUBLE, MAXDOUBLE);
			maximal_corner = Vector2d(-MAXDOUBLE, -MAXDOUBLE);
			for (int i = 0; i < points.size(); i++)
			{
				minimal_corner[0] = min(minimal_corner[0], points[i][0]);
				minimal_corner[1] = min(minimal_corner[1], points[i][1]);
				maximal_corner[0] = max(maximal_corner[0], points[i][0]);
				maximal_corner[1] = max(maximal_corner[1], points[i][1]);
			}
		}

		static void GetBoundingBox(const std::vector<std::vector<Vector2d>>& points, Vector2d& minimal_corner, Vector2d& maximal_corner)
		{
			minimal_corner = Vector2d(MAXDOUBLE, MAXDOUBLE);
			maximal_corner = Vector2d(-MAXDOUBLE, -MAXDOUBLE);
			for (int i = 0; i < points.size(); i++)
			{
				for (int j = 0; j < points[i].size(); j++)
				{
					minimal_corner[0] = min(minimal_corner[0], points[i][j][0]);
					minimal_corner[1] = min(minimal_corner[1], points[i][j][1]);
					maximal_corner[0] = max(maximal_corner[0], points[i][j][0]);
					maximal_corner[1] = max(maximal_corner[1], points[i][j][1]);
				}
			}
		}

		static double CircumCircleRaidius(const Vector2d& v0, const Vector2d& v1, const Vector2d& v2)
		{
			double a = GetLength(v0 - v1);
			double b = GetLength(v0 - v2);
			double c = GetLength(v1 - v2);
			double p = (a + b + c) / 2.0;
			double area = (4.0 * pow(p * (p - a) * (p - b) * (p - c), 0.5));
			double radius;

			if (IsAlmostZero(area))
			{
				double max_l = a;
				max_l = max(max_l, b);
				max_l = max(max_l, c);
				radius = 10 * max_l;
			}
			else
				radius = a * b * c / area;
			return radius;
		}

		static double GetTriangleArea(const Vector3d & v0, const Vector3d & v1, const Vector3d & v2)
		{
			double a = GetDistance(v0, v1);
			double b = GetDistance(v2, v1);
			double c = GetDistance(v0, v2);
			double p = (a + b + c) / 2.0;
			return sqrt(p * (p - a) * (p - b) * (p - c));
		}

		template <class Type>
		static double GetLength(const Type& v0, const Type& v1) {
			return GetLength(v0 - v1);
		}

		template <class Type>
		static double GetDistance(const Type& v0, const Type& v1) {
			return GetLength(v0 - v1);
		}

		template <class Type>
		static double GetDistance(const Type& v, const std::vector<Type>& vs)
		{
			double min_d = MAXDOUBLE;
			for (const auto& iter : vs)
				min_d = min(min_d, GetDistance(v, iter));
			return min_d;
		}

		template <class Type>
		static double GetDistance(const Type& v, const std::vector<std::vector<Type>>& vs)
		{
			double min_d = MAXDOUBLE;
			for (const auto& iter : vs)
				min_d = min(min_d, GetDistance(v, iter));
			return min_d;
		}

		template <class Type>
		static double GetDistance(const std::vector<std::vector<Type>>& vs1, const std::vector<std::vector<Type>>& vs2)
		{
			double total_d = 0.0;
			int nb = 0;
			for (auto& vec : vs1)
			{
				for (auto& v : vec)
				{
					total_d = total_d + GetDistance(v, vs2);
					nb++;
				}
			}
			if (nb != 0) total_d = total_d / (double)nb;
			return total_d;
		}

		template <class Type>
		static double GetDistanceNONE(
			const std::vector<std::vector<Type>>& vs1,
			const std::vector<std::vector<Type>>& vs2)
		{
			return (GetDistance(vs1, vs2) + GetDistance(vs2, vs1)) / 2.0;
		}

		template <class Type>
		static int GetNearestPointIndex(const Type& v, const std::vector<Type>& vs)
		{
			int index = -1;
			double min_d = MAXDOUBLE;
			for (int i = 0; i < vs.size(); i++)
			{
				double cur_d = GetDistance(v, vs[i]);
				if (cur_d < min_d)
				{
					min_d = cur_d;
					index = 0;
				}
			}
			return index;
		}

		static double GetLength(const std::vector<Vector2d>& points)
		{
			double length = 0.0;

			for (int i = 0; i < points.size(); i++)
				length += GetLength(points[i], points[(static_cast<int64_t>(i) + 1) % points.size()]);

			return length;
		}

		static double GetLength(const Vector3d1& points)
		{
			double length = 0.0;

			for (int i = 0; i < points.size(); i++)
				length += GetLength(points[i], points[(static_cast<int64_t>(i) + 1) % points.size()]);

			return length;
		}

		static Vector2d1 RemoveClosePoints(const Vector2d1& xys_)
		{
			Vector2d1 xys = xys_;
			std::vector<int> remove_int;
			if (xys.size() > 2)
			{
				for (int i = 0; i < xys.size() - 1; i++)
				{
					double d = GetDistance(xys[i], xys[(int)(i+1)]);
					if (d < 0.00001) remove_int.push_back(i + 1);
				}

				for (int i = (int)(remove_int.size() - 1); i >= 0; i--)
				{
					xys.erase(xys.begin() + remove_int[i]);
				}
			}
			return xys;
		}

		static Vector3d PlaneProject(const Vector3d& planar_location, Vector3d& planar_direction, const Vector3d& p)
		{
			if (IsAlmostZero(GetLength(planar_location, p)))
				return planar_location;

			double angle = GetAngleBetween(planar_direction, p - planar_location);
			double length = GetLength(planar_location, p);

			if (angle <= Math_PI / 2.0)
				return p - SetVectorLength(planar_direction, length * sin(Math_PI / 2.0 - angle));
			else
				return p + SetVectorLength(planar_direction, length * sin(angle - Math_PI / 2.0));

		}

		static Vector3d GetRandomDirection(const double& direction_length = 1.0)
		{
			double alpha_angle = rand() / double(RAND_MAX) * 2.0 * PGL::Math_PI;
			double alpha_beta = rand() / double(RAND_MAX) * 2.0 * PGL::Math_PI;
			auto direction_0 = PGL::Functs::RotationAxis(Vector3d(direction_length, 0.0, 0.0), alpha_angle, Vector3d(0.0, 1.0, 0.0));
			auto direction_axis = PGL::Functs::GetCrossproduct(direction_0, Vector3d(0.0, 1.0, 0.0));
			return PGL::Functs::RotationAxis(direction_0, alpha_beta, direction_axis);
		}



		//https://medium.com/@all2one/generating-uniformly-distributed-points-on-sphere-1f7125978c4c

		//method: NormalDeviate; TrigDeviate; CoordinateApproach;MinimalDistance;RegularDistribution;
		static Vector3d1 GetRandomDirections(const int& count, const string& method)
		{
			if (method == "NormalDeviate") return GetRandomDirections_Normal_Deviate(count);
			if (method == "TrigDeviate") return GetRandomDirections_Trig_method(count);
			if (method == "CoordinateApproach") return GetRandomDirections_Coordinate_Approach(count);
			if (method == "MinimalDistance") return GetRandomDirections_Minimal_Distance(count);
			if (method == "RegularDistribution") return GetRandomDirections_Regular_Distribution(count);

			MAssert("Input method does not be implemented: "+ method);
			return Vector3d1();
		}

		static Vector3d1 GetRandomDirections_Regular_Distribution(const int& count)
		{
			Vector3d1 directions;
			double a = 4.0 * Math_PI * 1.0 / static_cast<double>(count);
			double d = sqrt(a);
			size_t num_phi = (size_t)round(Math_PI / d);
			double d_phi = Math_PI / static_cast<double>(num_phi);
			double d_theta = a / d_phi;
			for (int m = 0; m < num_phi; ++m) {
				double phi = Math_PI * (m + 0.5) / num_phi;
				size_t num_theta = (size_t)round(2 * Math_PI * sin(phi) / d_theta);
				for (int n = 0; n < num_theta; ++n) {
					double theta = 2 * Math_PI * n / static_cast<double>(num_theta);
					Vector3d p;
					p.x = sin(phi) * cos(theta);
					p.y = sin(phi) * sin(theta);
					p.z = cos(phi);
					directions.push_back(p);
				}
			}
			return directions;
		}


		static Vector3d1 GetRandomDirections_Normal_Deviate(const int& count)
		{
			std::mt19937 rnd;
			std::normal_distribution<double> dist(0.0, 1.0);

			Vector3d1 directions;
			for (int i = 0; i < count; ++i)
			{
				bool bad_luck = false;
				do
				{
					double x = dist(rnd);
					double y = dist(rnd);
					double z = dist(rnd);
					double r2 = x * x + y * y + z * z;
					if (r2 == 0)
						bad_luck = true;
					else
					{
						bad_luck = false;
						double r = sqrt(r2);
						directions.push_back(Vector3d(x / r, y / r, z / r));
					}
				} while (bad_luck);
			}

			return directions;
		}

		static Vector3d1 GetRandomDirections_Trig_method(const int& count)
		{
			std::mt19937 rnd;
			std::uniform_real_distribution<double> dist(0.0, 1.0);

			Vector3d1 directions;
			for (int i = 0; i < count; ++i)
			{
				double z = 2.0 * dist(rnd) - 1.0;
				double t = 2.0 * Math_PI * dist(rnd);
				double r = sqrt(1.0 - z * z);
				directions.push_back(Vector3d(r * cos(t), r * sin(t), z));
			}
			return directions;
		};

		static Vector3d1 GetRandomDirections_Coordinate_Approach(const int& count)
		{
			std::mt19937 rnd;
			std::uniform_real_distribution<double> dist(-1.0, 1.0);

			Vector3d1 directions;
			for (int i = 0; i < count; ++i)
			{
				bool rejected = false;
				do
				{
					double u = dist(rnd);
					double v = dist(rnd);
					double s = u * u + v * v;
					if (s > 1.0)
						rejected = true;
					else
					{
						rejected = false;
						double a = 2.0 * sqrt(1.0 - s);
						directions.push_back(Vector3d(a * u, a * v, 2.0 * s - 1.0));
					}
				} while (rejected);
			}
			return directions;
		}

		//random sample a set of directions on the Gaussian Sphere
		static Vector3d1 GetRandomDirections_Minimal_Distance(const int& dns, const int dis_iters =100)
		{
			double gaussion_sphere_radius = 1.0;
			double idea_distance = 2 * gaussion_sphere_radius / sqrt(dns);

			Vector3d1 directions;
			for (int i = 0; i < dns; i++)
			{
				OutputIterInfo("Random Directions", dns, i, 10);

				for (int j = 0; j < dis_iters; j++)
				{
					Vector3d center(0.0, 0.0, 0.0);

					//glm::ballRand(gaussion_sphere_radius) is slower than my solution
					Vector3d random_direction = GetRandomDirection(gaussion_sphere_radius);

					auto dis = Functs::GetDistance(random_direction, directions);
					if (dis > idea_distance)
					{
						directions.push_back(random_direction);
						break;
					}
				}

			}
			return directions;
		}

		static void Connecting_Segments(const Vector3d2& segments, Vector3d2& lines)
		{
			//save connecting relations
			std::vector<bool> used(segments.size(), false);
			std::vector<int> relations;
#pragma region get_relations
			for (int i = 0; i < segments.size(); i++)
			{
				for (int j = i + 1; j < segments.size(); j++)
				{
					if (i != j && !used[i] && !used[j])
					{
						double l_0_0 = GetLength(segments[i][0], segments[j][0]);
						double l_0_1 = GetLength(segments[i][0], segments[j][1]);
						double l_1_0 = GetLength(segments[i][1], segments[j][0]);
						double l_1_1 = GetLength(segments[i][1], segments[j][1]);

						bool b_0_0 = IsAlmostZero_Double(l_0_0, DOUBLE_EPSILON);
						bool b_0_1 = IsAlmostZero_Double(l_0_1, DOUBLE_EPSILON);
						bool b_1_0 = IsAlmostZero_Double(l_1_0, DOUBLE_EPSILON);
						bool b_1_1 = IsAlmostZero_Double(l_1_1, DOUBLE_EPSILON);

						if ((b_0_0 && b_1_1) || (b_0_1 && b_1_0))
						{
							used[j] = true;
							continue;
						}

						if (b_0_0)
						{
							relations.push_back(i);
							relations.push_back(0);
							relations.push_back(j);
							relations.push_back(0);
							continue;
						}
						if (b_0_1)
						{
							relations.push_back(i);
							relations.push_back(0);
							relations.push_back(j);
							relations.push_back(1);
							continue;
						}
						if (b_1_0)
						{
							relations.push_back(i);
							relations.push_back(1);
							relations.push_back(j);
							relations.push_back(0);
							continue;
						}
						if (b_1_1)
						{
							relations.push_back(i);
							relations.push_back(1);
							relations.push_back(j);
							relations.push_back(1);
							continue;
						}
					}
				}
			}
#pragma endregion

			std::vector<std::vector<int>> ones;


			while (true)
			{
				int index = -1;
				int end = -1;

				for (int i = 0; i < segments.size(); i++)
				{
					if (!used[i]) {
						index = i;
						end = 0;
						used[i] = true;
						break;
					}
				}

				if (index < 0)break;

				Vector3d1 line(1, segments[index][end]);

				std::vector<int> one(1, index);

				while (true)
				{
					end = 1 - end;
					bool search = false;
					for (int i = 0; i < relations.size(); i = i + 4)
					{
						if (relations[i] == index && relations[static_cast<int64_t>(i) + 1] == end && !used[relations[static_cast<int64_t>(i) + 2]])
						{
							line.push_back(segments[relations[static_cast<int64_t>(i) + 2]][relations[static_cast<int64_t>(i) + 3]]);
							one.push_back(relations[static_cast<int64_t>(i) + 2]);
							index = relations[static_cast<int64_t>(i) + 2];
							end = relations[static_cast<int64_t>(i) + 3];
							used[index] = true;
							search = true;
							break;
						}
						if (relations[static_cast<int64_t>(i) + 2] == index && relations[static_cast<int64_t>(i) + 3] == end && !used[relations[i]])
						{
							line.push_back(segments[relations[i]][relations[static_cast<int64_t>(i) + 1]]);
							one.push_back(relations[i]);
							index = relations[i];
							end = relations[static_cast<int64_t>(i) + 1];
							used[index] = true;
							search = true;
							break;
						}
					}
					if (!search) { break; }
				}

				ones.push_back(one);
				lines.push_back(line);
			}
		}

		static Vector3d IntersectPointPlane2Ray(const Vector3d& planar_location, Vector3d& planar_direction, 
			const Vector3d& ray_location, Vector3d& ray_vector)
		{
			Vector3d project_point = PlaneProject(planar_location, planar_direction, ray_location);
			double distance = GetDistance(ray_location, project_point);
			if (IsAlmostZero(GetLength(project_point, ray_location)))
				return ray_location;
			double angle = GetAngleBetween(ray_vector, project_point - ray_location);
			double length = distance / cos(angle);
			return ray_location + SetVectorLength(ray_vector, length);
		}
		
		static Vector3d ComputeNormalFromPolyline(const Vector3d1& points)
		{
			Vector3d planar_direction;
			planar_direction = GetCrossproduct(points[0] - points[1], points[2] - points[1]);
			SetVectorLength(planar_direction, 1.0);
			return planar_direction;
		}

		static void  ComputePlanarFromPolyline(Vector3d& planar_location, Vector3d& planar_direction, const Vector3d1& points)
		{
			planar_location = points[0];
			planar_direction = GetCrossproduct(points[0] - points[1], points[2] - points[1]);
			SetVectorLength(planar_direction, 1.0);
		}



		// Compute barycentric coordinates (u, v, w) for
		// point p with respect to triangle (a, b, c)
		static void Barycentric(const Vector3d& p, const Vector3d& a, const Vector3d& b, const Vector3d& c, double& u, double& v, double& w)
		{
			Vector3d v0 = GetMinus(b, a), v1 = GetMinus(c, a), v2 = GetMinus(p, a);

			double d00 = GetDotproduct(v0, v0);
			double d01 = GetDotproduct(v0, v1);
			double d11 = GetDotproduct(v1, v1);
			double d20 = GetDotproduct(v2, v0);
			double d21 = GetDotproduct(v2, v1);
			double denom = d00 * d11 - d01 * d01;
			v = (d11 * d20 - d01 * d21) / denom;
			w = (d00 * d21 - d01 * d20) / denom;
			u = 1.0f - v - w;
		}

		static Vector3d Barycentric(const Vector3d& p, const Vector3d& a, const Vector3d& b, const Vector3d& c)
		{
			Vector3d v0 = GetMinus(b, a), v1 = GetMinus(c, a), v2 = GetMinus(p, a);
			double d00 = GetDotproduct(v0, v0);
			double d01 = GetDotproduct(v0, v1);
			double d11 = GetDotproduct(v1, v1);
			double d20 = GetDotproduct(v2, v0);
			double d21 = GetDotproduct(v2, v1);
			double denom = d00 * d11 - d01 * d01;
			double v = (d11 * d20 - d01 * d21) / denom;
			double w = (d00 * d21 - d01 * d20) / denom;
			double u = 1.0f - v - w;
			return Vector3d(u,v,w);
		}

		//===============================================================
		template <class Type>
		static bool DetectColinear(const Type& v, const Type& s, const Type& e, const double& angle_match_error, const double& dis_match_error)
		{
			if (DetectCoincident(v, s, dis_match_error) || DetectCoincident(v, e, dis_match_error)) return true;
			double angle = GetAngleBetween(v - s, e - s);
			if (IsAlmostZero_Double(angle, angle_match_error))return true;
			if (IsAlmostZero_Double(angle - Math_PI, angle_match_error))return true;
			return false;
		}

		template <class Type>
		static bool DetectVertical(const Type& direction_0, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
		{
			auto angle = GetAngleBetween(direction_0, direction_1);
			return IsAlmostZero_Double(angle - Math_PI / 2.0, angle_match_error);
		};

		template <class Type>
		static bool DetectVertical(const Type& seg_0_s, const Type& seg_0_e, const Type& seg_1_s, const Type& seg_1_e, const double& angle_match_error, const double& dis_match_error)
		{
			return DetectVertical(seg_0_e - seg_0_s, seg_0_e - seg_0_s, angle_match_error, dis_match_error);
		};

		template <class Type>
		static bool DetectVertical(const std::pair<Type, Type>& seg_0, const std::pair<Type, Type>& seg_1, 
			const double& angle_match_error, const double& dis_match_error)
		{
			return DetectVertical(seg_0.first, seg_0.second, seg_1.first, seg_1.second, angle_match_error, dis_match_error);
		};

		template <class Type>
		static bool DetectCoincident(const Type& v0, const Type& v1, const double& EPSILON = DOUBLE_EPSILON)
		{
			return IsAlmostZero_Double(GetDistance(v0, v1), EPSILON);
		}

		static bool DetectCoplanar(const Vector3d& planar_location_0, Vector3d& planar_direction_0,
			const Vector3d& planar_location_1, const Vector3d& planar_direction_1,
			const double& angle_match_error, const double& dis_match_error)
		{
			auto angle = GetAngleBetween(planar_direction_0, planar_direction_1);
			if (IsAlmostZero_Double(angle - Math_PI, angle_match_error) || IsAlmostZero_Double(angle, angle_match_error))
			{
				double dis = GetLength(PlaneProject(planar_location_0, planar_direction_0, planar_location_1), planar_location_1);
				return IsAlmostZero_Double(dis, dis_match_error);
			}
			else
				return false;
		}

		template <class Type>
		static bool DetectParallel(const Type& direction_0, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
		{
			auto angle = GetAngleBetween(direction_0, direction_1);
			return (IsAlmostZero_Double(angle - Math_PI, angle_match_error) || IsAlmostZero_Double(angle, angle_match_error));
		};

		template <class Type>
		static bool DetectParallel(const std::pair<Type, Type>& seg_0, const std::pair<Type, Type>& seg_1, 
			const double& angle_match_error, const double& dis_match_error)
		{
			return DetectParallel(seg_0.second - seg_0.first, seg_1.second - seg_1.first, angle_match_error, dis_match_error);
		};
		template <class Type>
		static bool DetectCoDirection(const Type& direction_0, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
		{
			auto angle = GetAngleBetween(direction_0, direction_1);
			return (IsAlmostZero_Double(angle, angle_match_error));
		};

		template <class Type>
		static bool DetectColinear_Direction(const Type& location_0, const Type& direction_0, const Type& location_1, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
		{
			if (DetectParallel(direction_0, direction_1, angle_match_error, dis_match_error))
			{
				if (IsAlmostZero_Double(GetLength(location_0 - location_1), dis_match_error))
					return true;
				return DetectParallel(direction_0, location_1 - location_0, angle_match_error, dis_match_error);
			}
			else
				return false;
		};

		static bool DetectAlign2D(const Vector2d& s_0, const Vector2d& e_0, const Vector2d& s_1, const Vector2d& e_1, const double& angle_match_error, const double& dis_match_error)
		{
			double d0 = GetDistance(s_0, s_1);
			double d1 = GetDistance(s_0, e_1);
			double d2 = GetDistance(e_0, s_1);
			double d3 = GetDistance(e_0, e_1);
			if (IsAlmostZero_Double(d0, dis_match_error) && IsAlmostZero_Double(d3, dis_match_error))
				return true;
			if (IsAlmostZero_Double(d1, dis_match_error) && IsAlmostZero_Double(d2, dis_match_error))
				return true;
			return false;
		};

		static bool DetectAlign3D(const Vector3d& s_0, const Vector3d& e_0, const Vector3d& s_1, const Vector3d& e_1, const double& angle_match_error, const double& dis_match_error)
		{
			double d0 = GetDistance(s_0, s_1);
			double d1 = GetDistance(s_0, e_1);
			double d2 = GetDistance(e_0, s_1);
			double d3 = GetDistance(e_0, e_1);
			if (IsAlmostZero_Double(d0, dis_match_error) && IsAlmostZero_Double(d3, dis_match_error))
				return true;
			if (IsAlmostZero_Double(d1, dis_match_error) && IsAlmostZero_Double(d2, dis_match_error))
				return true;
			return false;
		};

		//this function has bug ;
		//Do not use it
		template <class Type>
		static bool DetectColinear_Segment(const Type& s_0, const Type& e_0, const Type& s_1, const Type& e_1, const double& angle_match_error, const double& dis_match_error)
		{
			return DetectColinear_Direction(s_0, e_0 - s_0, s_1, e_1 - s_1, angle_match_error, dis_match_error);
		};

		static std::vector<Vector2d> Polygon_Clear(const std::vector<Vector2d>& vecs,
			const double& angle_match_error, const double& dis_match_error)
		{
			//remove duplicate points
			std::vector<Vector2d> vecs_0;
			for (int i = 0; i < vecs.size(); i++)
			{
				if (i != vecs.size() - 1)
				{
					if (vecs_0.empty()) vecs_0.emplace_back(vecs[i]);
					else
					{
						if (dis_match_error > 0)
						{
							if (!IsAlmostZero_Double(GetLength(vecs_0.back(), vecs[i]), dis_match_error))
								vecs_0.emplace_back(vecs[i]);
						}
						else
						{
							if (!IsAlmostZero(GetLength(vecs_0.back(), vecs[i])))
								vecs_0.emplace_back(vecs[i]);
						}
					}
				}
				else
				{
					if (vecs_0.empty())
						vecs_0.emplace_back(vecs[i]);
					else
					{
						if (dis_match_error > 0)
						{
							if (!IsAlmostZero_Double(GetLength(vecs_0.back(), vecs[i]), dis_match_error) &&
								!IsAlmostZero_Double(GetLength(vecs_0.front(), vecs[i]), dis_match_error))
								vecs_0.emplace_back(vecs[i]);
						}
						else
						{
							if (!IsAlmostZero(GetLength(vecs_0.back(), vecs[i])) &&
								!IsAlmostZero(GetLength(vecs_0.front(), vecs[i])))
								vecs_0.emplace_back(vecs[i]);
						}

					}
				}
			}
			//remove collinear points
			std::vector<Vector2d> vecs_1;
			for (int i = 0; i < vecs_0.size(); i++)
			{
				auto pre_v = vecs_0[(i + vecs_0.size() - 1) % vecs_0.size()];
				auto cur_v = vecs_0[i];
				auto next_v = vecs_0[(static_cast<int64_t>(i) + 1) % vecs_0.size()];
				double angle = GetAngleBetween(cur_v - pre_v, next_v - cur_v);
				if (angle_match_error > 0.0)
				{
					if (!IsAlmostZero_Double(angle, angle_match_error))
						vecs_1.emplace_back(cur_v);
				}
				else
				{
					if (!IsAlmostZero(angle))
						vecs_1.emplace_back(cur_v);
				}
			}

			return vecs_1;
		};


		static void GetUnitCube(Vector3d1& cube_vecs,
			Vector1i1& cube_face_id_0, Vector1i1& cube_face_id_1, Vector1i1& cube_face_id_2,
			const double& scale = 1.0)
		{
			cube_vecs.clear();
			cube_face_id_0.clear();
			cube_face_id_1.clear();
			cube_face_id_2.clear();

			cube_vecs.push_back(Vector3d(0.5, 0.5, 0.5));
			cube_vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
			cube_vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
			cube_vecs.push_back(Vector3d(0.5, 0.5, -0.5));

			cube_vecs.push_back(Vector3d(0.5, -0.5, 0.5));
			cube_vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
			cube_vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
			cube_vecs.push_back(Vector3d(0.5, -0.5, -0.5));

			auto sm = Functs::ScaleMatrix(Vector3d(scale, scale, scale));
			cube_vecs = Functs::PosApplyM(cube_vecs, sm);

			int face_index_0[4] = { 0, 1, 2, 3 };
			int face_index_1[4] = { 5, 1, 0, 4 };
			int face_index_2[4] = { 4, 0, 3, 7 };
			int face_index_3[4] = { 5, 4, 7, 6 };
			int face_index_4[4] = { 7, 3, 2, 6 };
			int face_index_5[4] = { 6, 2, 1, 5 };

			Vector1i2 quad_faces;
			quad_faces.push_back(Vector1i1(face_index_0, face_index_0 + 4));
			quad_faces.push_back(Vector1i1(face_index_1, face_index_1 + 4));
			quad_faces.push_back(Vector1i1(face_index_2, face_index_2 + 4));
			quad_faces.push_back(Vector1i1(face_index_3, face_index_3 + 4));
			quad_faces.push_back(Vector1i1(face_index_4, face_index_4 + 4));
			quad_faces.push_back(Vector1i1(face_index_5, face_index_5 + 4));

			for (auto qf : quad_faces)
			{
				cube_face_id_0.push_back(qf[2]);
				cube_face_id_1.push_back(qf[1]);
				cube_face_id_2.push_back(qf[0]);
				cube_face_id_0.push_back(qf[0]);
				cube_face_id_1.push_back(qf[3]);
				cube_face_id_2.push_back(qf[2]);
			}
		};


		static Vector3d1 EmumerateRotations()
		{
			Vector3d1 rotations;
			rotations.emplace_back(Vector3d(90.0, 90.0, 180.0));
			rotations.emplace_back(Vector3d(-90.0, 0.0, 90.0));
			rotations.emplace_back(Vector3d(0.0, 180.0, 0.0));

			rotations.emplace_back(Vector3d(90.0, 0.0, 180.0));
			rotations.emplace_back(Vector3d(0.0, 0.0, 90.0));
			rotations.emplace_back(Vector3d(0.0, -90.0, 0.0));

			rotations.emplace_back(Vector3d(0.0, 0.0, 0.0));
			rotations.emplace_back(Vector3d(90.0, 0.0, 90.0));
			rotations.emplace_back(Vector3d(90.0, -90.0, 180.0));

			rotations.emplace_back(Vector3d(0.0, 90.0, 0.0));
			rotations.emplace_back(Vector3d(180.0, 0.0, 90.0));
			rotations.emplace_back(Vector3d(90.0, -180.0, 180.0));

			rotations.emplace_back(Vector3d(0.0, 180.0, 90.0));
			rotations.emplace_back(Vector3d(-90.0, 90.0, 90.0));
			rotations.emplace_back(Vector3d(90.0, 180.0, 0.0));

			rotations.emplace_back(Vector3d(0.0, 0.0, 180.0));
			rotations.emplace_back(Vector3d(0.0, -90.0, 90.0));
			rotations.emplace_back(Vector3d(90.0, 0.0, -90.0));

			rotations.emplace_back(Vector3d(-90.0, 0.0, -90.0));
			rotations.emplace_back(Vector3d(180.0, 90.0, 90.0));
			rotations.emplace_back(Vector3d(0.0, -180.0, 180.0));

			rotations.emplace_back(Vector3d(90.0, 0.0, 0.0));
			rotations.emplace_back(Vector3d(90.0, -90.0, 90.0));
			rotations.emplace_back(Vector3d(0.0, 0.0, 270.0));

			for (auto& rotation : rotations)
			{
				rotation[0] = rotation[0] / 180.0 * Math_PI;
				rotation[1] = rotation[1] / 180.0 * Math_PI;
				rotation[2] = rotation[2] / 180.0 * Math_PI;
			}
			return rotations;
		};


#pragma endregion

#pragma region BasicMathFunctions

		static double RandomD(const double& min_d = 0.0, const double& max_d = 1.0)
		{
			std::uniform_real_distribution<> dis(min_d, max_d);
			return dis(MATHGEN);
		}

		//select a number among 0,1,2,...,s-1
		//s should be positive int
		static int RandomI(const int& s)
		{
			if (s == 1) return 0;
			double x = 1.0 / s;
			double d = RandomD();

			for (int i = 0; i < s; i++)
			{
				double min_d = i * x;
				double max_d = (static_cast<int64_t>(i) + 1) * x;
				if (d >= min_d && d <= max_d)
				{
					return i;
				}
			}

			return -1;
		}

		static double RandomDD(const double& min_d = 0.0, const double& max_d = 1.0)
		{
			std::uniform_real_distribution<> dis(min_d, max_d);
			return dis(MATHGEN);
		}

		//select a number among 0,1,2,...,s-1
		//s should be positive int
		static int RandomII(int s)
		{
			if (s == 1) return 0;
			double x = 1.0 / s;
			double d = RandomDD();

			for (int i = 0; i < s; i++)
			{
				double min_d = i * x;
				double max_d = (static_cast<int64_t>(i) + 1) * x;
				if (d >= min_d && d <= max_d)
				{
					return i;
				}
			}

			return -1;
		}

		static Vector3d GetCrossproduct(const Vector3d& v1, const Vector3d& v2) {
			return glm::cross(v1, v2);
		}

		template <class Type>
		static double GetDotproduct(const Type& v1, const Type& v2) {
			return glm::dot(v1, v2);
		}

		template <class Type>
		static Type GetMinus(const Type& a, const Type& b)
		{
			Type c=a;
			for (int i = 0; i < a.length(); i++)
				c[i] = a[i] - b[i];
			return c;
		}

		static bool AreAlmostEqual(const double& value1, const double& value2) {
			if (value1 == value2) {
				return true;
			}
			double eps = (glm::abs(value1) + glm::abs(value2) + 10.0) * DOUBLE_EPSILON;
			double delta = value1 - value2;
			return (-eps < delta) && (eps > delta);
		}

		static bool AreAlmostEqual_Double(const double& value1, const double& value2, const double& EPSILON) {
			return IsAlmostZero_Double(value1 - value2, EPSILON);
		}

		static bool IsAlmostZero(const double& value) {
			return (value < DOUBLE_EPSILON) && (value > -DOUBLE_EPSILON);
		}
		static bool IsAlmostZero_Double(const double& value, const double& EPSILON) {
			return (value < EPSILON) && (value > -EPSILON);
		}
		
		/// Returns true if two given floating point numbers are epsilon-equal.
		/// Method automatically adjust the epsilon to the absolute size of given numbers.
		static bool AreAlmostEqual(const float& value1, const float& value2) {
			// in case they are Infinities (then epsilon check does not work)
			if (value1 == value2) {
				return true;
			}
			// computes (|value1-value2| / (|value1| + |value2| + 10.0)) < SINGLE_EPSILON
			float eps = (float)((glm::abs(value1) + glm::abs(value2) + 10.0) * SINGLE_EPSILON);
			float delta = value1 - value2;
			return (-eps < delta) && (eps > delta);
		}

		static void ZeroVector(Vector3d& v)
		{
			if (IsAlmostZero(v[0]))v[0] = 0.0;
			if (IsAlmostZero(v[1]))v[1] = 0.0;
			if (IsAlmostZero(v[2]))v[2] = 0.0;
		}

		template<class Type>
		static double GetMax(const Type& vec)
		{
			double maxd = vec[0];
			for (int i = 0; i < vec.length(); i++)
				maxd = max(maxd,vec[i]);
			return maxd;
		}


		template<class Type>
		static bool VectorInsertNoDuplicate(std::vector<Type>& vecs, const Type& element)
		{
			if (CheckContain(vecs, element))
				return false;

			vecs.emplace_back(element);
			return true;
		}

		template<class Type>
		static void VectorInsertNoDuplicate(std::vector<Type>& vecs, const std::vector<Type>& elements)
		{
			for (auto& element : elements)
				VectorInsertNoDuplicate(vecs, element);
		}

		template<class Type>
		static std::vector<Type> VectorMerge(const std::vector<Type>& vecs_0, const std::vector<Type>& vecs_1)
		{
			std::vector<Type> result = vecs_0;
			result.insert(result.end(), vecs_1.begin(), vecs_1.end());
			return result;
		}

		template<class Type>
		static bool CheckContain(const std::vector <Type>& vecs, const Type& element)
		{
			return std::find(vecs.begin(), vecs.end(), element) != vecs.end();
		}

		template <class Type>
		static int VectorIndex(const std::vector <Type>& vecs, const Type& element)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i] == element)
				{
					return i;
				}
			}
			return -1;
		}

		template <class Type>
		static int VectorIndex(const std::vector <std::vector <Type>>& vecs, const Type& element)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (VectorIndex(vecs[i], element) >= 0)
					return i;
			}
			return -1;
		}

		template <class Type>
		static std::vector <Type> VectorAdd(const std::vector <Type>& vecs, const Type& element)
		{
			std::vector <Type> result = vecs;
			for (auto& r : result)
				r += element;
			return result;
		}

		static bool VectorContain(const std::vector<int>& vecs, const int& element)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i] == element)
					return true;
			}

			return false;
		}

		static int VectorContainReturnIndex(const std::vector<int>& vecs, const int& element)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i] == element)
					return i;
			}

			return -1;
		}

		static bool VectorContainForSpecialCase(const std::vector<std::vector<int>>& vecs, const std::vector<int>& element)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i][0] == element[0] && vecs[i][1] == element[1])
					return true;
			}

			return false;
		}
		static bool VectorContainForSpecialCase1(const std::vector<std::vector<int>>& vecs, const std::vector<int>& element)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i][0] == element[0] && vecs[i][1] == element[1]) return true;
				if (vecs[i][0] == element[1] && vecs[i][1] == element[0]) return true;
			}

			return false;
		}

		static int VectorContainForSpecialCase2(const std::vector<std::vector<int>>& vecs, const int& element_0, const int& element_1)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i][0] == element_0 && vecs[i][1] == element_1) return i;
				if (vecs[i][0] == element_1 && vecs[i][1] == element_0) return i;
			}
			return -1;
		}

		static int VectorContainForSpecialCase3(const std::vector<std::vector<int>>& vecs, const int& element_0, const int& element_1)
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				if (vecs[i][0] == element_0 && vecs[i][1] == element_1) return i;
			}
			return -1;
		}
#pragma endregion


#pragma region Transformation

		static Vector2d Vector3d2d(const Vector3d& v)
		{
			return Vector2d(v[0], v[1]);
		}

		static Vector2d1 Vector3d2d(const Vector3d1& vecs_3d)
		{
			Vector2d1 vecs_2d;
			for (auto& v : vecs_3d)
				vecs_2d.emplace_back(Vector3d2d(v));
			return vecs_2d;
		}

		static Vector2d2 Vector3d2d(const Vector3d2& vecs_3d)
		{
			Vector2d2 vecs_2d;
			for (auto v : vecs_3d)
				vecs_2d.emplace_back(Vector3d2d(v));
			return vecs_2d;
		}

		static Vector3d Vector2d3d(const Vector2d& v, const double& z = 0.0)
		{
			return Vector3d(v[0], v[1], z);
		}

		static Vector3d1 Vector2d3d(const Vector2d1& vecs_2d, double z = 0.0)
		{
			Vector3d1 vecs_3d;
			for (auto v : vecs_2d)
				vecs_3d.emplace_back(Vector2d3d(v, z));
			return vecs_3d;
		}

		static Vector3d2 Vector2d3d(const Vector2d2& vecs_2d, double z = 0.0)
		{
			Vector3d2 vecs_3d;
			for (auto v : vecs_2d)
				vecs_3d.emplace_back(Vector2d3d(v, z));
			return vecs_3d;
		}

		static Vector3d VecApplyM(const Vector3d& v, const  glm::dmat4& M)
		{
			return Vector3d(M * glm::vec4(v, 1.0)) - Vector3d(M * glm::vec4(Vector3d(0.0, 0.0, 0.0), 1.0));
		}

		static Vector3d1 VecApplyM(const Vector3d1& vecs, const glm::dmat4& M)
		{
			Vector3d1 ps;
			for (auto &p : vecs)
				ps.emplace_back(VecApplyM(p, M));
			return ps;
		}

		static Vector3d2 VecApplyM(const Vector3d2& veces, const glm::dmat4& M)
		{
			Vector3d2 pses;
			for (auto vecs : veces)
				pses.emplace_back(VecApplyM(vecs, M));
			return pses;
		}

		static Vector3d PosApplyM(const Vector3d& v, const glm::dmat4& M)
		{
			return Vector3d(M * glm::vec4(v, 1.0));
		}

		static Vector3d1 PosApplyM(const Vector3d1& vecs, const glm::dmat4& M)
		{
			Vector3d1 ps;
			for (auto &p : vecs)
				ps.emplace_back(PosApplyM(p, M));
			return ps;
		}

		static std::pair<Vector3d, Vector3d> PosApplyM(const std::pair<Vector3d, Vector3d>& vecs, const glm::dmat4& M)
		{
			return std::pair<Vector3d, Vector3d>(PosApplyM(vecs.first, M), PosApplyM(vecs.second, M));
		}

		static Vector3d2 PosApplyM(const Vector3d2& veces, const glm::dmat4& M)
		{
			Vector3d2 pses;
			for (auto vecs : veces)
				pses.emplace_back(PosApplyM(vecs, M));
			return pses;
		}

		static Vector3d3 PosApplyM(const Vector3d3& veces, const glm::dmat4& M)
		{
			Vector3d3 pses;
			for (auto vecs : veces)
				pses.emplace_back(PosApplyM(vecs, M));
			return pses;
		}
		
		static glm::dmat4 RotationMatrixXYZ(const Vector3d& xx, const Vector3d &yy, const Vector3d &zz)
		{
			auto x = xx;
			auto y = yy;
			auto z = zz;
		
			ZeroVector(x);
			ZeroVector(y);
			ZeroVector(z);
			x = x / (double)GetLength(x);
			y = y / (double)GetLength(y);
			z = z / (double)GetLength(z);
		
			glm::dmat4  rotationMatrix;
		
			rotationMatrix[0][0] = x[0];
			rotationMatrix[0][1] = y[0];
			rotationMatrix[0][2] = z[0];
			rotationMatrix[0][3] = 0.0;
		
			rotationMatrix[1][0] = x[1];
			rotationMatrix[1][1] = y[1];
			rotationMatrix[1][2] = z[1];
			rotationMatrix[1][3] = 0.0;
		
			rotationMatrix[2][0] = x[2];
			rotationMatrix[2][1] = y[2];
			rotationMatrix[2][2] = z[2];
			rotationMatrix[2][3] = 0.0;
		
			rotationMatrix[3][0] = 0.0;
			rotationMatrix[3][1] = 0.0;
			rotationMatrix[3][2] = 0.0;
			rotationMatrix[3][3] = 1.0;
		
			return rotationMatrix;
		
		}


		static glm::dmat4 RotationMatrix(const Vector3d& o, const Vector3d &t, const Vector3d &n)
		{
			double angle = GetAngleBetween(o, t);
	
			if (IsAlmostZero(angle))
			{
				glm::dmat4  rotationMatrix;
	
				rotationMatrix[0][0] = 1.0;
				rotationMatrix[0][1] = 0.0;
				rotationMatrix[0][2] = 0.0;
				rotationMatrix[0][3] = 0.0;
				rotationMatrix[1][0] = 0.0;
				rotationMatrix[1][1] = 1.0;
				rotationMatrix[1][2] = 0.0;
				rotationMatrix[1][3] = 0.0;
				rotationMatrix[2][0] = 0.0;
				rotationMatrix[2][1] = 0.0;
				rotationMatrix[2][2] = 1.0;
				rotationMatrix[2][3] = 0.0;
				rotationMatrix[3][0] = 0.0;
				rotationMatrix[3][1] = 0.0;
				rotationMatrix[3][2] = 0.0;
				rotationMatrix[3][3] = 1.0;
	
				return rotationMatrix;
			}
			else
			{
				return RotationMatrix(n, angle);
			}
		}
	
		static void AAA(const Vector3d& v,Vector3d &n)
		{
			auto a = v[0];
			auto b = v[1];
			auto c = v[2];
			bool bx = IsAlmostZero(a);
			bool by = IsAlmostZero(b);
			bool bz = IsAlmostZero(c);
		
			if (bx&&by&&bz)
			{
				std::cerr << "if (bx&&by&&bz)" << std::endl;
				system("pause");
			}
		
			if (bx&&by&&!bz)
			{
				n[0] = 1.0;
				n[1] = 1.0;
				n[2] = 0.0;
			}
			if (bx&&!by&&bz)
			{
				n[0] = 1.0;
				n[1] = 0.0;
				n[2] = 1.0;
			}
		
			if (bx&&!by&&!bz)
			{
				n[0] = 1.0;
				n[1] = 1.0;
				n[2] = -b/c;
			}
		
			if (!bx&&by&&bz)
			{
				n[0] = 0.0;
				n[1] = 1.0;
				n[2] = 1.0;
			}
		
			if (!bx&&by&&!bz)
			{
				n[0] = 1.0;
				n[1] = 1.0;
				n[2] = -a / c;
			}
			if (!bx&&!by&&bz)
			{
				n[0] = 1.0;
				n[1] = -a/b;
				n[2] = 1.0;
			}
		
			if (!bx&&!by&&!bz)
			{
				n[0] = 1.0;
				n[1] = 1.0;
				n[2] = -(a+b) / c;
			}
		
		}

		static glm::dmat4 RotationMatrix(const Vector3d& o, const Vector3d &t)
		{
			Vector3d n = GetCrossproduct(o, t);
			double angle = GetAngleBetween(o, t);
	
			if (IsAlmostZero(angle - Math_PI))
			{
				AAA(o,n);
			}
	
			if (IsAlmostZero(angle))
			{
				glm::dmat4  rotationMatrix;
	
				rotationMatrix[0][0] = 1.0;
				rotationMatrix[0][1] = 0.0;
				rotationMatrix[0][2] = 0.0;
				rotationMatrix[0][3] = 0.0;
				rotationMatrix[1][0] = 0.0;
				rotationMatrix[1][1] = 1.0;
				rotationMatrix[1][2] = 0.0;
				rotationMatrix[1][3] = 0.0;
				rotationMatrix[2][0] = 0.0;
				rotationMatrix[2][1] = 0.0;
				rotationMatrix[2][2] = 1.0;
				rotationMatrix[2][3] = 0.0;
				rotationMatrix[3][0] = 0.0;
				rotationMatrix[3][1] = 0.0;
				rotationMatrix[3][2] = 0.0;
				rotationMatrix[3][3] = 1.0;
	
				return rotationMatrix;
			}
			else
			{
				return RotationMatrix(n, angle);
			}
		}



		static glm::dmat4 RotationMatrix(const Vector3d& n, const double& angle)
		{
			//return glm::rotate(angle, n);
			double u = n[0];
			double v = n[1];
			double w = n[2];

			glm::dmat4  rotationMatrix;

			double L = (u * u + v * v + w * w);

			//angle = angle * M_PI / 180.0; //converting to radian value	
			double u2 = u * u;
			double v2 = v * v;
			double w2 = w * w;

			rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
			rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
			rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
			rotationMatrix[0][3] = 0.0;

			rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
			rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
			rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
			rotationMatrix[1][3] = 0.0;

			rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
			rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
			rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
			rotationMatrix[2][3] = 0.0;

			rotationMatrix[3][0] = 0.0;
			rotationMatrix[3][1] = 0.0;
			rotationMatrix[3][2] = 0.0;
			rotationMatrix[3][3] = 1.0;

			return rotationMatrix;
		}

		static Vector3d RotationAxis(const Vector3d& p, const double& angle, const Vector3d& n)
		{
			// auto rtv = glm::rotate(angle, n)* glm::dvec4(p, 1.0);
			// return Vector3d(rtv[0], rtv[1], rtv[2]);
			glm::dmat4 inputMatrix(0.0);
			inputMatrix[0][0] = p[0];
			inputMatrix[1][0] = p[1];
			inputMatrix[2][0] = p[2];
			inputMatrix[3][0] = 1.0;
			double u = n[0];
			double v = n[1];
			double w = n[2];

			glm::dmat4  rotationMatrix;

			double L = (u * u + v * v + w * w);

			//angle = angle * M_PI / 180.0; //converting to radian value
			double u2 = u * u;
			double v2 = v * v;
			double w2 = w * w;

			rotationMatrix[0][0] = (u2 + (v2 + w2) * glm::cos(angle)) / L;
			rotationMatrix[0][1] = (u * v * (1 - glm::cos(angle)) - w * glm::sqrt(L) * glm::sin(angle)) / L;
			rotationMatrix[0][2] = (u * w * (1 - glm::cos(angle)) + v * glm::sqrt(L) * glm::sin(angle)) / L;
			rotationMatrix[0][3] = 0.0;

			rotationMatrix[1][0] = (u * v * (1 - glm::cos(angle)) + w * glm::sqrt(L) * glm::sin(angle)) / L;
			rotationMatrix[1][1] = (v2 + (u2 + w2) * glm::cos(angle)) / L;
			rotationMatrix[1][2] = (v * w * (1 - glm::cos(angle)) - u * glm::sqrt(L) * glm::sin(angle)) / L;
			rotationMatrix[1][3] = 0.0;

			rotationMatrix[2][0] = (u * w * (1 - glm::cos(angle)) - v * glm::sqrt(L) * glm::sin(angle)) / L;
			rotationMatrix[2][1] = (v * w * (1 - glm::cos(angle)) + u * glm::sqrt(L) * glm::sin(angle)) / L;
			rotationMatrix[2][2] = (w2 + (u2 + v2) * glm::cos(angle)) / L;
			rotationMatrix[2][3] = 0.0;

			rotationMatrix[3][0] = 0.0;
			rotationMatrix[3][1] = 0.0;
			rotationMatrix[3][2] = 0.0;
			rotationMatrix[3][3] = 1.0;

			double outputMatrix[4][1] = { 0.0, 0.0, 0.0, 0.0 };

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 1; j++) {
					outputMatrix[i][j] = 0;
					for (int k = 0; k < 4; k++) {
						outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
					}
				}
			}
			return Vector3d(outputMatrix[0][0], outputMatrix[0][1], outputMatrix[0][2]);
		}

		static glm::dmat4 TranslationMatrix(const Vector3d& v)
		{
			glm::dmat4  translationMatrix;
			translationMatrix[0][0] = 1.0;
			translationMatrix[0][1] = 0.0;
			translationMatrix[0][2] = 0.0;
			translationMatrix[0][3] = 0.0;

			translationMatrix[1][0] = 0.0;
			translationMatrix[1][1] = 1.0;
			translationMatrix[1][2] = 0.0;
			translationMatrix[1][3] = 0.0;

			translationMatrix[2][0] = 0.0;
			translationMatrix[2][1] = 0.0;
			translationMatrix[2][2] = 1.0;
			translationMatrix[2][3] = 0.0;

			translationMatrix[3][0] = v[0];
			translationMatrix[3][1] = v[1];
			translationMatrix[3][2] = v[2];
			translationMatrix[3][3] = 1.0;

			return translationMatrix;
		}

		static glm::dmat4 ScaleMatrix(const Vector3d& v)
		{
			glm::dmat4  translationMatrix;
			translationMatrix[0][0] = v[0];
			translationMatrix[0][1] = 0.0;
			translationMatrix[0][2] = 0.0;
			translationMatrix[0][3] = 0.0;

			translationMatrix[1][0] = 0.0;
			translationMatrix[1][1] = v[1];
			translationMatrix[1][2] = 0.0;
			translationMatrix[1][3] = 0.0;

			translationMatrix[2][0] = 0.0;
			translationMatrix[2][1] = 0.0;
			translationMatrix[2][2] = v[2];
			translationMatrix[2][3] = 0.0;

			translationMatrix[3][0] = 0;
			translationMatrix[3][1] = 0;
			translationMatrix[3][2] = 0;
			translationMatrix[3][3] = 1.0;

			return translationMatrix;
		}

		static Vector2d RotationAxis2d(const Vector2d& p, const double& angle, const Vector2d& center)
		{
			Vector3d r = RotationAxis(Vector3d(p[0] - center[0], 0.0, p[1] - center[1]),
				angle, Vector3d(0.0, 1.0, 0.0)) + Vector3d(center[0], 0.0, center[1]);
			return Vector2d(r[0], r[2]);
		}

		static Vector3d RotationAxis(const Vector3d& p, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
		{
			return RotationAxis(p - ray_point, angle, ray_vector) + ray_point;
		}

		static void RotationAxis(Vector3d1& points, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
		{
			for (int i = 0; i < points.size(); i++)
				points[i] = RotationAxis(points[i], angle, ray_point, ray_vector);
		}

		static void RotationAxis(Vector3d2& points, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
		{
			for (int i = 0; i < points.size(); i++)
				RotationAxis(points[i], angle, ray_point, ray_vector);
		}
		static void RotationAxis(Vector3d3& pointses, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
		{
			for (int i = 0; i < pointses.size(); i++)
				RotationAxis(pointses[i], angle, ray_point, ray_vector);
		}

		static Vector3d Translate(const Vector3d& p, const Vector3d& v)
		{
			return p + v;
		}

		static void Translate(Vector3d1& points, const Vector3d& v)
		{
			for (int i = 0; i < points.size(); i++)
				points[i] = Translate(points[i], v);
		}

		static void Translate(Vector3d2& points, const Vector3d& v)
		{
			for (int i = 0; i < points.size(); i++)
				Translate(points[i], v);
		}
#pragma endregion


#pragma region IOFunctions

		static bool LoadExisting(const std::string& path)
		{	
			std::ifstream file(path, std::ios::in);
			if (!file) return false;
			return true;
		}
	
		static bool LoadVectors(const std::string& path, Vector3d3 &vec_3)
		{
			//zigzag_final_path
			int nb_0,nb_1,nb_2;
			std::ifstream file(path, std::ios::in);
		
			if (!file) return false;
			
			file >> nb_0;
			for (int i = 0; i < nb_0; i++)
			{
				file >> nb_1;
				Vector3d2 vec_2;
				for (int j = 0; j < nb_1; j++)
				{
					file >> nb_2;
					Vector3d1 vec_1(nb_2,Vector3d(0.0,0.0,0.0));
					for (int k = 0; k < nb_2; k++)
						file >> vec_1[k][0] >> vec_1[k][1] >> vec_1[k][2];
					vec_2.emplace_back(vec_1);
				}
				vec_3.emplace_back(vec_2);
			}
			file.clear();
			file.close();
		
			return true;
		}
		
		
		static bool LoadVectors(const std::string& path, Vector3d1 &vec_3)
		{
			//zigzag_final_path
			std::ifstream file(path, std::ios::in);
		
			if (!file) return false;
		
			int nb;
			file >> nb;
			for (int i = 0; i < nb; i++)
			{
				vec_3.emplace_back(Vector3d());
				file >> vec_3.back()[0] >> vec_3.back()[1] >> vec_3.back()[2];
			}
			file.clear();
			file.close();
		
			return true;
		}

		static void OutputVectors(const std::string& out_path, const Vector3d3 &vecs)
		{
			std::ofstream file(out_path);
			file << vecs.size() << std::endl;
		
			for (int i = 0; i < vecs.size(); i++){
				file << vecs[i].size() << std::endl;
				for (int j = 0; j < vecs[i].size(); j++){
					file << vecs[i][j].size() << std::endl;
					for (int k = 0; k < vecs[i][j].size(); k++)
						file<< vecs[i][j][k][0] << " " << vecs[i][j][k][1] << " " << vecs[i][j][k][2] << " ";
					file << "" << std::endl;
				}
			}
		
			file.clear();
			file.close();
		}
		
		static void OutputVectors(const std::string& out_path, const Vector3d1 &vecs)
		{
			std::ofstream file(out_path);
			file << vecs.size() << std::endl;
			for (int i = 0; i < vecs.size(); i++)
				file << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] <<std::endl;
			file.clear();
			file.close();
		}


		static HMODULE LoadHMODULE(const string& dll_path)
		{
			if(!DetectExisting(dll_path))
				MAssert("The dll does not exist: "+dll_path);

			HMODULE hModule = LoadLibrary(_T(dll_path.c_str()));
			if (!hModule)
			{
				DWORD dw = GetLastError(); // returns 0xc1 (193)
				MAssert("LoadLibrary failed with error code " + std::to_string(dw));
			}
			else
				std::cerr << "LoadLibrary success\n";

			return hModule;
		};

		static void LoadObj3d(const char* path, std::vector<double>&coords, std::vector<int>&tris)
		{
			auto get_first_integer = [](const char* v)
			{
				int ival;
				std::string s(v);
				std::replace(s.begin(), s.end(), '/', ' ');
				sscanf(s.c_str(), "%d", &ival);
				return ival;
			};

			double x, y, z;
			char line[1024], v0[1024], v1[1024], v2[1024];

			// open the file, return if open fails
			FILE* fp = fopen(path, "r");
			if (!Functs::DetectExisting(path))
			{
				Functs::MAssert("This file does not exist: " + std::string(path));
				return;
			};

			int i = 0;
			while (fgets(line, 1024, fp)) 
			{
				if (line[0] == 'v') 
				{
					sscanf(line, "%*s%lf%lf%lf", &x, &y, &z);
					coords.push_back(x);
					coords.push_back(y);
					coords.push_back(z);
				}
				else
				{
					if (line[0] == 'f')
					{
						sscanf(line, "%*s%s%s%s", v0, v1, v2);
						tris.push_back(get_first_integer(v0) - 1);
						tris.push_back(get_first_integer(v1) - 1);
						tris.push_back(get_first_integer(v2) - 1);
					}
				}
			}
			fclose(fp);
		};

		static void LoadObj3d(const char* path_, Vector3d1 & vecs, Vector1i1 & face_id_0, Vector1i1 & face_id_1, Vector1i1 & face_id_2)
		{
			std::string path = path_;
			if (path.substr(path.size() - 3, path.size()) == "obj")
			{
				std::vector<double> coords;
				Vector1i1 tris;

				LoadObj3d(path.c_str(), coords, tris);

				if (coords.size() == 0)
				{
					return;
				}

				for (int i = 0; i < (int)coords.size(); i += 3)
				{
					vecs.push_back(Vector3d(coords[i + 0], coords[i + 1], coords[i + 2]));
				}

				for (int i = 0; i < (int)tris.size(); i += 3)
				{
					face_id_0.push_back(tris[i + 0]);
					face_id_1.push_back(tris[i + 1]);
					face_id_2.push_back(tris[i + 2]);
				}
				/*********************************************************************************/
			}
		};


		static void OutputRectangle2d(const std::string& path, const std::vector<Vector2d>& points)
		{
			std::ofstream file(path);

			for (int i = 0; i < points.size(); i++)
			{
				file << "v " << points[i][0] << " " << points[i][1] << " " << 0.0 << std::endl;
			}

			int nb = 1;

			file << "f ";
			for (int i = 0; i < points.size(); i++)
			{
				file << IntString(nb) << " ";
				nb++;
			}
			file << "" << std::endl;

			file.clear();
			file.close();
		}

		static void OutputObj3d(const std::string& path, const Vector3d1& points)
		{
			std::ofstream file(path);

			for (auto& p : points)
			{
				file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			}

			int nb = 1;
			file << "f ";
			for (auto& p : points)
			{
				file << IntString(nb) << " ";
				nb++;
			}
			file << "" << std::endl;

			file.clear();
			file.close();
		};

		static void OutputObj3d(const std::string& path, const Vector3d2& points, const int output_index = 1, const string str = "")
		{
			std::ofstream file(path);

			for (auto& points_ : points)
			{
				for (auto& p : points_)
				{
					file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
				}
			}

			if (output_index == 2) file << "g " + str + "_p_0" << std::endl;

			int nb = 1;
			for (int i = 0; i < points.size(); i++)
			{
				if (output_index == 1)
					file << "g " + str + "_f_" << i << std::endl;
				auto points_ = points[i];
				file << "f ";
				for (auto& p : points_)
				{
					file << IntString(nb) << " ";
					nb++;
				}
				file << "" << std::endl;
			}

			file.clear();
			file.close();
		};

		static void OutputObj3d(const std::string& path, const Vector3d3& points, const int output_index = 1)
		{
			std::ofstream file(path);

			for (auto& points_ : points)
			{
				for (auto& points__ : points_)
				{
					for (auto& p : points__)
					{
						file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
					}
				}
			}

			int nb = 1;
			if (output_index == 3) file << "g ps_0" << std::endl;
			for (int i = 0; i < points.size(); i++)
			{
				if (output_index == 2) file << "g p_" << i << std::endl;

				for (int j = 0; j < points[i].size(); j++)
				{
					if (output_index == 1) file << "g f_" << i << "_" << j << std::endl;

					auto points_ = points[i][j];
					file << "f ";
					for (auto& p : points_)
					{
						file << IntString(nb) << " ";
						nb++;
					}
					file << "" << std::endl;
				}
			}
			file.clear();
			file.close();
		};

		static void OutputObj3d(const std::string& path, const Vector3d3& points, const Vector3d1& colors, const int& output_index = 1)
		{
			std::ofstream file(path);

			int index = 0;
			for (int i = 0; i < points.size(); i++)
			{
				auto color = colors[i];
				for (int j = 0; j < points[i].size(); j++)
				{
					for (int k = 0; k < points[i][j].size(); k++)
					{
						file << "v " << points[i][j][k][0] << " " << points[i][j][k][1] << " " << points[i][j][k][2] << " " << color[0] << " " << color[1] << " " << color[2] << std::endl;
					}
				}
			}

			int nb = 1;
			if (output_index == 3) file << "g ps_0" << std::endl;
			for (int i = 0; i < points.size(); i++)
			{
				if (output_index == 2) file << "g p_" << i << std::endl;

				for (int j = 0; j < points[i].size(); j++)
				{
					if (output_index == 1) file << "g f_" << i << "_" << j << std::endl;

					auto points_ = points[i][j];
					file << "f ";
					for (auto& p : points_)
					{
						file << IntString(nb) << " ";
						nb++;
					}
					file << "" << std::endl;
				}
			}
			file.clear();
			file.close();
		};

		static void OutputObj3d(const std::string& path, const Vector3d3& points, const Vector3d2& colors, int output_index = 1)
		{
			std::ofstream file(path);

			int index = 0;
			for (int i = 0; i < points.size(); i++)
			{
				for (int j = 0; j < points[i].size(); j++)
				{
					auto color = colors[i][j];
					for (int k = 0; k < points[i][j].size(); k++)
					{
						file << "v " << points[i][j][k][0] << " " << points[i][j][k][1] << " " << points[i][j][k][2] << " " << color[0] << " " << color[1] << " " << color[2] << std::endl;
					}
				}
			}

			int nb = 1;
			if (output_index == 3) file << "g ps_0" << std::endl;
			for (int i = 0; i < points.size(); i++)
			{
				if (output_index == 2) file << "g p_" << i << std::endl;

				for (int j = 0; j < points[i].size(); j++)
				{
					if (output_index == 1) file << "g f_" << i << "_" << j << std::endl;

					auto points_ = points[i][j];
					file << "f ";
					for (auto& p : points_)
					{
						file << IntString(nb) << " ";
						nb++;
					}
					file << "" << std::endl;
				}
			}
			file.clear();
			file.close();
		};

		static void Output_tree(const int& nodes_nb, const std::vector<int>& edges,
			const std::string& path, const std::vector<string> labels = std::vector<string>())
		{
			std::ofstream file(path);

			file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
			file << "graph" << std::endl;
			file << "[" << std::endl;
			file << "  directed 0" << std::endl;

			for (int i = 0; i < nodes_nb; i++)
			{
				file << "node" << std::endl;
				file << "[" << std::endl;
				file << "id " << i << std::endl;

				if (labels.size() == nodes_nb)
					file << "label " << labels[i] << std::endl;
				else
					file << "label " << i << std::endl;

				file << "]" << std::endl;
			}

			for (int i = 0; i < edges.size(); i = i + 2)
			{
				file << "edge" << std::endl;
				file << "[" << std::endl;

				file << "source " << edges[i] << std::endl;
				file << "target " << edges[static_cast<int64_t>(i) + 1] << std::endl;

				file << "]" << std::endl;
			}

			file << "]" << std::endl;

			file.clear();
			file.close();
		}

		static bool DetectExisting(const std::string& path)
		{
			return (_access(path.c_str(), 0) != -1);
		}

		static void ClearFolder(const std::string& path)
		{
			if (_access(path.c_str(), 0) == -1)
			{
				if (_mkdir(path.c_str())) {};
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

		static std::string EXP(const std::string& py_path)
		{
			std::string path = std::string(_pgmptr).substr(0, std::string(_pgmptr).find_last_of('\\')) + py_path;
			//std::string path = std::string(_pgmptr).substr(0, std::string(_pgmptr).find_last_of('\\')) + py_path;
			if (!Functs::DetectExisting(path))
				Functs::MAssert("std::string EXP(const std::string py_path="")");
			return path;
		}

		static VectorStr1 GetFilesInDirectory(const std::string& path)
		{
			vector<string> names;
			string search_path = path + "/*.*";
			WIN32_FIND_DATA fd;
			HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
			if (hFind != INVALID_HANDLE_VALUE) {
				do {
					// read all (real) files in current folder
					// , delete '!' read other 2 default folder . and ..
					if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
						names.push_back(fd.cFileName);
					}
				} while (::FindNextFile(hFind, &fd));
				::FindClose(hFind);
			}
			return names;
		}





#pragma endregion

#pragma region Graph

		static std::vector<int> MinimalSpanningTreeGeneral(
			const std::vector<int>& edges,
			const std::vector<double>& costs)
		{
			auto unique_nodes = UniqueSet(edges);

			int node_nb = static_cast<int>(unique_nodes.size());

			std::map<int, int> unique_map_0;
			std::map<int, int> unique_map_1;
			for (int i = 0; i < unique_nodes.size(); i++)
			{
				unique_map_0.insert(std::pair<int, int>(unique_nodes[i], i));
				unique_map_1.insert(std::pair<int, int>(i, unique_nodes[i]));
			}

			std::vector<int> map_edges;
			for (auto edge : edges) map_edges.emplace_back(unique_map_0.at(edge));

			auto mst = MinimalSpanningTree(node_nb, map_edges, costs);

			std::vector<int> map_mst;
			for (auto m : mst) map_mst.emplace_back(unique_map_1.at(m));

			return map_mst;
		}

		static std::vector<int> MinimalSpanningTree(
			const int& node_nb, const std::vector<int>& edges, const std::vector<double>& costs)
		{
			std::vector<int> mst;
			std::vector<int> nodes;

			for (int i = 0; i < node_nb; i++) nodes.push_back(i);

			std::vector<std::vector<int>> containers;
			for (int i = 0; i < nodes.size(); i++)
			{
				std::vector<int> container;
				container.push_back(nodes[i]);
				containers.push_back(container);
				std::vector<int>().swap(container);
			}

			std::vector<bool> edges_used;
			for (int i = 0; i < costs.size(); i++)
				edges_used.push_back(false);

			do
			{
				//find a minimal cost edge
				int minimal_cost_edge_index = -1;
				double minimal_cost = MAXDOUBLE;
#pragma region find_a_minimal_cost_edge

				for (int j = 0; j < costs.size(); j++)
				{
					if (!edges_used[j])
					{
						if (costs[j] < minimal_cost)
						{
							minimal_cost = costs[j];
							minimal_cost_edge_index = j;
						}
					}
				}
#pragma endregion

				if (minimal_cost_edge_index < 0)
					break;

				//check valid
				int edge_index_0 = static_cast<int>(2) * minimal_cost_edge_index;
				int edge_index_1 = static_cast<int>(2) * minimal_cost_edge_index + 1;
				int node_index_0 = edges[edge_index_0];
				int node_index_1 = edges[edge_index_1];

				if (node_index_0 == 16 && node_index_1 == 18)
				{
					int dsad = 0;
				}

				int container_0 = -1;
				int container_0_0 = -1;
				int container_1 = -1;
				int container_1_0 = -1;

				for (int j = 0; j < containers.size() && (container_0 < 0 || container_1 < 0); j++)
				{
					for (int k = 0; k < containers[j].size() && (container_0 < 0 || container_1 < 0); k++)
					{
						if (node_index_0 == containers[j][k])
						{
							container_0 = j;
							container_0_0 = k;
						}
						if (node_index_1 == containers[j][k])
						{
							container_1 = j;
							container_1_0 = k;
						}
					}
				}

				if (!(container_0 >= 0 && container_1 >= 0))
				{
					break;
				}

				if (container_0 == container_1)
				{
					edges_used[minimal_cost_edge_index] = true;
				}
				else
				{
					mst.push_back(node_index_0);
					mst.push_back(node_index_1);
					edges_used[minimal_cost_edge_index] = true;

					for (int i = 0; i < containers[container_1].size(); i++)
					{
						containers[container_0].push_back(containers[container_1][i]);
					}

					containers.erase(containers.begin() + container_1);
				}

			} while (containers.size() != 1);

			std::vector<bool>().swap(edges_used);
			std::vector<std::vector<int>>().swap(containers);

			std::vector<int>().swap(nodes);


			return mst;
		};

		static std::vector<std::vector<int>> ConnectedComponents(const int& node_nb, const std::vector<std::pair<int, int>>& tree)
		{
			std::vector<int> tree_;
			for (auto& o : tree)
			{
				tree_.emplace_back(o.first);
				tree_.emplace_back(o.second);
			}
			return ConnectedComponents(node_nb, tree_);
		}

		static std::vector<std::vector<int>> ConnectedComponentsGeneral(const Vector1i1& nodes, const std::vector<int>& tree)
		{
			int node_nb = static_cast<int>(nodes.size());

			std::map<int, int> unique_map_0;
			std::map<int, int> unique_map_1;
			for (int i = 0; i < nodes.size(); i++)
			{
				unique_map_0.insert(std::pair<int, int>(nodes[i], i));
				unique_map_1.insert(std::pair<int, int>(i, nodes[i]));
			}

			std::vector<int> map_edges;
			for (auto edge : tree) map_edges.emplace_back(unique_map_0.at(edge));

			auto components = ConnectedComponents(node_nb, map_edges);

			std::vector<std::vector<int>> map_components;
			for (auto component : components)
			{
				map_components.emplace_back(std::vector<int>());
				for (auto c : component)
					map_components.back().emplace_back(unique_map_1.at(c));
			}
			return map_components;
		}

		static std::vector<std::vector<int>> ConnectedComponents(const int& node_nb, const std::vector<int>& tree)
		{
			std::vector<std::vector<int>> components;
			std::vector<int> index(node_nb, -1);
			int nb = 0;
			for (int i = 0; i < tree.size(); i = i + 2)
			{
				int ii = i + 1;
				if (index[tree[i]] == -1 && index[tree[ii]] == -1)
				{
					index[tree[i]] = nb;
					index[tree[ii]] = nb;
					nb++;
				}
				if (index[tree[i]] == -1 && index[tree[ii]] != -1)
				{
					index[tree[i]] = index[tree[ii]];
				}
				if (index[tree[i]] != -1 && index[tree[ii]] == -1)
				{
					index[tree[ii]] = index[tree[i]];
				}

				if (index[tree[i]] != -1 && index[tree[ii]] != -1)
				{
					int min_index = std::min(index[tree[i]], index[tree[ii]]);
					int max_index = std::max(index[tree[i]], index[tree[ii]]);
					for (auto& index_ : index)
					{
						if (index_ == max_index)index_ = min_index;
					}
				}

			}

			for (int i = 0; i < nb; i++)
			{
				std::vector<int> one;
				for (int j = 0; j < index.size(); j++)
					if (index[j] == i)
						one.emplace_back(j);
				if (!one.empty())components.emplace_back(one);
			}

			for (int j = 0; j < index.size(); j++) if (index[j] == -1)components.emplace_back(std::vector<int>(1, j));

			for (auto& component : components)
				std::sort(component.begin(), component.end());
			return components;
		};

#pragma endregion

#pragma region DevelopmentRelated
		template <class Type>
		static bool CerrLine(const Type& line, const int level=0)
		{
			for (int i = 0; i < level; i++)
				std::cerr << CERR_ITER;
			std::cerr << std::to_string(line) << std::endl;
			return true;
		}

		template <class Type>
		static bool CerrLine(ofstream &file, const Type& line, const int level = 0)
		{
			for (int i = 0; i < level; i++)
			{	
				file << CERR_ITER;
				std::cerr << CERR_ITER;
			}
			file << std::to_string(line) << std::endl;
			std::cerr << std::to_string(line) << std::endl;
			return true;
		}


		//tn: total number of iterations
		//cn: current iteration
		//fn: output frequency number
		static void OutputIterInfo(const string& title, const int& tn, const int& cn, const int& fn, const int level = 0)
		{
			int delta = fn > tn ? 1 : tn / fn;

			if (cn % delta == 0)
			{
				if (cn == 0) 
				{ 
					for (int i = 0; i < level; i++)
						std::cerr << CERR_ITER;
					std::cerr << title << ": "; 
				}
				std::cerr << Functs::DoubleString((double)(100.0*cn/tn), 1) << "% ";
				if (cn + delta >= tn) std::cerr << std::endl;
			}
		}

		static void MAssert(const std::string& str, const double sleep_seconds=-1)
		{
			std::cerr <<"Bug: " << str << std::endl;
			if (sleep_seconds <= 0)
				system("pause");
			else
				MSleep(sleep_seconds);
		}

		static void MAssert(const char* str, const double sleep_seconds = -1)
		{
			std::cerr << "Bug: " << str << std::endl;
			if (sleep_seconds <= 0)
				system("pause");
			else
				MSleep(sleep_seconds);
		}

		static void MSleep(const double& second)
		{
			this_thread::sleep_for(chrono::milliseconds((int)(second*1000)));
		}

		static void RunPY(const std::string& py_path, const std::string& paras)
		{
			std::string cmd = "python " + Functs::EXP(py_path) + " " + paras;
			system(cmd.c_str());
		}

		static void RunCMD(const std::string& cmd_str)
		{
			std::cerr << "Command String: " << cmd_str << std::endl;;
			system(cmd_str.c_str());
		}
		
		static std::string WinGetUserName()
		{
			char* user = getenv("username");
			return std::string(user);
		}

		static std::string WinGetCurDirectory()
		{
			char tmp[256];
			if (_getcwd(tmp, 256)) {};
			return std::string(tmp);
		}

		static bool WinCopy(const std::string& source_file, const std::string& target_folder)
		{
			if (!Functs::DetectExisting(source_file))
			{
				MAssert("Source file does not exist: " + source_file);
				return false;
			}

			if (!Functs::DetectExisting(target_folder))
			{
				MAssert("Target folder does not exist: " + target_folder);
				return false;
			}

			std::string str = "copy " + source_file + " " + target_folder;
			std::cerr << "Command string: " << str << std::endl;
			system(str.c_str());
			return true;
		}

		static bool WinDel(const std::string& source_file)
		{
			if (!Functs::DetectExisting(source_file))
			{
				MAssert("Source file does not exist: " + source_file);
				return false;
			}
			std::string str = "del " + source_file;
			std::cerr << "Command string: " << str << std::endl;
			system(str.c_str());
			return true;
		}

		static bool WinRename(const std::string& source_file, const std::string& rename_file)
		{
			if (!Functs::DetectExisting(source_file))
			{
				MAssert("Source file does not exist: " + source_file);
				return false;
			}

			std::string str = "rename " + source_file + " " + rename_file;
			std::cerr << "Command string: " << str << std::endl;
			system(str.c_str());
			return true;
		}

		//template <class Type>
		//static void XMLP(Type& t, const tinyxml2::XMLElement* params, const std::string& name)
		//{
		//	string t_type = typeid(t).name();
		//	if (params->FirstChildElement(name.c_str()))
		//	{
		//		if(t_type== typeid(int).name())
		//			t = atoi(params->FirstChildElement(name.c_str())->GetText());
		//		if (t_type == typeid(double).name())
		//			t = atof(params->FirstChildElement(name.c_str())->GetText());
		//		if (t_type == typeid(string).name())
		//			t = params->FirstChildElement(name.c_str())->GetText();
		//	}
		//	else
		//		Functs::MAssert("Did not find the parameter: " + name);
		//};

#pragma endregion


		static Vector3d ColorMapping(const double& isolevel)
		{
			double output_c_0, output_c_1, output_c_2;
			ColorMapping(isolevel, output_c_0, output_c_1, output_c_2);
			return Vector3d(output_c_0, output_c_1, output_c_2);
		}

		static void ColorMapping(const double& isolevel, double& output_c_0, double& output_c_1, double& output_c_2)
		{
			Vector3d v;
			if (isolevel >= 0 && isolevel <= 0.25)
			{
				v[0] = 0;
				v[1] = isolevel / 0.25;
				v[2] = 1;
			}

			if (isolevel > 0.25 && isolevel <= 0.50)
			{
				v[0] = 0;
				v[1] = 1;
				v[2] = 1 - (isolevel - 0.25) / 0.25;
			}

			if (isolevel > 0.50 && isolevel <= 0.75)
			{
				v[0] = (isolevel - 0.50) / 0.25;
				v[1] = 1;
				v[2] = 0;
			}

			if (isolevel > 0.75 && isolevel <= 1.0)
			{
				v[0] = 1;
				v[1] = 1 - (isolevel - 0.75) / 0.25;
				v[2] = 0;
			}

			if (isolevel < 0.0)
			{
				v[0] = 0.0;
				v[1] = 0.0;
				v[2] = 0.0;
			}

			if (isolevel > 1.0)
			{
				v[0] = 0.5;
				v[1] = 0.0;
				v[2] = 0.0;
			}
			output_c_0 = v[0];
			output_c_1 = v[1];
			output_c_2 = v[2];
		}
	};

		//To debug a release build
		//Open the Property Pages dialog box for the project.
		//Click the C / C++ node. Set Debug Information Format to C7 compatible(/ Z7) or Program Database(/ Zi).
		//Expand Linker and click the General node.Set Enable Incremental Linking to No(/ INCREMENTAL:NO).
		//Select the Debugging node.Set Generate Debug Info to Yes(/ DEBUG).
		//Select the Optimization node.Set References to / OPT:REF and Enable COMDAT Folding to / OPT : ICF.

}
#endif 

