#pragma once
#include <fstream>
#include <cmath>
#include <vector>
class USmooth {
	int nx;
	int ny;
	double dx;
	double dy;
	double* data;
public:
	USmooth() {
		nx = 256;
		ny = 256;
		dx = 2.0 / nx;
		dy = 2.0 / ny;
		data = new double[(nx + 1) * (ny + 1)];
	}
	USmooth(int nx_, int ny_) {
		nx = nx_;
		ny = ny_;
		dx = 2.0 / nx;
		dy = 2.0 / ny;
		data = new double[(nx + 1) * (ny + 1)];
		for (int i = 0; i < ny + 1; ++i) {
			for (int j = 0; j < nx + 1; ++j) {
				data[i * (nx + 1) + j] = (-1.0 + dx * j - 0.5) * (-1.0 + dx * j - 0.5) + (1.0 - dy * i + 0.5) * (1.0 - dy * i +
					0.5);
			}
		}
	}
	void GetDiff(double*& res) const {
		std::ofstream out;
		out.open("diff.txt");
		if (out.is_open())
		{
			for (int i = 0; i < (ny + 1); ++i) {
				for (int j = 0; j < nx + 1; ++j) {
					out << fabs(data[i * (nx + 1) + j] - res[i * (nx + 1) + j]) << ' ';
				}
				out << std::endl;
			}
		}
		out.close();
	}
	~USmooth() {
		delete[] data;
	}
};
bool CheckRight(std::vector<double>& dots, std::vector < double >::iterator idx) {
	if (dots[0] < dots[1]) {
		return *idx == dots[1];
	}
	else {
		return *idx == dots[0];
	}
}
class Ellipse {
	double lx, ly;
	double centerx, centery;
	double angle;
	double eps_ = 0.00001;
public:
	Ellipse() {
		
			lx = 0.2;
		ly = 0.12;
		centerx = 0.5;
		centery = 0.5;
		angle = 0.0;
	}
	Ellipse(double a_, double b_, double x_, double y_, double angle_) {
		lx = a_;
		ly = b_;
		centerx = x_;
		centery = y_;
		angle = angle_;
	}
	std::vector<double> Get(const double& y, const double& eps) {
		std::vector<double> res;
		double y0 = y - centery;
		double a = (cos(angle) * cos(angle)) / (lx * lx) + (sin(angle) * sin(angle)) / (ly * ly);
		double b = 2.0 * y0 * sin(angle) * cos(angle) * ((1 / (lx * lx)) - (1 / (ly * ly)));
		double c = y0 * y0 * (sin(angle) * sin(angle) / (lx * lx) + cos(angle) * cos(angle) / (ly * ly)) - 1.0;
		double D = b * b - 4.0 * a * c;
		if (D >= 0.0) {
			double x1 = (-b + sqrt(D)) / (2 * a);
			double x2 = (-b - sqrt(D)) / (2 * a);
			x1 += centerx;
			x2 += centerx;
			res.push_back(x1);
			res.push_back(x2);
		}
		return res;
	}
	~Ellipse() = default;
};
class Shape {
	double* data;
	int nx, ny;
	double* borders;
public:
	Shape() {
		nx = 128;
		ny = 128;
		data = new double[(nx + 1) * (ny + 1)];
		double pi = 3.141592653;
		for (int i = 0; i < ny + 1; ++i) {
			for (int j = 0; j < nx + 1; ++j) {
				double xt = borders[0] + j * 2.0 / nx;
				double yt = borders[3] - i * 2.0 / ny;
				data[i * (nx + 1) + j] = 16.0 * pi * cos(2 * pi * xt * xt + 6.0 * pi * yt * yt + pi * xt * yt) -
					(4.0 * pi * xt + pi * yt) * (4.0 * pi * xt + pi * yt) * sin(2 * pi * xt * xt + 6.0 * pi * yt * yt + pi * xt * yt)
					-
					(12.0 * pi * yt + pi * xt) * (12.0 * pi * yt + pi * xt) * sin(2 * pi * xt * xt + 6.0 * pi * yt * yt + pi * xt * yt);
			}
		}
		borders[0] = -1.0;
		borders[1] = 1.0;
		borders[2] = -1.0;
		borders[3] = 1.0;
	}
	Shape(int& nx_, int& ny_) {
		nx = nx_;
		ny = ny_;
		data = new double[(nx + 1) * (ny + 1)];
		borders = new double[4];
		for (int i = 0; i < (nx + 1) * (ny + 1); ++i) data[i] = 4.0;
		borders[0] = -1.0;
		borders[1] = 1.0;
		borders[2] = -1.0;
		borders[3] = 1.0;
		Ellipse e(0.3, 0.2, 0.5, 0.5, 0.5236);
		std::vector<double> tmp;
		double min_ = 1000.0;
		double dx = (borders[1] - borders[0]) / nx;
		double dy = (borders[3] - borders[2]) / ny;
		double ksi = 0.0;
		std::ofstream out;
		std::ofstream ell;
		out.open("ksi.txt");
		ell.open("elly.txt");
		
			for (int i = 0; i < ny / 2 + 1; ++i) {
				tmp = e.Get(borders[3] - i * dy, 01e-16);
				if (tmp.size() != 0) {
					for (auto it = tmp.begin(); it != tmp.end(); ++it) {
						int idx = 0;
						for (int j = nx / 2; j < nx + 1; ++j) {
							if (fabs(*it - borders[0] - j * dx) < min_) {
								idx = j;
								min_ = fabs(*it - borders[0] - j * dx);
							}
						}
						ell << i << ' ';
						ksi = borders[0] + idx * dx - (*it);
						double inner, outer, edge;
						/*if (ksi <= 0.01 * dx) {
						outer = 2.0;
						inner = -2.0;
						edge = 0.0;
						}*/
						//else {
						double det = 3.0 * (ksi * ksi) / (dx * dx);
						inner = -(1.0 / det) * (2.0 * ksi / dx + 4.0);
						outer = (1.0 / det) * (4.0 - 2.0 * ksi / dx);
						edge = 4.0 * dx / (3.0 * ksi);
						//}
						if (CheckRight(tmp, it)) {
							data[i * (nx + 1) + idx - 1] = -outer;
							data[i * (nx + 1) + idx] = -edge;
							data[i * (nx + 1) + idx + 1] = -inner;
						}
						else {
							data[i * (nx + 1) + idx - 1] = outer;
							data[i * (nx + 1) + idx] = edge;
							data[i * (nx + 1) + idx + 1] = inner;
						}
						out << ksi << ' ';
						ksi = 0.0;
						min_ = 1000.0;
					}
				}
				tmp.clear();
			}
		out.close();
		ell.close();
	}
	void PrintRhs() {
		std::ofstream out;
		out.open("GapRhs.txt");
		if (out.is_open()) {
			for (int i = 0; i < ny + 1; ++i) {
				for (int j = 0; j < nx + 1; ++j) {
					out << data[i * (nx + 1) + j] << ' ';
				}
				out << std::endl;
			}
		}
		out.close();
	}
	void Copy(double*& rhs) {
		for (int i = 0; i < ny + 1; ++i) {
			for (int j = 0; j < nx + 1; ++j) {
				rhs[i * (nx + 1) + j] = data[i * (nx + 1) + j];
			}
		}
	}
	~Shape() {
		delete[] data;
		delete[] borders;
	}
};