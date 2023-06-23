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
				data[i * (nx + 1) + j] = (-1.0 + dx * j - 0.5) * (-1.0 + dx * j - 0.5) + (1.0 - dy * i + 0.5) * (1.0 - dy * i + 0.5);
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
		double a = (cos(angle) * cos(angle)) / (lx * lx) + (sin(angle) * sin(angle)) / (ly * ly);
		double b = 2.0 * (cos(angle) * (y * sin(angle) - centerx) / (lx * lx) - sin(angle) * (y * cos(angle) - centery) / (ly * ly));
		double c = (y * sin(angle) - centerx) * (y * sin(angle) - centerx) / (lx * lx) + (y * cos(angle) - centery) * (y * cos(angle) - centery) / (ly * ly) - 1.0;
		double D = b * b - 4.0 * a * c;
		if (D >= 0.0) {
			double x1 = (-b + sqrt(D)) / (2 * a);
			double x2 = (-b - sqrt(D)) / (2 * a);
			res.push_back(x1);
			if (fabs(x1 - x2) > eps) res.push_back(x2);
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
		nx = 256;
		ny = 256;
		data = new double[(nx + 1) * (ny + 1)];
		for (int i = 0; i < (nx + 1) * (ny + 1); ++i) data[i] = 0.0;
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
		for (int i = 0; i < (nx + 1) * (ny + 1); ++i) data[i] = 0.0;
		borders[0] = -1.0;
		borders[1] = 1.0;
		borders[2] = -1.0;
		borders[3] = 1.0;
		Ellipse e(0.2, 0.2, 0.5, 0.5, 0.5236);
		std::vector<double> tmp;
		double min_ = 1000.0;
		double dx = (borders[1] - borders[0]) / nx;
		double dy = (borders[3] - borders[2]) / ny;
		double ksi = 0.0;
		for (int i = 1; i < ny + 1; ++i) {
			tmp = e.Get(borders[3] - i * dy, 01e-10);
			for (auto it = tmp.begin(); it != tmp.end(); ++it) {
				int idx = 0;
				for (int j = 1; j < nx; ++j) {
					if (fabs(*it - borders[0] - j * dx) < min_) {
						idx = j;
						min_ = fabs(*it - borders[0] - j * dx);
					}
				}
				if (tmp.size() != 0) {
					ksi = borders[0] + idx * dx - (*it);
					std::cout << ksi << std::endl;
					double inner, outer, edge;
					/*if (ksi < 0.0001) {
						outer = 2.0;
						inner = -2.0;
						edge = 0.0;
					}
					else {*/
						double det = 3.0 * (ksi * ksi) / (dx * dx);
						inner = (1 / det) * (-2.0 * ksi / dx - 4.0);
						outer = (1 / det) * (-2.0 * ksi / dx + 4.0);
						edge = 4.0 * dx / (3.0 * ksi);
					//}
					/*if (CheckRight(tmp, it)) {
						data[i * (nx + 1) + idx - 1] = -outer;
						data[i * (nx + 1) + idx] = -edge;
						data[i * (nx + 1) + idx + 1] = -inner;
					}
					else {*/
						data[i * (nx + 1) + idx - 1] = outer;
						data[i * (nx + 1) + idx] = edge;
						data[i * (nx + 1) + idx + 1] = inner;
					//}
				}
				ksi = 0.0;
				min_ = 1000.0;

			}
			tmp.clear();
		}
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
