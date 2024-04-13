#include <mpi.h>
#include <iostream>
#include "Header.h"
int ipow(int a, int n) {
	int r = 1;
	for (int i = 0; i < n; i++) {
		r *= a;
	}
	return r;
}
void free(double*& Uright, double*& Uleft, double*& Utop, double*& Ubot, double*& rhs, double*& data) {
	delete[] Uright;
	delete[] Uleft;
	delete[] Utop;
	delete[] Ubot;
	delete[] rhs;
	delete[] data;
}
int GetH(int& size) {
	int res = 1;
	int tmp = 0;
	while (true) {
		tmp += ipow(4, res - 1);
		if (tmp >= size) break;
		res++;
	}
	return res;
}
void InterpolateY(double*& border, const double& dx, const double& dy, int nx, int ny, double*& right,
	double*& left, double*& top, double*& bot, double& f, double*& borders, double*& center, int size_) {
	double lft = 0.0, rgt = 0.0, cntr = 0.0, lft_bdr = borders[3], cntr_bdr = 0.0, rgt_bdr = 0.0;
	double between = (borders[3] - borders[2]) / (size_ - 1);
	double step = 2.0 * between;
	int idx = 0;
	for (int st = 0; st < size_ - 2; st += 2) {
		lft = center[st];
		cntr = center[st + 1];
		rgt = center[st + 2];
		rgt_bdr = lft_bdr - step;
		cntr_bdr = lft_bdr - between;
		for (int i = idx; i < idx + 2 * ny / (size_ - 1); ++i) {
			double yk = borders[3] - dy * i;
			border[i] = lft * (yk - cntr_bdr) * (yk - rgt_bdr) / ((lft_bdr - cntr_bdr) * (lft_bdr - rgt_bdr)) +
				cntr * (yk - lft_bdr) * (yk - rgt_bdr) / ((cntr_bdr - lft_bdr) * (cntr_bdr - rgt_bdr)) +
				rgt * (yk - lft_bdr) * (yk - cntr_bdr) / ((rgt_bdr - lft_bdr) * (rgt_bdr - cntr_bdr));
		}
		idx += 2 * ny / (size_ - 1);
		lft_bdr -= step;
	}
	
	if (rgt_bdr != borders[2]) {
		lft = center[size_ - 3];
		cntr = center[size_ - 2];
		rgt = center[size_ - 1];
		lft_bdr = borders[2] + step;
		rgt_bdr = borders[2];
		cntr_bdr = borders[2] + between;
		for (int i = ny - 2 * ny / (size_ - 1); i < ny; ++i) {
			double yk = borders[3] - dy * i;
			border[i] = lft * (yk - cntr_bdr) * (yk - rgt_bdr) / ((lft_bdr - cntr_bdr) * (lft_bdr - rgt_bdr)) +
				cntr * (yk - lft_bdr) * (yk - rgt_bdr) / ((cntr_bdr - lft_bdr) * (cntr_bdr - rgt_bdr)) +
				rgt * (yk - lft_bdr) * (yk - cntr_bdr) / ((rgt_bdr - lft_bdr) * (rgt_bdr - cntr_bdr));
		}
	}
	
	border[0] = top[nx / 2];
	border[ny] = bot[nx / 2];
}

void InterpolateX(double*& border, const double& dx, const double& dy, int nx, int ny, double*& right, double*& left,
	double*& top, double*& bot, double& f, double*& borders, double*& center, int size_) {
	double lft = 0.0, rgt = 0.0, cntr = 0.0, lft_bdr = borders[0], cntr_bdr = 0.0, rgt_bdr = 0.0;
	double between = (borders[1] - borders[0]) / (size_ - 1);
	double step = 2.0 * between;
	int idx = 0;
	for (int st = 0; st < size_ - 2; st += 2) {
		lft = center[st];
		cntr = center[st + 1];
		rgt = center[st + 2];
		rgt_bdr = lft_bdr + step;
		cntr_bdr = lft_bdr + between;
		for (int i = idx; i <= idx + nx / (size_ - 1); ++i) {
			double xk = borders[0] + dx * i;
			border[i] = lft * (xk - cntr_bdr) * (xk - rgt_bdr) / ((lft_bdr - cntr_bdr) * (lft_bdr - rgt_bdr)) +
				cntr * (xk - lft_bdr) * (xk - rgt_bdr) / ((cntr_bdr - lft_bdr) * (cntr_bdr - rgt_bdr)) +
				rgt * (xk - lft_bdr) * (xk - cntr_bdr) / ((rgt_bdr - lft_bdr) * (rgt_bdr - cntr_bdr));
		}
		idx += nx / (size_ - 1);
		lft_bdr += step;
	}

	if (rgt_bdr != borders[1]) {
		lft = center[size_ - 3];
		cntr = center[size_ - 2];
		rgt = center[size_ - 1];
		lft_bdr = borders[1] - step;
		rgt_bdr = borders[1];
		cntr_bdr = borders[1] - between;
		for (int i = nx / 2 - nx / (size_ - 1); i <= nx / 2; ++i) {
			double yk = borders[0] + dx * i;
			border[i] = lft * (yk - cntr_bdr) * (yk - rgt_bdr) / ((lft_bdr - cntr_bdr) * (lft_bdr - rgt_bdr)) +
				cntr * (yk - lft_bdr) * (yk - rgt_bdr) / ((cntr_bdr - lft_bdr) * (cntr_bdr - rgt_bdr)) +
				rgt * (yk - lft_bdr) * (yk - cntr_bdr) / ((rgt_bdr - lft_bdr) * (rgt_bdr - cntr_bdr));
		}
	}

	border[0] = left[ny / 2];
	border[nx / 2] = right[ny / 2];
}

void step(size_t rank, int& nx, int& ny, double& dx, double& dy, double*& rhs, double*& Utop, double*& Ubot, double*& Uleft, double*&
	Uright, const int* kids, const int parent, double*& borders) {
	MPI_Status status;
	if (rank != 0) {
		MPI_Recv(rhs, (ny + 1) * (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(Utop, (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(Ubot, (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(Uright, (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(Uleft, (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(borders, 4, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	}
	double* Eborder = new double[4];
	Eborder[0] = 0.5 * (borders[0] + borders[1]);
	Eborder[1] = borders[1];
	Eborder[2] = borders[2];
	Eborder[3] = borders[3];
	double* Oborder = new double[4];
	Oborder[0] = borders[0];
	Oborder[1] = Eborder[0];
	Oborder[2] = borders[2];
	Oborder[3] = borders[3];
	double* b00 = new double[4];
	b00[0] = Eborder[0];
	b00[1] = Eborder[1];
	b00[2] = 0.5 * (Eborder[2] + Eborder[3]);
	b00[3] = Eborder[3];
	double* b01 = new double[4];
	b01[0] = Eborder[0];
	b01[1] = Eborder[1];
	b01[2] = Eborder[2];
	b01[3] = b00[2];
	double* b10 = new double[4];
	b10[0] = Oborder[0];
	b10[1] = Oborder[1];
	b10[2] = b00[2];
	b10[3] = b00[3];
	double* b11 = new double[4];
	b11[0] = Oborder[0];
	b11[1] = Oborder[1];
	b11[2] = b01[2];
	b11[3] = b01[3];

	double* EvenRhs = new double[(ny + 1) * (nx / 2 + 1)];
	double* OddRhs = new double[(ny + 1) * (nx / 2 + 1)];
	for (int i = 0; i < ny + 1; ++i) {
		for (int j = nx / 2; j < nx + 1; ++j) EvenRhs[i * (nx / 2 + 1) + j - nx / 2] = 0.5 * (rhs[i * (nx + 1) + j] + rhs[i * (nx + 1) + nx - j]);
		for (int j = 0; j < (nx / 2) + 1; ++j) OddRhs[i * (nx / 2 + 1) + j] = 0.5 * (rhs[i * (nx + 1) + j] - rhs[i * (nx + 1) + nx - j]);
	}
	double* UtopE = new double[nx / 2 + 1];
	double* UbotE = new double[nx / 2 + 1];
	double* UtopO = new double[nx / 2 + 1];
	
		double* UbotO = new double[nx / 2 + 1];
	for (int i = (nx / 2); i < nx + 1; ++i) {
		UtopE[i - nx / 2] = 0.5 * (Utop[i] + Utop[nx - i]);
		UtopO[i - nx / 2] = 0.5 * (Utop[i - (nx / 2)] - Utop[nx - i + nx / 2]);
		UbotE[i - nx / 2] = 0.5 * (Ubot[i] + Ubot[nx - i]);
		UbotO[i - nx / 2] = 0.5 * (Ubot[i - (nx / 2)] - Ubot[nx - (i - nx / 2)]);
	}
	double* UrightE = new double[ny + 1];
	double* UrightO = new double[ny + 1];
	double* UleftE = new double[ny + 1];
	double* UleftO = new double[ny + 1];
	for (int i = 0; i < ny + 1; ++i) {
		UrightE[i] = 0.5 * (Uright[i] + Uleft[i]);
		UleftO[i] = 0.5 * (Uleft[i] - Uright[i]);
		UrightO[i] = 0.0;
	}

	int size_ = 33;
	double* tmp_rhs = new double[size_ * size_];
	double* gz_data = new double[size_ * size_];
	for (int i = 0; i < size_; ++i) {
		for (int j = 0; j < size_; ++j) {
			tmp_rhs[i * size_ + j] = rhs[i * (ny / (size_ - 1)) * (nx + 1) + j * (nx / (size_ - 1))];
		}
	}
	for (int i = 0; i < size_; ++i) {
		for (int j = 0; j < size_; ++j) {
			gz_data[i * size_ + j] = 0.0;
		}
	}
	
	for (int i = 0; i < size_; ++i) {
		gz_data[i] = Utop[i * (nx / (size_ - 1))];
		gz_data[size_ * (size_ - 1) + i] = Ubot[i * (nx / (size_ - 1))];
		gz_data[i * size_] = Uleft[i * nx / (size_ - 1)];
		gz_data[size_ * i + size_ - 1] = Uright[i * nx / (size_ - 1)];
	}
	double  max_ = 0.0;
	double eps = 1e-6;
	double h = (borders[1] - borders[0]) / (size_ - 1);
	do {
		max_ = 0.0;
		for (int i = 1; i < size_ - 1; i++) {
			for (int j = 1; j < size_ - 1; j++) {
				double u0 = gz_data[i * size_ + j];
				gz_data[i * size_ + j] = 0.5 * (h * h * (gz_data[(i - 1) * size_ + j] + gz_data[(i + 1) * size_ + j]) / (h * h + h * h) 
					+ h * h * (gz_data[i * size_ + j - 1] + gz_data[i * size_ + j + 1]) / (h * h + h * h) - 
					h * h * h * h * tmp_rhs[i * size_ + j] / (h * h + h * h));
				double d = fabs(gz_data[i * size_ + j] - u0);
				if (d > max_) max_ = d;
			}
		}
	} while (max_ > eps);
	double* cntr = new double[size_];
	for (int i = 0; i < size_; ++i) {
		cntr[i] = gz_data[size_ / 2 + size_ * i];
	}
	delete[] gz_data;
	delete[] tmp_rhs;

	InterpolateY(UleftE, dx, dy, nx, ny, Uright, Uleft, Utop, Ubot, rhs[(nx + 1) * ny / 2 + nx / 2], borders, cntr, size_);
	delete[] cntr;

	double* F00 = new double[(1 + nx / 2) * (1 + ny / 2)];
	double* F01 = new double[(1 + nx / 2) * (1 + ny / 2)];
	double* F10 = new double[(1 + nx / 2) * (1 + ny / 2)];
	double* F11 = new double[(1 + nx / 2) * (1 + ny / 2)];
	for (int i = 0; i < 1 + (ny / 2); ++i) {
		for (int j = 0; j < 1 + (nx / 2); ++j) {
			F00[i * (1 + nx / 2) + j] = 0.5 * (EvenRhs[i * (1 + nx / 2) + j] + EvenRhs[(ny - i) * (1 + nx / 2) + j]);
			F10[i * (1 + nx / 2) + j] = 0.5 * (OddRhs[i * (1 + nx / 2) + j] + OddRhs[(ny - i) * (1 + nx / 2) + j]);
		}
	}
	for (int i = ny / 2; i < ny + 1; ++i) {
		for (int j = 0; j < 1 + (nx / 2); ++j) {
			F01[(i - ny / 2) * (1 + nx / 2) + j] = 0.5 * (EvenRhs[i * (1 + nx / 2) + j] - EvenRhs[(ny - i) * (1 + nx / 2) + j]);
			F11[(i - ny / 2) * (1 + nx / 2) + j] = 0.5 * (OddRhs[i * (1 + nx / 2) + j] - OddRhs[(ny - i) * (1 + nx / 2) + j]);
		}
	}
	double* top00 = new double[1 + nx / 2];
	double* top01 = new double[1 + nx / 2];
	double* bot00 = new double[1 + nx / 2];
	double* bot01 = new double[1 + nx / 2];
	double* right00 = new double[1 + ny / 2];
	double* right01 = new double[1 + ny / 2];
	double* left00 = new double[1 + ny / 2];
	double* left01 = new double[1 + ny / 2];
	for (int i = 0; i < 1 + nx / 2; ++i) {
		top00[i] = 0.5 * (UtopE[i] + UbotE[i]);
		bot01[i] = 0.5 * (UbotE[i] - UtopE[i]);
		top01[i] = 0.0;
	}
	for (int i = 0; i < 1 + ny / 2; ++i) {
		right00[i] = 0.5 * (UrightE[i] + UrightE[ny - i]);
		right01[i] = 0.5 * (UrightE[i + ny / 2] - UrightE[(ny / 2) - i]);
	}
	for (int i = 0; i < 1 + ny / 2; ++i) {
		left00[i] = 0.5 * (UleftE[i] + UleftE[ny - i]);
		left01[i] = 0.5 * (UleftE[i + ny / 2] - UleftE[(ny / 2) - i]);
	}

	tmp_rhs = new double[size_ * size_];
	gz_data = new double[size_ * size_];
	for (int i = 0; i < size_; ++i) {
		for (int j = 0; j < size_; ++j) {
			tmp_rhs[i * size_ + j] = EvenRhs[i * (ny / (size_ - 1)) * (1 + nx / 2) + j * ((nx / 2) / (size_ - 1))];
		}
	}

	for (int i = 0; i < size_; ++i) {
		for (int j = 0; j < size_; ++j) {
			gz_data[i * size_ + j] = 0.0;
		}
	}

	for (int i = 0; i < size_; ++i) {
		gz_data[i] = UtopE[i * ((nx / 2) / (size_ - 1))];
		gz_data[size_ * (size_ - 1) + i] = UbotE[i * ((nx / 2) / (size_ - 1))];
		gz_data[i * size_] = UleftE[i * ny / (size_ - 1)];
		gz_data[size_ * i + size_ - 1] = UrightE[i * ny / (size_ - 1)];
	}

	max_ = 0.0;
	eps = 1e-6;
	int cnt = 0;
	double hx = (Eborder[1] - Eborder[0]) / (size_ - 1);
	double hy = (Eborder[3] - Eborder[2]) / (size_ - 1);
	do {
		max_ = 0.0;
		for (int i = 1; i < size_ - 1; i++) {
			for (int j = 1; j < size_ - 1; j++) {
				double u0 = gz_data[i * size_ + j];
				gz_data[i * size_ + j] = 0.5 * (hx * hx * (gz_data[(i - 1) * size_ + j] + gz_data[(i + 1) * size_ + j]) / (hx * hx + hy * hy) +
					hy * hy * (gz_data[i * size_ + j - 1] + gz_data[i * size_ + j + 1]) / (hx * hx + hy * hy) - 
					(hx * hx * hy * hy) * tmp_rhs[i * size_ + j] / (hx * hx + hy * hy));
				double d = fabs(gz_data[i * size_ + j] - u0);
				if (d > max_) max_ = d;
			}
		}
		cnt += 1;
	} while (max_ > eps);
	//std::cout << cnt << "\n";
	double* center = new double[size_];
	for (int i = 0; i < size_; ++i) center[i] = gz_data[size_ * (size_ / 2) + i];
	delete[] gz_data;
	delete[] tmp_rhs;
	//std::cout << "EVEN" << "\n";
	InterpolateX(bot00, dx, dy, nx, ny, UrightE, UleftE, UtopE, UbotE, EvenRhs[(1 + nx / 2) * (ny / 2) + (nx / 4)], Eborder, center, size_);
	delete[] center;

	double* top10 = new double[1 + nx / 2];
	double* top11 = new double[1 + nx / 2];
	double* bot10 = new double[1 + nx / 2];
	double* bot11 = new double[1 + nx / 2];
	double* right10 = new double[1 + ny / 2];
	double* right11 = new double[1 + ny / 2];
	double* left10 = new double[1 + ny / 2];
	double* left11 = new double[1 + ny / 2];
	for (int i = 0; i < 1 + nx / 2; ++i) {
		top10[i] = 0.5 * (UtopO[i] + UbotO[i]);
		bot11[i] = 0.5 * (UbotO[i] - UtopO[i]);
		top11[i] = 0.0;
	}
	for (int i = 0; i < 1 + ny / 2; ++i) {
		right10[i] = 0.5 * (UrightO[i] + UrightO[ny - i]);
		
			right11[i] = 0.5 * (UrightO[i + ny / 2] - UrightO[(ny / 2) - i]);
	}
	for (int i = 0; i < 1 + ny / 2; ++i) {
		left10[i] = 0.5 * (UleftO[i] + UleftO[ny - i]);
		left11[i] = 0.5 * (UleftO[i + ny / 2] - UleftO[(ny / 2) - i]);
	}

	tmp_rhs = new double[size_ * size_];
	gz_data = new double[size_ * size_];
	for (int i = 0; i < size_; ++i) {
		for (int j = 0; j < size_; ++j) {
			tmp_rhs[i * size_ + j] = OddRhs[i * (ny / (size_ - 1)) * (1 + nx / 2) + j * ((nx / 2) / (size_ - 1))];
		}
	}
	for (int i = 0; i < size_; ++i) {
		for (int j = 0; j < size_; ++j) {
			gz_data[i * size_ + j] = 0.0;
		}
	}

	for (int i = 0; i < size_; ++i) {
		gz_data[i] = UtopO[i * ((nx / 2) / (size_ - 1))];
		gz_data[size_ * (size_ - 1) + i] = UbotO[i * ((nx / 2) / (size_ - 1))];
		gz_data[i * size_] = UleftO[i * ny / (size_ - 1)];
		gz_data[size_ * i + size_ - 1] = UrightO[i * ny / (size_ - 1)];
	}

	max_ = 0.0;
	eps = 1e-6;
	cnt = 0;
	hx = (Oborder[1] - Oborder[0]) / (size_ - 1);
	hy = (Oborder[3] - Oborder[2]) / (size_ - 1);
	do {
		max_ = 0.0;
		for (int i = 1; i < size_ - 1; i++) {
			for (int j = 1; j < size_ - 1; j++) {
				double u0 = gz_data[i * size_ + j];
				gz_data[i * size_ + j] = 0.5 * (hx * hx * (gz_data[(i - 1) * size_ + j] + gz_data[(i + 1) * size_ + j]) / (hx * hx + hy * hy) +
					hy * hy * (gz_data[i * size_ + j - 1] + gz_data[i * size_ + j + 1]) / (hx * hx + hy * hy) -
					hx * hx * hy * hy * tmp_rhs[i * size_ + j] / (hx * hx + hy * hy));
				double d = fabs(gz_data[i * size_ + j] - u0);
				if (d > max_) max_ = d;
			}
		}
		cnt += 1;
	} while (max_ > eps);
	//std::cout << cnt << '\n';
	center = new double[size_];
	for (int i = 0; i < size_; ++i) center[i] = gz_data[size_ * (size_ / 2) + i];
	delete[] gz_data;
	delete[] tmp_rhs;
	//std::cout << "ODD" << "\n";
	InterpolateX(bot10, dx, dy, nx, ny, UrightO, UleftO, UtopO, UbotO, OddRhs[(1 + nx / 2) * ny / 2 + nx / 4], Oborder, center, size_);
	delete[] center;
	//to 00:
	MPI_Send(F00, (1 + nx / 2) * (1 + ny / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
	MPI_Send(top00, (1 + nx / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
	MPI_Send(bot00, (1 + nx / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
	MPI_Send(right00, (1 + ny / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
	MPI_Send(left00, (1 + ny / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
	MPI_Send(b00, 4, MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
	//to 01:
	MPI_Send(F01, (1 + nx / 2) * (1 + ny / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD);
	MPI_Send(top01, (1 + nx / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD);
	MPI_Send(bot01, (1 + nx / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD);
	MPI_Send(right01, (1 + ny / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD);
	MPI_Send(left01, (1 + ny / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD);
	MPI_Send(b01, 4, MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD);
	//to 10:
	MPI_Send(F10, (1 + nx / 2) * (1 + ny / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD);
	MPI_Send(top10, (1 + nx / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD);
	MPI_Send(bot10, (1 + nx / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD);
	MPI_Send(right10, (1 + ny / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD);
	MPI_Send(left10, (1 + ny / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD);
	MPI_Send(b10, 4, MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD);
	//to 11:
	MPI_Send(F11, (1 + nx / 2) * (1 + ny / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD);
	MPI_Send(top11, (1 + nx / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD);
	MPI_Send(bot11, (1 + nx / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD);
	MPI_Send(right11, (1 + ny / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD);
	MPI_Send(left11, (1 + ny / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD);
	MPI_Send(b11, 4, MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD);
	delete[] b00;
	delete[] b01;
	delete[] b10;
	delete[] b11;
}
int StageStart(int N) {
	int res = 0;
	for (int i = 0; i < N - 1; ++i) {
		res += ipow(4, i);
	}
	return res;
}
void prod(int rank, int& nx, int& ny, double& dx, double*& rhs,
	double*& Utop, double*& Ubot, double*& Uleft, double*& Uright, const int* kids, const int parent, double*& data) {
	double eps = 0.00001;
	for (int i = 0; i < ny + 1; ++i) {
		for (int j = 0; j < nx + 1; ++j) {
			data[i * (nx + 1) + j] = 0.0;
		}
	}
	for (int i = 0; i < nx + 1; ++i) {
		data[i] = Utop[i];
		data[(nx + 1) * ny + i] = Ubot[i];
	}
	for (int i = 0; i < ny + 1; ++i) {
		data[i * (nx + 1)] = Uleft[i];
		data[nx + i * (nx + 1)] = Uright[i];
	}
	double max_ = 0.0;
	do {
		max_ = 0.0;
		for (int i = 1; i < ny; i++)
			for (int j = 1; j < nx; j++) {
				double u0 = data[i * (nx + 1) + j];
				
					data[i * (nx + 1) + j] = 0.25 * (data[(i - 1) * (nx + 1) + j] + data[(i + 1) * (nx + 1) + j] + data[i * (nx + 1) + j - 1]
						+ data[i * (nx + 1) + j + 1] - dx * dx * rhs[i * (nx + 1) + j]);
				double d = fabs(data[i * (nx + 1) + j] - u0);
				if (d > max_) max_ = d;
			}
	} while (max_ > eps);
}
void Gauss(int rank, int& nx, int& ny, double& dx, double*& rhs,
	double*& Utop, double*& Ubot, double*& Uleft, double*& Uright, const int* kids, const int parent, double*& data) {
	MPI_Status status;
	MPI_Recv(rhs, (ny + 1) * (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(Utop, (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(Ubot, (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(Uright, (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(Uleft, (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	double* b = new double[4];
	int counter = 0;
	MPI_Recv(b, 4, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
	double eps = 0.00001;
	for (int i = 0; i < ny + 1; ++i) {
		for (int j = 0; j < nx + 1; ++j) {
			data[i * (nx + 1) + j] = 0.0;
		}
	}
	for (int i = 0; i < nx + 1; ++i) {
		data[i] = Utop[i];
		data[(nx + 1) * ny + i] = Ubot[i];
	}
	for (int i = 0; i < ny + 1; ++i) {
		data[i * (nx + 1)] = Uleft[i];
		data[nx + i * (nx + 1)] = Uright[i];
	}
	double max_ = 0.0;
	do {
		counter++;
		max_ = 0.0;
		for (int i = 1; i < ny; i++)
			for (int j = 1; j < nx; j++) {
				double u0 = data[i * (nx + 1) + j];
				data[i * (nx + 1) + j] = 0.25 * (data[(i - 1) * (nx + 1) + j] + data[(i + 1) * (nx + 1) + j] + data[i * (nx + 1) + j - 1]
					+ data[i * (nx + 1) + j + 1] - dx * dx * rhs[i * (nx + 1) + j]);
				double d = fabs(data[i * (nx + 1) + j] - u0);
				if (d > max_) max_ = d;
			}
	} while (max_ > eps);
	delete[] b;
}
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double dx = 0.0, dy = 0.0, x1 = -1.0, x2 = 1.0, y1 = -1.0, y2 = 1.0;
	int nx = 1024, ny = 1024;
	if (rank == 0) {
		std::cin >> nx;
		std::cin >> ny;
	}
	MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	dx = (x2 - x1) / nx;
	dy = (y2 - y1) / ny;
	//найти уровень каждого процесса
	int h = 1;
	int tmp = 0;
	while (tmp <= size) {
		if (rank > tmp) h++;
		tmp += ipow(4, h - 1);
	}
	nx /= ipow(2, h - 1);
	
		ny /= ipow(2, h - 1);
	int kids[4] = { 4 * rank + 1, 4 * rank + 2, 4 * rank + 3, 4 * rank + 4 };
	int parent = StageStart(h - 1) + (rank - StageStart(h)) / 4;
	if (rank > 0 && rank <= 4) parent = 0;
	double* rhs = new double[(ny + 1) * (nx + 1)];
	double* Uleft = new double[ny + 1];
	double* Uright = new double[ny + 1];
	double* Ubot = new double[nx + 1];
	double* Utop = new double[nx + 1];
	double* data = new double[(ny + 1) * (nx + 1)];
	for (int i = 0; i < (ny + 1) * (nx + 1); ++i) data[i] = 0.0;
	double* getres1 = new double[(1 + ny / 2) * (1 + nx / 2)];
	double* getres2 = new double[(1 + ny / 2) * (1 + nx / 2)];
	double* getres3 = new double[(1 + ny / 2) * (1 + nx / 2)];
	double* getres4 = new double[(1 + ny / 2) * (1 + nx / 2)];
	double* tmp_odd = new double[(ny + 1) * (1 + nx / 2)];
	for (int i = 0; i < (1 + ny / 2) * (1 + nx / 2); ++i) {
		getres1[i] = 0.0;
		getres2[i] = 0.0;
		getres3[i] = 0.0;
		getres4[i] = 0.0;
	}
	if (rank == 0) {
		h = 1;
		double pi = 3.141592653;
		for (int i = 0; i < ny + 1; ++i) {
			for (int j = 0; j < nx + 1; ++j) rhs[i * (nx + 1) + j] = 4.0;
		}
		
		for (int i = 0; i < ny + 1; ++i) {
			for (int j = 0; j < nx + 1; ++j) {
				double xt = -1.0 + j * 2.0 / nx;
				double yt = 1.0 - i * 2.0 / ny;
				rhs[i * (nx + 1) + j] = 16.0 * pi * cos(2 * pi * xt * xt + 6.0 * pi * yt * yt + pi * xt * yt) -
					(4.0 * pi * xt + pi * yt) * (4.0 * pi * xt + pi * yt) * sin(2 * pi * xt * xt + 6.0 * pi * yt * yt + pi * xt * yt)
					-
					(12.0 * pi * yt + pi * xt) * (12.0 * pi * yt + pi * xt) * sin(2 * pi * xt * xt + 6.0 * pi * yt * yt + pi * xt * yt);
			}
		}
		
		std::ofstream analyt;
		analyt.open("GapRhs.txt");
		if (analyt.is_open())
		{
			for (int i = 0; i < (ny + 1); ++i) {
				for (int j = 0; j < nx + 1; ++j) {
					double xt = -1.0 + j * 2.0 / nx;
					double yt = 1.0 - i * 2.0 / ny;
					analyt << sin(2 * pi * xt * xt + 6 * pi * yt * yt + pi * xt * yt) << ' '; //(xt - 0.5) * (xt - 0.5) + (yt + 0.5) * (yt + 0.5) << ' ';
				}
				analyt << std::endl;
			}
		}
		analyt.close();
		
		
		//Shape shape(nx, ny);
		//shape.PrintRhs();
		//shape.Copy(rhs);
		//for (int i = 0; i < ny + 1; ++i) Uleft[i] = 2.25 + (0.5 + (y2 - static_cast<double>(i * dy))) * (0.5 + (y2 - static_cast<double>(i * dy)));
		//for (int i = 0; i < ny + 1; ++i) Uright[i] = 0.25 + (0.5 + (y2 - static_cast<double>(i * dy))) * (0.5 + (y2 - static_cast<double>(i * dy)));
		//for (int i = 0; i < nx + 1; ++i) Ubot[i] = 0.25 + (-0.5 + (x1 + static_cast<double>(i * dx))) * (-0.5 + (x1 + static_cast<double>(i * dx)));
		//for (int i = 0; i < nx + 1; ++i) Utop[i] = 2.25 + (-0.5 + (x1 + static_cast<double>(i * dx))) * (-0.5 + (x1 + static_cast<double>(i * dx)));
		for (int i = 0; i < ny + 1; ++i) Uleft[i] = sin(2 * pi + 6.0 * pi * (y2 - dy * i) * (y2 - dy * i) - pi * (y2 - dy * i));
		for (int i = 0; i < ny + 1; ++i) Uright[i] = sin(2 * pi + 6.0 * pi * (y2 - dy * i) * (y2 - dy * i) + pi * (y2 - dy * i));
		for (int i = 0; i < nx + 1; ++i) Ubot[i] = sin(2 * pi * (x1 + dx * i) * (x1 + dx * i) + 6.0 * pi - pi * (x1 + dx * i));
		for (int i = 0; i < nx + 1; ++i) Utop[i] = sin(2 * pi * (x1 + dx * i) * (x1 + dx * i) + 6.0 * pi + pi * (x1 + dx * i));
	}
	int hight = GetH(size);
	//x1, x2, y1, y2 для каждого процесса
	double* borders = new double[4];
	if (rank == 0) {
		borders[0] = -1.0;
		borders[1] = 1.0;
		borders[2] = -1.0;
		borders[3] = 1.0;
	}
	double starttime = MPI_Wtime();
	/*if (rank == 0) {
	prod(rank, nx, ny, dx, rhs, Utop, Ubot, Uleft, Uright, kids, parent, data);
	std::cout << MPI_Wtime() - starttime << std::endl;
	}*/
	for (int st = 1; st <= hight; ++st) {
		if (h == st) {
			if (h == hight) {
				Gauss(rank, nx, ny, dx, rhs, Utop, Ubot, Uleft, Uright, kids, parent, data);
			}
			else step(rank, nx, ny, dx, dy, rhs, Utop, Ubot, Uleft, Uright, kids, parent, borders);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (h == hight) MPI_Send(data, (nx + 1) * (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
	for (int st = hight - 1; st >= 0; --st) {
		if (h == st) {
			MPI_Status status;
			MPI_Recv(getres1, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD, &status);
			MPI_Recv(getres2, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD, &status);
			MPI_Recv(getres3, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD, &status);
			MPI_Recv(getres4, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD, &status);
			for (int i = 0; i < ny / 2 + 1; ++i) {
				for (int j = 0; j < nx / 2 + 1; ++j) {
					data[i * (1 + nx) + j] += getres1[i * (1 + nx / 2) + j];
					data[i * (1 + nx) + nx / 2 + j] += getres2[i * (1 + nx / 2) + j];
					data[(i + ny / 2) * (1 + nx) + j] += getres3[i * (1 + nx / 2) + j];
					data[(i + ny / 2) * (1 + nx) + nx / 2 + j] += getres4[i * (1 + nx / 2) + j];
				}
				
			}
			for (int i = 0; i < nx + 1; ++i) {
				for (int j = 0; j < ny / 2 + 1; ++j) {
					data[j * (nx + 1) + i] -= data[(ny - j) * (nx + 1) + i];
				}
			}
			for (int i = 0; i < nx / 2 + 1; ++i) {
				for (int j = 0; j < ny / 2 + 1; ++j) {
					data[(j + ny / 2) * (nx + 1) + i] += getres1[(ny / 2 - j) * (nx / 2 + 1) + i];
					data[(j + ny / 2) * (nx + 1) + nx / 2 + i] += getres2[(ny / 2 - j) * (nx / 2 + 1) + i];
				}
			}
			for (int k = 0; k < ny + 1; ++k) {
				for (int q = 0; q < 1 + nx / 2; ++q) {
					tmp_odd[k * (1 + nx / 2) + q] = data[k * (1 + nx) + q];
				}
			}
			for (int k = 0; k < ny + 1; ++k) {
				for (int q = 0; q < 1 + nx / 2; ++q) {
					data[k * (1 + nx) + q] += data[k * (1 + nx) + nx - q];
				}
			}
			for (int k = 0; k < ny + 1; ++k) {
				for (int q = 0; q < 1 + nx / 2; ++q) {
					data[k * (1 + nx) + q + nx / 2] -= tmp_odd[k * (1 + nx / 2) + nx / 2 - q];
				}
			}
			if (h != 1) {
				MPI_Send(data, (nx + 1) * (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
			}
		}
	}
	double endtime = MPI_Wtime();
	if (rank == 0) std::cout << "time : " << endtime - starttime << std::endl;
	if (h == 1) {
		for (int k = 1; k < ny; ++k) {
			for (int q = 1; q < nx; ++q) {
				data[k * (nx + 1) + q] = (dx * dx * dy * dy / (2.0 * (dx * dx + dy * dy))) *
					((data[(k + 1) * (nx + 1) + q] + data[(k - 1) * (nx + 1) + q]) / (dx * dx) +
						(data[k * (nx + 1) + q + 1] + data[k * (nx + 1) + q - 1]) / (dy * dy) - rhs[k * (nx + 1) + q]);
			}
		}
		std::ofstream out;
		out.open("hello.txt");
		if (out.is_open())
		{
			for (int i = 0; i < (ny + 1); ++i) {
				for (int j = 0; j < nx + 1; ++j) {
					out << data[i * (nx + 1) + j] << ' ';
				}
				out << std::endl;
			}
		}
		out.close();

		if (rank == 0) std::cout << "fullfilled" << std::endl;
		USmooth res = USmooth(nx, ny);
		res.GetDiff(data);
	}
	delete[] Uright;
	delete[] Uleft;
	delete[] Utop;
	delete[] Ubot;
	delete[] rhs;
	delete[] data;
	delete[] getres1;
	delete[] getres2;
	delete[] getres3;
	delete[] getres4;
	delete[] tmp_odd;
	MPI_Finalize();
	return 0;
}

