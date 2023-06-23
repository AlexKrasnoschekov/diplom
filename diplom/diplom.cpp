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
    double*& left, double*& top, double*& bot, double& f, double*& borders) {
    double center = 0.2 * (right[ny / 2] + left[ny / 2] + bot[nx / 2] + top[nx / 2]) 
        + 0.05 * (top[nx] + top[0] + bot[nx] + bot[0]) - 0.3 * f * dy * dy;
    border[0] = top[nx / 2];
    border[ny] = bot[nx / 2];
    for (int i = 1; i < ny; ++i) {
        double yk = borders[3] - dy * i;
        double yc = (borders[3] + borders[2]) / 2;
        border[i] = top[nx / 2] * (yk - borders[2]) * (yk - yc) / ((borders[3] - yc) * (borders[3] - borders[2]))
            + center * (yk - borders[3]) * (yk - borders[2]) / ((yc - borders[3]) * (yc - borders[2])) 
            + bot[nx / 2] * (yk  - borders[3]) * (yk - yc) / ((borders[2] - borders[3]) * (borders[2] - yc));
    }
}

void InterpolateX(double*& border, const double& dx, const double& dy, int nx, int ny, double*& right, double*& left,
    double*& top, double*& bot, double& f, double*& borders) {
    double center = 0.2 * (right[ny / 2] + left[ny / 2] + bot[nx / 4] + top[nx / 4]) 
        + 0.05 * (top[nx / 2] + top[0] + bot[nx / 2] + bot[0]) - 0.3 * f * dx * dx;
    border[0] = left[ny / 2];
    border[nx / 2] = right[ny / 2];
    for (int i = 1; i < nx / 2; ++i) {
        double xk = borders[0] + i * dx;
        double xc = (borders[0] + borders[1]) / 2;
        border[i] = left[ny / 2] * (xk - xc) * (xk - borders[1]) / ((borders[0] - borders[1]) * (borders[0] - xc)) 
            + center * (xk - borders[0]) * (xk - borders[1]) / ((xc - borders[0]) * (xc - borders[1])) 
            + right[ny / 2] * (xk - borders[0]) * (xk - xc) / ((borders[1] - xc) * (borders[1] - borders[0]));
    }
}


void step(size_t rank, int& nx, int& ny, double& dx, double& dy, double*& rhs, double*& Utop, double*& Ubot, double*& Uleft, double*& Uright, const int* kids, const int parent, double*& borders) {
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
    InterpolateY(UleftE, dx, dy, nx, ny, Uright, Uleft, Utop, Ubot, rhs[(nx + 1) * ny / 2 + nx / 2], borders);




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
    // ЗДЕСЬ ПОЛУЧАЕТСЯ bot00
    InterpolateX(bot00, dx, dy, nx, ny, UrightE, UleftE, UtopE, UbotE, EvenRhs[(1 + nx / 2) * (ny / 2) + (nx / 4)], Eborder);
    // ДАЛЬШЕ РАСКЛАДЫВАЕТСЯ ODD НА 10 И 01
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
    InterpolateX(bot10, dx, dy, nx, ny, UrightO, UleftO, UtopO, UbotO, OddRhs[(1 + nx / 2) * ny / 2 + nx / 4], Oborder);

    //to 00:
    MPI_Send(F00, (1 + nx / 2)* (1 + ny / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD);
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
void Gauss(int rank, int& nx, int& ny, double& dx, double*& rhs, double*& Utop, double*& Ubot, double*& Uleft, double*& Uright, const int* kids, const int parent, double*& data) {
    MPI_Status status;
    MPI_Recv(rhs, (ny + 1) * (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(Utop, (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(Ubot, (nx + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(Uright, (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(Uleft, (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, &status);
    double* b = new double[4];
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
    double starttime = MPI_Wtime();
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
    int kids[4] = {4 * rank + 1, 4 * rank + 2, 4 * rank + 3, 4 * rank + 4};
    int parent = StageStart(h - 1) + (rank - StageStart(h)) / 4;
    if (rank > 0 && rank <= 4) parent = 0;
    double* rhs = new double [(ny + 1) * (nx + 1)];
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
        /*for (int i = 0; i < ny + 1; ++i) {
            for (int j = 0; j < nx + 1; ++j) rhs[i * (nx + 1) + j] = 4.0;
        }*/
        Shape shape(nx, ny);
        shape.PrintRhs();
        shape.Copy(rhs);
        for (int i = 0; i < ny + 1; ++i) Uleft[i] = 2.25 + (0.5 + (y2 - static_cast<double>(i * dy))) * (0.5 + (y2 - static_cast<double>(i * dy)));
        for (int i = 0; i < ny + 1; ++i) Uright[i] = 0.25 + (0.5 + (y2 - static_cast<double>(i * dy))) * (0.5 + (y2 - static_cast<double>(i * dy)));
        for (int i = 0; i < nx + 1; ++i) Ubot[i] = 0.25 + (-0.5 + (x1 + static_cast<double>(i * dx))) * (-0.5 + (x1 + static_cast<double>(i * dx)));
        for (int i = 0; i < nx + 1; ++i) Utop[i] = 2.25 + (-0.5 + (x1 + static_cast<double>(i * dx))) * (-0.5 + (x1 + static_cast<double>(i * dx)));
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
    for (int st = hight - 1; st >= 1; --st) {
            if (h == st) {
                MPI_Status status;
                MPI_Recv(getres1, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[0], 0, MPI_COMM_WORLD, &status);
                MPI_Recv(getres2, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[1], 0, MPI_COMM_WORLD, &status);
                MPI_Recv(getres3, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[2], 0, MPI_COMM_WORLD, &status);
                MPI_Recv(getres4, (1 + ny / 2) * (1 + nx / 2), MPI_DOUBLE, kids[3], 0, MPI_COMM_WORLD, &status);
                for (int i = 0; i < ny / 2 + 1; ++i) {
                    for (int j = 0; j < nx / 2 + 1; ++j) {
                        data[i * (1 + nx) + j] += getres1[i * (1 + nx / 2) + j];
                        data[i * (1 + nx) + j] -= getres3[((ny / 2) - i) * (1 + nx / 2) + j];
                        data[(i + ny / 2) * (1 + nx) + j] += getres1[((ny / 2) - i) * (1 + nx / 2) + j];
                        data[(i + ny / 2) * (1 + nx) + j] += getres3[i * (1 + nx / 2) + j];
                    }
                }
                
                for (int i = 0; i < ny / 2 + 1; ++i) {
                    for (int j = 0; j < nx / 2 + 1; ++j) {
                        data[i * (1 + nx) + nx / 2 + j] += getres2[i * (1 + nx / 2) + j];
                        data[i * (1 + nx) + nx / 2 + j] -= getres4[((ny / 2) - i) * (1 + nx / 2) + j];
                        data[(i + ny / 2) * (1 + nx) + nx / 2 + j] += getres2[(ny / 2 - i) * (1 + nx / 2) + j];
                        data[(i + ny / 2) * (1 + nx) + nx / 2 + j] += getres4[i * (1 + nx / 2) + j];
                    }
                }
                
                for (int k = 0; k < ny + 1; ++k) {
                    for (int q = 0; q < 1 + nx / 2; ++q) {
                        tmp_odd[k * (1 + nx / 2) + q] = - data[k * (1 + nx) + nx / 2 - q];
                    }
                }
                
                for (int k = 0; k < ny + 1; ++k) {
                    for (int q = 0; q < 1 + nx / 2; ++q) {
                        data[k * (1 + nx) + q] += data[k * (1 + nx) + nx - q];
                    }
                }
                
                for (int k = 0; k < ny + 1; ++k) {
                    for (int q = 0; q < 1 + nx / 2; ++q) {
                        data[k * (1 + nx) + q + nx / 2] += tmp_odd[k * (1 + nx / 2) + q];
                    }
                }
                for (int k = 0; k < nx + 1; ++k) data[(ny / 2) * (nx + 1) + k] /= 2.0;
                if (h != 1) {
                    MPI_Send(data, (nx + 1) * (ny + 1), MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
                }
            }
    }
    double endtime = MPI_Wtime();
    if (rank == 0) std::cout << "time : " << endtime - starttime << std::endl;
    if (h == 1) {
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

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
