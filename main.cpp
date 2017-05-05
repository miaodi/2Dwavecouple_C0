#define EIGEN_NO_DEBUG

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iomanip>
#include "spline.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/KroneckerProduct>


using namespace std;
using namespace spline;
using namespace Eigen;
const double pi = 3.141592653589793238462643383;

typedef complex<double> comp;

void GenerateKnot(int order, int refine, int insert, double *insert_knot,
                  double *&knot, int &m);

void CombineKnot(double *knot1, int m1, double *knot2, int m2, double *&knot,
                 int &m);

void Geometry(double xi, double eta, double &pxpxi, double &pxpeta,
              double &pypxi, double &pypeta);

void CompToPhy_patch1(double xi, double eta, double &x, double &y);

void CompToPhy_patch2(double xi, double eta, double &x, double &y);

MatrixXd weightFunction(int p, int m, double *knot) {
    MatrixXd C = MatrixXd::Zero((m - 2 * p), m - p);
    double gaussian[] = {-sqrt(3.0 / 5), 0, sqrt(3.0 / 5)};
    double weight[] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    VectorXd area = VectorXd::Zero(m - p);

    for (int j = 0; j < m - p; j++)
        for (int i = 0; i < m - 2 * p; i++) {
            double J = (knot[i + p + 1] - knot[i + p]) / 2;
            double Middle = (knot[i + p + 1] + knot[i + p]) / 2;
            for (int k = 0; k < 3; k++) {
                C(i, j) += weight[k]
                           * pow(DersOneBasisFUn(p, m, knot, j,
                                                 J * gaussian[k] + Middle, 0), 1) * J;
            }
            area(j) += C(i, j);
        }

    for (int j = 0; j < m - p; j++)
        for (int i = 0; i < m - 2 * p; i++)
            C(i, j) /= area(j);

    MatrixXd MM = MatrixXd::Zero((m - 2 * p), m - p);
    for (int j = 0; j < m - p; j++) {
        double test = 0;
        int location = 0;
        for (int i = 0; i < m - 2 * p; i++) {
            if (test <= C(i, j)) {
                test = C(i, j);
                location = i;
            }
        }
        MM(location, j) = 1;
    }

    MatrixXd CC = MatrixXd::Zero(m - p, (p + 1) * (m - 2 * p));
    for (int j = 0; j < m - p; j++)
        for (int i = 0; i < m - 2 * p; i++) {
            CC(j, p * i + j) = MM(i, j);
        }
    return CC;
}

MatrixXd localProjection(double *knot, int p, int m) {
    vector<double> knots;
    MatrixXd CC = MatrixXd::Identity(m - p, m - p);
    for (int i = 0; i < m + 1; i++)
        knots.push_back(knot[i]);
    for (int i = 0; i < m + 1 - 2 * (p + 1); i++) {
        int k = (i + 1) * (p + 1);
        for (int j = 0; j < p; j++) {
            MatrixXd C = MatrixXd::Zero(m - p + i * p + j + 1,
                                        m - p + i * p + j);
            for (int ii = 0; ii < k - p + 1; ii++)
                C(ii, ii) = 1;
            for (int ii = k - p + 1; ii <= k; ii++) {
                C(ii, ii) = (knot[i + p + 1] - knots[ii])
                            / (knots[ii + p] - knots[ii]);
                C(ii, ii - 1) = 1
                                - (knot[i + p + 1] - knots[ii])
                                  / (knots[ii + p] - knots[ii]);
            }
            for (int ii = k + 1; ii < m - p + i * p + j + 1; ii++)
                C(ii, ii - 1) = 1;
            CC *= C.transpose();
            knots.insert(knots.begin() + i * (p + 1) + p + 1, knot[i + p + 1]);
        }
    }
    MatrixXd B = MatrixXd::Zero((p + 1) * (m - 2 * p), (p + 1) * (m - 2 * p));
    for (int i = 0; i < m - 2 * p; i++)
        B.block(i * (p + 1), i * (p + 1), p + 1, p + 1) = CC.block(i,
                                                                   i * (p + 1), p + 1, p + 1).transpose();
    return B;
}

void Bezier_Knot(double order, double *knot, int m, double *&knot_projection,
                 int &m_projection) {
    vector<double> Knot;
    for (int i = 0; i < m; i++) {
        if (knot[i] != knot[i + 1]) {
            for (int j = 0; j <= order; j++)
                Knot.push_back(knot[i]);
        }
    }
    for (int j = 0; j <= order; j++)
        Knot.push_back(1);
    m_projection = Knot.end() - Knot.begin() - 1;
    knot_projection = new double[m_projection + 1];
    for (int i = 0; i < m_projection + 1; i++)
        knot_projection[i] = Knot[i];
}

comp Analytical(double x, double y, double omega, double theta) {
    comp result(cos(omega * x * cos(theta) + omega * y * sin(theta)),
                sin(omega * x * cos(theta) + omega * y * sin(theta)));
    return result;
}

int main() {
    double theta = pi/4;
    double omega = 10;
    double *gaussian = x7;
    double *weight = w7;
    int gaussian_points = 7;
    int order;
    int refine;
    cin >> order >> refine;
    int m_x_patch1, m_y_patch1;
    int m_x_patch2, m_y_patch2, m_y_patch2_projection;
    int p_x = order, p_y = order;
    double *knots_x_patch1, *knots_y_patch1;
    double *knots_x_patch2, *knots_y_patch2;
    double *knots_y_patch2_projection;
    double insertion_patch1[] = {.5};
    double insertion_patch2[] = {1.0 / 3, 2.0 / 3};
    GenerateKnot(p_x, refine, 1, insertion_patch1, knots_x_patch1, m_x_patch1);
    GenerateKnot(p_y, refine, 1, insertion_patch1, knots_y_patch1, m_y_patch1);
    GenerateKnot(p_x, refine, 2, insertion_patch2, knots_x_patch2, m_x_patch2);
    GenerateKnot(p_y, refine, 2, insertion_patch2, knots_y_patch2, m_y_patch2);
    Bezier_Knot(p_y, knots_y_patch2, m_y_patch2, knots_y_patch2_projection,
                m_y_patch2_projection);
    double *knots_y_coupling;
    int m_y_coupling;
    CombineKnot(knots_y_patch1, m_y_patch1, knots_y_patch2, m_y_patch2,
                knots_y_coupling,
                m_y_coupling);
    int dof_x_patch1 = m_x_patch1 - p_x, dof_y_patch1 = m_y_patch1 - p_y,
            dof_x_patch2 = m_x_patch2 - p_x,
            dof_y_patch2 = m_y_patch2 - p_y,
            dof_y_patch2_projection = m_y_patch2_projection - p_y;
    int dof_patch1 = dof_x_patch1 * dof_y_patch1,
            dof_patch2 = dof_x_patch2 * dof_y_patch2;
    int elements_x_patch1 = m_x_patch1 - 2 * p_x,
            elements_y_patch1 = m_y_patch1 - 2 * p_y,
            elements_x_patch2 = m_x_patch2 - 2 * p_x,
            elements_y_patch2 = m_y_patch2 - 2 * p_y;
    int elements_y_coupling = m_y_coupling - 2 * p_y;

        MatrixXd M = MatrixXd::Zero(dof_y_patch2, dof_y_patch2),
                  N2N1 = MatrixXd::Zero(dof_y_patch2, dof_y_patch1);
        for (int ii_y = 0; ii_y < elements_y_coupling; ii_y++) {
            double J_y = (knots_y_coupling[ii_y + p_y + 1] - knots_y_coupling[ii_y + p_y]) /
                              2;
            double Middle_y = (knots_y_coupling[ii_y + p_y + 1] + knots_y_coupling[ii_y +
                                    p_y]) / 2;
            int i_y_patch1 = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
            int i_y_patch2 = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
            for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                double eta = Middle_y + J_y * gaussian[jj_y];
                double **ders_y_patch1, **ders_y_patch2;
                DersBasisFuns(i_y_patch1, eta, p_y, knots_y_patch1, 0, ders_y_patch1);
                DersBasisFuns(i_y_patch2, eta, p_y, knots_y_patch2, 0, ders_y_patch2);
                VectorXd Neta_patch1 = VectorXd::Zero(dof_y_patch1),
                          Neta_patch2 = VectorXd::Zero(dof_y_patch2);
                for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                    Neta_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[0][kk_y];
                    Neta_patch2(i_y_patch2 - p_y + kk_y) = ders_y_patch2[0][kk_y];
                }
                for (int k = 0; k < 1; k++)
                    delete ders_y_patch1[k];
                delete[] ders_y_patch1;
                for (int k = 0; k < 1; k++)
                    delete ders_y_patch2[k];
                delete[] ders_y_patch2;
                M += weight[jj_y] * Neta_patch2 * Neta_patch2.transpose() * J_y;;
                N2N1 += weight[jj_y] * Neta_patch2 * Neta_patch1.transpose() * J_y;
            }
        }
        N2N1.leftCols(1) -= M.leftCols(1);
        N2N1.rightCols(1) -= M.rightCols(1);
        MatrixXd MN2N1 = M.block(1, 1, dof_y_patch2 - 2,
                                  dof_y_patch2 - 2).fullPivHouseholderQr().solve(N2N1.block(1, 0,
                                          dof_y_patch2 - 2,
                                          dof_y_patch1));
        MatrixXd MN2N1_new = MatrixXd::Zero(dof_y_patch2, dof_y_patch1);
        MN2N1_new(0, 0) = 1;
        MN2N1_new(dof_y_patch2 - 1, dof_y_patch1 - 1) = 1;
        MN2N1_new.block(1, 0, dof_y_patch2 - 2, dof_y_patch1) = MN2N1;


    /*
    MatrixXd M = MatrixXd::Zero(dof_y_patch2_projection, dof_y_patch2_projection),
            N2N1 = MatrixXd::Zero(dof_y_patch2_projection, dof_y_patch1);
    for (int ii_y = 0; ii_y < elements_y_coupling; ii_y++) {
        double J_y = (knots_y_coupling[ii_y + p_y + 1] - knots_y_coupling[ii_y + p_y]) /
                     2;
        double Middle_y = (knots_y_coupling[ii_y + p_y + 1] + knots_y_coupling[ii_y +
                                                                               p_y]) / 2;
        int i_y_patch1 = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
        int i_y_patch2_projection = Findspan(m_y_patch2_projection, p_y,
                                             knots_y_patch2_projection, Middle_y);
        for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
            double eta = Middle_y + J_y * gaussian[jj_y];
            double **ders_y_patch1, **ders_y_patch2_projection;
            DersBasisFuns(i_y_patch1, eta, p_y, knots_y_patch1, 0, ders_y_patch1);
            DersBasisFuns(i_y_patch2_projection, eta, p_y, knots_y_patch2_projection, 0,
                          ders_y_patch2_projection);
            VectorXd Neta_patch1 = VectorXd::Zero(dof_y_patch1),
                    Neta_patch2 = VectorXd::Zero(dof_y_patch2_projection);
            for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                Neta_patch1(i_y_patch1 - p_y + kk_y) = ders_y_patch1[0][kk_y];
                Neta_patch2(i_y_patch2_projection - p_y + kk_y) =
                        ders_y_patch2_projection[0][kk_y];
            }
            for (int k = 0; k < 1; k++)
                delete ders_y_patch1[k];
            delete[] ders_y_patch1;
            for (int k = 0; k < 1; k++)
                delete ders_y_patch2_projection[k];
            delete[] ders_y_patch2_projection;
            M += weight[jj_y] * Neta_patch2 * Neta_patch2.transpose() * J_y;
            N2N1 += weight[jj_y] * Neta_patch2 * Neta_patch1.transpose() * J_y;
        }
    }
    MatrixXd MN2N1 = M.fullPivHouseholderQr().solve(N2N1);
    MatrixXd weightness = weightFunction(p_y, m_y_patch2, knots_y_patch2);
    MatrixXd CC = localProjection(knots_y_patch2, p_y, m_y_patch2);
    MatrixXd MN2N1_projection = weightness * CC.fullPivHouseholderQr().solve(MN2N1);
     */
    MatrixXd couplingMatrix = MatrixXd::Zero(dof_patch2,
                                             (dof_patch2 - dof_y_patch2 + dof_y_patch1));

    couplingMatrix.block(0, 0, dof_y_patch2, dof_y_patch1) = MN2N1_new;
    couplingMatrix.block(dof_y_patch2, dof_y_patch1,
                         (dof_patch2 - dof_y_patch2),
                         (dof_patch2 - dof_y_patch2)) = MatrixXd::Identity(
            (dof_patch2 - dof_y_patch2), (dof_patch2 - dof_y_patch2));
    MatrixXd K_patch1 = MatrixXd::Zero(dof_patch1, dof_patch1);
    MatrixXd M_patch1 = MatrixXd::Zero(dof_patch1, dof_patch1);
    for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
        double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
        for (int ii_y = 0; ii_y < elements_y_patch1; ii_y++) {
            double J_y = (knots_y_patch1[ii_y + p_y + 1] - knots_y_patch1[ii_y + p_y]) / 2;
            double Middle_y = (knots_y_patch1[ii_y + p_y + 1] + knots_y_patch1[ii_y +
                                                                               p_y]) / 2;
            int i_y = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 1, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y_patch1, 1, ders_y);
                    VectorXd Nxi(p_x + 1), Nxi_xi(p_x + 1), Neta(p_y + 1), Neta_eta(p_y + 1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi(kk_x) = ders_x[0][kk_x];
                        Nxi_xi(kk_x) = ders_x[1][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta(kk_y) = ders_y[0][kk_y];
                        Neta_eta(kk_y) = ders_y[1][kk_y];
                    }
                    for (int k = 0; k < 2; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 2; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    VectorXd Nxi_xiNeta, NxiNeta_eta, NxiNeta;
                    Nxi_xiNeta = kroneckerProduct(Nxi_xi, Neta);
                    NxiNeta_eta = kroneckerProduct(Nxi, Neta_eta);
                    NxiNeta = kroneckerProduct(Nxi, Neta);
                    double pxpxi, pxpeta, pypxi, pypeta;
                    Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    VectorXd Nx_xNy, NxNy_y;
                    Nx_xNy = 1.0 / Jacobian * (Nxi_xiNeta * pypeta - NxiNeta_eta * pypxi);
                    NxNy_y = 1.0 / Jacobian * (-Nxi_xiNeta * pxpeta + NxiNeta_eta * pxpxi);
                    for (int kkx = 0; kkx < (p_x + 1) * (p_y + 1); kkx++) {
                        for (int kky = 0; kky < (p_x + 1) * (p_y + 1); kky++) {
                            MatrixXd Bx(2, 1);
                            MatrixXd By(2, 1);
                            Bx << Nx_xNy(kkx), NxNy_y(kkx);
                            By << Nx_xNy(kky), NxNy_y(kky);
                            K_patch1.block(((m_y_patch1 - p_y) * (kkx / (p_y + 1) + i_x - p_x) + kkx %
                                                                                                 (p_y + 1) + i_y - p_y),
                                           ((m_y_patch1 - p_y) * (kky /
                                                                  (p_y + 1) + i_x - p_x) + kky
                                                                                           % (p_y + 1) + i_y - p_y), 1,
                                           1) += weight[jj_x] * weight[jj_y] * Jacobian * Bx.transpose() * By * J_x *
                                                 J_y;
                            M_patch1(((m_y_patch1 - p_y) * (kkx / (p_y + 1) + i_x - p_x) + kkx %
                                                                                           (p_y + 1) + i_y - p_y),
                                     ((m_y_patch1 - p_y) * (kky /
                                                            (p_y + 1) + i_x - p_x) + kky
                                                                                     % (p_y + 1) + i_y - p_y)) +=
                                    weight[jj_x] * weight[jj_y] * Jacobian * NxiNeta(kkx) * NxiNeta(kky) * J_x *
                                    J_y;
                        }
                    }
                }
            }
        }
    }
    K_patch1 -= pow(omega, 2) * M_patch1;
    MatrixXcd K_compl_patch1 = MatrixXcd::Zero(dof_patch1, dof_patch1);
    K_compl_patch1.real() = K_patch1;
    K_patch1.resize(1, 1);
    M_patch1.resize(1, 1);
    MatrixXd K_patch2 = MatrixXd::Zero(dof_patch2, dof_patch2);
    MatrixXd M_patch2 = MatrixXd::Zero(dof_patch2, dof_patch2);
    for (int ii_x = 0; ii_x < elements_x_patch2; ii_x++) {
        double J_x = (knots_x_patch2[ii_x + p_x + 1] - knots_x_patch2[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch2[ii_x + p_x + 1] + knots_x_patch2[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch2, p_x, knots_x_patch2, Middle_x);
        for (int ii_y = 0; ii_y < elements_y_patch2; ii_y++) {
            double J_y = (knots_y_patch2[ii_y + p_y + 1] - knots_y_patch2[ii_y + p_y]) / 2;
            double Middle_y = (knots_y_patch2[ii_y + p_y + 1] + knots_y_patch2[ii_y +
                                                                               p_y]) / 2;
            int i_y = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x_patch2, 1, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y_patch2, 1, ders_y);
                    VectorXd Nxi(p_x + 1), Nxi_xi(p_x + 1), Neta(p_y + 1), Neta_eta(p_y + 1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi(kk_x) = ders_x[0][kk_x];
                        Nxi_xi(kk_x) = ders_x[1][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta(kk_y) = ders_y[0][kk_y];
                        Neta_eta(kk_y) = ders_y[1][kk_y];
                    }
                    for (int k = 0; k < 2; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 2; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    VectorXd Nxi_xiNeta, NxiNeta_eta, NxiNeta;
                    Nxi_xiNeta = kroneckerProduct(Nxi_xi, Neta);
                    NxiNeta_eta = kroneckerProduct(Nxi, Neta_eta);
                    NxiNeta = kroneckerProduct(Nxi, Neta);
                    double pxpxi, pxpeta, pypxi, pypeta;
                    Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    VectorXd Nx_xNy, NxNy_y;
                    Nx_xNy = 1.0 / Jacobian * (Nxi_xiNeta * pypeta - NxiNeta_eta * pypxi);
                    NxNy_y = 1.0 / Jacobian * (-Nxi_xiNeta * pxpeta + NxiNeta_eta * pxpxi);
                    for (int kkx = 0; kkx < (p_x + 1) * (p_y + 1); kkx++) {
                        for (int kky = 0; kky < (p_x + 1) * (p_y + 1); kky++) {
                            MatrixXd Bx(2, 1);
                            MatrixXd By(2, 1);
                            Bx << Nx_xNy(kkx), NxNy_y(kkx);
                            By << Nx_xNy(kky), NxNy_y(kky);
                            K_patch2.block(((m_y_patch2 - p_y) * (kkx / (p_y + 1) + i_x - p_x) + kkx %
                                                                                                 (p_y + 1) + i_y - p_y),
                                           ((m_y_patch2 - p_y) * (kky /
                                                                  (p_y + 1) + i_x - p_x) + kky
                                                                                           % (p_y + 1) + i_y - p_y), 1,
                                           1) += weight[jj_x] * weight[jj_y] * Jacobian * Bx.transpose() * By * J_x *
                                                 J_y;
                            M_patch2(((m_y_patch2 - p_y) * (kkx / (p_y + 1) + i_x - p_x) + kkx %
                                                                                           (p_y + 1) + i_y - p_y),
                                     ((m_y_patch2 - p_y) * (kky /
                                                            (p_y + 1) + i_x - p_x) + kky
                                                                                     % (p_y + 1) + i_y - p_y)) +=
                                    weight[jj_x] * weight[jj_y] * Jacobian * NxiNeta(kkx) * NxiNeta(kky) * J_x *
                                    J_y;
                        }
                    }
                }
            }
        }
    }
    K_patch2 -= pow(omega, 2) * M_patch2;
    MatrixXcd K_compl_patch2 = MatrixXcd::Zero(dof_patch2, dof_patch2);
    K_compl_patch2.real() = K_patch2;
    K_patch2.resize(1, 1);
    M_patch2.resize(1, 1);
    MatrixXcd K = MatrixXcd::Zero((dof_patch1 + dof_patch2 -
                                   dof_y_patch2), (dof_patch1 + dof_patch2 - dof_y_patch2));
    K.block(0, 0, dof_patch1, dof_patch1) = K_compl_patch1;
    K.block((dof_patch1 - dof_y_patch1), (dof_patch1 - dof_y_patch1),
            (dof_patch2 - dof_y_patch2 + dof_y_patch1),
            (dof_patch2 - dof_y_patch2 + dof_y_patch1)) += couplingMatrix.transpose()
                                                           * K_compl_patch2 * couplingMatrix;

    MatrixXcd x_basis_patch1 = MatrixXcd::Zero(dof_x_patch1, dof_x_patch1);
    VectorXcd U_x_bottom_rhs_patch1 = VectorXcd::Zero(dof_x_patch1);
    for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
        double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
        for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
            double eta = 0;
            double xi = Middle_x + J_x * gaussian[jj_x];
            double **ders_x;
            DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 0, ders_x);
            VectorXd Nxi = VectorXd::Zero(dof_x_patch1);
            for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                Nxi(i_x - p_x + kk_x) = ders_x[0][kk_x];
            }
            for (int k = 0; k < 1; k++)
                delete ders_x[k];
            delete[] ders_x;
            double pxpxi, pxpeta, pypxi, pypeta;
            Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
            double x, y;
            CompToPhy_patch1(xi, eta, x, y);
            auto u_analytical = Analytical(x, y, omega, theta);
            x_basis_patch1.real() += weight[jj_x] * Nxi * MatrixXd(Nxi.transpose()) * J_x * pxpxi;
            U_x_bottom_rhs_patch1 += weight[jj_x] * Nxi * u_analytical * J_x * pxpxi;
        }
    }
    VectorXcd U_x_bottom_dirichlet_patch1 = x_basis_patch1.fullPivHouseholderQr().solve(U_x_bottom_rhs_patch1);
    VectorXcd U_x_top_rhs_patch1 = VectorXd::Zero(dof_x_patch1);
    for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
        double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
        for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
            double eta = .999999999999;
            double xi = Middle_x + J_x * gaussian[jj_x];
            double **ders_x;
            DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 0, ders_x);
            VectorXd Nxi = VectorXd::Zero(dof_x_patch1);
            for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                Nxi(i_x - p_x + kk_x) = ders_x[0][kk_x];
            }
            for (int k = 0; k < 1; k++)
                delete ders_x[k];
            delete[] ders_x;
            double pxpxi, pxpeta, pypxi, pypeta;
            Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
            double x, y;
            CompToPhy_patch1(xi, eta, x, y);
            auto u_analytical = Analytical(x, y, omega, theta);
            U_x_top_rhs_patch1 += weight[jj_x] * Nxi * u_analytical * J_x * pxpxi;
        }
    }
    VectorXcd U_x_top_dirichlet_patch1 = x_basis_patch1.fullPivHouseholderQr().solve(U_x_top_rhs_patch1);

    MatrixXcd y_basis_patch1 = MatrixXcd::Zero(dof_y_patch1, dof_y_patch1);
    VectorXcd U_y_rhs_patch1 = VectorXcd::Zero(dof_y_patch1);
    for (int ii_y = 0; ii_y < elements_y_patch1; ii_y++) {
        double J_y = (knots_y_patch1[ii_y + p_y + 1] - knots_y_patch1[ii_y + p_y]) / 2;
        double Middle_y = (knots_y_patch1[ii_y + p_y + 1] + knots_y_patch1[ii_y + p_y]) / 2;
        int i_y = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
        for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
            double eta = Middle_y + J_y * gaussian[jj_y];
            double xi = 0;
            double **ders_y;
            DersBasisFuns(i_y, eta, p_y, knots_y_patch1, 0, ders_y);
            VectorXd Neta = VectorXd::Zero(dof_y_patch1);
            for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                Neta(i_y - p_y + kk_y) = ders_y[0][kk_y];
            }
            for (int k = 0; k < 1; k++)
                delete ders_y[k];
            delete[] ders_y;
            double pxpxi, pxpeta, pypxi, pypeta;
            Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
            double x, y;
            CompToPhy_patch1(xi, eta, x, y);
            auto u_analytical = Analytical(x, y, omega, theta);
            y_basis_patch1 += weight[jj_y] * Neta * Neta.transpose() * J_y * pypeta;
            U_y_rhs_patch1 += weight[jj_y] * Neta * u_analytical * J_y * pypeta;
        }
    }
    VectorXcd U_y_dirichlet_patch1 = y_basis_patch1.fullPivHouseholderQr().solve(U_y_rhs_patch1);
    cout << U_y_dirichlet_patch1 << endl;

    MatrixXcd x_basis_patch2 = MatrixXcd::Zero(dof_x_patch2, dof_x_patch2);
    VectorXcd U_x_bottom_rhs_patch2 = VectorXcd::Zero(dof_x_patch2);
    for (int ii_x = 0; ii_x < elements_x_patch2; ii_x++) {
        double J_x = (knots_x_patch2[ii_x + p_x + 1] - knots_x_patch2[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch2[ii_x + p_x + 1] + knots_x_patch2[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch2, p_x, knots_x_patch2, Middle_x);
        for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
            double eta = 0;
            double xi = Middle_x + J_x * gaussian[jj_x];
            double **ders_x;
            DersBasisFuns(i_x, xi, p_x, knots_x_patch2, 0, ders_x);
            VectorXd Nxi = VectorXd::Zero(dof_x_patch2);
            for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                Nxi(i_x - p_x + kk_x) = ders_x[0][kk_x];
            }
            for (int k = 0; k < 1; k++)
                delete ders_x[k];
            delete[] ders_x;
            double pxpxi, pxpeta, pypxi, pypeta;
            Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
            double x, y;
            CompToPhy_patch2(xi, eta, x, y);
            auto u_analytical = Analytical(x, y, omega, theta);
            x_basis_patch2.real() += weight[jj_x] * Nxi * MatrixXd(Nxi.transpose()) * J_x * pxpxi;
            U_x_bottom_rhs_patch2 += weight[jj_x] * Nxi * u_analytical * J_x * pxpxi;
        }
    }
    VectorXcd U_x_bottom_dirichlet_patch2 = x_basis_patch2.fullPivHouseholderQr().solve(U_x_bottom_rhs_patch2);
    VectorXcd U_x_top_rhs_patch2 = VectorXcd::Zero(dof_x_patch2);
    for (int ii_x = 0; ii_x < elements_x_patch2; ii_x++) {
        double J_x = (knots_x_patch2[ii_x + p_x + 1] - knots_x_patch2[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch2[ii_x + p_x + 1] + knots_x_patch2[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch2, p_x, knots_x_patch2, Middle_x);
        for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
            double eta = .9999999999999;
            double xi = Middle_x + J_x * gaussian[jj_x];
            double **ders_x;
            DersBasisFuns(i_x, xi, p_x, knots_x_patch2, 0, ders_x);
            VectorXd Nxi = VectorXd::Zero(dof_x_patch2);
            for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                Nxi(i_x - p_x + kk_x) = ders_x[0][kk_x];
            }
            for (int k = 0; k < 1; k++)
                delete ders_x[k];
            delete[] ders_x;
            double pxpxi, pxpeta, pypxi, pypeta;
            Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
            double x, y;
            CompToPhy_patch2(xi, eta, x, y);
            auto u_analytical = Analytical(x, y, omega, theta);
            U_x_top_rhs_patch2 += weight[jj_x] * Nxi * u_analytical * J_x * pxpxi;
        }
    }
    VectorXcd U_x_top_dirichlet_patch2 = x_basis_patch2.fullPivHouseholderQr().solve(U_x_top_rhs_patch2);
    cout << U_x_top_dirichlet_patch2 << endl;

    MatrixXcd y_basis_patch2 = MatrixXcd::Zero(dof_y_patch2, dof_y_patch2);
    VectorXcd U_y_rhs_patch2 = VectorXcd::Zero(dof_y_patch2);
    for (int ii_y = 0; ii_y < elements_y_patch2; ii_y++) {
        double J_y = (knots_y_patch2[ii_y + p_y + 1] - knots_y_patch2[ii_y + p_y]) / 2;
        double Middle_y = (knots_y_patch2[ii_y + p_y + 1] + knots_y_patch2[ii_y + p_y]) / 2;
        int i_y = Findspan(m_y_patch2, p_y, knots_y_patch2, Middle_y);
        for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
            double eta = Middle_y + J_y * gaussian[jj_y];
            double xi = .99999999999999;
            double **ders_y;
            DersBasisFuns(i_y, eta, p_y, knots_y_patch2, 0, ders_y);
            VectorXd Neta = VectorXd::Zero(dof_y_patch2);
            for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                Neta(i_y - p_y + kk_y) = ders_y[0][kk_y];
            }
            for (int k = 0; k < 1; k++)
                delete ders_y[k];
            delete[] ders_y;
            double pxpxi, pxpeta, pypxi, pypeta;
            Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
            double x, y;
            CompToPhy_patch2(xi, eta, x, y);
            auto u_analytical = Analytical(x, y, omega, theta);
            y_basis_patch2.real() += weight[jj_y] * Neta * Neta.transpose() * J_y * pypeta;
            U_y_rhs_patch2 += weight[jj_y] * Neta * u_analytical * J_y * pypeta;
        }
    }
    VectorXcd U_y_dirichlet_patch2 = y_basis_patch2.fullPivHouseholderQr().solve(U_y_rhs_patch2);
    cout << U_y_dirichlet_patch2 << endl;

    VectorXcd Dirichlet = VectorXcd::Zero((dof_patch1 + dof_patch2 - dof_y_patch2));
    for (int i = 0; i < dof_y_patch1; i++)
        Dirichlet[i] = U_y_dirichlet_patch1[i];
    for (int i = 0; i < dof_x_patch1; i++)
        Dirichlet[dof_y_patch1 * (i + 1) - 1] = U_x_top_dirichlet_patch1[i];
    for (int i = 0; i < dof_x_patch1; i++)
        Dirichlet[dof_y_patch1 * i] = U_x_bottom_dirichlet_patch1[i];
    for (int i = 1; i < dof_x_patch2; i++)
        Dirichlet[dof_patch1 + dof_y_patch2 * i - 1] = U_x_top_dirichlet_patch2[i];
    for (int i = 1; i < dof_x_patch2; i++)
        Dirichlet[dof_patch1 + dof_y_patch2 * (i - 1)] = U_x_bottom_dirichlet_patch2[i];
    for (int i = 0; i < dof_y_patch2; i++)
        Dirichlet[dof_patch1 + dof_patch2 - 2 * dof_y_patch2 + i] = U_y_dirichlet_patch2[i];

    VectorXcd F = -K * Dirichlet;
    MatrixXd transform = MatrixXd::Zero(dof_patch1 + dof_patch2 -
                                        dof_y_patch2 - dof_y_patch1 - dof_y_patch2 - 2 * (dof_x_patch1 + dof_x_patch2) +
                                        6, dof_patch1 + dof_patch2 -
                                           dof_y_patch2);

    int x_it = 0, y_it = 0;
    for (int i = 0; i < dof_x_patch1; i++) {
        for (int j = 0; j < dof_y_patch1; j++) {
            if ((i == 0) || (j == 0) || (j == dof_y_patch1 - 1)) {
                y_it++;
            } else {
                transform(x_it, y_it) = 1;
                x_it++;
                y_it++;
            }
        }
    }
    for (int i = 0; i < dof_x_patch2 - 1; i++) {
        for (int j = 0; j < dof_y_patch2; j++) {
            if ((i == dof_x_patch2 - 2) || (j == 0) || (j == dof_y_patch2 - 1)) {
                y_it++;
            } else {
                transform(x_it, y_it) = 1;
                x_it++;
                y_it++;
            }
        }
    }
    MatrixXcd K_cal = transform * K * transform.transpose();
    VectorXcd F_cal = transform * F;
    VectorXcd U_cal = K_cal.fullPivHouseholderQr().solve(F_cal);
    cout << U_cal;

    VectorXcd U_patch1 = VectorXcd::Zero(dof_patch1);
    VectorXcd U_patch2 = VectorXcd::Zero(dof_patch2);
    x_it = 0;
    for (int i = 0; i < dof_x_patch1; i++) {
        for (int j = 0; j < dof_y_patch1; j++) {
            if (j == 0) {
                U_patch1[j + i * dof_y_patch1] = Dirichlet[j + i * dof_y_patch1];
            } else if (i == 0) {
                U_patch1[j + i * dof_y_patch1] = Dirichlet[j + i * dof_y_patch1];
            } else if (j == dof_y_patch1 - 1) {
                U_patch1[j + i * dof_y_patch1] = Dirichlet[j + i * dof_y_patch1];
            } else {
                U_patch1[j + i * dof_y_patch1] = U_cal[x_it];
                x_it++;
            }
        }
    }
    double L2_norm = 0, L2_norm_error = 0;
    for (int ii_x = 0; ii_x < elements_x_patch1; ii_x++) {
        double J_x = (knots_x_patch1[ii_x + p_x + 1] - knots_x_patch1[ii_x + p_x]) / 2;
        double Middle_x = (knots_x_patch1[ii_x + p_x + 1] + knots_x_patch1[ii_x +
                                                                           p_x]) / 2;
        int i_x = Findspan(m_x_patch1, p_x, knots_x_patch1, Middle_x);
        for (int ii_y = 0; ii_y < elements_y_patch1; ii_y++) {
            double J_y = (knots_y_patch1[ii_y + p_y + 1] - knots_y_patch1[ii_y + p_y]) / 2;
            double Middle_y = (knots_y_patch1[ii_y + p_y + 1] + knots_y_patch1[ii_y +
                                                                               p_y]) / 2;
            int i_y = Findspan(m_y_patch1, p_y, knots_y_patch1, Middle_y);
            for (int jj_x = 0; jj_x < gaussian_points; jj_x++) {
                for (int jj_y = 0; jj_y < gaussian_points; jj_y++) {
                    double xi = Middle_x + J_x * gaussian[jj_x];
                    double eta = Middle_y + J_y * gaussian[jj_y];
                    double **ders_x, **ders_y;
                    DersBasisFuns(i_x, xi, p_x, knots_x_patch1, 1, ders_x);
                    DersBasisFuns(i_y, eta, p_y, knots_y_patch1, 1, ders_y);
                    VectorXd Nxi = VectorXd::Zero(dof_x_patch1),
                            Neta = VectorXd::Zero(dof_y_patch1);
                    for (int kk_x = 0; kk_x < p_x + 1; kk_x++) {
                        Nxi(i_x - p_x + kk_x) = ders_x[0][kk_x];
                    }
                    for (int kk_y = 0; kk_y < p_y + 1; kk_y++) {
                        Neta(i_y - p_y + kk_y) = ders_y[0][kk_y];
                    }
                    for (int k = 0; k < 2; k++)
                        delete ders_x[k];
                    delete[] ders_x;
                    for (int k = 0; k < 2; k++)
                        delete ders_y[k];
                    delete[] ders_y;
                    VectorXd NxiNeta;
                    NxiNeta = kroneckerProduct(Nxi, Neta);
                    double pxpxi, pxpeta, pypxi, pypeta;
                    Geometry(xi, eta, pxpxi, pxpeta, pypxi, pypeta);
                    double Jacobian = pxpxi * pypeta - pxpeta * pypxi;
                    MatrixXcd U = (NxiNeta.transpose()) * U_patch1;
                    double x, y;
                    CompToPhy_patch1(xi, eta, x, y);
                    auto u_analytical = Analytical(x, y, omega, theta);
                    L2_norm_error += weight[jj_x] * weight[jj_y] * (pow(abs(U(0, 0) - u_analytical), 2)) * J_x * J_y
                                     * Jacobian;
                    L2_norm += weight[jj_x] * weight[jj_y] * (pow(abs(u_analytical), 2)) * J_x * J_y * Jacobian;
                }
            }
        }
    }
    cout << sqrt(L2_norm_error) / sqrt(L2_norm) << endl;

    return 0;
}


void GenerateKnot(int order, int refine, int insert, double *insert_knot,
                  double *&knot, int &m) {
    vector<double> Knot;
    for (int i = 0; i <= order; i++) {
        Knot.push_back(0);
    }
    for (int i = 0; i < insert; i++) {
        Knot.push_back(insert_knot[i]);
    }
    for (int i = 0; i <= order; i++) {
        Knot.push_back(1);
    }
    for (int i = 0; i < refine; i++) {
        for (int j = 0; j < Knot.end() - Knot.begin() - 1; j++) {
            double insertion;
            if (Knot[j] != Knot[j + 1]) {
                insertion = (Knot[j] + Knot[j + 1]) / 2;
                Knot.insert(Knot.begin() + j + 1, insertion);
                j++;
            }
        }
    }
    m = Knot.end() - Knot.begin() - 1;
    knot = new double[m + 1];
    for (int i = 0; i < m + 1; i++)
        knot[i] = Knot[i];
}

void CombineKnot(double *knot1, int m1, double *knot2, int m2, double *&knot,
                 int &m) {
    vector<double> Knot;
    int i1 = 0, i2 = 0;
    while (i1 <= m1) {
        if (knot1[i1] == knot2[i2]) {
            Knot.push_back(knot1[i1]);
            i1++, i2++;
        } else if (knot1[i1] < knot2[i2]) {
            Knot.push_back(knot1[i1]);
            i1++;
        } else {
            Knot.push_back(knot2[i2]);
            i2++;
        }
    }

    m = Knot.end() - Knot.begin() - 1;
    knot = new double[m + 1];
    for (int i = 0; i < m + 1; i++)
        knot[i] = Knot[i];
}

void Geometry(double xi, double eta, double &pxpxi, double &pxpeta,
              double &pypxi, double &pypeta) {
    pxpxi = .5;
    pxpeta = 0;
    pypxi = 0;
    pypeta = 1;
}

void CompToPhy_patch1(double xi, double eta, double &x, double &y) {
    y = eta;
    x = .5 * xi;
}

void CompToPhy_patch2(double xi, double eta, double &x, double &y) {
    y = eta;
    x = .5 * xi + .5;
}
