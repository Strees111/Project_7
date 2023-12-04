#include <iostream>
#include <string>
#include <fstream>

using namespace std;

inline void printMatrix(double **a, const int &n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << *(*(a + i) + j) << ' ';
        }
        cout << endl;
    }
}

inline void newmatrix(double **&a, double **&b, double **&c, const int &n) {
    a = new double *[n];
    b = new double *[n];
    c = new double *[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n]{};
        b[i] = new double[n]{};
        c[i] = new double[n]{};
    }
}

inline void Del_ete(double **a, double **b, double **c, const int &n) {
    for (int i = 0; i < n; ++i) {
        delete[] a[i];
        delete[] b[i];
        delete[] c[i];
    }
    delete[] a;
    delete[] b;
    delete[] c;
}

inline void determinant(double **a, const int &n) {
    double d = 1;
    for (int k = 0; k < n - 1; ++k) {
        double c = a[k][k];
        if (c == 0) {
            for (int i = k + 1; i < n; ++i) {
                if (a[i][k] != 0) {
                    for (int j = 0; j < n; ++j) {
                        swap(a[k][j], a[i][j]);
                    }
                    d *= -1;
                    c = a[k][k];
                    break;
                }
            }
            if (c == 0) {
                cout << c << endl;
                return;
            }
        }
        for (int i = k + 1; i < n; ++i) {
            double tmp = -a[i][k] / c;
            for (int j = k; j < n; ++j) {
                a[i][j] += a[k][j] * tmp;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        d *= a[i][i];
    }
    cout << d << endl;
}

void multiply(double **result, double **M, double **K, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = 0.0;
            for (int k = 0; k < n; ++k) {
                result[i][j] += M[i][k] * K[k][j];
            }
        }
    }
}

void CharPolCoefs(double *coefs, const double *traces, int N) {
    coefs[1] = -traces[1];
    for (int i = 2; i <= N; ++i) {
        coefs[i] = -traces[i];
        for (int j = 1; j <= i - 1; ++j) {
            coefs[i] -= coefs[j] * traces[i - j];
        }
        coefs[i] /= i;
    }
}

void characteristicPolynomial(double **matrix, double *coefs, int N) {
    double **tmp = new double *[N];
    for (int i = 0; i < N; ++i) {
        tmp[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            tmp[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    double *trace = new double[N + 1]{1.0};

    double **prev = new double *[N];
    for (int i = 0; i < N; ++i) {
        prev[i] = new double[N]{};
    }
    // Вычисление следа матрицы в степени
    for (int i = 1; i <= N; ++i) {
        multiply(prev, tmp, matrix, N);
        double tr = 0.0;
        for (int j = 0; j < N; ++j) {
            tr += prev[j][j];
        }
        trace[i] = tr;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                tmp[j][k] = prev[j][k];
                prev[j][k] = 0.0;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        delete [] prev[i];
    }
    delete [] prev;
    CharPolCoefs(coefs, trace, N);

    delete[] trace;
    for (int i = 0; i < N; ++i) {
        delete [] tmp[i];
    }
    delete [] tmp;
}

inline void inverse_matrix(double **a, double **b, const int &n) {
    for (int k = 0; k < n - 1; ++k) {
        double c = a[k][k];
        if (c == 0) {
            for (int i = k + 1; i < n; ++i) {
                if (a[i][k] != 0) {
                    swap(a[k], a[i]);
                    swap(b[k], b[i]);
                    c = a[k][k];
                    break;
                }
            }
            if (c == 0) {
                throw "Matrix is singular. Cannot compute inverse.";
            }
        }
        for (int i = k + 1; i < n; ++i) {
            double tmp = -a[i][k] / c;
            for (int j = 0; j < n; ++j) {
                a[i][j] += a[k][j] * tmp;
                b[i][j] += b[k][j] * tmp;
            }
        }
    }
    for (int k = n - 1; k > 0; k--) {
        double c = a[k][k];
        if (c == 0) {
            for (int i = k + 1; i < n; ++i) {
                if (a[i][k] != 0) {
                    swap(a[k], a[i]);
                    swap(b[k], b[i]);
                    c = a[k][k];
                    break;
                }
            }
            if (c == 0) {
                throw "Matrix is singular. Cannot compute inverse.";
            }
        }
        for (int i = k - 1; i >= 0; i--) {
            double d = -a[i][k] / c;
            for (int j = 0; j < n; j++) {
                a[i][j] += a[k][j] * d;
                b[i][j] += b[k][j] * d;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        double c = a[i][i];
        if (c == 0) {
            throw "Matrix is singular. Cannot compute inverse.";
        }
        for (int j = 0; j < n; ++j) {
            a[i][j] /= c;
            b[i][j] /= c;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << b[i][j] << ' ';
        }
        cout << endl;
    }
}


int main() {
    system("chcp 65001 > nul");
    int q = 0;
    cout << "Введите 1, если хотите вывести определитель" << endl;
    cout << "Введите 2, если хотите вывести характеристический полином" << endl;
    cout << "Введите 3, если хотите вывести обратную матрицу" << endl;
    cin >> q;
    switch (q) {
        case 1:
            for (int u = 1; u <= 10; ++u) {
                string path = ".\\tests\\" + to_string(u) + "matrix.txt";
                ifstream fin(path);
                if (!fin.is_open()) {
                    cout << "Can not open to file " << path << endl;
                    return 1;
                }
                int rows = 0, cols = 0;
                fin >> rows >> cols;
                if (rows != cols) {
                    cout << u << ")" << " determinant" << endl << "Error" << endl << endl;
                    fin.close();
                    continue;
                }
                int n = rows;
                double **a, **b, **c;
                newmatrix(a, b, c, n);
                for (int i = 0; i < rows; ++i) {
                    b[i][i] = 1;
                    for (int j = 0; j < cols; ++j) {
                        fin >> a[i][j];
                    }
                }
                cout << u << ")" << " determinant" << endl;
                determinant(a, n);
                cout << endl;
                Del_ete(a, b, c, n);
                fin.close();
            }
            break;
        case 2:
            for (int u = 1; u <= 10; ++u) {
                string path = ".\\tests\\" + to_string(u) + "matrix.txt";
                ifstream fin(path);
                if (!fin.is_open()) {
                    cout << "Can not open to file " << path << endl;
                    return 1;
                }
                int rows = 0, cols = 0;
                fin >> rows >> cols;
                if (rows != cols) {
                    cout << u << ")" << " characteristic_polynomial" << endl << "Error" << endl << endl;
                    fin.close();
                    continue;
                }
                int n = rows;
                double **a, **b, **c;
                newmatrix(a, b, c, n);
                for (int i = 0; i < rows; ++i) {
                    b[i][i] = 1;
                    for (int j = 0; j < cols; ++j) {
                        fin >> a[i][j];
                        c[i][j] = a[i][j];
                    }
                }
                cout << u << ")" << " characteristic_polynomial" << endl;
                double *coefs = new double[cols + 1]{1.0};
                cout << 1 << ' ' << cols + 1 << std::endl;
                characteristicPolynomial(a, coefs, rows);
                for (int i = 0; i < cols + 1; ++i) {
                    cout << coefs[i] << ' ';
                }
                cout << endl;
                delete[] coefs;
                Del_ete(a, b, c, n);
                fin.close();
            }
            break;
        case 3:
            for (int u = 1; u <= 10; ++u) {
                string path = ".\\tests\\" + to_string(u) + "matrix.txt";
                ifstream fin(path);
                if (!fin.is_open()) {
                    cout << "Can not open to file " << path << endl;
                    return 1;
                }
                int rows = 0, cols = 0;
                fin >> rows >> cols;
                if (rows != cols) {
                    cout << u << ")" << " inverse_matrix" << endl << "Error" << endl << endl;
                    fin.close();
                    continue;
                }
                int n = rows;
                double **a, **b, **c;
                newmatrix(a, b, c, n);
                for (int i = 0; i < rows; ++i) {
                    b[i][i] = 1;
                    for (int j = 0; j < cols; ++j) {
                        fin >> a[i][j];
                    }
                }
                cout << u << ")" << " inverse_matrix" << endl;
                try {
                    inverse_matrix(a, b, n);
                    cout << endl;
                }
                catch (const char *) {
                    cout << "Singular" << endl << endl;
                    Del_ete(a, b, c, n);
                    fin.close();
                    continue;
                }
                Del_ete(a, b, c, n);
                fin.close();
            }
            break;
        default:
            cout << "Ты совсем?" << endl;
    }
    return 0;
}