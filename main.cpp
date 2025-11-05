#include <iostream>
#include "cmath"
#include <exception>
#include <chrono>
#include <vector>
// #include <omp.h>

class Rational{
protected:
    void eq_denom(Rational& fr, Rational& sc){
        if(fr.denom != sc.denom){
            int temp_d = fr.denom;
            fr.num *= sc.denom;
            fr.denom *= sc.denom;
            sc.num *= temp_d;
            sc.denom *= temp_d;
        }
    }

    void normalize(Rational& fr){
        if(fr.num == 0)
            fr.denom = 1;
        else if(fr.num % fr.denom == 0){
            fr.num /= fr.denom;
            fr.denom = 1;
        }
        else if(fr.denom & fr.num == 0){
            fr.denom /= fr.num;
            fr.num = 1;
        }
    }
public:
    int num;
    int denom;

    Rational(int n, int d){
        if(d == 0)
            throw "Denominator can't be zero\n";
        num = n;
        denom = d;
    }

    Rational(int n) : num{n}, denom{1} { };

    Rational() { };

    static Rational from_int(int i){
        return Rational{i, 1};
    }

    // operator int() const{
    //     if(denom != 1)
    //         throw "Error";

    //     return num;
    // }

    Rational operator + (Rational& other){
        eq_denom(*this, other);
        Rational answer(num + other.num, denom);
        normalize(answer);

        return answer;
    }

    Rational operator - (Rational& other){
        eq_denom(*this, other);
        Rational answer(num - other.num, denom);
        normalize(answer);

        return answer;
    }

    Rational operator * (const Rational& other){
        Rational answer(num * other.num, denom * other.denom);
        normalize(answer);

        return answer;
    }

    Rational operator / (const Rational& other){
        Rational answer(num * other.denom, denom * other.num);
        normalize(answer);

        return answer;
    }

    Rational& operator = (const Rational& other){
        if(other.num != num || other.denom != denom){
            num = other.num;
            denom = other.denom;
        }
        return *this;
    }

    Rational& operator = (const int& other){
        num = other;
        denom = 1;

        return *this;
    }

    Rational& operator *= (const Rational& other){
        num *= other.num;
        denom *= other.denom;
        normalize(*this);

        return *this;
    }

    Rational& operator += (Rational& other){
        eq_denom(*this, other);
        num += other.num;
        normalize(*this);
        
        return *this;
    }

    Rational& operator -= (Rational& other){
        eq_denom(*this, other);
        num -= other.num;
        normalize(*this);
        
        return *this;
    }

    bool operator < (Rational& other){
        eq_denom(*this, other);
        bool ans = num < other.num;
        normalize(*this);
        normalize(other);
    
        return ans;
    }

    friend std::ostream& operator << (std::ostream& os, const Rational& rat){
        if(rat.denom == 1)
            os << rat.num;
        else if(rat.num == 0)
            os << 0;
        else
            os << rat.num << "/" << rat.denom;
        return os;
    }
};

long int det(float** A, int n);

long int M(float** A, int n, int J, int I = 0);

Rational det_Gausses(std::vector<std::vector<Rational>> A, int n);

int main() {
    srand(time(0));

    const int N = 8;

    std::vector<std::vector<Rational>> A;

    for(size_t i = 0; i < N; i++){
        std::vector<Rational> temp;
        for(size_t j = 0; j < N; j++){
            Rational rat(Rational(rand() % 10));
            temp.push_back(rat);
            std::cout << rat << ' ';
        }
        A.push_back(temp);
        std::cout << '\n';
    }

    std::cout << '\n' << "====TRIANG====\n";
    det_Gausses(A, N);

    return 0;
}

long int det(float** A, int n){

    long int det = 0;

    if(n < 1){
        throw "Unable to find a det";
    }
    else if(n == 1){
        return A[0][0];
    }
    else if(n == 2){
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }
    else{
        for(int j = 0; j < n; j++){
            long int Min = M(A, n, j);
            det += pow(-1, j) * A[0][j] * Min;
        }
    }

    return det;
}

long int M(float** A, int n, int J, int I){

    int ki = 0, kj = 0;

    float** B = new float*[n-1]{};
    for(size_t i = 0; i < n-1; i++)
        B[i] = new float[n-1];

    for(size_t i = 0; i < n-1; i++){
        kj = 0;
        if(i == I)
            ki = 1;
        for(size_t j = 0; j < n-1; j++){
            if(j == J)
                kj = 1;
            B[i][j] = A[i+ki][j+kj];
        }
    }

    long int M = det(B, n-1);

    for(size_t i = 0; i < n-1; i++){
        delete[] B[i];
    }
    delete[] B;

    return M;
}

Rational det_Gausses(std::vector<std::vector<Rational>> A, int n){
    Rational det = 1;

    Rational max_val = 0;
    int row_max_val = 0;
    for(size_t i = 0; i < n; i++){
        if(max_val < A[i][0]){
            max_val = A[i][0];
            row_max_val = i;
        }
    }

    // std::cout << '\n' << '\n' << "====AFTERMAX====\n";
    if(row_max_val != 0){
        std::swap(A[0], A[row_max_val]);
        det *= Rational(-1);
    }

    // for(size_t i = 0; i < n; i++){
    //     for(size_t j = 0; j < n; j++){
    //         std::cout << A[i][j] << ' ';
    //     }
    //     std::cout << '\n';
    // }

    det *= A[0][0];

    // std::cout << '\n' << '\n' << "====AFTERSHIT====\n";
    // if(std::signbit(A[1][0]))
    //     coef *= -1;
    try{
        for(size_t t = 1; t < n; t++){
            Rational coef = A[t][0] / A[0][0];
            for(size_t k = 0; k < n; k++){
                auto at = coef * A[0][k];
                // std::cout << A[t][k] << " coef: " << coef << " at: " << at << '\n';
                A[t][k] -= at;
            }
        }
    }
    catch(const std::exception& e){
        std::cout << e.what() << '\n';
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            std::cout << A[i][j] << ' ';
        }
        std::cout << '\n';
    }

    return det;
}