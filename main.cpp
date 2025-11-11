#include <iostream>
#include "cmath"
#include <exception>
#include <chrono>
#include <vector>
// #include <omp.h>

class Rational{
protected:
    void eq_denom(Rational& fr, Rational& sc){
        if(fr.denom > 0 && sc.denom < 0 || fr.denom < 0 && sc.denom > 0){
            sc.denom *= -1;
            sc.num *= -1;
        }

        if(fr.denom != sc.denom){
            int temp_d = fr.denom;
            fr.num *= sc.denom;
            fr.denom *= sc.denom;
            sc.num *= temp_d;
            sc.denom *= temp_d;
        }
    }

    void normalize(Rational& fr){
        if (fr.denom == 0)
            throw "Denominator can't be zero";
        if (fr.num == 0) {
            fr.denom = 1;
            return;
        }
        int gcd_num = gcd(fr.num, fr.denom);
        fr.num /= gcd_num;
        fr.denom /= gcd_num;

        if (fr.denom < 0) {
            fr.num = -fr.num;
            fr.denom = -fr.denom;
        }
    }

    int gcd(int a, int b){
        a = std::abs(a);
        b = std::abs(b);

        while(b != 0){
            int temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }
public:
    long int num;
    long int denom;

    Rational(long int n) : num{n}, denom{1} { };

    Rational(long int n, long int d){
        if(d == 0)
            throw "Denominator can't be zero\n";
        num = n;
        denom = d;
    }

    Rational() { };

    Rational operator + (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);
        Rational answer(temp_this.num + temp_other.num, temp_this.denom);
        normalize(answer);

        return answer;
    }

    Rational operator - (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);        
        Rational answer(temp_this.num - temp_other.num, temp_this.denom);
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
        normalize(*this);

        return *this;
    }

    Rational& operator = (const long int& other){
        num = other;
        denom = 1;
        normalize(*this);

        return *this;
    }

    Rational& operator *= (const Rational& other){
        num *= other.num;
        denom *= other.denom;
        normalize(*this);

        return *this;
    }

    Rational& operator += (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);        
        temp_this.num += temp_other.num;
        normalize(temp_this);
        *this = temp_this;

        return *this;
    }

    Rational& operator -= (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);        
        temp_this.num -= temp_other.num;
        normalize(temp_this);
        *this = temp_this;
        
        return *this;
    }

    bool operator < (Rational& other){
        eq_denom(*this, other);
        bool ans = num < other.num;
        normalize(*this);
        normalize(other);
    
        return ans;
    }

    bool operator == (Rational other){
        eq_denom(*this, other);
        bool ans = num == other.num;
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

void new_main();

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

    auto start_time = std::chrono::high_resolution_clock::now();

    Rational ans = det_Gausses(A, N);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur = end_time - start_time;

    std::cout << "Answer is: " << ans << " Time is: " << dur.count() << '\n';

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

    for(size_t i = 0; i < n; i++){
        Rational max_val = A[i][i];
        int row_max_val = i;
        for(size_t k = i; k < n; k++){
            if(max_val < A[k][i]){
                max_val = A[k][i];
                row_max_val = k;
            }
        }

        if(row_max_val != i){
            std::swap(A[i], A[row_max_val]);
            det *= Rational(-1);
        }

        det *= A[i][i];

        try{
            // std::cout << "Coef: ";
            for(size_t t = i+1; t < n; t++){
                // if(A[t][i] == (Rational(0)))
                //     continue;
                Rational coef = A[t][i] / A[i][i];
                std::cout << coef << ' ';
                for(size_t k = i; k < n; k++){
                    auto at = coef * A[i][k];
                    A[t][k] -= at;
                }
            }
            // std::cout << '\n';
        }
        catch(const std::exception& e){
            std::cout << e.what() << '\n';
        }

        std::cout << "Stage: " << i << '\n';
        for(size_t k = 0; k < n; k++){
            for(size_t j = 0; j < n; j++){
                std::cout << A[k][j] << ' ';
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }

    std::cout << "\n\n" << "=====FINAL=====\n";

    for(int k = 0; k < n; k++){
            for(int j = 0; j < n; j++){
                std::cout << A[k][j] << ' ';
            }
            std::cout << '\n';
        }

    return det;
}