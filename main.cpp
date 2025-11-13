#include <iostream>
#include "cmath"
#include <exception>
#include <chrono>
#include <vector>
// #include <omp.h>

class Rational{
protected:
    void eq_denom(Rational& fr, Rational& sc){
        if (fr.denom < 0) {
            fr.num = -fr.num;
            fr.denom = -fr.denom;
        }
        if (sc.denom < 0) {
            sc.num = -sc.num;
            sc.denom = -sc.denom;
        }

        if(fr.denom != sc.denom){
            auto lcm_num = fr.denom * sc.denom / gcd(fr.denom, sc.denom);

            fr.num *= lcm_num / fr.denom;
            fr.denom *= lcm_num / fr.denom;
            sc.num *= lcm_num / sc.denom;
            sc.denom *= lcm_num / sc.denom;
        }
    }

    void normalize(){
        if (this->denom < 0) {
            this->num = -this->num;
            this->denom = -this->denom;
        }
        if (this->denom == 0)
            throw std::invalid_argument("Denominator can't be zero");
        if (this->num == 0) {
            this->denom = 1;
            return;
        }
        auto gcd_num = gcd(this->num, this->denom);
        this->num /= gcd_num;
        this->denom /= gcd_num;
    }

    long long int gcd(long long int a, long long int b){
        a = std::abs(a); b = std::abs(b);

        while(b != 0){
            int temp = a % b;
            a = b;
            b = temp;
        }

        return a;
    }

    
public:
    long long int num;
    long long int denom;

    Rational(long long int n) : num{n}, denom{1} { };

    Rational(long long int n, long long int d){
        if(d == 0)
            throw std::invalid_argument("Denominator can't be zero");
        if(d < 0){
            d *= -1;
            n *= -1;
        }
        num = n;
        denom = d;
    }

    Rational() { };

    Rational operator + (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);
        Rational answer(temp_this.num + temp_other.num, temp_this.denom);
        answer.normalize();

        return answer;
    }

    Rational operator - (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);        
        Rational answer(temp_this.num - temp_other.num, temp_this.denom);
        answer.normalize();

        return answer;
    }

    Rational operator * (const Rational& other){
        Rational answer(num * other.num, denom * other.denom);
        answer.normalize();

        return answer;
    }

    Rational operator / (const Rational& other){
        Rational answer(num * other.denom, denom * other.num);
        answer.normalize();

        return answer;
    }

    Rational& operator = (const Rational& other){
        if(other.num != num || other.denom != denom){
            num = other.num;
            denom = other.denom;
        }
        this->normalize();

        return *this;
    }

    Rational& operator = (const long int& other){
        num = other;
        denom = 1;
        this->normalize();

        return *this;
    }

    Rational& operator *= (const Rational& other){
        num *= other.num;
        denom *= other.denom;
        this->normalize();

        return *this;
    }

    Rational& operator += (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);        
        temp_this.num += temp_other.num;
        temp_this.normalize();
        *this = temp_this;

        return *this;
    }

    Rational& operator -= (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);        
        temp_this.num -= temp_other.num;
        temp_this.normalize();
        *this = temp_this;
        
        return *this;
    }

    bool operator < (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);
        bool ans = temp_this.num < temp_other.num;
    
        return ans;
    }

    bool operator == (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        eq_denom(temp_this, temp_other);
        bool ans = temp_this.num == temp_other.num;
    
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

    auto start_time = std::chrono::high_resolution_clock::now();

    Rational ans = det_Gausses(A, N);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur = end_time - start_time;

    std::cout << "Answer is: " << ans << " Time is: " << dur.count() << '\n';

    return 0;
}

Rational det_Gausses(std::vector<std::vector<Rational>> A, int n){
    Rational det = 1;

    for(size_t i = 0; i < n; i++){
        Rational max_val = A[i][i];
        auto row_max_val = i;
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
            for(size_t t = i+1; t < n; t++){
                Rational coef = A[t][i] / A[i][i];
                for(size_t k = i; k < n; k++){
                    auto at = coef * A[i][k];
                    A[t][k] -= at;
                }
            }
        }
        catch(const std::exception& e){
            std::cout << e.what() << '\n';
        }

        std::cout << "\n====STAGE: " << i << " ====\n";
        for(size_t k = 0; k < n; k++){
            for(size_t l = 0; l < n; l++){
                std::cout << A[k][l] << ' ';
            }
            std::cout << '\n';
        }
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