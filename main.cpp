#include <iostream>
#include "cmath"
#include <exception>
#include <chrono>
#include <vector>
#include <algorithm>
// #include <omp.h>

class Rational{
protected:
    void eq_denom(Rational& oth){
        if (this->denom < 0) {
            this->num *= -1;
            this->denom *= -1;
        }
        if (oth.denom < 0) {
            oth.num *= -1;
            oth.denom *= -1;
        }

        if(this->denom != oth.denom){
            auto lcm_num = this->denom / gcd(this->denom, oth.denom) * oth.denom;

            this->num = (lcm_num / this->denom) * this->num;
            this->denom = lcm_num;
            oth.num = (lcm_num / oth.denom) * oth.num;
            oth.denom = lcm_num;
        }
    }

    static __int128 gcd(__int128 a, __int128 b){
        if(a < 0) a = -a;
        if(b < 0) b = -b;

        while(b != 0){
            int temp = a % b;
            a = b;
            b = temp;
        }

        return a;
    }

    
public:
    __int128 num;
    __int128 denom;

    Rational(__int128 n) : num{n}, denom{1} { };

    Rational(__int128 n, __int128 d){
        if(d == 0)
            throw std::invalid_argument("Division by zero");
        if(d < 0){
            d *= -1;
            n *= -1;
        }
        num = n;
        denom = d;
    }

    Rational() { };

    void normalize(){
        if (this->denom < 0) {
            this->num = -this->num;
            this->denom = -this->denom;
        }
        if (this->denom == 0)
            throw std::invalid_argument("Division by zero");
        if (this->num == 0) {
            this->denom = 1;
            return;
        }
        auto gcd_num = gcd(this->num, this->denom);
        this->num /= gcd_num;
        this->denom /= gcd_num;
    }

    Rational operator + (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        temp_this.eq_denom(temp_other);
        Rational answer(temp_this.num + temp_other.num, temp_this.denom);
        answer.normalize();

        return answer;
    }

    Rational operator - (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        temp_this.eq_denom(temp_other);
        Rational answer(temp_this.num - temp_other.num, temp_this.denom);
        answer.normalize();

        return answer;
    }

    Rational operator * (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        auto g1 = gcd(temp_this.num, temp_other.denom);
        auto g2 = gcd(temp_other.num, temp_this.denom);

        Rational answer((temp_this.num / g1) * (temp_other.num / g2), (temp_this.denom / g2) * (temp_other.denom / g1));
        answer.normalize();

        return answer;
    }

    Rational operator / (const Rational& other){
        if(other.num == 0)
            throw std::invalid_argument("Division by zero");
        Rational temp_other;

        temp_other.num = other.denom;
        temp_other.denom = other.num;

        if(temp_other.denom < 0){
            temp_other.num = -temp_other.num;
            temp_other.denom = -temp_other.denom;
        }

        return *this * temp_other;
    }

    Rational& operator = (const Rational& other){
        if(other.num != num || other.denom != denom){
            num = other.num;
            denom = other.denom;
        }
        this->normalize();

        return *this;
    }

    Rational& operator = (const __int128& other){
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

        temp_this.eq_denom(temp_other);
        temp_this.num += temp_other.num;
        temp_this.normalize();
        *this = temp_this;

        return *this;
    }

    Rational& operator -= (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        temp_this.eq_denom(temp_other);
        temp_this.num -= temp_other.num;
        temp_this.normalize();
        *this = temp_this;
        
        return *this;
    }

    bool operator < (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        temp_this.eq_denom(temp_other);
        bool ans = temp_this.num < temp_other.num;
    
        return ans;
    }

    bool operator == (const Rational& other){
        Rational temp_this = *this;
        Rational temp_other = other;

        temp_this.eq_denom(temp_other);
        bool ans = temp_this.num == temp_other.num;
    
        return ans;
    }

    friend std::ostream& operator << (std::ostream& os, Rational rat){
        std::string num;

        if(rat.num == 0)
            num = "0";
        else{
            bool negative;
            if(rat.num < 0){
                negative = true;
                rat.num = -rat.num;
            }
            else
                negative = false;

            while(rat.num > 0){
                num.push_back('0' + (rat.num % 10));
                rat.num /= 10;
            }
            if(negative)
                num.push_back('-');

            std::reverse(num.begin(), num.end());
        }

        if(rat.denom == 1){
            os << num;
        }
        else if(num == "0")
            os << 0;
        else{
            std::string denom;
            while(rat.denom > 0){
                denom.push_back('0' + (rat.denom % 10));
                rat.denom /= 10;
            }

            std::reverse(denom.begin(), denom.end());

            os << num << "/" << denom;
        }
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

        std::cout << "\n\n====AFTER MAXING====\n";
        for(size_t l = 0; l < n; l++){
            for(size_t p = 0; p < n; p++){
                std::cout << A[l][p] << ' ';
            }
            std::cout << '\n';
        }
        std::cout << '\n';

        det *= A[i][i];

        try{
            for(size_t t = i+1; t < n; t++){
                Rational coef = A[t][i] / A[i][i];
                std::cout << "cur coef: " << coef << '\n';
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