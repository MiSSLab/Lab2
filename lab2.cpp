#include <iostream>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include <iomanip>
#include <vector>
#include <math.h>

using namespace std;

#define N 16777216
#define k_chi 10

void print(mpz_t dividend, mpz_t divisor,  int d) {

    mpz_t quotient, quotient_to_print;
    mpz_t remainder;
    mpz_init(quotient);
    mpz_init(quotient_to_print);
    mpz_init(remainder);

    mpz_tdiv_qr(quotient, remainder, dividend, divisor);
    gmp_printf("%Zd", quotient);
    if(mpz_cmp_si(remainder,0) != 0){
        gmp_printf(".");

        for(int i = 0; i < d && mpz_cmp_si(remainder,0) != 0; i++){
            mpz_mul_ui(remainder, remainder, 10);  
            mpz_tdiv_qr(quotient, remainder, remainder, divisor);
            if( i == d-1 && mpz_cmp_si(quotient,0) == 0){
                break;
            }else{
                gmp_printf("%Zd", quotient);
            }
        }
    }
    gmp_printf("\n");
    return;
}

void chi_parameters(mpq_class ranges_chi[], mpq_class probabilities_chi[]){
    mpq_class range_length;

    for (int i = 0; i <= k_chi; i++)
    {
        //count numerator i^2
        ranges_chi[i].get_num() = pow(i, 2);

        //count range length
        if(i!=0){
            range_length = ranges_chi[i].get_num() - ranges_chi[i-1].get_num();
            probabilities_chi[i-1].get_num() = range_length;
            probabilities_chi[i-1].get_den() = 100;
        }
        ranges_chi[i].get_den() = 100;
    }
}

void chi_check_range(mpq_class input[], mpq_class Y_chi[], mpq_class ranges_chi[], unsigned long int n){
    //check range for given point
    for(int i = 0; i < n; i++){
        for (int a = 1; a <= k_chi; a++){
            if (input[i] < ranges_chi[a]){
                Y_chi[a-1] ++;
                break;
            }
        }
    }
}

void chi_calculate(mpq_class Y_chi[], mpq_class probabilities_chi[], unsigned long int n, int d){
    mpq_class result, tmp_q, npi;

    //calculate result
    for(int i = 0; i < k_chi; i++){
        npi = n * probabilities_chi[i];
        tmp_q = Y_chi[i] - npi;
        tmp_q *= tmp_q;
        tmp_q /= npi;
        result += tmp_q;
    }

    print(result.get_num().get_mpz_t(), result.get_den().get_mpz_t(), d);
}

void sort(mpq_class input[], int left, int right) {
      int i = left, j = right;
      mpq_class tmp;
      mpq_class pivot = input[(left + right) / 2];

      while (i <= j) {
            while (input[i] < pivot)
                  i++;
            while (input[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = input[i];
                  input[i] = input[j];
                  input[j] = tmp;
                  i++;
                  j--;
            }
      };

      if (left < j)  sort(input, left, j);
      if (i < right) sort(input, i, right);
}

int main(int argc, char* argv[]) {

    if(argc!= 2){
        printf("Wrong number of arguments.\n");
        return 1;
    }

    int d = atoi(argv[1]);

    //set precision
    mp_bitcnt_t prec = ceil(d * log(10)/log(2));
    mpf_set_default_prec(prec + 128);

    //read rational numbers from input
    unsigned long int n = 0;

    string new_line;
    string input_number;
    mpq_class* input = new mpq_class[N];

    getline(cin, new_line);
    stringstream input_stream(new_line);

    while(true){
        if(!(input_stream >> input_number)) {
            break;
        }
        input[n] = input_number;
        n++;
    }

    //chi square
    mpq_class ranges_chi[k_chi + 1];
    mpq_class probabilities_chi[k_chi];
    mpq_class Y_chi[k_chi];

    for (int i = 0; i < k_chi; i++){
        Y_chi[i] = 0;
    }

    chi_parameters(ranges_chi, probabilities_chi);
    chi_check_range(input, Y_chi, ranges_chi, n);
    chi_calculate(Y_chi, probabilities_chi, n, d);

    //kolmogorov
    mpq_class k_plus = 0;
    mpq_class k_minus = 0;
    mpq_class k_plus_tmp = 0;
    mpq_class k_minus_tmp = 0;
    mpq_class n_mpq = n;

    sort(input, 0, n);

    for (int i = 0; i < n; i++)
    {
        k_plus_tmp  = (i/n_mpq - input[i]);
        k_minus_tmp = (input[i] - (i-1)/n_mpq);

        if (k_plus_tmp > k_plus){
            k_plus = k_plus_tmp;
        }
        if (k_minus_tmp > k_minus){
            k_minus = k_minus_tmp;
        }
    }

    print(k_plus.get_num().get_mpz_t(), k_plus.get_den().get_mpz_t(), d);
    print(k_minus.get_num().get_mpz_t(), k_minus.get_den().get_mpz_t(), d);

    return 0;
}