#include <nfl.hpp>
#include <iostream>
#include <gmpxx.h>

int main() {
    constexpr int degree = 4; // degree of polynomial
    constexpr int nModul = 2; // number of moduli to present Z_t
    using poly_t = nfl::poly_p<uint16_t, degree, nModul>;
    poly_t ply;
    ply.set({1,2,3, 4});

    poly_t ply2, result;
    ply2.set({4, 3, 2, 1});

    std::cout << ply << std::endl;
    std::cout << ply2 << std::endl;
    // convert to NTT form take O(nlogn)
    ply.ntt_pow_phi();
    ply2.ntt_pow_phi();
    // element-wise multiplication take O(n)
    result = ply * ply2;

    // convert back to poly form, take O(nlogn)
    result.invntt_pow_invphi();
    std::cout << result << std::endl; 

    // The product to moduli is return as a mpz object
    gmp_printf ("P = %Zd\n", poly_t::moduli_product());
    return 0;
}

