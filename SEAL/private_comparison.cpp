#include <iostream>
#include <numeric>
#include <iomanip>

#include <seal/seal.h>
#include <NTL/ZZ.h>

struct Comparator {
    seal::Plaintext test;
    seal::Plaintext half;
};

void init(Comparator& comp, seal::SEALContext const &context)
{
    long t = context.plain_modulus().value();
    uint64_t half = t - static_cast<uint64_t >(NTL::InvMod(2, t));

    const int n = context.poly_modulus().coeff_count() - 1;
    comp.test.resize(n);
    for (long i = 0; i < n; ++i)
        comp.test[i] = half;
    comp.half.resize(1);
    comp.half[0] = t - half;
}

seal::Ciphertext private_compare(
        seal::Ciphertext const& c1,
        seal::Ciphertext c2,
        Comparator const& comparator,
        seal::Evaluator & evl,
        seal::GaloisKeys const& gkeys,
        seal::EvaluationKeys const& evks,
        seal::SEALContext const& context)
{
    const uint64_t n = static_cast<uint64_t >(context.poly_modulus().coeff_count() - 1);
    evl.apply_galois(c2, (n << 1) - 1, gkeys); // X^b --> X^{2n - b} which is X^{-b} over modulo X^n + 1.
    evl.multiply(c2, c1);
    evl.multiply_plain(c2, comparator.test);
    evl.add_plain(c2, comparator.half);
    evl.relinearize(c2, evks); // convert back to size-2 ciphertext
    return c2;
}


int main() {
    seal::EncryptionParameters params;
    params.set_poly_modulus("1x^4096 + 1");
    params.set_coeff_modulus(seal::coeff_modulus_128(4096));
    params.set_plain_modulus(1013);

    seal::SEALContext context(params);
    seal::KeyGenerator keygen(context);

    seal::PublicKey const& pk = keygen.public_key();
    seal::SecretKey const& sk = keygen.secret_key();

    const int bit_decomp = 30;
    seal::GaloisKeys gal_keys;
    seal::EvaluationKeys evl_keys;
    const uint64_t m_1 = (static_cast<uint64_t>(context.poly_modulus().coeff_count() - 1) * 2) - 1;

    keygen.generate_galois_keys(bit_decomp, {m_1}, gal_keys);
    keygen.generate_evaluation_keys(bit_decomp, evl_keys); // by default s^2 -> s

    seal::Encryptor encryptor(context, pk);
    seal::Decryptor decryptor(context, sk);
    seal::Evaluator evaluator(context);

    seal::Plaintext p("1x^34");
    seal::Plaintext p2("1x^33");
    seal::Ciphertext cp, cp2;
    encryptor.encrypt(p, cp);
    encryptor.encrypt(p2, cp2);


    Comparator comparator;
    init(comparator, context);

    {
        auto _st = std::clock();
        auto ans = private_compare(cp, cp2, comparator, evaluator,
                                   gal_keys, evl_keys, context);
        auto _end = std::clock();
        std::cout << "compare " << (_end - _st) / (double) CLOCKS_PER_SEC * 1000. << "\n";
        decryptor.decrypt(ans, p);
        std::cout << p[0] << "\n"; // 34 > 33, so this should be 1
    }

    {
        auto _st = std::clock();
        auto ans = private_compare(cp2, cp, comparator, evaluator,
                                   gal_keys, evl_keys, context);
        auto _end = std::clock();
        std::cout << "compare " << (_end - _st) / (double) CLOCKS_PER_SEC * 1000. << "\n";
        decryptor.decrypt(ans, p);
        std::cout << p[0] << "\n"; // the opposite side, this should be 0
    }
    return 0;
}
