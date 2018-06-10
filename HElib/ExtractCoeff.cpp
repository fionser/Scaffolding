#include "HElib/intraSlot.h"
#include "HElib/FHEContext.h"
#include "HElib/EncryptedArray.h"
#include "HElib/NumbTh.h"

#include <iostream>
#include <algorithm>
#include <ctime>
NTL::ZZX pack_polys_as_poly(const EncryptedArray *ea) {
    const long l = ea->size();
    const long d = ea->getDegree();
    // pack `l` polys each with `d` degree
    std::vector<NTL::ZZX> polys(l);
    for (auto &poly : polys) {
        poly.SetLength(1);
        poly[0] = 1;
        //for (int _d = 0; _d < d; ++_d)
        //    NTL::SetCoeff(poly, _d, _d + 1);
    }

    std::cout << "packed: " << polys << "\n";
    NTL::ZZX ans;
    ea->encode(ans, polys);
    return ans;
}

std::vector<long> gen_vector(long length, long range) {
    std::vector<long> vec(length);
    for (int i = 0; i < length; ++i)
        vec[i] = NTL::RandomBnd(range);
    return vec;
}

NTL::ZZX gen_packed_vectors(std::vector<std::vector<long>> const& vecs,
                            const EncryptedArray *ea) {
    const long l = ea->size();
    const long d = ea->getDegree();

    std::vector<NTL::ZZX> slots(l);
    for (long i = 0; i < l; ++i) {
        slots[i].SetLength(vecs.at(i).size());
        for (long k = 0; k < d; ++k)
            NTL::SetCoeff(slots[i], k, vecs[i][k]);
    }

    NTL::ZZX packed;
    ea->encode(packed, slots);
    return packed;
}

long inner_product(std::vector<long> const& v1, 
                   std::vector<long> const& v2,
                   long range) {
    long ip = 0;
    size_t len = v1.size();
    assert(v2.size() == len);
    for (size_t i = 0; i < len; ++i) {
        ip += (v1[i]) * (v2[i]);
        ip %= range;
    }
    return ip;
}


template<class Poly>
void myApplyLinPolyLL(Ctxt& ctxt, 
                      const vector<Poly>& encodedC, 
                      long d)
{
    assert(d == lsize(encodedC));

    ctxt.cleanUp();  // not sure, but this may be a good idea
    ctxt.multByConstant(encodedC[0]);

    Ctxt origin(ctxt);

    for (long j = 0; j < d - 1; j++) {
        Ctxt power_j(origin);
        power_j.frobeniusAutomorph(j + 1);
        ctxt += power_j;
    }
}

int main(int argc, char *argv[]) {
    long m = 64;
    long p = 11;
    ArgMapping amap;
    amap.arg("m", m, "m");
    amap.arg("p", p, "p");
    amap.parse(argc, argv);

    FHEcontext context(m, p, 1);
    context.bitsPerLevel = 40;
    buildModChain(context, 3);
    std::cout << "kappa = " << context.securityLevel() << std::endl;
    FHESecKey sk(context);
    sk.GenSecKey(64);
    addFrbMatrices(sk);
    //addSomeFrbMatrices(sk);
    //addBSGSFrbMatrices(sk);
    std::cout << "|ks| = " << sk.keySWlist().size() << std::endl;
    {
        std::ofstream ss("sk.sk");
        if (ss.is_open()) {
            ss << sk;
            ss.close();
        }
    }
    const FHEPubKey &pk = sk;

    const EncryptedArray *ea = context.ea;
    long l = ea->size();
    long d = ea->getDegree();
    std::cout << "l = " << l << " d = " << d << std::endl;
    auto _start = std::clock();
    std::vector<NTL::ZZX> L(d, NTL::to_ZZX(0));
    L[d - 1] = NTL::to_ZZX(1); // extract the coefficent of X^{d-1} as the coefficent of X^0
    std::vector<NTL::ZZX> coeff;
    ea->buildLinPolyCoeffs(coeff, L);
    std::vector<NTL::ZZX> encodedC(d);
    for (int i = 0; i < d; ++i) {
        std::vector<NTL::ZZX> tmp(l, coeff[i]); // l copies of coeff[i];
        NTL::ZZX tmp2;
        ea->encode(tmp2, tmp);
        encodedC[i] = tmp2;
    }
    auto _end = std::clock();
    std::cout << "prepare coeffs " << (_end - _start) / (double) CLOCKS_PER_SEC << std::endl;

    std::vector<std::vector<long>> vecs_A(l);
    std::vector<std::vector<long>> vecs_B(l);
    std::vector<long> excepted(l);
    for (long i = 0; i < l; ++i) {
        vecs_A[i] = gen_vector(d, p);
        vecs_B[i] = gen_vector(d, p);
        excepted[i] = inner_product(vecs_A[i], vecs_B[i], p);
        std::reverse(vecs_B[i].begin(), vecs_B[i].end());
    }

    NTL::ZZX packed_vecs_A = gen_packed_vectors(vecs_A, ea);
    NTL::ZZX packed_vecs_B = gen_packed_vectors(vecs_B, ea);

    Ctxt enc_vecs_A(pk), enc_vecs_B(pk);
    pk.Encrypt(enc_vecs_A, packed_vecs_A);
    {
        std::ofstream css("ctx.ctx");
        if (css.is_open()) {
            css << enc_vecs_A;
            css.close();
        }
    }
    //std::cout << "threads " << NTL::AvailableThreads() << std::endl;
    _start = std::clock();
    Ctxt computed(enc_vecs_A);
    computed.multByConstant(packed_vecs_B);
    _end = std::clock();
    std::cout << "computed inner_product: " << (_end - _start) / (double) CLOCKS_PER_SEC << std::endl;
    //computed.multiplyBy(enc_vecs_B);

    //_start = std::clock();
    FHE_NTIMER_START(LINEAR_MAP);
    myApplyLinPolyLL(computed, encodedC, d);
    FHE_NTIMER_STOP(LINEAR_MAP);

    printNamedTimer(std::cout, "LINEAR_MAP");
    //_end = std::clock();
    //std::cout << "apply linear map: " << (_end - _start) / (double) CLOCKS_PER_SEC << std::endl;


    _start = std::clock();
    std::vector<NTL::ZZX> slots;
    ea->decrypt(computed, sk, slots);
    _end = std::clock();
    std::cout << "decryption: " << (_end - _start) / (double) CLOCKS_PER_SEC << std::endl;

    std::vector<long> _computed;
    for (auto &slot : slots) {
        _computed.emplace_back(NTL::to_long(NTL::coeff(slot, 0)));
    }

    std::cout << excepted.size() << "\n";
    std::cout << excepted << "\n";
    std::cout << _computed << "\n";
    printNamedTimer(std::cout, "MyFrobenius");
    printNamedTimer(std::cout, "MyMultByConstant");

    return 0;
}
