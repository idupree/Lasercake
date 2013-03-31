/*

    Copyright Eli Dupree and Isaac Dupree, 2011, 2012, 2013

    This file is part of Lasercake.

    Lasercake is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Lasercake is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Lasercake.  If not, see <http://www.gnu.org/licenses/>.

*/

// When you add a new tests file, define a new name here and with
// DECLARE_TESTS_FILE near the top of test_header.hpp, and put at
// the bottom of your tests file:
// REGISTER_TESTS // This must come last in the file.
#define TESTS_FILE bignum_tests
#include "test_header.hpp"


#include "../data_structures/bignum.hpp"


#include "../data_structures/bounds_checked_int.hpp"
#include "../units.hpp"

/*
extern void g(bignum<6> const&);
extern unsigned int q,p,r,s,t,u,v,w;
void f() { g(long_multiply_signed<6>(bignum<4>{{q,p,t,u}},bignum<4>{{r,s,v,w}})); }

template class bigint<192>;
template class biguint<192>;
template class bigint<0>;
template class biguint<0>;
*/

static void bignum_compile_test() {
  bigint<128> aa = bigint<64>{};
  bigint<128> ab(bigint<192>{});
  biguint<128> ac = biguint<64>{};
  biguint<128> ad(biguint<192>{});
  +aa;+ac;
  -aa;-ac;
  aa+aa;ac+ac;
  aa-aa;ac-ac;
  aa*aa;ac*ac;
  multiply_to<64>(aa,aa);
  multiply_to<64>(ac,ac);
  lossless_multiply(aa,aa);
  lossless_multiply(ac,ac);
  aa&aa;ac&ac;
  aa|aa;ac|ac;
  aa^aa;ac^ac;
  ~aa;~ac;
  aa>>90;ac>>90;
  aa<<90;ac<<90;
  aa==aa;ac==ac;
  aa!=aa;ac!=ac;
  aa<aa;ac<ac;
  aa<=aa;ac<=ac;
  aa>aa;ac>ac;
  aa>=aa;ac>=ac;
  
  aa+=aa;ac+=ac;
  aa-=aa;ac-=ac;
  aa*=aa;ac*=ac;
  aa&=aa;ac&=ac;
  aa|=aa;ac|=ac;
  aa^=aa;ac^=ac;
  aa>>=90;ac>>=90;
  aa<<=90;ac<<=90;
  
  (++aa)++; (++ac)++;
  (--aa)--; (--ac)--;
  
  !aa; !ac;
  (bool)(aa); (bool)(ac);
  
  //aa = ac;
  // explicit sign conversion
  aa = decltype(aa)(ac);
  ac = decltype(ac)(aa);
  
  std::numeric_limits<decltype(aa)>::max();
  std::numeric_limits<decltype(ac)>::max();
  
  (float)(aa);
  (float)(ac);
  (double)(aa);
  (double)(ac);
  //(long double)(aa);
  //(long double)(ac);
  
  abs(aa);
  abs(ac);
  hash_value(aa);
  hash_value(ac);
  ilog2(aa);
  ilog2(ac);
  isqrt(aa);
  isqrt(ac);
  popcount(aa);
  popcount(ac);
  count_trailing_zeroes(aa);
  count_trailing_zeroes(ac);
  
  //bounds_checked_int<bigint<256> > fsddfs;
  physical_quantity<bigint<256>, units<dim::meter<1>> > dist;
  dist+dist;
  dist*dist*dist;
}

template<typename BigInt>
static void nonnegative_tests() {
  //LOG << __PRETTY_FUNCTION__ << '\n';
  BOOST_CHECK_EQUAL(BigInt(15) >> 2, BigInt(3));
  BOOST_CHECK_EQUAL(BigInt(16) >> 2, BigInt(4));
  BOOST_CHECK_EQUAL(BigInt(17) >> 2, BigInt(4));
  BOOST_CHECK_EQUAL(BigInt(13) << 2, BigInt(52));

  BOOST_CHECK_EQUAL(BigInt(15) >> 0, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 0, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) >> 32, BigInt(0));
  BOOST_CHECK_EQUAL(BigInt(15) >> 64, BigInt(0));
  BOOST_CHECK_GT(BigInt(15) << 32, BigInt(999999));
  BOOST_CHECK_GT(BigInt(15) << 64, BigInt(999999));

  BOOST_CHECK_EQUAL(BigInt(15) << 32 >> 32, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 64 >> 64, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 64 >> 32 >> 32, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 64 >> 21 >> 21 >> 22, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 64 >> 40 >> 24, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 64 >> 24 >> 40, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 21 << 21 << 22 >> 64, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 40 << 24 >> 64, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 24 << 40 >> 64, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 66 >> 22 >> 22 >> 22, BigInt(15));
  BOOST_CHECK_EQUAL(BigInt(15) << 22 << 22 << 22 >> 66, BigInt(15));
  
  BOOST_CHECK_EQUAL(BigInt(15) * 3, BigInt(45));
  BOOST_CHECK_EQUAL(BigInt(0x100000000) * 0x100000000, BigInt(0x100000000) << 32);
  BOOST_CHECK_EQUAL(BigInt(0x100000001) * 0x100000001, (BigInt(0x100000000) << 32) + BigInt(0x200000001));

  BOOST_CHECK_GT(BigInt(1e20), BigInt(1) << 66);
  BOOST_CHECK_LT(BigInt(1e20), BigInt(1) << 67);
  BOOST_CHECK_GT(BigInt(1e20f), BigInt(1) << 66);
  BOOST_CHECK_LT(BigInt(1e20f), BigInt(1) << 67);
}



BOOST_AUTO_TEST_CASE(bignum_runtests) {
  //bigint<128> aa = bigint<64>(1000000000LL);
  //aa = aa*aa*aa*(aa*(aa*aa));
  bigint<512> a1 = bigint<64>(int64_t(1000000000LL));
  bigint<512> aa = a1*a1*a1*(a1*(a1*a1));
  //LOG << aa << "\n";
  //LOG << decltype(aa)(2) << "\n";
  
  nonnegative_tests<biguint<128>>();
  nonnegative_tests<bigint<128>>();
  nonnegative_tests<biguint<192>>();
  nonnegative_tests<bigint<192>>();

  BOOST_CHECK_EQUAL(bigint<128>(-1) >> 1, bigint<128>(-1));
  BOOST_CHECK_EQUAL(bigint<128>(-2) >> 1, bigint<128>(-1));
  BOOST_CHECK_EQUAL(bigint<128>(-3) >> 1, bigint<128>(-2));
  BOOST_CHECK_EQUAL(bigint<128>(-4) >> 1, bigint<128>(-2));
  BOOST_CHECK_EQUAL(bigint<128>(-5) >> 1, bigint<128>(-3));
  BOOST_CHECK_EQUAL(bigint<128>(-5) << 66 >> 22 >> 22 >> 22, bigint<128>(-5));
  BOOST_CHECK_EQUAL(bigint<128>(-5) << 64 >> 21 >> 21 >> 22, bigint<128>(-5));
  BOOST_CHECK_EQUAL(bigint<128>(-5) << 22 << 22 << 22 >> 66, bigint<128>(-5));
  BOOST_CHECK_EQUAL(bigint<128>(-5) << 22 << 21 << 21 >> 64, bigint<128>(-5));
  
  BOOST_CHECK_GT(double(aa), 1e53);
  BOOST_CHECK_LT(double(aa), 1e55);
  // it exceeds range of float, so +Inf
  BOOST_CHECK_EQUAL(float(aa), 1.0f/0.0f);
  
  BOOST_CHECK_LT(double(-aa), -1e53);
  BOOST_CHECK_GT(double(-aa), -1e55);
  BOOST_CHECK_EQUAL(float(-aa), -1.0f/0.0f);
  
  BOOST_CHECK_EQUAL(biguint<128>(-1e20), biguint<128>(0));
  BOOST_CHECK_EQUAL(biguint<128>(-1.0), biguint<128>(0));
  BOOST_CHECK_EQUAL(biguint<128>(-0.0), biguint<128>(0));
  BOOST_CHECK_LT(bigint<128>(-1e20), bigint<128>(-1) << 66);
  BOOST_CHECK_GT(bigint<128>(-1e20), bigint<128>(-1) << 67);
  BOOST_CHECK_LT(bigint<128>(-1e20f), bigint<128>(-1) << 66);
  BOOST_CHECK_GT(bigint<128>(-1e20f), bigint<128>(-1) << 67);
  
  BOOST_CHECK_EQUAL(isqrt(aa), a1*a1*a1);
  BOOST_CHECK_EQUAL(isqrt(aa+1), a1*a1*a1);
  BOOST_CHECK_EQUAL(isqrt(aa-1), a1*a1*a1 - 1);
  BOOST_CHECK_EQUAL(isqrt(aa + a1*a1*a1*2), a1*a1*a1);
  BOOST_CHECK_EQUAL(isqrt(aa + a1*a1*a1*2 + 1), a1*a1*a1 + 1);
  BOOST_CHECK_EQUAL(isqrt(aa + a1*a1*a1*2 + 2), a1*a1*a1 + 1);
  BOOST_CHECK_EQUAL(isqrt(biguint<128>(0)), biguint<128>(0));
  BOOST_CHECK_EQUAL(isqrt(biguint<128>(1)), biguint<128>(1));
  BOOST_CHECK_EQUAL(isqrt(biguint<128>(2)), biguint<128>(1));
  BOOST_CHECK_EQUAL(isqrt(biguint<128>(3)), biguint<128>(1));
  BOOST_CHECK_EQUAL(isqrt(biguint<128>(4)), biguint<128>(2));
  BOOST_CHECK_EQUAL(isqrt(biguint<128>(5)), biguint<128>(2));

  BOOST_CHECK_EQUAL(ilog2(biguint<128>(1)), biguint<128>(0));
  BOOST_CHECK_EQUAL(ilog2(biguint<128>(2)), biguint<128>(1));
  BOOST_CHECK_EQUAL(ilog2(biguint<128>(3)), biguint<128>(1));
  BOOST_CHECK_EQUAL(ilog2(biguint<128>(4)), biguint<128>(2));
  BOOST_CHECK_EQUAL(ilog2(biguint<128>(5)), biguint<128>(2));
  BOOST_CHECK_EQUAL(ilog2(aa), 179);
  
  //BOOST_CHECK(biguint<128>(reciprocal_unsigned(biguint<128>(13)).reciprocal));
  BOOST_CHECK_EQUAL(biguint<128>(37) / biguint<128>(13), 2);
  BOOST_CHECK_EQUAL(biguint<128>(38) / biguint<128>(13), 2);
  BOOST_CHECK_EQUAL(biguint<128>(39) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(40) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(41) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(42) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(43) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(44) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(45) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(46) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(47) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(48) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(49) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(50) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(51) / biguint<128>(13), 3);
  BOOST_CHECK_EQUAL(biguint<128>(52) / biguint<128>(13), 4);
  BOOST_CHECK_EQUAL(biguint<128>(53) / biguint<128>(13), 4);
  BOOST_CHECK_EQUAL(biguint<128>(54) / biguint<128>(13), 4);

  BOOST_CHECK_EQUAL(biguint<128>(0) / biguint<128>(3), 0);
  BOOST_CHECK_EQUAL(biguint<128>(1) / biguint<128>(3), 0);
  BOOST_CHECK_EQUAL(biguint<128>(2) / biguint<128>(3), 0);
  BOOST_CHECK_EQUAL(biguint<128>(3) / biguint<128>(3), 1);
  BOOST_CHECK_EQUAL(biguint<128>(4) / biguint<128>(3), 1);
  BOOST_CHECK_EQUAL(biguint<128>(5) / biguint<128>(3), 1);
  BOOST_CHECK_EQUAL(biguint<128>(6) / biguint<128>(3), 2);
  BOOST_CHECK_EQUAL(biguint<128>(7) / biguint<128>(3), 2);

  BOOST_CHECK_EQUAL(biguint<128>(0) / biguint<128>(6), 0);
  BOOST_CHECK_EQUAL(biguint<128>(1) / biguint<128>(6), 0);
  BOOST_CHECK_EQUAL(biguint<128>(2) / biguint<128>(6), 0);
  BOOST_CHECK_EQUAL(biguint<128>(3) / biguint<128>(6), 0);
  BOOST_CHECK_EQUAL(biguint<128>(4) / biguint<128>(6), 0);
  BOOST_CHECK_EQUAL(biguint<128>(5) / biguint<128>(6), 0);
  BOOST_CHECK_EQUAL(biguint<128>(6) / biguint<128>(6), 1);
  BOOST_CHECK_EQUAL(biguint<128>(7) / biguint<128>(6), 1);
  BOOST_CHECK_EQUAL(biguint<128>(8) / biguint<128>(6), 1);
  BOOST_CHECK_EQUAL(biguint<128>(9) / biguint<128>(6), 1);
  BOOST_CHECK_EQUAL(biguint<128>(10) / biguint<128>(6), 1);
  BOOST_CHECK_EQUAL(biguint<128>(11) / biguint<128>(6), 1);
  BOOST_CHECK_EQUAL(biguint<128>(12) / biguint<128>(6), 2);
  BOOST_CHECK_EQUAL(biguint<128>(13) / biguint<128>(6), 2);

  BOOST_CHECK_EQUAL((biguint<128>(1) << 32) / biguint<128>(3), 0x55555555);
  BOOST_CHECK_EQUAL((biguint<128>(1) << 31) / biguint<128>(3), 0x2aaaaaaa);
  BOOST_CHECK_EQUAL((biguint<128>(1) << 64) / biguint<128>(3), 0x5555555555555555);
  BOOST_CHECK_EQUAL((biguint<128>(1) << 124) / biguint<128>(3), (biguint<128>(0x0555555555555555)<<64) + 0x5555555555555555);
  BOOST_CHECK_EQUAL((biguint<128>(1) << 125) / biguint<128>(3), (biguint<128>(0x0aaaaaaaaaaaaaaa)<<64) + 0xaaaaaaaaaaaaaaaa);
  BOOST_CHECK_EQUAL((biguint<128>(1) << 126) / biguint<128>(3), (biguint<128>(0x1555555555555555)<<64) + 0x5555555555555555);
  BOOST_CHECK_EQUAL((biguint<128>(1) << 127) / biguint<128>(3), (biguint<128>(0x2aaaaaaaaaaaaaaa)<<64) + 0xaaaaaaaaaaaaaaaa);//

  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) - 0) / ((biguint<128>(1) << 32) - 0), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) - 1) / ((biguint<128>(1) << 32) - 0), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) - 1) / ((biguint<128>(1) << 32) - 1), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) - 2) / ((biguint<128>(1) << 32) - 1), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) + 0) / ((biguint<128>(1) << 32) + 0), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) + 0) / ((biguint<128>(1) << 32) + 1), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) + 1) / ((biguint<128>(1) << 32) + 1), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 32) + 1) / ((biguint<128>(1) << 32) + 2), 0);

  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) - 0) / ((biguint<128>(1) << 63) - 0), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) - 1) / ((biguint<128>(1) << 63) - 0), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) - 1) / ((biguint<128>(1) << 63) - 1), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) - 2) / ((biguint<128>(1) << 63) - 1), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) + 0) / ((biguint<128>(1) << 63) + 0), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) + 0) / ((biguint<128>(1) << 63) + 1), 0);//
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) + 1) / ((biguint<128>(1) << 63) + 1), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 63) + 1) / ((biguint<128>(1) << 63) + 2), 0);

  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) - 0) / ((biguint<128>(1) << 64) - 0), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) - 1) / ((biguint<128>(1) << 64) - 0), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) - 1) / ((biguint<128>(1) << 64) - 1), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) - 2) / ((biguint<128>(1) << 64) - 1), 0);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) + 0) / ((biguint<128>(1) << 64) + 0), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) + 0) / ((biguint<128>(1) << 64) + 1), 0);//
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) + 1) / ((biguint<128>(1) << 64) + 1), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 64) + 1) / ((biguint<128>(1) << 64) + 2), 0);

  BOOST_CHECK_EQUAL(biguint<128>(0) / biguint<128>(2), 0);
  BOOST_CHECK_EQUAL(biguint<128>(1) / biguint<128>(2), 0);
  BOOST_CHECK_EQUAL(biguint<128>(2) / biguint<128>(2), 1);
  BOOST_CHECK_EQUAL(biguint<128>(3) / biguint<128>(2), 1);
  BOOST_CHECK_EQUAL(biguint<128>(4) / biguint<128>(2), 2);
  BOOST_CHECK_EQUAL(biguint<128>(5) / biguint<128>(2), 2);
  BOOST_CHECK_EQUAL(biguint<128>(6) / biguint<128>(2), 3);

  BOOST_CHECK_THROW(biguint<128>(0) / biguint<128>(0), std::logic_error);
  BOOST_CHECK_THROW(biguint<128>(1) / biguint<128>(0), std::logic_error);

  BOOST_CHECK_EQUAL(biguint<128>(0) / biguint<128>(1), 0);
  BOOST_CHECK_EQUAL(biguint<128>(1) / biguint<128>(1), 1);
  BOOST_CHECK_EQUAL(biguint<128>(2) / biguint<128>(1), 2);
  BOOST_CHECK_EQUAL(biguint<128>(3) / biguint<128>(1), 3);
  BOOST_CHECK_EQUAL(biguint<128>(4) / biguint<128>(1), 4);
  BOOST_CHECK_EQUAL(biguint<128>(5) / biguint<128>(1), 5);
  BOOST_CHECK_EQUAL(biguint<128>(6) / biguint<128>(1), 6);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 128) - 1) / biguint<128>(1), ((biguint<128>(1) << 128) - 1));
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 128) - 2) / biguint<128>(1), ((biguint<128>(1) << 128) - 2));

  BOOST_CHECK_EQUAL(((biguint<128>(1) << 120) - 1) / ((biguint<128>(1) << 120) - 3), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 120) - 3) / ((biguint<128>(1) << 120) - 3), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 120) - 4) / ((biguint<128>(1) << 120) - 3), 0);

  BOOST_CHECK_EQUAL(((biguint<128>(1) << 127) - 1) / ((biguint<128>(1) << 127) - 3), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 127) - 3) / ((biguint<128>(1) << 127) - 3), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 127) - 4) / ((biguint<128>(1) << 127) - 3), 0);

  BOOST_CHECK_EQUAL(((biguint<128>(1) << 128) - 1) / ((biguint<128>(1) << 128) - 3), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 128) - 3) / ((biguint<128>(1) << 128) - 3), 1);
  BOOST_CHECK_EQUAL(((biguint<128>(1) << 128) - 4) / ((biguint<128>(1) << 128) - 3), 0);

  BOOST_CHECK(+aa);
  BOOST_CHECK(-aa);
  BOOST_CHECK_EQUAL(-(-aa), aa);
  BOOST_CHECK_NE(-aa, aa);
  BOOST_CHECK_EQUAL(aa, aa);
  BOOST_CHECK_EQUAL(aa-aa, decltype(aa)(0));
  BOOST_CHECK_NE(aa+aa, aa);
  BOOST_CHECK_NE(aa+aa, decltype(aa)(0));
  BOOST_CHECK_LT(decltype(aa)(-1), decltype(aa)(0));
  BOOST_CHECK(!(aa-aa));
  BOOST_CHECK(!-(aa-aa));
  BOOST_CHECK(~(aa-aa));
  BOOST_CHECK(!(aa+(-aa)));
  BOOST_CHECK(!((-aa)+aa));
  BOOST_CHECK_EQUAL(-1*aa, -aa);
  BOOST_CHECK_EQUAL(aa*(aa*aa), (aa*aa)*aa);
  BOOST_CHECK_EQUAL(aa*(-aa*aa), (-aa*-aa)*-aa);
  BOOST_CHECK_EQUAL(aa*(aa+aa), aa*aa+aa*aa);
  BOOST_CHECK_EQUAL(aa*(aa-aa), aa*aa-aa*aa);
  BOOST_CHECK_EQUAL(aa*(aa+(-aa)), aa*aa-aa*aa);
  BOOST_CHECK_EQUAL(aa*(aa+(-aa)), aa*aa+(-aa*aa));
  BOOST_CHECK_EQUAL(aa+aa, 2*aa);
  BOOST_CHECK_EQUAL(aa+aa+aa, 3*aa);

  BOOST_CHECK_EQUAL(lossless_multiply(aa, lossless_multiply(aa,aa)), lossless_multiply(lossless_multiply(aa, aa), aa));
  BOOST_CHECK_EQUAL(lossless_multiply(aa, lossless_multiply(-aa,aa)), lossless_multiply(lossless_multiply(-aa, -aa), -aa));
  BOOST_CHECK_EQUAL(lossless_multiply(aa, aa+aa), lossless_multiply(aa, aa) + lossless_multiply(aa, aa));
  BOOST_CHECK_EQUAL(lossless_multiply(aa, aa-aa), lossless_multiply(aa, aa) - lossless_multiply(aa, aa));
/*
  aa&aa;
  aa|aa;
  aa^aa;
  ~aa;~ac;
  aa>>90;ac>>90;
  aa<<90;ac<<90;
  aa==aa;ac==ac;
  aa!=aa;ac!=ac;
  aa<aa;ac<ac;
  aa<=aa;ac<=ac;
  aa>aa;ac>ac;
  aa>=aa;ac>=ac;
  
  aa+=aa;ac+=ac;
  aa-=aa;ac-=ac;
  aa*=aa;ac*=ac;
  aa&=aa;ac&=ac;
  aa|=aa;ac|=ac;
  aa^=aa;ac^=ac;
  aa>>=90;ac>>=90;
  aa<<=90;ac<<=90;*/
}


REGISTER_TESTS // This must come last in the file.
