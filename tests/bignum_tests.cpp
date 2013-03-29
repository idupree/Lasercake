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
}

BOOST_AUTO_TEST_CASE(bignum_runtests) {
  //bigint<128> aa = bigint<64>(1000000000LL);
  //aa = aa*aa*aa*(aa*(aa*aa));
  bigint<512> aa = bigint<64>(1000000000LL);
  aa = aa*aa*aa*(aa*(aa*aa));
  //LOG << aa << "\n";
  //LOG << decltype(aa)(2) << "\n";
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
