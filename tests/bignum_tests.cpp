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
  multiply_to_fit(aa,aa);
  multiply_to_fit(ac,ac);
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
}


REGISTER_TESTS // This must come last in the file.
