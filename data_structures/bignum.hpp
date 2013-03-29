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

#ifndef LASERCAKE_BIGNUM_HPP__
#define LASERCAKE_BIGNUM_HPP__

// Arbitrary fixed-size integral types.
// We need bit depths that are too large for mere uint64_t but which
// are small enough that register/stack allocation and inlining are
// more efficient than heap allocation and function calls.
// With -O3, most operations unroll to straight-line code with no
// loops or branches.

#include "../config.hpp"
#include <boost/utility/enable_if.hpp>
#include <ostream>
#include <cmath>
#include "../cxx11/hash.hpp"
#include <boost/functional/hash.hpp>
#include "numbers.hpp"

/*
template<size_t Bits>
struct fixed_size_bignum {
  static_assert(Bits % 64 == 0, "er.");
  uint64_t limbs[Bits / 64];
};*/

// can be 32ey later for small platforms

namespace bignum {

// 'limb' following GMP lib terminology
typedef uint64_t limb_type;
typedef int64_t signed_limb_type;
static const int limb_bits = 64;
typedef DETECTED_uint128_t twice_limb_type;
template<size_t Limbs>
struct bignum {
  // This array's order is little-endian, regardless of the bit order
  // within each word.
  limb_type limbs[Limbs];
};
template<size_t Limbs>
struct bignum_with_overflow {
  bignum<Limbs> num;
  bool overflow;
};
template<size_t LimbsOut, size_t LimbsA, size_t LimbsB>
inline bignum<LimbsOut> long_multiply_unsigned(bignum<LimbsA> a, bignum<LimbsB> b) {
  bignum<LimbsOut> result = {{}};
  for(size_t ia = 0; ia != LimbsA; ++ia) {
    uint32_t carry = 0;
    for(size_t ib = 0; ib != LimbsB; ++ib) {
      //short circuit if a1[i] or a2[j] (or their product?) is zero?
      const size_t ir = ia+ib;
      if(ir+1 < LimbsOut) {
        //full_t z = width_doubling_multiply(a1[i], a2[j]);
        //z[0] z[1]
        const twice_limb_type subproduct = (twice_limb_type)a.limbs[ia] * b.limbs[ib];
        const limb_type low = (limb_type)(subproduct);
        limb_type high = (limb_type)(subproduct >> limb_bits) + carry;
        const limb_type oldlimb0 = result.limbs[ir];
        const limb_type newlimb0 = oldlimb0 + low;
        // Carry.
        // This can't overflow high because e.g. 9*9 == 81 < 99, 0b11 * 0b11 == 0b1001 < 0b1111,
        // max uint32 * max uint32 < max uint64.
        high += (newlimb0 < oldlimb0);
        result.limbs[ir] = newlimb0;
        const limb_type oldlimb1 = result.limbs[ir+1];
        const limb_type newlimb1 = oldlimb1 + high;
        carry = (newlimb1 < oldlimb1);
        result.limbs[ir+1] = newlimb1;
      }
      else if(ir < LimbsOut) {
        limb_type low = a.limbs[ia] * b.limbs[ib];
        result.limbs[ir] += low;
      }
    }
    const size_t ir = ia+LimbsB+1;
    if(ir < LimbsOut) {
      // hasn't been touched yet
      assert(result.limbs[ir] == 0);
      result.limbs[ir] = carry;
    }
  }
  return result;
}

//cast_unsigned
template<size_t LimbsOut, size_t LimbsA>
inline bignum<LimbsOut> zero_extend(bignum<LimbsA> a) {
  bignum<LimbsOut> result;
  for(size_t ir = 0; ir < LimbsA && ir < LimbsOut; ++ir) {
    result.limbs[ir] = a.limbs[ir];
  }
  for(size_t ir = LimbsA; ir < LimbsOut; ++ir) {
    result.limbs[ir] = 0;
  }
  return result;
}
template<size_t LimbsOut>
inline bignum<LimbsOut> zero_extend_from_limb(limb_type a) {
  bignum<LimbsOut> result = {{}};
  if(LimbsOut > 0) {
    result.limbs[0] = a;
  }
  return result;
}

template<size_t LimbsOut>
inline bignum<LimbsOut> zero_extend_from_uint64(uint64_t a) {
  static_assert(limb_bits >= 64, "bug");
  return zero_extend_from_limb<LimbsOut>(a);
}
template<size_t LimbsOut>
inline bignum<LimbsOut> zero_extend_from_uint32(uint32_t a) {
  static_assert(limb_bits >= 32, "bug");
  return zero_extend_from_limb<LimbsOut>(a);
}


template<size_t LimbsOut, size_t LimbsA>
inline bignum<LimbsOut> sign_extend(bignum<LimbsA> a) {
  bignum<LimbsOut> result;
  const bool neg = is_negative(a);
  const limb_type sign_ext = (neg ? (limb_type)(-1) : (limb_type)(0));
  for(size_t ir = 0; ir != LimbsA && ir != LimbsOut; ++ir) {
    result.limbs[ir] = a.limbs[ir];
  }
  if(LimbsOut > LimbsA) {
    for(size_t ir = LimbsA; ir != LimbsOut; ++ir) {
      result.limbs[ir] = sign_ext;
    }
  }
  return result;
}

template<size_t LimbsOut>
inline bignum<LimbsOut> sign_extend_from_limb(signed_limb_type a) {
  return sign_extend<LimbsOut>(bignum<1>{{(limb_type)(a)}});
}

template<size_t LimbsOut>
inline bignum<LimbsOut> sign_extend_from_int64(int64_t a) {
  static_assert(limb_bits >= 64, "bug");
  return sign_extend_from_limb<LimbsOut>((limb_type)(signed_limb_type)(a));
}
template<size_t LimbsOut>
inline bignum<LimbsOut> sign_extend_from_int32(int32_t a) {
  static_assert(limb_bits >= 32, "bug");
  return sign_extend_from_limb<LimbsOut>((limb_type)(signed_limb_type)(a));
}

template<size_t LimbsOut>
inline bignum<LimbsOut> zero() {
  bignum<LimbsOut> result = {{}};
  return result;
}
template<size_t LimbsOut>
inline bignum<LimbsOut> max_unsigned() {
  bignum<LimbsOut> result;
  for(size_t ir = 0; ir != LimbsOut; ++ir) {
    result.limbs[ir] = (limb_type)(-1);
  }
  return result;
}
template<size_t LimbsOut>
inline bignum<LimbsOut> max_signed() {
  bignum<LimbsOut> result;
  for(size_t ir = 0; ir != LimbsOut-1; ++ir) {
    result.limbs[ir] = (limb_type)(-1);
  }
  result.limbs[LimbsOut-1] = ((limb_type)(-1) >> 1);
  return result;
}
template<size_t LimbsOut>
inline bignum<LimbsOut> min_signed() {
  bignum<LimbsOut> result = {{}};
  result.limbs[LimbsOut-1] = ((limb_type)(1) << (limb_bits-1));;
  return result;
}


inline bool is_negative_limb(limb_type a) {
  return a & ((limb_type)(1) << (limb_bits-1));
}
template<size_t Limbs>
inline bool is_negative(bignum<Limbs> a) {
  return a.limbs[Limbs-1] & ((limb_type)(1) << (limb_bits-1));
}
// This *could* allow Limbs+1 for carry overflow...

template<size_t Limbs>
inline bignum_with_overflow<Limbs> negate_overflow(bignum<Limbs> a) {
  bignum_with_overflow<Limbs> result;
  result.overflow = false;
  bool carry = true;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    const limb_type newlimb = ~a.limbs[ir] + (limb_type)(carry);
    //LOG << std::hex << "?<" << a.limbs[ir] << ":->" << newlimb << ">";
    result.num.limbs[ir] = newlimb;
    if(ir == Limbs-1 && carry && newlimb == ((limb_type)(1) << (limb_bits-1))) {
      result.overflow = true;
    }
    carry = (carry && (newlimb == 0));
  }
  //LOG << "\n";
  return result;
}
template<size_t Limbs>
inline bignum<Limbs> negate(bignum<Limbs> a) {
  bignum<Limbs> result = negate_overflow(a).num;
  return result;
}
template<size_t LimbsOut, size_t LimbsA, size_t LimbsB>
inline bignum<LimbsOut> long_multiply_signed(bignum<LimbsA> a, bignum<LimbsB> b) {
  static_assert(LimbsA > 0 && LimbsB > 0, "signed ints need a sign bit");
  // If a is negative
  const bool isnega = is_negative(a);
  const bool isnegb = is_negative(b);
  bignum<LimbsOut> result = {{}};
  const bignum_with_overflow<LimbsA> nega = negate_overflow(a);
  const bignum_with_overflow<LimbsB> negb = negate_overflow(b);
  // max negative two's complement int requires special treatment.
  // unspecified zeroes here are implicit from "= {}" of result above.
  if(isnega && nega.overflow) {
    for(size_t ib = 0; ib+LimbsA < LimbsOut && ib < LimbsB; ++ib) {
      const size_t ir = ib+LimbsA;
      result.limbs[ir] = negb.num.limbs[ib];
    }
  }
  else if(isnegb && negb.overflow) {
    for(size_t ia = 0; ia+LimbsB < LimbsOut && ia < LimbsA; ++ia) {
      const size_t ir = ia+LimbsB;
      result.limbs[ir] = nega.num.limbs[ia];
    }
  }
  else {
    const bignum<LimbsA> absa = (isnega ? nega.num : a);
    const bignum<LimbsB> absb = (isnegb ? negb.num : b);
    result = long_multiply_unsigned<LimbsOut>(absa, absb);
    //LOG << "Well?" << (isnega) << ','<<isnegb<<".\n";
    if(isnega != isnegb) {
      result = negate(result);
    }
  }
  return result;
}

template<size_t Limbs>
inline bignum<Limbs> add(bignum<Limbs> a, bignum<Limbs> b) {
  bignum<Limbs> result;
  bool carry = false;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    result.limbs[ir] = (limb_type)(carry);
    carry = false;
    const limb_type oldlimba = result.limbs[ir];
    const limb_type newlimba = oldlimba + a.limbs[ir];
    carry = (carry || newlimba < oldlimba);
    result.limbs[ir] = newlimba;
    const limb_type oldlimbb = result.limbs[ir];
    const limb_type newlimbb = oldlimbb + b.limbs[ir];
    carry = (carry || newlimbb < oldlimbb);
    result.limbs[ir] = newlimbb;
  }
  return result;
}
template<size_t Limbs>
inline bignum<Limbs> subtract(bignum<Limbs> a, bignum<Limbs> b) {
  bignum<Limbs> result;
  bool carry = true;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    result.limbs[ir] = (limb_type)(carry);
    carry = false;
    const limb_type oldlimba = result.limbs[ir];
    const limb_type newlimba = oldlimba + a.limbs[ir];
    carry = (carry || newlimba < oldlimba);
    result.limbs[ir] = newlimba;
    const limb_type oldlimbb = result.limbs[ir];
    const limb_type newlimbb = oldlimbb + ~b.limbs[ir];
    carry = (carry || newlimbb < oldlimbb);
    result.limbs[ir] = newlimbb;
  }
  return result;
}
template<size_t Limbs>
inline bignum<Limbs> bitwise_and(bignum<Limbs> a, bignum<Limbs> b) {
  bignum<Limbs> result;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    result.limbs[ir] = a.limbs[ir] & b.limbs[ir];
  }
  return result;
}
template<size_t Limbs>
inline bignum<Limbs> bitwise_or(bignum<Limbs> a, bignum<Limbs> b) {
  bignum<Limbs> result;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    result.limbs[ir] = a.limbs[ir] | b.limbs[ir];
  }
  return result;
}
template<size_t Limbs>
inline bignum<Limbs> bitwise_xor(bignum<Limbs> a, bignum<Limbs> b) {
  bignum<Limbs> result;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    result.limbs[ir] = a.limbs[ir] ^ b.limbs[ir];
  }
  return result;
}
template<size_t Limbs>
inline bignum<Limbs> bitwise_complement(bignum<Limbs> a) {
  bignum<Limbs> result;
  for(size_t ir = 0; ir != Limbs; ++ir) {
    result.limbs[ir] = ~a.limbs[ir];
  }
  return result;
}
template<size_t Limbs>
inline bool nonzero(bignum<Limbs> a) {
  for(size_t ir = 0; ir != Limbs; ++ir) {
    if(a.limbs[ir]) {
      return true;
    }
  }
  return false;
}
template<size_t Limbs>
inline bool equal(bignum<Limbs> a, bignum<Limbs> b) {
  for(size_t ir = 0; ir != Limbs; ++ir) {
    if(a.limbs[ir] != b.limbs[ir]) {
      return false;
    }
  }
  return true;
}
template<size_t Limbs>
inline bool less_than_unsigned(bignum<Limbs> a, bignum<Limbs> b) {
  for(size_t ir = 0; ir != Limbs; ++ir) {
    const limb_type av = a.limbs[Limbs-1-ir];
    const limb_type bv = b.limbs[Limbs-1-ir];
    if(av < bv) {
      return true;
    }
    if(av > bv) {
      return false;
    }
  }
  return false;
}
template<size_t Limbs>
inline bool less_than_signed(bignum<Limbs> a, bignum<Limbs> b) {
  for(size_t ir = 0; ir != Limbs; ++ir) {
    limb_type av = a.limbs[Limbs-1-ir];
    limb_type bv = b.limbs[Limbs-1-ir];
    if(ir == 0) {
      const limb_type adj = (limb_type)(1) << (limb_bits-1); 
      av ^= adj;
      bv ^= adj;
    }
    if(av < bv) {
      return true;
    }
    if(av > bv) {
      return false;
    }
  }
  return false;
}


template<size_t Limbs>
inline int32_t popcount(bignum<Limbs> a) {
  int32_t result = 0;
  for(size_t ia = 0; ia != Limbs; ++ia) {
    result += ::popcount(a.limbs[ia]);
  }
  return result;
}

template<size_t Limbs>
inline int32_t log2_unsigned(bignum<Limbs> a) {
  for(int ir = Limbs-1; ir >= 0; --ir) {
    if(a.limbs[ir] != 0) {
      return ilog2(a.limbs[ir]) + ir*limb_bits;
    }
  }
  caller_error("the logarithm of zero is undefined");
}
template<size_t Limbs>
inline int32_t log2_signed(bignum<Limbs> a) {
  caller_error_if(is_negative(a), "logarithm is only defined on positive numbers");
  return log2_unsigned(a);
}

template<size_t Limbs>
inline int32_t count_trailing_zeroes(bignum<Limbs> a) {
  for(int ir = 0; ir != Limbs; ++ir) {
    if(a.limbs[ir] != 0) {
      return ::count_trailing_zeroes(a.limbs[ir]) + ir*limb_bits;
    }
  }
  caller_error("the number of trailing zeroes of zero is undefined");
}

#if 0
// TODO is that a good shift argument type?
template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> shift_right_zero_extend(bignum<Limbs> a, uint32_t shift) {
  bignum<LimbsOut> result;
  size_t limb_offset = shift / limb_bits;
  size_t sub_limb_offset_0 = shift % limb_bits;
  size_t sub_limb_offset_1 = limb_bits - sub_limb_offset_0;
  for(size_t ir = 0; ir != LimbsOut; ++ir) {
    result.limbs[ir] = 0;
    if(ir+limb_offset < Limbs) {
      result.limbs[ir] |= (a.limbs[ir+limb_offset+1] >> sub_limb_offset_0);
      if(ir+limb_offset+1 < Limbs && sub_limb_offset_0 != 0) {
        result.limbs[ir] |= (a.limbs[ir+limb_offset+1] << sub_limb_offset_1);
      }
    }
  }
  return result;
}
#endif
template<bool SignExtend, size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> shift_impl(bignum<Limbs> a, ptrdiff_t limb_offset, size_t sub_limb_offset_0) {
  bignum<LimbsOut> result;
  const size_t sub_limb_offset_1 = limb_bits - sub_limb_offset_0;
  const bool neg = is_negative(a);
  const limb_type sign_ext = ((SignExtend && neg) ? (limb_type)(-1) : (limb_type)(0));
  if(sub_limb_offset_0 == 0) {
    for(ptrdiff_t ir = 0; ir != LimbsOut; ++ir) {
      const size_t ia = (size_t)(ir-limb_offset);
      // size_t is modulo, so this also checks for < 0.
      if(ia < Limbs) {
        result.limbs[ir] = a.limbs[ia];
      }
      else if(ir < limb_offset) {
        result.limbs[ir] = 0;
      }
      else {
        result.limbs[ir] = sign_ext;
      }
    }
  }
  else {
    for(ptrdiff_t ir = 0; ir != LimbsOut; ++ir) {
      const size_t ia0 = ir-limb_offset;
      const size_t ia1 = ir-limb_offset+1;
      if(ia0 < Limbs && ia1 < Limbs) {
        result.limbs[ir] =
          (a.limbs[ia0] >> sub_limb_offset_0) |
          (a.limbs[ia1] << sub_limb_offset_1);
      }
      else if(ia0 < Limbs) {
        if(SignExtend) {
          result.limbs[ir] = ((signed_limb_type)(a.limbs[ia0]) >> sub_limb_offset_0);
        }
        else {
          result.limbs[ir] = ((a.limbs[ia0]) >> sub_limb_offset_0);
        }
      }
      else if(ia1 < Limbs) {
        result.limbs[ir] = (a.limbs[ia1] << sub_limb_offset_1);
      }
      else if(ir < limb_offset) {
        result.limbs[ir] = 0;
      }
      else {
        result.limbs[ir] = sign_ext;
      }
    }
  }
  return result;
}
// TODO is that a good shift argument type?
template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> shift_left_zero_extend(bignum<Limbs> a, uint32_t shift) {
  const size_t limb_offset = shift / limb_bits;
  const size_t sub_limb_offset_0 = shift % limb_bits;
  return shift_impl<false, LimbsOut>(a, limb_offset, sub_limb_offset_0);
}
template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> shift_right_zero_extend(bignum<Limbs> a, uint32_t shift) {
  ptrdiff_t limb_offset = shift / limb_bits;
  ptrdiff_t sub_limb_offset_0 = shift % limb_bits;
  if(sub_limb_offset_0 != 0) {
    ++limb_offset;
    sub_limb_offset_0 = limb_bits - sub_limb_offset_0;
  }
  return shift_impl<false, LimbsOut>(a, -limb_offset, sub_limb_offset_0);
}
template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> shift_left_sign_extend(bignum<Limbs> a, uint32_t shift) {
  const size_t limb_offset = shift / limb_bits;
  const size_t sub_limb_offset_0 = shift % limb_bits;
  return shift_impl<false, LimbsOut>(a, limb_offset, sub_limb_offset_0);
}
template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> shift_right_sign_extend(bignum<Limbs> a, uint32_t shift) {
  ptrdiff_t limb_offset = shift / limb_bits;
  ptrdiff_t sub_limb_offset_0 = shift % limb_bits;
  if(sub_limb_offset_0 != 0) {
    ++limb_offset;
    sub_limb_offset_0 = limb_bits - sub_limb_offset_0;
  }
  return shift_impl<true, LimbsOut>(a, -limb_offset, sub_limb_offset_0);
}



template<size_t LimbsOut, size_t Limbs, size_t OperationLimbs = Limbs/2+1>
inline bignum<LimbsOut> sqrt_unsigned(bignum<Limbs> radicand) {
  typedef bignum<OperationLimbs> operation_t;
  // in fact, lower_bound and mid stay within (Limbs+1)/2 ... hmm.

  bignum<LimbsOut> result = {{}};

  // log2(0) doesn't exist, but sqrt(0) does, so we have to check for it here.
  if(!nonzero(radicand)) {
    return result;
  }

  //shift is the log base 2 of radicand, divided by two, rounded down.
  const int32_t shift = log2_unsigned(radicand) >> 1;

  //bounds are [lower_bound, upper_bound), a half-open range.
  //lower_bound is guaranteed to be less than or equal to the answer.
  //upper_bound is guaranteed to be greater than the answer.
  operation_t lower_bound = {{}};  // = half_t(1) << shift;
  lower_bound.limbs[shift / limb_bits] = ((limb_type)(1) << (shift % limb_bits));

  //upper_bound is twice the original lower_bound;
  //upper_bound is    2**(floor(log2(radicand) / 2)+1)
  //which is equal to 2**ceil((log2(radicand)+1) / 2)
  operation_t upper_bound = {{}};
  upper_bound.limbs[(shift+1) / limb_bits] = ((limb_type)(1) << ((shift+1) % limb_bits));

  /* while(lower_bound < upper_bound - 1) */
  while(less_than_unsigned(lower_bound, add(upper_bound, max_unsigned<OperationLimbs>())))
  {
    const operation_t mid = shift_impl<false, OperationLimbs>(add(upper_bound, lower_bound), -1, limb_bits-1) /* >> 1 */;
    if(less_than_unsigned(radicand, long_multiply_unsigned<Limbs>(mid, mid))) {
      upper_bound = mid;
    }
    else {
      lower_bound = mid;
    }
  }
  
  result = zero_extend<LimbsOut>(lower_bound);
  return result;
}

template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> sqrt_signed(bignum<Limbs> radicand) {
  caller_error_if(is_negative(radicand), "sqrt of a negative number");
  // We know there's one bit free now (the sign bit), which
  // is enough to fit that awkward initial upper_bound in.
  // Optimization!
  return sqrt_unsigned<LimbsOut, Limbs, ((Limbs+1)/2)>(radicand);
}


#if 0
// problem: sign extension whetherness.
template<size_t LimbsOut, size_t LimbsA, size_t LimbsB>
inline bignum<LimbsOut> add(bignum<LimbsA> a, bignum<LimbsB> b) {
  bignum<LimbsOut> result = {{}};
  bool carry = false;
  for(size_t ir = 0; ir != LimbsOut; ++ir) {
    result.limbs[ir] = (limb_type)(carry);
    carry = false;
    if(ir < LimbsA) {
      const limb_type oldlimb = result.limbs[ir];
      const limb_type newlimb = oldlimb + a.limbs[ir];
      carry = (carry || newlimb < oldlimb);
    }
    if(ir < LimbsB) {
      const limb_type oldlimb = result.limbs[ir];
      const limb_type newlimb = oldlimb + b.limbs[ir];
      carry = (carry || newlimb < oldlimb);
    }
  }
  return result;
}
#endif
//TODO subtract, bitwise


template<size_t Limbs>
float to_float_unsigned(bignum<Limbs> a) {
  static_assert(limb_bits >= 32, "bug");
  int ir;
  for(ir = Limbs-1; ir > 0; --ir) {
    if(a.limbs[ir] != 0) {
      break;
    }
  }
  if(ir == 0) {
    return (float)(a.limbs[0]);
  }
  else {
    return
      ldexpf((float)(a.limbs[ir]), ir*limb_bits) +
      ldexpf((float)(a.limbs[ir-1]), (ir-1)*limb_bits);
  }
}
template<size_t Limbs>
double to_double_unsigned(bignum<Limbs> a) {
  static_assert(limb_bits >= 64, "bug");
  int ir;
  for(ir = Limbs-1; ir > 0; --ir) {
    if(a.limbs[ir] != 0) {
      break;
    }
  }
  if(ir == 0) {
    return (double)(a.limbs[0]);
  }
  else {
    return
      ldexp((double)(a.limbs[ir]), ir*limb_bits) +
      ldexp((double)(a.limbs[ir-1]), (ir-1)*limb_bits);
  }
}
template<size_t Limbs>
long double to_long_double_unsigned(bignum<Limbs> a) {
  static_assert(Limbs-Limbs, "unimplemented");
  return to_double_unsigned(a);
}

template<size_t Limbs>
float to_float_signed(bignum<Limbs> a) {
  if(is_negative(a)) { return -to_float_unsigned(negate(a)); }
  else { return to_float_unsigned(a); }
}
template<size_t Limbs>
double to_double_signed(bignum<Limbs> a) {
  if(is_negative(a)) { return -to_double_unsigned(negate(a)); }
  else { return to_double_unsigned(a); }
}
template<size_t Limbs>
long double to_long_double_signed(bignum<Limbs> a) {
  if(is_negative(a)) { return -to_long_double_unsigned(negate(a)); }
  else { return to_long_double_unsigned(a); }
}

template<size_t LimbsOut, size_t Limbs>
inline bignum<LimbsOut> divide_by_limb_unsigned(bignum<Limbs> a, limb_type b) {
  
}
template<size_t Limbs>
inline void show_limbs_hex_bigendian(bignum<Limbs> a, char out[Limbs*(limb_bits/4 + 1)]) {
  size_t i = 0;
  for(size_t ia = 0; ia != Limbs; ++ia) {
    limb_type val = a.limbs[Limbs-1-ia];
    //LOG <<std::hex<< "?" << val << "?";
    for(int off = limb_bits - 4; off >= 0; off -= 4) {
      const char nybble = (char)((val >> off) & 0xf);
      out[i++] = ((nybble < 0xa) ? ('0'+nybble) : ('a'+nybble-0xa));
    }
    out[i++] = ' ';
  }
  out[i-1] = '\0';
}


//bounds checking could be stuck in these fairly runtime-efficiently.
//these are fairly macro-able.

//TODO increment could be made faster

template<size_t Bits> struct bigint;
template<size_t Bits> struct biguint;


template<size_t Limbs> inline size_t hash_value(bignum<Limbs> a) {
  size_t seed = 0;
  for(size_t ia = 0; ia != Limbs; ++ia) {
    boost::hash_combine(seed, a.limbs[ia]);
  }
  return seed;
}

template<size_t Bits>
struct biguint : bignum<(Bits/limb_bits)> {
  static_assert(Bits % limb_bits == 0, "Bigint only supports bit sizes that are a multiple of limb size.");
  static const size_t bignum_limbs = Bits/limb_bits;
  typedef bignum<bignum_limbs> bignum_type;

  biguint() {}

  biguint(uint64_t i) : bignum_type(zero_extend_from_uint64<bignum_limbs>(i)) {}

  explicit biguint(bignum_type i) : bignum_type(i) {}

  template<size_t OtherBits>
  biguint(biguint<OtherBits> i,
    typename boost::enable_if_c<(OtherBits < Bits)>::type* = 0
  ) : bignum_type(zero_extend<bignum_limbs>(i)) {}
  template<size_t OtherBits>
  explicit biguint(biguint<OtherBits> i,
    typename boost::enable_if_c<(OtherBits > Bits)>::type* = 0
  ) : bignum_type(zero_extend<bignum_limbs>(i)) {}

  //template<size_t OtherBits>
  //explicit biguint(bigint<OtherBits> i) : bignum_type(i) {}

  explicit operator bool()const { return nonzero(*this); }

  explicit operator float()const { return to_float_unsigned(*this); }
  explicit operator double()const { return to_double_unsigned(*this); }
  explicit operator long double()const { return to_long_double_unsigned(*this); }

friend inline biguint<Bits> operator*(biguint<Bits> a, biguint<Bits> b)
{ return biguint<Bits>(long_multiply_unsigned<biguint<Bits>::bignum_limbs>(a, b)); }

friend inline biguint<Bits> operator+(biguint<Bits> a) { return a; }
friend inline biguint<Bits> operator-(biguint<Bits> a) { return biguint<Bits>(negate(a)); }
friend inline biguint<Bits> abs(biguint<Bits> a) { return a; }
friend inline biguint<Bits> operator+(biguint<Bits> a, biguint<Bits> b) { return biguint<Bits>(add(a, b)); }
friend inline biguint<Bits> operator-(biguint<Bits> a, biguint<Bits> b) { return biguint<Bits>(subtract(a, b)); }

friend inline biguint<Bits> operator&(biguint<Bits> a, biguint<Bits> b) { return biguint<Bits>(bitwise_and(a, b)); }
friend inline biguint<Bits> operator|(biguint<Bits> a, biguint<Bits> b) { return biguint<Bits>(bitwise_or(a, b)); }
friend inline biguint<Bits> operator^(biguint<Bits> a, biguint<Bits> b) { return biguint<Bits>(bitwise_xor(a, b)); }
friend inline biguint<Bits> operator~(biguint<Bits> a) { return biguint<Bits>(bitwise_complement(a)); }

friend inline biguint<Bits> operator<<(biguint<Bits> a, uint32_t shift)
{ return biguint<Bits>(shift_left_zero_extend<biguint<Bits>::bignum_limbs>(a, shift)); }
friend inline biguint<Bits> operator>>(biguint<Bits> a, uint32_t shift)
{ return biguint<Bits>(shift_right_zero_extend<biguint<Bits>::bignum_limbs>(a, shift)); }

friend inline bool operator==(biguint<Bits> a, biguint<Bits> b) { return equal(a, b); }
friend inline bool operator!=(biguint<Bits> a, biguint<Bits> b) { return !equal(a, b); }
friend inline bool operator<(biguint<Bits> a, biguint<Bits> b) { return less_than_unsigned(a, b); }
friend inline bool operator>(biguint<Bits> a, biguint<Bits> b) { return less_than_unsigned(b, a); }
friend inline bool operator>=(biguint<Bits> a, biguint<Bits> b) { return !less_than_unsigned(a, b); }
friend inline bool operator<=(biguint<Bits> a, biguint<Bits> b) { return !less_than_unsigned(b, a); }


friend inline biguint<Bits>& operator+=(biguint<Bits>& a, biguint<Bits> b) { a = a + b; return a; }
friend inline biguint<Bits>& operator-=(biguint<Bits>& a, biguint<Bits> b) { a = a - b; return a; }
friend inline biguint<Bits>& operator*=(biguint<Bits>& a, biguint<Bits> b) { a = a * b; return a; }
friend inline biguint<Bits>& operator&=(biguint<Bits>& a, biguint<Bits> b) { a = a & b; return a; }
friend inline biguint<Bits>& operator|=(biguint<Bits>& a, biguint<Bits> b) { a = a | b; return a; }
friend inline biguint<Bits>& operator^=(biguint<Bits>& a, biguint<Bits> b) { a = a ^ b; return a; }
friend inline biguint<Bits>& operator<<=(biguint<Bits>& a, uint32_t shift) { a = a << shift; return a; }
friend inline biguint<Bits>& operator>>=(biguint<Bits>& a, uint32_t shift) { a = a >> shift; return a; }

friend inline biguint<Bits>& operator++(biguint<Bits>& a) { a = a + 1; return a; }
friend inline biguint<Bits> operator++(biguint<Bits>& a, int) { biguint<Bits> old = a; a = a + 1; return old; }
friend inline biguint<Bits>& operator--(biguint<Bits>& a) { a = a - 1; return a; }
friend inline biguint<Bits> operator--(biguint<Bits>& a, int) { biguint<Bits> old = a; a = a - 1; return old; }

friend inline biguint<((bigint<Bits>::bignum_limbs+1)/2)*limb_bits> isqrt(biguint<Bits> a)
{ return biguint<((bigint<Bits>::bignum_limbs+1)/2)*limb_bits>(sqrt_unsigned<((bigint<Bits>::bignum_limbs+1)/2)>(a)); }
friend inline int32_t ilog2(biguint<Bits> a) { return log2_unsigned(a); }

friend std::ostream& operator<<(std::ostream& os, biguint<Bits> a) { return os; }
};



template<size_t Bits>
struct bigint : bignum<(Bits/limb_bits)> {
  static_assert(Bits % limb_bits == 0, "Bigint only supports bit sizes that are a multiple of limb size.");
  static const size_t bignum_limbs = Bits/limb_bits;
  typedef bignum<bignum_limbs> bignum_type;

  bigint() {}

  bigint(int64_t i) : bignum_type(sign_extend_from_int64<bignum_limbs>(i)) {}

  explicit bigint(bignum_type i) : bignum_type(i) {}

  template<size_t OtherBits>
  bigint(bigint<OtherBits> i,
    typename boost::enable_if_c<(OtherBits < Bits)>::type* = 0
  ) : bignum_type(sign_extend<bignum_limbs>(i)) {}
  template<size_t OtherBits>
  explicit bigint(bigint<OtherBits> i,
    typename boost::enable_if_c<(OtherBits > Bits)>::type* = 0
  ) : bignum_type(sign_extend<bignum_limbs>(i)) {}

  explicit operator bool()const { return nonzero(*this); }

  explicit operator float()const { return to_float_signed(*this); }
  explicit operator double()const { return to_double_signed(*this); }
  explicit operator long double()const { return to_long_double_signed(*this); }

friend inline bigint<Bits> operator*(bigint<Bits> a, bigint<Bits> b)
{ return bigint<Bits>(long_multiply_signed<bigint<Bits>::bignum_limbs>(a, b)); }
friend inline bigint<Bits> operator+(bigint<Bits> a) { return a; }
friend inline bigint<Bits> operator-(bigint<Bits> a) { return bigint<Bits>(negate(a)); }
friend inline bigint<Bits> abs(bigint<Bits> a) { return (a < 0 ? -a : a); }
friend inline bigint<Bits> operator+(bigint<Bits> a, bigint<Bits> b) { return bigint<Bits>(add(a, b)); }
friend inline bigint<Bits> operator-(bigint<Bits> a, bigint<Bits> b) { return bigint<Bits>(subtract(a, b)); }

friend inline bigint<Bits> operator&(bigint<Bits> a, bigint<Bits> b) { return bigint<Bits>(bitwise_and(a, b)); }
friend inline bigint<Bits> operator|(bigint<Bits> a, bigint<Bits> b) { return bigint<Bits>(bitwise_or(a, b)); }
friend inline bigint<Bits> operator^(bigint<Bits> a, bigint<Bits> b) { return bigint<Bits>(bitwise_xor(a, b)); }
friend inline bigint<Bits> operator~(bigint<Bits> a) { return bigint<Bits>(bitwise_complement(a)); }

friend inline bool operator==(bigint<Bits> a, bigint<Bits> b) { return equal(a, b); }
friend inline bool operator!=(bigint<Bits> a, bigint<Bits> b) { return !equal(a, b); }
friend inline bool operator<(bigint<Bits> a, bigint<Bits> b) { return less_than_signed(a, b); }
friend inline bool operator>(bigint<Bits> a, bigint<Bits> b) { return less_than_signed(b, a); }
friend inline bool operator>=(bigint<Bits> a, bigint<Bits> b) { return !less_than_signed(a, b); }
friend inline bool operator<=(bigint<Bits> a, bigint<Bits> b) { return !less_than_signed(b, a); }

friend inline bigint<Bits> operator<<(bigint<Bits> a, uint32_t shift)
{ return bigint<Bits>(shift_left_sign_extend<bigint<Bits>::bignum_limbs>(a, shift)); }
friend inline bigint<Bits> operator>>(bigint<Bits> a, uint32_t shift)
{ return bigint<Bits>(shift_right_sign_extend<bigint<Bits>::bignum_limbs>(a, shift)); }

friend inline bigint<Bits>& operator+=(bigint<Bits>& a, bigint<Bits> b) { a = a + b; return a; }
friend inline bigint<Bits>& operator-=(bigint<Bits>& a, bigint<Bits> b) { a = a - b; return a; }
friend inline bigint<Bits>& operator*=(bigint<Bits>& a, bigint<Bits> b) { a = a * b; return a; }
friend inline bigint<Bits>& operator&=(bigint<Bits>& a, bigint<Bits> b) { a = a & b; return a; }
friend inline bigint<Bits>& operator|=(bigint<Bits>& a, bigint<Bits> b) { a = a | b; return a; }
friend inline bigint<Bits>& operator^=(bigint<Bits>& a, bigint<Bits> b) { a = a ^ b; return a; }
friend inline bigint<Bits>& operator<<=(bigint<Bits>& a, uint32_t shift) { a = a << shift; return a; }
friend inline bigint<Bits>& operator>>=(bigint<Bits>& a, uint32_t shift) { a = a >> shift; return a; }

friend inline bigint<Bits>& operator++(bigint<Bits>& a) { a = a + 1; return a; }
friend inline bigint<Bits> operator++(bigint<Bits>& a, int) { bigint<Bits> old = a; a = a + 1; return old; }
friend inline bigint<Bits>& operator--(bigint<Bits>& a) { a = a - 1; return a; }
friend inline bigint<Bits> operator--(bigint<Bits>& a, int) { bigint<Bits> old = a; a = a - 1; return old; }

//Imagine taking the sqrt of a signed 8 bit value.  Can it fit into 4 bits?
// sqrt(127) is 11.something, which is less than 15 (unsigned 4bit max) but
// greater than 7 (signed 4bit max).  So we conservatively leave enough space here.
friend inline bigint<(bigint<Bits>::bignum_limbs/2+1)*limb_bits> isqrt(bigint<Bits> a)
{ return bigint<(bigint<Bits>::bignum_limbs/2+1)*limb_bits>(sqrt_signed<(bigint<Bits>::bignum_limbs/2+1)>(a)); }
friend inline int32_t ilog2(bigint<Bits> a) { return log2_signed(a); }

friend std::ostream& operator<<(std::ostream& os, bigint<Bits> a) {
  char out[bigint<Bits>::bignum_limbs*(limb_bits/4 + 1)];
  show_limbs_hex_bigendian(a, out);
  os << out;
  return os;
}
};

// Or templatize on both bits and take the max?
// Might be faster depending on inlining, but requires mixed operators with non-bigint
// types to be explicitly defined.
/*
template<size_t Bits> inline biguint<Bits> operator*(biguint<Bits> a, biguint<Bits> b)
{ return long_multiply_unsigned<biguint<Bits>::bignum_limbs>(a, b); }
template<size_t Bits> inline bigint<Bits> operator*(bigint<Bits> a, bigint<Bits> b)
{ return long_multiply_signed<biguint<Bits>::bignum_limbs>(a, b); }
*/
// For explicit multiplication-precision control:
template<size_t Bits, size_t BitsA, size_t BitsB> inline biguint<Bits>
multiply_to(biguint<BitsA> a, biguint<BitsB> b) { return biguint<Bits>(long_multiply_unsigned<biguint<Bits>::bignum_limbs>(a, b)); }
template<size_t Bits, size_t BitsA, size_t BitsB> inline bigint<Bits>
multiply_to(bigint<BitsA> a, bigint<BitsB> b) { return bigint<Bits>(long_multiply_signed<biguint<Bits>::bignum_limbs>(a, b)); }
template<size_t BitsA, size_t BitsB> inline biguint<(BitsA+BitsB)>
lossless_multiply(biguint<BitsA> a, biguint<BitsB> b) { return biguint<(BitsA+BitsB)>(long_multiply_unsigned<biguint<(BitsA+BitsB)>::bignum_limbs>(a, b)); }
template<size_t BitsA, size_t BitsB> inline bigint<(BitsA+BitsB)>
lossless_multiply(bigint<BitsA> a, bigint<BitsB> b) { return bigint<(BitsA+BitsB)>(long_multiply_signed<biguint<(BitsA+BitsB)>::bignum_limbs>(a, b)); }

} /* end namespace bignum */

using bignum::biguint;
using bignum::bigint;
using bignum::multiply_to;
using bignum::lossless_multiply;

namespace std {
template<size_t Bits>
class numeric_limits< biguint<Bits> >
{
public:
  static constexpr bool is_specialized = true;
  static biguint<Bits> min() noexcept { return biguint<Bits>(bignum::zero<biguint<Bits>::bignum_limbs>()); }
  static biguint<Bits> max() noexcept { return biguint<Bits>(bignum::max_unsigned<biguint<Bits>::bignum_limbs>()); }
  static biguint<Bits> lowest() noexcept { return biguint<Bits>(bignum::zero<biguint<Bits>::bignum_limbs>()); }

  static constexpr int  digits = Bits;
  static constexpr int  digits10 = Bits * std::log10(2);
  static constexpr bool is_signed = false;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr int  radix = 2;

  static constexpr bool is_iec559 = false;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = true;

  static constexpr bool traps = false;
};
template<size_t Bits>
class numeric_limits< bigint<Bits> >
{
public:
  static constexpr bool is_specialized = true;
  static bigint<Bits> min() noexcept { return bigint<Bits>(bignum::min_signed<bigint<Bits>::bignum_limbs>()); }
  static bigint<Bits> max() noexcept { return bigint<Bits>(bignum::max_signed<bigint<Bits>::bignum_limbs>()); }
  static bigint<Bits> lowest() noexcept { return bigint<Bits>(bignum::min_signed<bigint<Bits>::bignum_limbs>()); }

  static constexpr int  digits = Bits-1;
  static constexpr int  digits10 = (Bits-1) * std::log10(2);
  static constexpr bool is_signed = true;
  static constexpr bool is_integer = true;
  static constexpr bool is_exact = true;
  static constexpr int  radix = 2;

  static constexpr bool is_iec559 = false;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = true;

  static constexpr bool traps = false;
};
}

namespace boost {
  template<size_t Bits> struct   make_signed<  bigint<Bits> > { typedef  bigint<Bits> type; };
  template<size_t Bits> struct   make_signed< biguint<Bits> > { typedef  bigint<Bits> type; };
  template<size_t Bits> struct make_unsigned<  bigint<Bits> > { typedef biguint<Bits> type; };
  template<size_t Bits> struct make_unsigned< biguint<Bits> > { typedef biguint<Bits> type; };
}

namespace std {
template<size_t Bits> inline bigint<Bits> abs(bigint<Bits> a) { return (a < 0) ? -a : a; }
template<size_t Bits> inline biguint<Bits> abs(biguint<Bits> a) { return a; }
}

namespace HASH_NAMESPACE {
  template<size_t Bits>
  struct hash< bigint<Bits> > {
    inline size_t operator()(bigint<Bits> a)const {
      return hash_value(a);
    }
  };
  template<size_t Bits>
  struct hash< biguint<Bits> > {
    inline size_t operator()(bigint<Bits> a)const {
      return hash_value(a);
    }
  };
}

#if 0
// gcc explorer tests
#include <stdint.h>
#include <stdlib.h>
using namespace std;
typedef __uint128_t DETECTED_uint128_t;

extern void g(bignum<6> const&);
extern int q,p,r,s,t,u,v,w;
void f() { g(long_multiply<6>(bignum<4>{q,p,t,u},bignum<4>{r,s,v,w})); }
#endif

//lossless_multiply
//multiply_to_fit

//division: FORCEINLINE to see likely shared reciprocal?

#endif

