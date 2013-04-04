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

#ifndef LASERCAKE_COMPILE_TIME_BIGNUM_HPP__
#define LASERCAKE_COMPILE_TIME_BIGNUM_HPP__

#include <boost/type_traits/conditional.hpp>

// Compile-time arbitrary-precision naturals, integers, and rational numbers.

// ct = compile-time
namespace ct {

// Data types:

// or million since that's simpler for multiplication and also cross-language-dialect
// little-endian:/, base 1 billion (10^9), negatives have all components negated,
// high zeroes non-canonical
// ct = compile-time
// 10-based to make compile errors somewhat readable
// zero: ?

typedef uint64_t milliodigit;
static const milliodigit base = 1000000;
typedef int32_t shift_type;

// for the sake of operator overloading, derive
// these classes from a common base using CRTP
// (Curiously Recurring Template Pattern: base takes derived as a
// template argument).
template<typename Num>
struct number {};

// domain/universe for operators?
//hmm.

//ct: compile-time
//nat: natural number (integer >= 0)
template<milliodigit... Milliodigits>
struct nat : number<nat<Milliodigits...>> {};
// < 0
template<typename Nonnegative>
struct negative : number<negative<Nonnegative>> {};
//rational number
//normalized
// integers are nat
//positive: rational<nat<..>, nat<..>>
//zero: rational<nat<>, nat<1>>
//negative: negative<rational<nat<..>, nat<..>>>
template<typename Num, typename Den>
struct rational : number<rational<Num,Den>> {};

struct divide_by_zero {};

/*
// represents Nat * base^milliodigits_shifted_left
template<int32_t milliodigits_shifted_left, typename Nat>
struct floating {};
*/

//template<intmax_t Int> struct integer_literal;
template<typename... Num> struct add;
template<typename... Num> struct multiply;
template<typename NumA, typename NumB> struct subtract;
// Natural number subtraction (negative results become 'below_zero')
template<typename NatA, typename NatB> struct subtract_nat;
template<typename Num> struct negate;
template<typename Num, intmax_t Exponent> struct power;
// Compare: ::value is -1 if lhs < rhs; 0 if equal; 1 if lhs > rhs
template<typename NatA, typename NatB> struct compare;
template<typename Integer> struct even;
template<typename Integer> struct odd;
// Rounds down; ::type and ::quot are quotient, and ::rem is remainder.
template<typename NatA, typename NatB> struct divide_nat;
// Like divide_nat; rounds towards zero, like int operator/ does.
template<typename NatA, typename NatB> struct divide_integer;
// Exact
template<typename RatA, typename RatB> struct divide_rational;
template<typename Rat> struct reciprocal;

template<typename Num> struct round_down_to_nat;


// add and multiply are associative
template<typename NatA, typename NatB, typename NatC, typename...Nat>
struct add<NatA, NatB, NatC, Nat...> : add<typename add<NatA, NatB>::type, NatC, Nat...> {};
template<typename NatA, typename NatB, typename NatC, typename...Nat>
struct multiply<NatA, NatB, NatC, Nat...> : multiply<typename multiply<NatA, NatB>::type, NatC, Nat...> {};




///////////////////////////
//// Basic natural-number operations
///////////////////////////


template<milliodigit... Milliodigits> struct round_down_to_nat<nat<Milliodigits...>> {
  typedef nat<Milliodigits...> type;
};


template<typename Rest, milliodigit Least>
struct make_nat_from_rest_and_least;
template<milliodigit... Milliodigits, milliodigit Least>
struct make_nat_from_rest_and_least<nat<Milliodigits...>, Least> {
  typedef nat<Least, Milliodigits...> type;
};
template<>
struct make_nat_from_rest_and_least<nat<>, 0> {
  typedef nat<> type;
};
// useful for subtract
struct below_zero {};
template<milliodigit Least>
struct make_nat_from_rest_and_least<below_zero, Least> {
  typedef below_zero type;
};

template<milliodigit Milliodigit>
struct make_nat_from_milliodigit {
  static_assert(Milliodigit < base, "'milliodigit' too high");
  typedef nat<Milliodigit> type;
};
template<>
struct make_nat_from_milliodigit<0> {
  typedef nat<> type;
};

template<uintmax_t Int> struct uinteger_literal {
  // TODO do this more generically
  // This supports 64 bit literals
  typedef
    typename make_nat_from_rest_and_least<
      typename make_nat_from_rest_and_least<
        typename make_nat_from_rest_and_least<
          typename make_nat_from_rest_and_least<
            nat<>,
            (Int/base/base/base%base)>::type,
          (Int/base/base%base)>::type,
        (Int/base%base)>::type,
      (Int%base)>::type type;
};
/*
template<intmax_t Int> struct integer_literal
  : boost::conditional<(Int < 0),
      negative<typename uinteger_literal<(-Int)>::type>,
      typename uinteger_literal<Int>::type> {};
      */
template<typename Int, Int Value, bool Negative = (Value < 0)> struct literal;
template<typename Int, Int Value> struct literal<Int, Value, true>
  : uinteger_literal<static_cast<uintmax_t>(-Value)> {};
template<typename Int, Int Value> struct literal<Int, Value, false>
  : uinteger_literal<static_cast<uintmax_t>(Value)> {};

// Natural number addition
// implemented by recursive, low-to-high addition on the digits
template<milliodigit... MilliodigitsA, milliodigit MilliodigitA0,
         milliodigit... MilliodigitsB, milliodigit MilliodigitB0>
struct add<nat<MilliodigitA0, MilliodigitsA...>,
           nat<MilliodigitB0, MilliodigitsB...>>
  : make_nat_from_rest_and_least<
    typename add<nat<MilliodigitsA...>,
      typename add<nat<MilliodigitsB...>,
        typename make_nat_from_milliodigit<
          ((MilliodigitA0 + MilliodigitB0) / base)>::type>::type>::type,
    ((MilliodigitA0 + MilliodigitB0) % base)> {};
// recursion base cases
template<typename NatA> struct add<NatA, nat<>> { typedef NatA type; };
template<typename NatB> struct add<nat<>, NatB> { typedef NatB type; };
template<> struct add<nat<>, nat<>> { typedef nat<> type; };

// Natural number subtraction (negative results become 'below_zero')
// implemented by recursive, low-to-high subtraction on the digits
template<typename NatA, typename NatB, bool AIsBigger> struct subtract_nat_aux;
template<milliodigit... MilliodigitsA, milliodigit MilliodigitA0,
         milliodigit... MilliodigitsB, milliodigit MilliodigitB0>
struct subtract_nat<nat<MilliodigitA0, MilliodigitsA...>,
                    nat<MilliodigitB0, MilliodigitsB...>>
  : subtract_nat_aux<
      nat<MilliodigitA0, MilliodigitsA...>,
      nat<MilliodigitB0, MilliodigitsB...>,
      (MilliodigitA0 >= MilliodigitB0)> {};
template<milliodigit... MilliodigitsA, milliodigit MilliodigitA0,
         milliodigit... MilliodigitsB, milliodigit MilliodigitB0>
struct subtract_nat_aux<nat<MilliodigitA0, MilliodigitsA...>,
                        nat<MilliodigitB0, MilliodigitsB...>,
                        true>
  : make_nat_from_rest_and_least<
        typename subtract_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...> >::type,
        (MilliodigitA0 - MilliodigitB0)> {};
template<milliodigit... MilliodigitsA, milliodigit MilliodigitA0,
         milliodigit... MilliodigitsB, milliodigit MilliodigitB0>
struct subtract_nat_aux<nat<MilliodigitA0, MilliodigitsA...>,
                        nat<MilliodigitB0, MilliodigitsB...>,
                        false>
  : make_nat_from_rest_and_least<
        typename subtract_nat<nat<MilliodigitsA...>,
          typename add<nat<MilliodigitsB...>, nat<1> >::type>::type,
        (base - (MilliodigitB0 - MilliodigitA0))> {};
// recursion base cases
template<typename NatA> struct subtract_nat<NatA, nat<>> { typedef NatA type; };
template<typename NatB> struct subtract_nat<nat<>, NatB> { typedef below_zero type; };
template<> struct subtract_nat<nat<>, nat<>> { typedef nat<> type; };

// Natural number multiplication
// implemented by recursively calling a
// recursive multiply-by-one-digit 'multiply_by_milliodigit'
template<typename Nat, milliodigit MilliodigitB>
struct multiply_by_milliodigit;
// implemented by recursive, low-to-high multiplication of the digits
template<milliodigit... MilliodigitsA, milliodigit MilliodigitA0,
         milliodigit MilliodigitB>
struct multiply_by_milliodigit<nat<MilliodigitA0, MilliodigitsA...>,
                               MilliodigitB>
  : make_nat_from_rest_and_least<
      typename add<
        typename multiply_by_milliodigit<nat<MilliodigitsA...>, MilliodigitB>::type,
        typename make_nat_from_milliodigit<
          ((MilliodigitA0 * MilliodigitB) / base)>::type>::type,
      ((MilliodigitA0 * MilliodigitB) % base)> {};
// recursion base case
template<milliodigit MilliodigitB> struct multiply_by_milliodigit<nat<>, MilliodigitB> {
  typedef nat<> type;
};
// full multiplication
// implemented by recursive, low-to-high addition on the digits
template<milliodigit... MilliodigitsA,
         milliodigit... MilliodigitsB, milliodigit MilliodigitB0>
struct multiply<nat<MilliodigitsA...>,
                nat<MilliodigitB0, MilliodigitsB...>>
  : add<
      typename multiply_by_milliodigit<nat<MilliodigitsA...>, MilliodigitB0>::type,
      typename make_nat_from_rest_and_least<
        typename multiply<nat<MilliodigitsA...>, nat<MilliodigitsB...>>::type,
        0>::type
    > {};
// recursion base case
template<typename NatA> struct multiply<NatA, nat<>> { typedef nat<> type; };

// Natural number comparison
// -1 if lhs < rhs; 0 if equal; 1 if lhs > rhs
// lazy-programmer implementation: subtract and see what the result is like
template<typename SubtractResult> struct compare_impl : boost::integral_constant<int, 1> {};
template<> struct compare_impl<nat<>> : boost::integral_constant<int, 0> {};
template<> struct compare_impl<below_zero> : boost::integral_constant<int, (-1)> {};
template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct compare<nat<MilliodigitsA...>, nat<MilliodigitsB...>>
  : compare_impl<typename subtract_nat<nat<MilliodigitsA...>,
                                       nat<MilliodigitsB...>>::type> {};

// Exponentiation based on multiplication
template<typename Base, intmax_t Exponent, bool Negative = (Exponent < 0)>
struct power_impl1;
template<typename Base, uintmax_t Exponent, bool Odd = (Exponent % 2 == 1)>
struct power_impl2;
template<typename B> struct power_impl2<B, 0, false> { typedef nat<1> type; };
template<typename B> struct power_impl2<B, 1, true> { typedef B type; };
template<typename B> struct power_impl2<B, 2, false> : multiply<B, B> {};
template<typename B, uintmax_t EvenExponent> struct power_impl2<B, EvenExponent, false> {
  typedef typename power_impl2<B, EvenExponent/2>::type half;
  typedef typename multiply<half, half>::type type;
};
template<typename B, uintmax_t OddExponent> struct power_impl2<B, OddExponent, true> {
  typedef typename power_impl2<B, OddExponent/2>::type half;
  typedef typename multiply<B, typename multiply<half, half>::type>::type type;
};
template<typename Base, intmax_t Exponent> struct power_impl1<Base, Exponent, false>
  : power_impl2<Base, (uintmax_t(Exponent))> {};
template<typename Base, intmax_t Exponent>
struct power : power_impl1<Base, Exponent> {};


///////////////////////////
//// Natural-number truncating division
///////////////////////////

// Natural number division
template<typename Nat, shift_type Shift> struct shift_left_by_milliodigits;
template<milliodigit... Milliodigits, shift_type Shift>
struct shift_left_by_milliodigits<nat<Milliodigits...>, Shift>
  : shift_left_by_milliodigits<nat<0, Milliodigits...>, (Shift-1)> {};
template<milliodigit... Milliodigits>
struct shift_left_by_milliodigits<nat<Milliodigits...>, 0> { typedef nat<Milliodigits...> type; };
template<shift_type Shift>
struct shift_left_by_milliodigits<nat<>, Shift> { typedef nat<> type; };
template<>
struct shift_left_by_milliodigits<nat<>, 0> { typedef nat<> type; };
template<milliodigit... Milliodigits, shift_type Shift>
struct shift_left_by_milliodigits<negative<nat<Milliodigits...>>, Shift> {
  typedef negative<typename shift_left_by_milliodigits<nat<Milliodigits...>, Shift>::type> type;
};

template<typename Nat, shift_type Shift> struct shift_right_by_milliodigits;
template<milliodigit... Milliodigits, milliodigit Milliodigit0, shift_type Shift>
struct shift_right_by_milliodigits<nat<Milliodigit0, Milliodigits...>, Shift>
  : shift_right_by_milliodigits<nat<Milliodigits...>, (Shift-1)> {};
template<shift_type Shift> struct shift_right_by_milliodigits<nat<>, Shift> { typedef nat<> type; };
template<milliodigit... Milliodigits, milliodigit Milliodigit0>
struct shift_right_by_milliodigits<nat<Milliodigit0, Milliodigits...>, 0> { typedef nat<Milliodigit0, Milliodigits...> type; };
template<milliodigit... Milliodigits>
struct shift_right_by_milliodigits<nat<Milliodigits...>, 0> { typedef nat<Milliodigits...> type; };
template<> struct shift_right_by_milliodigits<nat<>, 0> { typedef nat<> type; };
// TODO this rounds towards zero not -Inf.. is that what we want?
// Now changed to round towards -Inf
template<milliodigit... Milliodigits, shift_type Shift>
struct shift_right_by_milliodigits<negative<nat<Milliodigits...>>, Shift> {
  typedef negative<typename shift_right_by_milliodigits<typename add<nat<Milliodigits...>, nat<(base-1)>>::type, Shift>::type> type;
};

template<typename Nat, shift_type Shift, bool Left = (Shift >= 0)> struct shift_by_milliodigits;
template<typename Nat, shift_type Shift> struct shift_by_milliodigits<Nat, Shift, true>
  : shift_left_by_milliodigits<Nat, Shift> {};
template<typename Nat, shift_type Shift> struct shift_by_milliodigits<Nat, Shift, false>
  : shift_right_by_milliodigits<Nat, (-Shift)> {};

template<typename Nat> struct milliodigit_count;
template<milliodigit... Milliodigits, milliodigit Milliodigit0>
struct milliodigit_count<nat<Milliodigit0, Milliodigits...>>
  : boost::integral_constant<shift_type, (1 + milliodigit_count<nat<Milliodigits...>>::value)> {};
template<> struct milliodigit_count<nat<>> : boost::integral_constant<shift_type, 0> {};


// represents Nat * base^^(-milliodigits_shifted_right)
// represents Nat * base^^milliodigits_shifted_left
// only used internally, because it's redundant with 'rational'
// always used with negative exponent
template<shift_type milliodigits_shifted_left, typename Nat>
struct floating {};
//template<shift_type new_precision_in_milliodigits, typename Floating> struct trim_floating_precision;
// trims precision or adds excess precision:
template<shift_type new_base_left_shift, typename Floating> struct set_floating_milliodigit_exponent_rounding_down;
//template<shift_type new_base_right_shift, typename Floating> struct set_floating_milliodigit_exponent;

template<shift_type new_base_left_shift, typename Nat, shift_type Shift>
struct set_floating_milliodigit_exponent_rounding_down<new_base_left_shift, floating<Shift, Nat>> {
  typedef floating<new_base_left_shift,
    typename shift_by_milliodigits<Nat, (Shift - new_base_left_shift)>::type> type;
};
/*
//_towards_zero?
template<shift_type new_base_left_shift, typename Nat, shift_type Shift>
struct set_floating_milliodigit_exponent_rounding_down<new_base_left_shift, floating<Shift, negative<Nat>>> {
  typedef floating<new_base_left_shift,
    negative<typename shift_by_milliodigits<Nat, (Shift - new_base_left_shift)>::type>> type;
};*/

// implement the operations needed for below
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct multiply<floating<ShiftA, NatA>, floating<ShiftB, NatB>> {
  typedef floating<(ShiftA + ShiftB),
                   typename multiply<NatA, NatB>::type> type;
};
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB, bool LT = (ShiftA < ShiftB)>
struct add_floating;
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct add_floating<ShiftA, NatA, ShiftB, NatB, true> {
  typedef floating<ShiftA,
    typename add<
      NatA,
      typename shift_left_by_milliodigits<NatB, (ShiftB - ShiftA)>::type
      >::type> type;
};
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct add_floating<ShiftA, NatA, ShiftB, NatB, false> {
  typedef floating<ShiftB,
    typename add<
      typename shift_left_by_milliodigits<NatA, (ShiftA - ShiftB)>::type,
      NatB
      >::type> type;
};
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB, bool LT = (ShiftA < ShiftB)>
struct subtract_floating;
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct subtract_floating<ShiftA, NatA, ShiftB, NatB, true> {
  typedef floating<ShiftA,
    typename subtract<
      NatA,
      typename shift_left_by_milliodigits<NatB, (ShiftB - ShiftA)>::type
      >::type> type;
};
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct subtract_floating<ShiftA, NatA, ShiftB, NatB, false> {
  typedef floating<ShiftB,
    typename subtract<
      typename shift_left_by_milliodigits<NatA, (ShiftA - ShiftB)>::type,
      NatB
      >::type> type;
};

template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct add<floating<ShiftA, NatA>, floating<ShiftB, NatB>>
  : add_floating<ShiftA, NatA, ShiftB, NatB> {};
template<shift_type ShiftA, typename NatA, shift_type ShiftB, typename NatB>
struct subtract<floating<ShiftA, NatA>, floating<ShiftB, NatB>>
  : subtract_floating<ShiftA, NatA, ShiftB, NatB> {};

template<shift_type Shift, typename Nat> struct round_down_to_nat<floating<Shift, Nat> >
 : shift_by_milliodigits<Nat, Shift> {};

template<typename Nat> struct initial_recip_estimate;// { static_assert(compare<Nat, nat<>>::value, "divide by zero"); };
/*wrong endianness
template<milliodigit MilliodigitM0, milliodigit MilliodigitM1, milliodigit... Milliodigits>
struct initial_recip_estimate<nat<MilliodigitM0, MilliodigitM1, Milliodigits...>> {
  typedef nat<(base*base / (MilliodigitM0*base + MilliodigitM1))> type;
  static_assert(value < base, "bug");
};*/
template<milliodigit Milliodigit0, milliodigit...Milliodigits>
struct initial_recip_estimate<nat<Milliodigit0, Milliodigits...>>
  : initial_recip_estimate<nat<Milliodigits...>> {};
#if 0
template<milliodigit MilliodigitM0, milliodigit MilliodigitM1>
struct initial_recip_estimate<nat<MilliodigitM1, MilliodigitM0>> {
  typedef nat<(base*base / (MilliodigitM0*base + MilliodigitM1))> type;
};
template<milliodigit Milliodigit>
struct initial_recip_estimate<nat<Milliodigit>> {
  typedef nat<(base / Milliodigit)> type;
};
template<>
struct initial_recip_estimate<nat<1>> {
  typedef nat<0, 1> type;
};
template<>
struct initial_recip_estimate<nat<0, 1>> {
  typedef nat<0, 1> type;
};
#endif
// Always round down, so that this code never need deal with a negative adjustment,
// and so that base / Milliodigit < base even when Milliodigit is 1.
template<milliodigit MilliodigitM0, milliodigit MilliodigitM1>
struct initial_recip_estimate<nat<MilliodigitM1, MilliodigitM0>> {
  typedef nat<(base*base / (MilliodigitM0*base + MilliodigitM1 + 1))> type;
};
template<milliodigit Milliodigit>
struct initial_recip_estimate<nat<Milliodigit>> {
  typedef nat<(base / (Milliodigit + 1))> type;
};

template<typename Nat, shift_type milliodigits_below_1, typename EstimateFloating>
struct reciprocal_nat_as_floating_impl1;
template<typename Nat, shift_type milliodigits_below_1, typename EstimateFloating, typename Adjust>
struct reciprocal_nat_as_floating_impl2
  : reciprocal_nat_as_floating_impl1<Nat, milliodigits_below_1, typename add<Adjust, EstimateFloating>::type> {};
template<typename EstimateFloating, typename TestMultiply> struct reciprocal_nat_as_floating_impl3;
template<typename Nat, shift_type milliodigits_below_1, typename EstimateFloating, shift_type Shift>
struct reciprocal_nat_as_floating_impl2<Nat, milliodigits_below_1, EstimateFloating, floating<Shift, nat<>>>
  : reciprocal_nat_as_floating_impl3<EstimateFloating,
      typename round_down_to_nat<typename multiply<floating<0, Nat>, EstimateFloating>::type>::type> {};
template<shift_type Shift, typename Nat>
struct reciprocal_nat_as_floating_impl3<floating<Shift, Nat>, nat<>> {
  typedef floating<Shift, typename add<Nat, nat<1>>::type> type;
};
template<shift_type Shift, typename Nat>
struct reciprocal_nat_as_floating_impl3<floating<Shift, Nat>, nat<1>> {
  typedef floating<Shift, Nat> type;
};
template<typename Nat, shift_type milliodigits_below_1, typename EstimateFloating>
struct reciprocal_nat_as_floating_impl1
#if 0
{
  typedef typename set_floating_milliodigit_exponent_rounding_down<
        (-milliodigits_below_1),
        typename multiply<
          EstimateFloating,
          typename subtract<
            floating<0, nat<1>>,
            typename multiply<
              floating<0, Nat>,
              EstimateFloating
              >::type>::type>::type>::type type;
};
#else
  : reciprocal_nat_as_floating_impl2<Nat, milliodigits_below_1, EstimateFloating,
      typename set_floating_milliodigit_exponent_rounding_down<
        (-milliodigits_below_1),
        typename multiply<
          EstimateFloating,
          typename subtract<
            floating<0, nat<1>>,
            typename multiply<
              floating<0, Nat>,
              EstimateFloating
              >::type>::type>::type>::type> {};
#endif
// use (milliodigits_of_precision+1) because there is no attempt for the
// top milliodigit to have a full digit's worth of information.
template<typename Nat, shift_type milliodigits_of_precision>
struct reciprocal_nat_as_floating {
  static const shift_type actual_precision = milliodigits_of_precision+1;
  static const shift_type milliodigits_below_1 = actual_precision + (milliodigit_count<Nat>::value - 1);
  typedef typename reciprocal_nat_as_floating_impl1<
      Nat,
      milliodigits_below_1,
      floating<
        (-milliodigits_below_1),
        //(-actual_precision-(milliodigit_count<Nat>::value-1)),
        typename shift_left_by_milliodigits<
          typename initial_recip_estimate<Nat>::type,
          (actual_precision-1)
        >::type
      >
    >::type type;
};

template<typename NatA, typename NatB> struct divide_nat_optimize {
#if 0
  typedef typename multiply<
      floating<0, NatA>,
      typename reciprocal_nat_as_floating<
        NatB,
        milliodigit_count<NatA>::value>::type>::type quotz;
  typedef typename reciprocal_nat_as_floating<
        NatB,
        milliodigit_count<NatA>::value>::type quot;
  typedef quot type;
  typedef nat<> rem;
#else
  typedef typename round_down_to_nat<typename multiply<
      floating<0, NatA>,
      typename reciprocal_nat_as_floating<
        NatB,
        milliodigit_count<NatA>::value>::type>::type>::type quot;
  typedef quot type;
  typedef typename subtract_nat<NatA, typename multiply<quot, NatB>::type>::type rem;
#endif
};
#if 0
template<typename NatA> struct divide_nat_optimize<NatA, nat<1>> {
  typedef NatA quot;
  typedef quot type;
  typedef nat<> rem;
};
template<typename Nat> struct divide_nat_optimize<Nat, Nat> {
  typedef nat<1> quot;
  typedef quot type;
  typedef nat<> rem;
};
template<typename NatA> struct divide_nat_optimize<nat<1>, nat<1>> {
  typedef nat<1> quot;
  typedef quot type;
  typedef nat<> rem;
};
template<typename NatA, milliodigit MilliodigitB> struct divide_nat_optimize<NatA, nat<MilliodigitB>>
  : divide_by_milliodigit<NatA, MilliodigitB> {};
#endif
template<typename NatA, typename NatB> struct divide_nat
  : divide_nat_optimize<NatA, NatB> {};
template<typename NatA> struct divide_nat<NatA, nat<>> {
  typedef divide_by_zero quot;
  typedef quot type;
  typedef divide_by_zero rem;
};

///////////////////////////
//// Integer operations as extensions of natural number operations
///////////////////////////

// Integer subtraction
template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct subtract<nat<MilliodigitsA...>, nat<MilliodigitsB...>>
  : boost::conditional<
      (compare<nat<MilliodigitsA...>, nat<MilliodigitsB...>>::value == -1),
        negative<typename subtract_nat<nat<MilliodigitsB...>,
                                       nat<MilliodigitsA...>>::type>,
        typename subtract_nat<nat<MilliodigitsA...>,
                              nat<MilliodigitsB...>>::type> {};

// extend add and subtract to operations involving one or more negatives
template<typename T1, typename T2> struct add<negative<T1>, negative<T2>> {
  typedef negative<typename add<T1, T2>::type> type;
};
template<typename T1, typename T2> struct add<negative<T1>, T2> : subtract<T2, T1> {};
template<typename T1, typename T2> struct add<T1, negative<T2>> : subtract<T1, T2> {};
template<typename T1, typename T2> struct subtract<negative<T1>, negative<T2>> : subtract<T2, T1> {};
template<typename T1, typename T2> struct subtract<T1, negative<T2>> : add<T1, T2> {};
template<typename T1, typename T2> struct subtract<negative<T1>, T2> {
  typedef negative<typename add<T1, T2>::type> type;
};

// Integer multiplication
template<typename T1, typename T2> struct multiply<negative<T1>, negative<T2>> : multiply<T1, T2> {};
template<typename T1, typename T2> struct multiply<negative<T1>, T2> {
  typedef negative<typename multiply<T1, T2>::type> type;
};
template<typename T1, typename T2> struct multiply<T1, negative<T2>> {
  typedef negative<typename multiply<T1, T2>::type> type;
};

// Integer comparison
template<typename T1, typename T2> struct compare<negative<T1>, negative<T2>> : compare<T2, T1> {};
template<typename T1, typename T2> struct compare<negative<T1>, T2> : boost::integral_constant<int, (-1)> {};
template<typename T1, typename T2> struct compare<T1, negative<T2>> : boost::integral_constant<int, 1> {};

// Integer negation
template<typename T> struct negate {
  typedef negative<T> type;
};
template<typename T> struct negate<negative<T>> {
  typedef T type;
};

// Integer division following int operator/ semantics
template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<nat<MilliodigitsA...>, nat<MilliodigitsB...>>
  : divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> {};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<negative<nat<MilliodigitsA...>>, nat<MilliodigitsB...>> {
  typedef divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> unsigned_result_;
  typedef negative<typename unsigned_result_::quot> type;
  typedef negative<typename unsigned_result_::quot> quot;
  typedef negative<typename unsigned_result_::rem> rem;
};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<nat<MilliodigitsA...>, negative<nat<MilliodigitsB...>>> {
  typedef divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> unsigned_result_;
  typedef negative<typename unsigned_result_::quot> type;
  typedef negative<typename unsigned_result_::quot> quot;
  typedef typename unsigned_result_::rem rem;
};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<negative<nat<MilliodigitsA...>>, negative<nat<MilliodigitsB...>>> {
  typedef divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> unsigned_result_;
  typedef typename unsigned_result_::quot type;
  typedef typename unsigned_result_::quot quot;
  typedef negative<typename unsigned_result_::rem> rem;
};


///////////////////////////
//// GCD (Greatest Common Denominator)
///////////////////////////

template<typename Integer> struct even;
template<> struct even<nat<>> : boost::true_type {};
template<milliodigit... Milliodigits, milliodigit Milliodigit0>
struct even<nat<Milliodigit0, Milliodigits...>>
  : boost::integral_constant<bool, ((Milliodigit0 & 1) == 0)> {};
template<milliodigit... Milliodigits>
struct even<negative<nat<Milliodigits...>>> : even<nat<Milliodigits...>> {};
template<typename Integer> struct odd
  : boost::integral_constant<bool, (!even<Integer>::value)> {};

//This implementation works because 2 (half: divide-by-2) is a factor of 'base'.
template<typename Nat> struct half_nat_rounding_down;
// implemented by recursive, low-to-high right-shifting
template<milliodigit... Milliodigits, milliodigit Milliodigit1, milliodigit Milliodigit0>
struct half_nat_rounding_down<nat<Milliodigit0, Milliodigit1, Milliodigits...>>
  : make_nat_from_rest_and_least<
      typename half_nat_rounding_down<nat<Milliodigit1, Milliodigits...>>::type,
      ((Milliodigit0 >> 1) + ((Milliodigit1 & 1) * (base >> 1)))> {};
// recursion base cases
template<milliodigit Milliodigit0>
struct half_nat_rounding_down<nat<Milliodigit0>> { typedef nat<(Milliodigit0 >> 1)> type; };
template<> struct half_nat_rounding_down<nat<1>> { typedef nat<> type; };
template<> struct half_nat_rounding_down<nat<>> { typedef nat<> type; };

template<typename Nat> struct half_nat_exactly : half_nat_rounding_down<Nat> {
  static_assert(even<Nat>::value, "half_nat_exactly used on an odd number");
};

// Lazy-programmer implementation of 'twice':
template<typename Num> struct twice : add<Num, Num> {};

// Binary GCD algorithm: https://en.wikipedia.org/wiki/Binary_GCD_algorithm
template<typename NatA, typename NatB, bool AIsEven, bool BIsEven> struct gcd_aux;
template<typename NatA, typename NatB, int Compare> struct gcd_aux_2;
template<typename NatA, typename NatB> struct gcd
  : gcd_aux<NatA, NatB, even<NatA>::value, even<NatB>::value> {};
template<typename T> struct gcd<T, T> {
  typedef T type;
};
template<typename NatA> struct gcd<NatA, nat<>> {
  typedef NatA type;
};
template<typename NatB> struct gcd<nat<>, NatB> {
  typedef NatB type;
};
template<typename NatA, typename NatB> struct gcd_aux<NatA, NatB, true, true>
  : twice<typename gcd<
      typename half_nat_exactly<NatA>::type,
      typename half_nat_exactly<NatB>::type>::type> {};
template<typename NatA, typename NatB> struct gcd_aux<NatA, NatB, true, false>
  : gcd<typename half_nat_exactly<NatA>::type, NatB> {};
template<typename NatA, typename NatB> struct gcd_aux<NatA, NatB, false, true>
  : gcd<NatA, typename half_nat_exactly<NatB>::type> {};
template<typename NatA, typename NatB> struct gcd_aux<NatA, NatB, false, false>
  : gcd_aux_2<NatA, NatB, compare<NatA, NatB>::value> {};
template<typename NatA, typename NatB> struct gcd_aux_2<NatA, NatB, 1>
  : gcd<typename half_nat_exactly<typename subtract<NatA, NatB>::type>::type, NatB> {};
template<typename NatA, typename NatB> struct gcd_aux_2<NatA, NatB, (-1)>
  : gcd<typename half_nat_exactly<typename subtract<NatB, NatA>::type>::type, NatA> {};

///////////////////////////
//// Rational numbers
///////////////////////////

template<typename IntNum, typename IntDen> struct make_rational;
template<typename Num, typename Den>
struct make_rational_impl1 {
  typedef rational<Num, Den> type;
};
template<typename Num>
struct make_rational_impl1<Num, nat<1>> {
  typedef Num type;
};
template<typename Num, typename Den>
struct make_rational {
  typedef typename gcd<Num, Den>::type gcd_;
  typedef typename make_rational_impl1<
    typename divide_nat<Num, gcd_>::type,
    typename divide_nat<Den, gcd_>::type>::type type;
};
template<typename Num, typename Den>
struct make_rational<negative<Num>, negative<Den>>
  : make_rational<Num, Den> {};
template<typename Num, typename Den> struct make_rational<negative<Num>, Den> {
  typedef negative<typename make_rational<Num, Den>::type> type;
};
template<typename Num, typename Den> struct make_rational<Num, negative<Den>> {
  typedef negative<typename make_rational<Num, Den>::type> type;
};
template<typename Num>
struct make_rational<Num, nat<>> {
  // error: divide by zero
  typedef divide_by_zero type;
};
template<>
struct make_rational<nat<>, nat<>> {
  // error: divide by zero
  typedef divide_by_zero type;
};
template<typename Den>
struct make_rational<nat<>, Den> {
  typedef nat<> type;
};

template<typename Num> struct as_rational {
  typedef Num type; //assume rational already
};
template<milliodigit...Milliodigits> struct as_rational<nat<Milliodigits...>> {
  typedef rational<nat<Milliodigits...>, nat<1>> type;
};
template<milliodigit...Milliodigits> struct as_rational<negative<nat<Milliodigits...>>> {
  typedef negative<rational<nat<Milliodigits...>, nat<1>>> type;
};

// use negative<> overloads for dealing with negative rationals
template<typename NumA, typename DenA, typename NumB, typename DenB>
struct add<rational<NumA, DenA>, rational<NumB, DenB>>
  : make_rational<
      typename add<
        typename multiply<NumA, DenB>::type,
        typename multiply<DenA, NumB>::type>::type,
      typename multiply<DenA, DenB>::type> {};

template<typename NumA, typename DenA, typename NumB, typename DenB>
struct subtract<rational<NumA, DenA>, rational<NumB, DenB>>
  : make_rational<
      typename subtract<
        typename multiply<NumA, DenB>::type,
        typename multiply<DenA, NumB>::type>::type,
      typename multiply<DenA, DenB>::type> {};

template<typename NumA, typename DenA, typename NumB, typename DenB>
struct multiply<rational<NumA, DenA>, rational<NumB, DenB>>
  : make_rational<
      typename multiply<NumA, NumB>::type,
      typename multiply<DenA, DenB>::type> {};

template<typename NumA, typename DenA, typename NumB, typename DenB>
struct divide_rational<rational<NumA, DenA>, rational<NumB, DenB>>
  : make_rational<
      typename multiply<NumA, DenB>::type,
      typename multiply<DenA, NumB>::type> {};

template<typename Num, typename Den>
struct reciprocal<rational<Num, Den>> {
  typedef rational<Den, Num> type;
};
template<milliodigit...Milliodigits>
struct reciprocal<nat<Milliodigits...>> {
  typedef rational<nat<1>, nat<Milliodigits...>> type;
};
template<>
struct reciprocal<nat<>> {
  // error: divide by zero
  typedef divide_by_zero type;
};
template<>
struct reciprocal<nat<1>> {
  typedef nat<1> type;
};
template<milliodigit...Milliodigits>
struct reciprocal<negative<nat<Milliodigits...>>> {
  typedef negative<typename reciprocal<nat<Milliodigits...>>::type> type;
};


// dull conversions
template<typename NumA, typename DenA, milliodigit...MilliodigitsB>
struct add<rational<NumA, DenA>, nat<MilliodigitsB...>>
  : add<rational<NumA, DenA>, typename as_rational<nat<MilliodigitsB...>>::type> {};
template<typename NumA, typename DenA, milliodigit...MilliodigitsB>
struct subtract<rational<NumA, DenA>, nat<MilliodigitsB...>>
  : subtract<rational<NumA, DenA>, typename as_rational<nat<MilliodigitsB...>>::type> {};
template<typename NumA, typename DenA, milliodigit...MilliodigitsB>
struct multiply<rational<NumA, DenA>, nat<MilliodigitsB...>>
  : multiply<rational<NumA, DenA>, typename as_rational<nat<MilliodigitsB...>>::type> {};
template<typename NumA, typename DenA, milliodigit...MilliodigitsB>
struct divide_rational<rational<NumA, DenA>, nat<MilliodigitsB...>>
  : divide_rational<rational<NumA, DenA>, typename as_rational<nat<MilliodigitsB...>>::type> {};

template<milliodigit...MilliodigitsA, typename NumB, typename DenB>
struct add<nat<MilliodigitsA...>, rational<NumB, DenB>>
  : add<typename as_rational<nat<MilliodigitsA...>>::type, rational<NumB, DenB>> {};
template<milliodigit...MilliodigitsA, typename NumB, typename DenB>
struct subtract<nat<MilliodigitsA...>, rational<NumB, DenB>>
  : subtract<typename as_rational<nat<MilliodigitsA...>>::type, rational<NumB, DenB>> {};
template<milliodigit...MilliodigitsA, typename NumB, typename DenB>
struct multiply<nat<MilliodigitsA...>, rational<NumB, DenB>>
  : multiply<typename as_rational<nat<MilliodigitsA...>>::type, rational<NumB, DenB>> {};
template<milliodigit...MilliodigitsA, typename NumB, typename DenB>
struct divide_rational<nat<MilliodigitsA...>, rational<NumB, DenB>>
  : divide_rational<typename as_rational<nat<MilliodigitsA...>>::type, rational<NumB, DenB>> {};

template<milliodigit...MilliodigitsA, milliodigit...MilliodigitsB>
struct divide_rational<nat<MilliodigitsA...>, nat<MilliodigitsB...>>
  : divide_rational<
      typename as_rational<nat<MilliodigitsA...>>::type,
      typename as_rational<nat<MilliodigitsB...>>::type
    > {};

template<typename T1, typename T2> struct divide_rational<negative<T1>, negative<T2>>
  : divide_rational<T1, T2> {};
template<typename T1, typename T2> struct divide_rational<negative<T1>, T2> {
  typedef negative<typename divide_rational<T1, T2>::type> type;
};
template<typename T1, typename T2> struct divide_rational<T1, negative<T2>> {
  typedef negative<typename divide_rational<T1, T2>::type> type;
};

//template<typename Base, intmax_t Exponent> struct power_impl1<Base, Exponent, true>
//  : reciprocal<typename power_impl2<Base, (-uintmax_t(Exponent))>::type> {};

// numerator<>
// denominator<>
// digit<base, exp>
namespace impl {
struct make_any_number {
  template<typename T> constexpr inline operator T()const { return T(); }
};
}

// ADL (Argument Dependent Lookup) will find these templates when appropriate...
// but wait, will they - ok added number<> base
// div vs / ??

template<typename A> constexpr inline A
operator+(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline typename negate<A>::type
operator-(number<A>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline typename add<A, B>::type
operator+(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline typename subtract<A, B>::type
operator-(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline typename multiply<A, B>::type
operator*(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline typename divide_rational<A, B>::type
operator/(number<A>, number<B>) { return impl::make_any_number(); }

template<typename Quot, typename Rem>
struct ctdiv_t {
  Quot quot;
  Rem rem;
};
template<typename A, typename B> constexpr inline
ctdiv_t<typename divide_integer<A, B>::quot, typename divide_integer<A, B>::rem>
div(number<A>, number<B>) { return impl::make_any_number(); }


template<typename A, typename B> constexpr inline bool
operator==(number<A>, number<B>) { return compare<A, B>::value == 0; }
template<typename A, typename B> constexpr inline bool
operator!=(number<A>, number<B>) { return compare<A, B>::value != 0; }
template<typename A, typename B> constexpr inline bool
operator<(number<A>, number<B>) { return compare<A, B>::value == -1; }
template<typename A, typename B> constexpr inline bool
operator>(number<A>, number<B>) { return compare<A, B>::value == 1; }
template<typename A, typename B> constexpr inline bool
operator<=(number<A>, number<B>) { return compare<A, B>::value != 1; }
template<typename A, typename B> constexpr inline bool
operator>=(number<A>, number<B>) { return compare<A, B>::value != -1; }
/*
template<typename A> constexpr inline typename
round<A, rounding_strategy<round_down, negative_continuous_with_positive>>::type
floor(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline typename
round<A, rounding_strategy<round_up, negative_continuous_with_positive>>::type
ceil(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline typename
round<A, rounding_strategy<round_to_nearest_with_ties_rounding_to_even>>::type
nearbyint(number<A>) { return impl::make_any_number(); }
*/

// and operator bool, maybe explicit operator float, etc.
// I suppose the user-visible types should *anyway* be different from the intermediates...

//#define NATt(integer) typename ::ct::integer_literal<(integer)>::type
//#define NATv(integer) (typename ::ct::integer_literal<(integer)>::type())

//pow
//abs

// NATlong or NAT using user integer literals
// usually, and fixed-limit slower-compile stringify constant exprs
// for gcc 4.6

template<char...Digits> struct parse_nat;

// based on open-std document N3531 2013-03-08
// http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2013/n3531.pdf
template <unsigned Base, typename Accumulator, char... Digits> struct parse_digits;
template <char... Digits> struct parse_nat;
template <char... Digits> struct parse_nat<'0','x',Digits...> : parse_digits<16, nat<>, Digits...> {};
template <char... Digits> struct parse_nat<'0','X',Digits...> : parse_digits<16, nat<>, Digits...> {};
template <char... Digits> struct parse_nat<'0',Digits...> : parse_digits<8, nat<>, Digits...> {};
template <char... Digits> struct parse_nat : parse_digits<10, nat<>, Digits...> {};

template <unsigned Base, typename Accumulator>
struct parse_digits<Base, Accumulator> { typedef Accumulator type; };

template <unsigned Base, typename Accumulator, char... Digits>
struct parse_digits<Base, Accumulator, '\0', Digits...> { typedef Accumulator type; };

template <unsigned Base, typename Accumulator, char... Digits>
struct parse_digits<Base, Accumulator, '0', Digits...>
  : parse_digits<Base, typename multiply_by_milliodigit<Accumulator, Base>::type, Digits...> {};
// could do it 6-at-a-time for decimal..
#define PARSE_DIGIT(ch, n) \
template <unsigned Base, typename Accumulator, char... Digits> \
struct parse_digits<Base, Accumulator, ch, Digits...> \
  : parse_digits<\
      Base, \
      typename add< \
        typename multiply_by_milliodigit<Accumulator, Base>::type, \
        nat<n> >::type, \
      Digits...> {};
PARSE_DIGIT('1', 1)
PARSE_DIGIT('2', 2)
PARSE_DIGIT('3', 3)
PARSE_DIGIT('4', 4)
PARSE_DIGIT('5', 5)
PARSE_DIGIT('6', 6)
PARSE_DIGIT('7', 7)
PARSE_DIGIT('8', 8)
PARSE_DIGIT('9', 9)
PARSE_DIGIT('a', 0xa)
PARSE_DIGIT('b', 0xb)
PARSE_DIGIT('c', 0xc)
PARSE_DIGIT('d', 0xd)
PARSE_DIGIT('e', 0xe)
PARSE_DIGIT('f', 0xf)
PARSE_DIGIT('A', 0xA)
PARSE_DIGIT('B', 0xB)
PARSE_DIGIT('C', 0xC)
PARSE_DIGIT('D', 0xD)
PARSE_DIGIT('E', 0xE)
PARSE_DIGIT('F', 0xF)
#undef PARSE_DIGIT
/*
template <unsigned base> struct parse_digits<base> { typedef nat<> type; };
template <unsigned base, char... Digits>
struct parse_digits<base,'0',Digits...> {
  typedef typename parse_digits<base, Digits>::type type;
};
// could do it 6-at-a-time for decimal..
#define PARSE_DIGIT(ch, n) \
template <unsigned base, char... Digits> \
struct parse_digits<base,ch,Digits...> { \
static_assert(base > n,"invalid digit"); \
typedef typename add< \
    typename multiply_by_milliodigit< \
      typename power<base, sizeof...(Digits)>, \
      1>, \
    parse_digits<base,Digits...>::type \
  >::type type; \
};
PARSE_DIGIT('1', 1)
PARSE_DIGIT('2', 2)
PARSE_DIGIT('3', 3)
PARSE_DIGIT('4', 4)
PARSE_DIGIT('5', 5)
PARSE_DIGIT('6', 6)
PARSE_DIGIT('7', 7)
PARSE_DIGIT('8', 8)
PARSE_DIGIT('9', 9)
PARSE_DIGIT('a', 0xa)
PARSE_DIGIT('b', 0xb)
PARSE_DIGIT('c', 0xc)
PARSE_DIGIT('d', 0xd)
PARSE_DIGIT('e', 0xe)
PARSE_DIGIT('f', 0xf)
PARSE_DIGIT('A', 0xA)
PARSE_DIGIT('B', 0xB)
PARSE_DIGIT('C', 0xC)
PARSE_DIGIT('D', 0xD)
PARSE_DIGIT('E', 0xE)
PARSE_DIGIT('F', 0xF)
#undef PARSE_DIGIT
#endif
*/

template<char...Digits>
constexpr inline typename parse_nat<Digits...>::type
_NAT() { return impl::make_any_number(); }

} /* end namespace ct */

#if LASERCAKE_HAVE_USER_DEFINED_LITERALS

template<char...Digits>
constexpr inline typename ct::parse_nat<Digits...>::type
operator "" _NAT() { return ct::impl::make_any_number(); }

#define NAT(n) (BOOST_PP_CAT(n,_NAT))

#else

// This might not be the fastest possible; our GCC 4.6 support is grudging.
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#define CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS_AUX(z, n, str) \
  (sizeof(str) > n ? str[n] : '\0')
#define CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS(ident) \
  BOOST_PP_ENUM(BOOST_PP_LIMIT_REPEAT, CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS_AUX, #ident)
#define NAT(n) (::ct::_NAT<CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS(n)>())
#endif

// NATt NATv
//#define NATtype(integer) typename ::ct::literal<decltype((integer)), (integer)>::type
//#define NAT(integer) (make_literal<decltype((integer)), (integer)>())
//#define NATtype(integer) typename ::ct::literal<unsigned long long, (BOOST_PP_CAT(integer,ull))>::type
//#define NATtype(integer) typename ::ct::literal<uintmax_t, (UINTMAX_C(integer))>::type
//#define NAT(integer) (make_literal<uintmax_t, (UINTMAX_C(integer))>())

//domain
struct N;
struct Z;
struct Q;
struct Qplus; //?
struct NaN;// division by zero, subtraction in N
//Should we have +/-Inf? Pros? cons?
//current opinion: the above is unnecessary
//Oh hm we might need constexpr sqrt.  Then it needs to specify
//a rounding mode and/or precision?  Or we need to make real<> HAR HAR HAR
//(which contains a compile-time-function exponent->bit or similar)





#endif
