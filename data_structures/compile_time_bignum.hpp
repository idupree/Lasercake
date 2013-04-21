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
//TODO once used in units, try changing the base and see the different speeds.
//then reversey to bigendian, & 'int' for nicer error msgs (;base 1000)
//Ooh because we only use multiply for high stuff and because
//it can be emulated with 4 smaller multiplies, we could actually have
//the impl word size be 10^18!
//Need to get rid of that last user-facing static_assert so that they don't
//have to see the intermediate types.
//Rounding to a rational precision...well there is 'height'...

#include <boost/type_traits/conditional.hpp>
#include "rounding_strategy.hpp"

// Compile-time arbitrary-precision naturals, integers, and rational numbers
// whose value is represented in their type (specifically, in their template
// arguments).  Equal values are guaranteed to have identical types.
//
//  Example:
//      INTEGER(10000000000000000000000000000000000000000) * RATIONAL(2, 3)
//  represents exactly 20000000000000000000000000000000000000000/3
//
//  Example:
//    // yocto- is the SI prefix for 10^-24.
//    constexpr auto yocto_value = pow(INTEGER(10), INTEGER(-24));
//    typedef decltype(atto_value) yocto_type; //or
//    typedef decltype(pow(INTEGER(10), INTEGER(-24))) yocto_type;
//
// === TYPES (you will see these in compiler error messages) ===
// In compile error messages, natural numbers are represented in little-endian
// base-1-million digits, each digit (probably) shown in decimal.
// For example, a billion is ct::nat<0, 1000>.
//
// (TODO: make it big-endian because that's more readable, at the cost of
// compile time reversing the order to compute e.g. addition and multiplication.)
//
// Negative integers are negative<nat<..>>.
// Non-integer rationals are in canonical form and
//   rational<nat<..>, nat<..>>   or   negative<rational<nat<..>, nat<..>>>.
// All compile-time-number types T are wrapped in ct::number<T> to make
// operator/function overloading on them easier.
//
//
// === CONSTRUCTION ===
// Literals are macros
//   INTEGER(integer-literal)
//   RATIONAL(integer-literal-numerator, integer-literal-denominator)
//   INTEGERtype(integer-literal)
//   RATIONALtype(integer-literal-numerator, integer-literal-denominator)
// INTEGER(n) returns a value, and INTEGERtype(n) returns
// decltype(INTEGER(n)); likewise for RATIONAL.  Literals support
// any number of digits[1].
//
// Non numeric-literal integral-constant-expressions can be converted to
// this family of arbitrary-precision types thus:
//   ct::from_int<integral-constant-expression>::value
//   ct::from_uint<integral-constant-expression>::value
//   typename ct::from_int<integral-constant-expression>::type
//   typename ct::from_uint<integral-constant-expression>::type
//
// Each type is default-constructible and usable in arithmetic
// (TODO add conversions to other numeric types). It stores no runtime data.
//
// === OPERATIONS ===
// Operations are defined on the value level because it has nice syntax.
//
//   decltype(a)     //use one of these numbers as a type
//   a == b, a != b
//   a < b, a > b
//   a <= b, a >= b  // comparisons return constexpr 'bool' rather than
//                   // an arbitrary-precision-numeric type.
//   +a
//   -a
//   abs(a)
//   sign(a)         //signum (-1, 0 or 1); sign(a)*abs(a) == a
//   a + b
//   a - b
//   a * b
//   a / b           //rational division
//   reciprocal(a)   //rational reciprocal
//   div(a, b).quot  //integral division rounding towards zero
//   div(a, b).rem   //integral remainder rounding towards zero
//   pow(a, b)       //exact rational results permitted (i.e. all integer b
//                   //  and some rational (a, b)).
//   TODO LATER pow(a, b, rounding_strategy)
//   TODO isqrt(a)        //sqrt rounding down to the nearest integer
//   TODO LATER sqrt(a, rounding_strategy)
//   numerator(a)    //negative if a is negative
//   denominator(a)  //always nonnegative; '1' for integers
//   floor(a)        //round towards negative infinity to nearest integer
//   ceil(a)         //round towards positive infinity to nearest integer
//   nearbyint(a)    //round towards nearest integer (ties go to even)
//   round(a, rounding_strategy)
//                   // logarithms, rounded since we
//                   // can't represent irrational quantities:
//   log(base, number, rounding_strategy)
//   log2(number, rounding_strategy)
//   log10(number, rounding_strategy)
//                   // (TODO fix round-to-nearest-* rounding
//                   //  and allow different rounding precisions)
//
//   is_integer(a)   //these return constexpr bool
//   is_nonnegative_integer(a)
//   is_positive_integer(a)
//   is_negative(a)
//   is_nonnegative(a)
//   is_positive(a)
//   is_zero(a)
//
//   ct::to_int<IntegralType>(a)
//                   // returns constexpr IntegralType constructed by
//                   // casts from uintmax_t, addition, multiplication, and
//                   // negation.
//                   // "ct::" required because argument-dependent lookup
//                   // doesn't work with explicit function template arguments.
//
//TODO:
//   a % b  or  fmod(a, b) or remainder(a, b)   //rational modulus
//   sqrt(a, rounding_strategy)
//   nth root (also extends pow) (but requires rounding strategy,
//                       perhaps defaulting to require_exact_answer)
//   general: https://en.wikipedia.org/wiki/Newton-Raphson#Square_root_of_a_number
//   can iterate using rationals!
//   hmm.. simplify<rational<a,b>>? right shifts both a and b by the same amt?
// Should this support bit-shifts (same as value*pow(INT(2), shift) and
// value/pow(INT(2), shift); or should they be floor() of those?)?
// What about bitwise operations |&^~ ? Should they operate on the
// twos-complement representation?  Do we know how to compute this
// with rationals?
//
// [1] Potentially limited by compiler limits on number of template arguments
// (each character/digit is a template argument during our numeric-literal-parsing).
// For limited compilers with variadic templates but not
// user-defined-numeric-literals (GCC 4.6), it is also limited by
// Boost.Preprocessor's limit of 256 (BOOST_PP_LIMIT_REPEAT) digits.


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
struct nat {};//: number<nat<Milliodigits...>> {};
// < 0
template<typename Nonnegative>
struct negative {};//: number<negative<Nonnegative>> {};
//rational number
//normalized
// integers are nat
//positive: rational<nat<..>, nat<..>>
//zero: rational<nat<>, nat<1>>
//negative: negative<rational<nat<..>, nat<..>>>
template<typename Num, typename Den>
struct rational {};//: number<rational<Num,Den>> {};
// inheritance didn't give as strong function-overloading properties,
// may slow compile-time computations a bit,
// and we don't need it anymore.

struct divide_by_zero {};
struct even_root_of_negative {};

template<typename UInt, UInt...Values>
struct array {
  typedef UInt values_type[sizeof...(Values)];
  static values_type const& values() {
    static const values_type values_ = { Values... };
    return values_;
  }
};
template<typename...Array> struct concat;
template<typename Array> struct concat<Array> {
  typedef Array type;
};
template<typename UInt, UInt...ValuesA, UInt...ValuesB, typename...Arrays>
struct concat<array<UInt, ValuesA...>, array<UInt, ValuesB...>, Arrays...>
: concat<array<UInt, ValuesA..., ValuesB...>, Arrays...> {};

template<char...Values>
inline const char* to_string(array<char, Values...>) {
  return array<char, Values..., '\0'>::values();
}

//convert_to_little_endian_unsigned<uint32_t>
//TODO
//struct in_base
//struct integer_to_base

//big/little

// can be a wrapper fn:
// except then it needs to know whether twos-complement..
//auto-sized/num-limbs-and-??-if-overflow
//  error-if-overflow/modulo-if-overflow

//twos-complement/error-if-negative
//...I suppose this can be a wrapper:
//if negative, -x-1 (i.e. ~x), then
//compl all bits of the result

//zero-has-a-digit

//uint type
//base (default to max uint-type + 1)

namespace impl {


// Operations that have the same name as a function on numbers
// (e.g. abs) have an underscore appended, because otherwise the
// name-lookup rules find the struct and think there is an ambiguity.
template<typename... Num> struct add;
template<typename... Num> struct multiply;
template<typename NumA, typename NumB> struct subtract;
// Natural number subtraction (negative results become 'below_zero')
template<typename NatA, typename NatB> struct subtract_nat;
template<typename Num> struct negate;
template<typename Num> struct abs_;
template<typename Num> struct sign_;
template<typename Num, typename Exponent,
  typename RoundingStrategy = rounding_strategy<require_exact_answer> > struct power;
template<typename Num, typename Root,
  typename RoundingStrategy = rounding_strategy<require_exact_answer> > struct root;
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
template<typename Rat> struct reciprocal_;
// Since we can't represent irrational values,
// you must pass a RoundingStrategy to round this to a nearby
// integer (or half-integer).
template<typename Base, typename Number, typename RoundingStrategy> struct logarithm;

template<typename Num> struct is_integer_;
template<typename Num> struct is_nonnegative_integer_;
template<typename Num> struct is_positive_integer_;
template<typename Num> struct is_nonnegative_;
template<typename Num> struct is_positive_;
template<typename Num> struct is_negative_;
template<typename Num> struct is_zero_;

template<typename Num, typename RoundingStrategy> struct round_;


template<typename Num> struct round_down_to_nat;

// Works on any integer/rational; returns the canonical rational representation.
// For negative numbers, numerator is negative and denominator is positive.
template<typename Rat> struct numerator_;
template<typename Rat> struct denominator_;

// add and multiply are associative
template<typename NatA, typename NatB, typename NatC, typename...Nat>
struct add<NatA, NatB, NatC, Nat...> : add<typename add<NatA, NatB>::type, NatC, Nat...> {};
template<typename NatA, typename NatB, typename NatC, typename...Nat>
struct multiply<NatA, NatB, NatC, Nat...> : multiply<typename multiply<NatA, NatB>::type, NatC, Nat...> {};

// classification
template<typename Num>
struct is_integer_ : boost::false_type {};
template<milliodigit... Milliodigits>
struct is_integer_<nat<Milliodigits...>> : boost::true_type {};
template<milliodigit... Milliodigits>
struct is_integer_<negative<nat<Milliodigits...>>> : boost::true_type {};

template<typename Num>
struct is_nonnegative_integer_ : boost::false_type {};
template<milliodigit... Milliodigits>
struct is_nonnegative_integer_<nat<Milliodigits...>> : boost::true_type {};

template<typename Num>
struct is_positive_integer_ : boost::false_type {};
template<milliodigit... Milliodigits, milliodigit Milliodigit0>
struct is_positive_integer_<nat<Milliodigit0, Milliodigits...>> : boost::true_type {};

template<typename Num>
struct is_nonnegative_ : boost::true_type {};
template<typename T>
struct is_nonnegative_<negative<T>> : boost::false_type {};

template<typename Num>
struct is_negative_ : boost::false_type {};
template<typename T>
struct is_negative_<negative<T>> : boost::true_type {};

template<typename Num>
struct is_positive_ : boost::true_type {};
template<typename T>
struct is_positive_<negative<T>> : boost::false_type {};
template<>
struct is_positive_<nat<>> : boost::false_type {};

template<typename Num>
struct is_zero_ : boost::false_type {};
template<>
struct is_zero_<nat<>> : boost::true_type {};

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
  typedef typename make_nat_from_rest_and_least<
    typename uinteger_literal<(Int/base)>::type,
    (Int%base)>::type type;
};
template<> struct uinteger_literal<0> {
  typedef nat<> type;
};

template<typename Int, Int Value, bool Negative = (Value < 0)> struct literal;
template<typename Int, Int Value> struct literal<Int, Value, true> {
  typedef negative<typename uinteger_literal<static_cast<uintmax_t>(-Value)>::type> type;
};
template<typename Int, Int Value> struct literal<Int, Value, false>
  : uinteger_literal<static_cast<uintmax_t>(Value)> {};

template<typename NumericType, typename CTNum> struct to_int_;
template<typename NumericType> struct to_int_<NumericType, nat<>> {
  static constexpr NumericType value() { return NumericType(0); }
};
template<typename NumericType, milliodigit...Milliodigits, milliodigit Milliodigit0> struct to_int_<NumericType, nat<Milliodigit0, Milliodigits...>> {
  static constexpr NumericType value() { return
    NumericType(NumericType(Milliodigit0) +
      NumericType(base) * to_int_<NumericType, nat<Milliodigits...>>::value()); }
};
template<typename NumericType, typename CTNum> struct to_int_<NumericType, negative<CTNum>> {
  static constexpr NumericType value() { return NumericType(-to_int_<NumericType, CTNum>::value()); }
};

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
template<milliodigit... MilliodigitsA> struct add<nat<MilliodigitsA...>, nat<>> {
  typedef nat<MilliodigitsA...> type; };
template<milliodigit... MilliodigitsB> struct add<nat<>, nat<MilliodigitsB...>> {
  typedef nat<MilliodigitsB...> type; };
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
template<milliodigit... MilliodigitsA> struct subtract_nat<nat<MilliodigitsA...>, nat<>> {
  typedef nat<MilliodigitsA...> type; };
template<milliodigit... MilliodigitsB> struct subtract_nat<nat<>, nat<MilliodigitsB...>> {
  typedef below_zero type; };
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
template<milliodigit... MilliodigitsA> struct multiply<nat<MilliodigitsA...>, nat<>> { typedef nat<> type; };

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
template<typename Base, typename Exponent, typename RoundingStrategy>
struct power_impl1;
template<typename Base, intmax_t Exponent, bool Negative = (Exponent < 0)>
struct power_impl2;
template<typename Base, uintmax_t Exponent, bool Odd = (Exponent % 2 == 1)>
struct power_impl3;
template<typename B> struct power_impl3<B, 0, false> { typedef nat<1> type; };
template<typename B> struct power_impl3<B, 1, true> { typedef B type; };
template<typename B> struct power_impl3<B, 2, false> : multiply<B, B> {};
template<typename B, uintmax_t EvenExponent> struct power_impl3<B, EvenExponent, false> {
  typedef typename power_impl3<B, EvenExponent/2>::type half;
  typedef typename multiply<half, half>::type type;
};
template<typename B, uintmax_t OddExponent> struct power_impl3<B, OddExponent, true> {
  typedef typename power_impl3<B, OddExponent/2>::type half;
  typedef typename multiply<B, typename multiply<half, half>::type>::type type;
};
template<typename Base, intmax_t Exponent> struct power_impl2<Base, Exponent, false>
  : power_impl3<Base, (uintmax_t(Exponent))> {};

// Large exponents are large.  Even with a small base, 2,
// 2^(1000000*1000000) is sure not to fit in memory.  It has 10^12 bits:
// a bit more than 10^11 bytes, a.k.a. 100 GB.
// In practice, template-instantiation depth is likely to make the limit
// much lower.
template<typename Base, milliodigit ExponentMilliodigit, typename RS>
struct power_impl1<Base, nat<ExponentMilliodigit>, RS> : power_impl2<Base, ExponentMilliodigit> {};
template<typename Base, milliodigit ExponentMilliodigit0, milliodigit ExponentMilliodigit1, typename RS>
struct power_impl1<Base, nat<ExponentMilliodigit0, ExponentMilliodigit1>, RS>
  : power_impl2<Base, ExponentMilliodigit1*base + ExponentMilliodigit0> {};
template<typename Base, typename Exponent, typename RS>
struct power_impl1<Base, negative<Exponent>, RS>
  : reciprocal_<typename power_impl1<Base, Exponent, RS>::type> {};
template<typename Base, typename Exponent, typename RS>
struct power_impl1 {
  static_assert(sizeof(Base) && false, "Compile-time-exponentiation overflow");
};

// Certain combinations work with even huge exponents.
template<typename Exponent, typename RS> struct power<nat<>, Exponent, RS> { typedef nat<> type; };
template<typename Exponent, typename RS> struct power<nat<1>, Exponent, RS> { typedef nat<1> type; };
template<typename Base, typename RS> struct power<Base, nat<>, RS> { typedef nat<1> type; };
template<typename Exponent, typename RS> struct power<negative<nat<1>>, Exponent, RS>
  : boost::conditional<even<Exponent>::value, nat<1>, negative<nat<1>>> {};
// 0 to the 0 is more often 1 than 0.
template<typename RS> struct power<nat<>, nat<>, RS> { typedef nat<1> type; };
template<typename Base, typename Exponent, typename RoundingStrategy>
struct power : power_impl1<Base, Exponent, RoundingStrategy> {};

template<milliodigit... Milliodigits> struct numerator_<nat<Milliodigits...>> {
  typedef nat<Milliodigits...> type;
};
template<milliodigit... Milliodigits> struct denominator_<nat<Milliodigits...>> {
  typedef nat<1> type;
};

///////////////////////////
//// Natural-number truncating division
///////////////////////////

#if 0
// could be a compilation-time-optimization
template<milliodigit... Milliodigits, milliodigit Divisor>
struct divide_by_milliodigit<nat<Milliodigits...>, Divisor> {
  
};
#endif

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
// This rounds towards -Inf as ordinary twos-complement right-shifts do
template<milliodigit... Milliodigits, shift_type Shift>
struct shift_right_by_milliodigits<negative<nat<Milliodigits...>>, Shift> {
  typedef negative<typename shift_right_by_milliodigits<typename add<nat<Milliodigits...>, nat<(base-1)>>::type, Shift>::type> type;
  //static_assert(!boost::is_same<type, negative<nat<>>>::value, "bug");
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
#if 1
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
#else
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
#endif

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

template<typename NatA, typename NatB> struct divide_nat_optimize2 {
  typedef typename round_down_to_nat<typename multiply<
      floating<0, NatA>,
      typename reciprocal_nat_as_floating<
        NatB,
        milliodigit_count<NatA>::value>::type>::type>::type quot;
  typedef quot type;
  typedef typename subtract_nat<NatA, typename multiply<quot, NatB>::type>::type rem;
};
template<milliodigit MilliodigitA, milliodigit MilliodigitB>
struct divide_nat_optimize2<nat<MilliodigitA>, nat<MilliodigitB>> {
  typedef typename make_nat_from_milliodigit<(MilliodigitA / MilliodigitB)>::type quot;
  typedef quot type;
  typedef typename make_nat_from_milliodigit<(MilliodigitA % MilliodigitB)>::type rem;
};
template<typename NatA, typename NatB> struct divide_nat_optimize1
  : divide_nat_optimize2<NatA, NatB> {};
template<typename Nat> struct divide_nat_optimize1<Nat, Nat> {
  typedef nat<1> quot;
  typedef quot type;
  typedef nat<> rem;
};
template<typename NatA> struct divide_nat_optimize1<NatA, nat<1>> {
  typedef NatA quot;
  typedef quot type;
  typedef nat<> rem;
};
template<> struct divide_nat_optimize1<nat<1>, nat<1>> {
  typedef nat<1> quot;
  typedef quot type;
  typedef nat<> rem;
};
#if 0
template<typename NatA, milliodigit MilliodigitB> struct divide_nat_optimize<NatA, nat<MilliodigitB>>
  : divide_by_milliodigit<NatA, MilliodigitB> {};
#endif
template<typename NatA, typename NatB> struct divide_nat
  : divide_nat_optimize1<NatA, NatB> {};
template<typename NatA> struct divide_nat<NatA, nat<>> {
  typedef divide_by_zero quot;
  typedef quot type;
  typedef divide_by_zero rem;
};

///////////////////////////
//// Integer operations as extensions of natural number operations
///////////////////////////

// Integer subtraction
template<typename NatA, typename NatB,
  bool AIsLess = (compare<NatA, NatB>::value == -1)>
struct subtract_impl;

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct subtract_impl<nat<MilliodigitsA...>, nat<MilliodigitsB...>, true> {
  typedef negative<typename subtract_nat<nat<MilliodigitsB...>,
                                         nat<MilliodigitsA...>>::type> type;
};
template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct subtract_impl<nat<MilliodigitsA...>, nat<MilliodigitsB...>, false> {
  typedef typename subtract_nat<nat<MilliodigitsA...>,
                                nat<MilliodigitsB...>>::type type;
};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct subtract<nat<MilliodigitsA...>, nat<MilliodigitsB...>>
  : subtract_impl<nat<MilliodigitsA...>, nat<MilliodigitsB...>> {};

// Signed add/subtract
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

// Signed multiplication
template<typename T1, typename T2> struct multiply<negative<T1>, negative<T2>> : multiply<T1, T2> {};
template<typename T1, typename T2> struct multiply<negative<T1>, T2> {
  typedef negative<typename multiply<T1, T2>::type> type;
};
template<typename T1, typename T2> struct multiply<T1, negative<T2>> {
  typedef negative<typename multiply<T1, T2>::type> type;
};

// Signed comparison
template<typename T1, typename T2> struct compare<negative<T1>, negative<T2>> : compare<T2, T1> {};
template<typename T1, typename T2> struct compare<negative<T1>, T2> : boost::integral_constant<int, (-1)> {};
template<typename T1, typename T2> struct compare<T1, negative<T2>> : boost::integral_constant<int, 1> {};

// (Signed) negation
// (TODO if we get complex numbers then fix these.)
template<typename T> struct negate {
  typedef negative<T> type;
};
template<> struct negate<nat<>> {
  typedef nat<> type;
};
template<typename T> struct negate<negative<T>> {
  typedef T type;
};
template<typename T> struct abs_ {
  typedef T type;
};
template<typename T> struct abs_<negative<T>> {
  typedef T type;
};
template<typename T> struct sign_ {
  typedef nat<1> type;
};
template<> struct sign_<nat<>> {
  typedef nat<> type;
};
template<typename T> struct sign_<negative<T>> {
  typedef negative<nat<1>> type;
};

// Signed numerator/denominator
template<typename T> struct numerator_<negative<T>> {
  typedef negative<typename numerator_<T>::type> type;
};
template<typename T> struct denominator_<negative<T>> {
  typedef typename denominator_<T>::type type;
};

// Integer division following int operator/ semantics
template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<nat<MilliodigitsA...>, nat<MilliodigitsB...>>
  : divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> {};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<negative<nat<MilliodigitsA...>>, nat<MilliodigitsB...>> {
  typedef divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> unsigned_result_;
  typedef typename negate<typename unsigned_result_::quot>::type type;
  typedef type quot;
  typedef typename negate<typename unsigned_result_::rem>::type rem;
};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<nat<MilliodigitsA...>, negative<nat<MilliodigitsB...>>> {
  typedef divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> unsigned_result_;
  typedef typename negate<typename unsigned_result_::quot>::type type;
  typedef type quot;
  typedef typename unsigned_result_::rem rem;
};

template<milliodigit... MilliodigitsA, milliodigit... MilliodigitsB>
struct divide_integer<negative<nat<MilliodigitsA...>>, negative<nat<MilliodigitsB...>>> {
  typedef divide_nat<nat<MilliodigitsA...>, nat<MilliodigitsB...>> unsigned_result_;
  typedef typename unsigned_result_::quot type;
  typedef type quot;
  typedef typename negate<typename unsigned_result_::rem>::type rem;
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
struct make_rational_from_coprime_positive_nats {
  typedef rational<Num, Den> type;
};
template<typename Num>
struct make_rational_from_coprime_positive_nats<Num, nat<1>> {
  typedef Num type;
};
template<typename Num, typename Den>
struct make_rational {
  typedef typename gcd<Num, Den>::type gcd_;
  typedef typename make_rational_from_coprime_positive_nats<
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

template<typename T>
struct reciprocal_<negative<T>> {
  typedef negative<typename reciprocal_<T>::type> type;
};
template<typename Num, typename Den>
struct reciprocal_<rational<Num, Den>> {
  typedef rational<Den, Num> type;
};
template<typename Den>
struct reciprocal_<rational<nat<1>, Den>> {
  typedef Den type;
};
template<milliodigit...Milliodigits>
struct reciprocal_<nat<Milliodigits...>> {
  typedef rational<nat<1>, nat<Milliodigits...>> type;
};
template<>
struct reciprocal_<nat<1>> {
  typedef nat<1> type;
};
template<>
struct reciprocal_<nat<>> {
  // error: divide by zero
  typedef divide_by_zero type;
};

template<typename NumA, typename DenA, typename NumB, typename DenB>
struct compare<rational<NumA, DenA>, rational<NumB, DenB>>
  : compare<
      typename multiply<NumA, DenB>::type,
      typename multiply<NumB, DenA>::type> {};

template<typename Num, typename Den>
struct numerator_<rational<Num, Den>> {
  typedef Num type;
};
template<typename Num, typename Den>
struct denominator_<rational<Num, Den>> {
  typedef Den type;
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
template<typename NumA, typename DenA, milliodigit...MilliodigitsB>
struct compare<rational<NumA, DenA>, nat<MilliodigitsB...>>
  : compare<rational<NumA, DenA>, typename as_rational<nat<MilliodigitsB...>>::type> {};

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
template<milliodigit...MilliodigitsA, typename NumB, typename DenB>
struct compare<nat<MilliodigitsA...>, rational<NumB, DenB>>
  : compare<typename as_rational<nat<MilliodigitsA...>>::type, rational<NumB, DenB>> {};

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

template<typename Base, intmax_t Exponent> struct power_impl2<Base, Exponent, true>
  : reciprocal_<typename power_impl3<Base, (-uintmax_t(Exponent))>::type> {};

///////////////////////////
//// Miscellaneous
///////////////////////////

// Rounding to integers

template<typename Num, rounding_strategy_for_positive_numbers Strategy> struct round_nonnegative;

#if 0
template<milliodigit...Milliodigits, typename RoundingStrategy>
struct round_<nat<Milliodigits...>, RoundingStrategy> {
  typedef nat<Milliodigits...> type;
};
template<milliodigit...Milliodigits, typename RoundingStrategy>
struct round_<negative<nat<Milliodigits...>>, RoundingStrategy> {
  typedef negative<nat<Milliodigits...>> type;
};
#endif

template<bool Increment, typename Num> struct increment_if;
template<typename Num> struct increment_if<true, Num> : add<Num, nat<1>> {};
template<typename Num> struct increment_if<false, Num> {
  typedef Num type;
};
template<bool Increment, typename Num> struct decrement_if;
template<typename Num> struct decrement_if<true, Num> : subtract<Num, nat<1>> {};
template<typename Num> struct decrement_if<false, Num> {
  typedef Num type;
};

template<milliodigit...Milliodigits, rounding_strategy_for_positive_numbers Strategy>
struct round_nonnegative<nat<Milliodigits...>, Strategy> {
  typedef nat<Milliodigits...> type;
};
template<typename Num, typename Den>
struct round_nonnegative<rational<Num, Den>, require_exact_answer> {
  static_assert(sizeof(Num) && false, "Rounding a number that was required to be exact.");
};
template<typename Num>
struct round_nonnegative<rational<Num, nat<2>>, round_to_nearest_with_ties_rounding_down> {
  typedef typename divide_nat<Num, nat<2>>::type type;
};
template<typename Num>
struct round_nonnegative<rational<Num, nat<2>>, round_to_nearest_with_ties_rounding_up> {
  typedef typename add<typename divide_nat<Num, nat<2>>::type, nat<1>>::type type;
};
template<typename Num>
struct round_nonnegative<rational<Num, nat<2>>, round_to_nearest_with_ties_rounding_to_even> {
  typedef typename divide_nat<Num, nat<2>>::type div_result_;
  typedef typename increment_if<odd<div_result_>::value, div_result_>::type type;
};
template<typename Num>
struct round_nonnegative<rational<Num, nat<2>>, round_to_nearest_with_ties_rounding_to_odd> {
  typedef typename divide_nat<Num, nat<2>>::type div_result_;
  typedef typename increment_if<even<div_result_>::value, div_result_>::type type;
};
template<typename Num>
struct round_nonnegative<rational<Num, nat<2>>, round_to_nearest_with_ties_rounding_to_halfway> {
  typedef rational<Num, nat<2>> type;
};

template<typename Num, typename Den>
struct round_nonnegative<rational<Num, Den>, round_inexact_to_halfway> {
  typedef typename add<typename divide_nat<Num, Den>::type, rational<nat<1>, nat<2>>>::type type;
};

template<typename Num, typename Den>
struct round_nonnegative<rational<Num, Den>, round_down> {
  typedef typename divide_nat<Num, Den>::type type;
};
template<typename Num, typename Den>
struct round_nonnegative<rational<Num, Den>, round_up> {
  typedef typename add<typename divide_nat<Num, Den>::type, nat<1>>::type type;
};
template<typename Num, typename Den, rounding_strategy_for_positive_numbers RoundToNearestSomething>
struct round_nonnegative<rational<Num, Den>, RoundToNearestSomething>
  : round_nonnegative<typename add<rational<Num, Den>, rational<nat<1>, nat<2>>>::type, round_down> {};


template<typename Num, rounding_strategy_for_positive_numbers Strategy>
struct round_<Num, rounding_strategy<Strategy, negative_is_forbidden>> {
  static_assert(is_nonnegative_<Num>::value, "Negative number appeared where explicitly forbidden");
  typedef typename round_nonnegative<Num, Strategy>::type type;
};
template<typename Num, rounding_strategy_for_positive_numbers Strategy>
struct round_<Num, rounding_strategy<Strategy, negative_variant_doesnt_make_a_difference>> {
  typedef typename round_nonnegative<Num, Strategy>::type type;
};
template<typename Num, rounding_strategy_for_positive_numbers Strategy>
struct round_<negative<Num>, rounding_strategy<Strategy, negative_variant_doesnt_make_a_difference>> {
  static_assert(Strategy >= round_to_nearest_with_ties_rounding_to_even,
                "You lied! The negative variant does make a difference.");
  typedef typename negate<typename round_nonnegative<Num, Strategy>::type>::type type;
};
template<typename Num, rounding_strategy_for_positive_numbers Strategy>
struct round_<negative<Num>, rounding_strategy<Strategy, negative_mirrors_positive>>
  : negate<typename round_nonnegative<Num, Strategy>::type> {};
template<typename Num>
struct round_<negative<Num>, rounding_strategy<round_up, negative_continuous_with_positive>>
  : negate<typename round_nonnegative<Num, round_down>::type> {};
template<typename Num>
struct round_<negative<Num>, rounding_strategy<round_down, negative_continuous_with_positive>>
  : negate<typename round_nonnegative<Num, round_up>::type> {};
template<typename Num>
struct round_<negative<Num>, rounding_strategy<round_to_nearest_with_ties_rounding_up, negative_continuous_with_positive>>
  : negate<typename round_nonnegative<Num, round_to_nearest_with_ties_rounding_down>::type> {};
template<typename Num>
struct round_<negative<Num>, rounding_strategy<round_to_nearest_with_ties_rounding_down, negative_continuous_with_positive>>
  : negate<typename round_nonnegative<Num, round_to_nearest_with_ties_rounding_up>::type> {};
template<typename Num, rounding_strategy_for_positive_numbers PositiveStrategy,
  rounding_strategy_for_negative_numbers NegativeStrategy>
struct round_<Num, rounding_strategy<PositiveStrategy, NegativeStrategy>> {
  typedef typename round_nonnegative<Num, PositiveStrategy>::type type;
};

///////////////////////////
//// Logarithm
///////////////////////////

struct undefined_logarithm {};

template<typename Base, typename Number>
struct log_impl1;

template<typename Base, typename Number, typename RoundingStrategy>
struct logarithm {
  typedef typename round_<typename log_impl1<Base, Number>::type, RoundingStrategy>::type type;
};

// This function's domain is {x | x > 0 and x != 1}
template<typename Argument>
struct log_argument_is_less_than_1;
template<milliodigit...Milliodigits>
struct log_argument_is_less_than_1<nat<Milliodigits...>> : boost::false_type {};
template<typename Num, typename Den>
struct log_argument_is_less_than_1<rational<Num, Den>>
  : boost::integral_constant<bool, (compare<Num, Den>::value == -1)> {};

template<typename Base, typename Number>
struct log_impl2;
template<typename Base, typename Number,
  bool BaseIsLessThanOne = log_argument_is_less_than_1<Base>::value,
  bool NumberIsLessThanOne = log_argument_is_less_than_1<Number>::value
  >
struct log_impl3;
template<typename Base, typename Number>
struct log_impl4;

// Neither the base nor the number may be negative.
template<typename Base, typename Number>
struct log_impl1 : log_impl2<Base, Number> {};
template<typename Base, typename Number>
struct log_impl1<negative<Base>, Number> {
  typedef undefined_logarithm type;
};
template<typename Base, typename Number>
struct log_impl1<negative<Base>, negative<Number>> {
  typedef undefined_logarithm type;
};
template<typename Base, typename Number>
struct log_impl1<Base, negative<Number>> {
  typedef undefined_logarithm type;
};

// Neither the base nor the number may zero.
// If the number is 1 and the base is valid,
// the exponent is 0.
template<typename Base, typename Number>
struct log_impl2 : log_impl3<Base, Number> {};
template<typename Number>
struct log_impl2<nat<>, Number> {
  typedef undefined_logarithm type;
};
template<typename Number>
struct log_impl2<nat<1>, Number> {
  typedef undefined_logarithm type;
};
template<typename Base>
struct log_impl2<Base, nat<1>> {
  typedef nat<> type;
};
template<typename Base>
struct log_impl2<Base, nat<>> {
  typedef undefined_logarithm type;
};
template<>
struct log_impl2<nat<1>, nat<>> {
  typedef undefined_logarithm type;
};
// Should these be nat<> or undefined_logarithm?
template<>
struct log_impl2<nat<1>, nat<1>> {
  typedef undefined_logarithm type;
};
template<>
struct log_impl2<nat<>, nat<1>> {
  typedef undefined_logarithm type;
};
// nat<1> or undefined_logarithm?
template<>
struct log_impl2<nat<>, nat<>> {
  typedef undefined_logarithm type;
};

template<typename Base, typename Number>
struct log_impl3<Base, Number, false, false>
  : log_impl4<Base, Number> {};
template<typename Base, typename Number>
struct log_impl3<Base, Number, false, true> {
  typedef typename negate<
    typename log_impl4<Base, typename reciprocal_<Number>::type>::type
    >::type type;
};
template<typename Base, typename Number>
struct log_impl3<Base, Number, true, false> {
  typedef typename negate<
    typename log_impl4<typename reciprocal_<Base>::type, Number>::type
    >::type type;
};
template<typename Base, typename Number>
struct log_impl3<Base, Number, true, true> {
  typedef
    typename log_impl4<typename reciprocal_<Base>::type,
                       typename reciprocal_<Number>::type>::type
    type;
};

// Now Base and Number are both > 1.

// If Number < Base... or ==...

// 0 -> result
// -1 -> backtrack
// 1 -> recur
template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename TestAccum = typename multiply<TestAccumPrev, TestAccumPrev>::type,
  typename ExponentAccum = typename twice<ExponentAccumPrev>::type,
  int Compare = compare<Number, TestAccum>::value>
struct log_impl5;

template<typename Base, typename Number>
struct log_impl4 : log_impl5<Base, Number, nat<1>, nat<>, Base, nat<1>> {};

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename TestAccum,
  typename ExponentAccum
> struct log_impl5<
  Base, Number, TestAccumPrev, ExponentAccumPrev,
  TestAccum, ExponentAccum, 1
> : log_impl5<Base, Number, TestAccum, ExponentAccum> {};

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename TestAccum,
  typename ExponentAccum
> struct log_impl5<
  Base, Number, TestAccumPrev, ExponentAccumPrev,
  TestAccum, ExponentAccum, 0
> {
  typedef ExponentAccum type;
};

// 0 -> result
// -1 -> recur without keep
// 1 -> recur with keep
// also TestAccum 0 means stop
template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename ExponentAccumTestSizePrev,
  typename ExponentAccumTestSize = typename half_nat_rounding_down<ExponentAccumTestSizePrev>::type,
  typename TestAccum = typename multiply<TestAccumPrev, typename power<Base, ExponentAccumTestSize>::type>::type,
  typename ExponentAccum = typename add<ExponentAccumPrev, ExponentAccumTestSize>::type,
  int Compare = compare<Number, TestAccum>::value,
  bool NoMoreExponentsToTest = is_zero_<ExponentAccumTestSize>::value>
struct log_impl6;

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename TestAccum,
  typename ExponentAccum
> struct log_impl5<
  Base, Number, TestAccumPrev, ExponentAccumPrev,
  TestAccum, ExponentAccum, (-1)
> : log_impl6<Base, Number, TestAccumPrev, ExponentAccumPrev,
      ExponentAccumPrev> {};

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename ExponentAccumTestSizePrev,
  typename ExponentAccumTestSize,
  typename TestAccum,
  typename ExponentAccum
> struct log_impl6<
  Base, Number, TestAccumPrev, ExponentAccumPrev, ExponentAccumTestSizePrev,
  ExponentAccumTestSize, TestAccum, ExponentAccum, 0, false
> {
  typedef ExponentAccum type;
};

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename ExponentAccumTestSizePrev,
  typename ExponentAccumTestSize,
  typename TestAccum,
  typename ExponentAccum
> struct log_impl6<
  Base, Number, TestAccumPrev, ExponentAccumPrev, ExponentAccumTestSizePrev,
  ExponentAccumTestSize, TestAccum, ExponentAccum, 1, false
> : log_impl6<Base, Number, TestAccum, ExponentAccum, ExponentAccumTestSize> {};

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename ExponentAccumTestSizePrev,
  typename ExponentAccumTestSize,
  typename TestAccum,
  typename ExponentAccum
> struct log_impl6<
  Base, Number, TestAccumPrev, ExponentAccumPrev, ExponentAccumTestSizePrev,
  ExponentAccumTestSize, TestAccum, ExponentAccum, (-1), false
> : log_impl6<Base, Number, TestAccumPrev, ExponentAccumPrev, ExponentAccumTestSize> {};

template<
  typename Base,
  typename Number,
  typename TestAccumPrev,
  typename ExponentAccumPrev,
  typename ExponentAccumTestSizePrev,
  typename ExponentAccumTestSize,
  typename TestAccum,
  typename ExponentAccum,
  int Compare
> struct log_impl6<
  Base, Number, TestAccumPrev, ExponentAccumPrev, ExponentAccumTestSizePrev,
  ExponentAccumTestSize, TestAccum, ExponentAccum, Compare, true
> {
  typedef typename add<ExponentAccumPrev, rational<nat<1>, nat<2>>>::type type;
};


template<typename Base, typename NatRoot> struct root_impl2;
template<typename Base, typename NatRoot, typename Estimate> struct root_impl3;
template<typename Base, typename NatRoot, typename OldEstimate, typename Delta> struct root_impl4;
template<typename Estimate, int EstimateVsActual> struct root_impl5;

template<typename Base, typename NatRoot, typename RoundingStrategy, bool OddRoot = odd<NatRoot>::value>
struct root_impl1 {
  typedef typename round_<
    typename root_impl2<Base, NatRoot>::type
  , RoundingStrategy>::type type;
};
template<typename Base, typename NatRoot, typename RoundingStrategy>
struct root_impl1<negative<Base>, NatRoot, RoundingStrategy, true> {
  typedef typename round_<
    negative<typename root_impl2<Base, NatRoot>::type>
  , RoundingStrategy>::type type;
};
template<typename Base, typename NatRoot, typename RoundingStrategy>
struct root_impl1<negative<Base>, NatRoot, RoundingStrategy, false> {
  typedef even_root_of_negative type;
};
template<typename NatRoot, typename RoundingStrategy, bool OddRoot>
struct root_impl1<nat<>, NatRoot, RoundingStrategy, OddRoot> {
  typedef nat<> type;
};

template<typename BaseNum, typename BaseDen, typename NatRoot, bool OddRoot>
struct root_impl1<rational<BaseNum, BaseDen>, NatRoot, rounding_strategy<require_exact_answer>, OddRoot> {
  typedef typename root_impl2<BaseNum, NatRoot>::type num_root_;
  typedef typename root_impl2<BaseDen, NatRoot>::type den_root_;
  static_assert(compare<BaseNum, typename power<num_root_, NatRoot>::type>::value
    == 0, "inexact root that was required to be exact");
  //static_assert(den_root_() == 5789345789, "dfsdfshk");
  static_assert(compare<BaseDen, typename power<den_root_, NatRoot>::type>::value
    == 0, "inexact root that was required to be exact");
  typedef rational<num_root_, den_root_> type;
};
template<typename BaseNum, typename BaseDen, typename NatRoot>
struct root_impl1<negative<rational<BaseNum, BaseDen>>, NatRoot, rounding_strategy<require_exact_answer>, true> {
  typedef typename root_impl2<BaseNum, NatRoot>::type num_root_;
  typedef typename root_impl2<BaseDen, NatRoot>::type den_root_;
  static_assert(compare<BaseNum, typename power<num_root_, NatRoot>::type>::value
    == 0, "inexact root that was required to be exact");
  static_assert(compare<BaseDen, typename power<den_root_, NatRoot>::type>::value
    == 0, "inexact root that was required to be exact");
  typedef negative<rational<num_root_, den_root_>> type;
};

template<typename Base, typename NatRoot>
struct root_impl2 {
  typedef typename logarithm<nat<2>, Base, rounding_strategy<round_down>>::type log2_base_;
  typedef typename power<nat<2>,
    //typename round<
      typename divide_nat<log2_base_, NatRoot>::type
    //, rounding_strategy<round_down>>::type
      >::type initial_estimate_;
  typedef typename root_impl3<Base, NatRoot, initial_estimate_>::type type;
};
template<typename Base, typename NatRoot, typename Estimate>
struct root_impl3 {
  typedef typename round_<
    typename divide_rational<
      typename subtract<typename power<Estimate, NatRoot>::type, Base>::type,
      typename multiply<NatRoot,
        typename power<Estimate, typename subtract<NatRoot, nat<1>>::type>::type>::type
    >::type,
    rounding_strategy<round_up, negative_continuous_with_positive> >::type new_delta_;
  typedef typename root_impl4<Base, NatRoot, Estimate, new_delta_>::type type;
};
template<typename Base, typename NatRoot, typename OldEstimate, typename Delta>
struct root_impl4
: root_impl3<Base, NatRoot, typename subtract<OldEstimate, Delta>::type> {};
template<typename Base, typename NatRoot, typename OldEstimate>
struct root_impl4<Base, NatRoot, OldEstimate, nat<>>
: root_impl5<OldEstimate, compare<typename power<OldEstimate, NatRoot>::type, Base>::value> {};
template<typename Estimate>
struct root_impl5<Estimate, 0> {
  typedef Estimate type;
};
template<typename Estimate>
struct root_impl5<Estimate, (-1)> {
  //typedef typename add<Estimate, nat<1>>::type may_be_exact_;
  //typename power<may_be_exact_, NatRoot>::type
  typedef typename add<Estimate, rational<nat<1>, nat<2>>>::type type;
};
template<typename Estimate>
struct root_impl5<Estimate, 1> {
  typedef typename subtract<Estimate, rational<nat<1>, nat<2>>>::type type;
};

template<typename Num, typename Root, typename RoundingStrategy> struct root
  : power<Num, typename reciprocal_<Root>::type, RoundingStrategy> {};

template<typename Base, typename ExpNum, typename ExpDen, typename RoundingStrategy>
struct power_impl1<Base, rational<ExpNum, ExpDen>, RoundingStrategy>
  : root_impl1<typename power_impl1<Base, ExpNum, RoundingStrategy>::type, ExpDen, RoundingStrategy> {};

template<typename Factor, typename Factoree,
  typename AccumExponent = nat<>,
  typename Quot = typename divide_nat<Factoree, Factor>::quot,
  typename Rem = typename divide_nat<Factoree, Factor>::rem>
struct extract_factor_impl {
  typedef Factoree rest_of_factoree;
  typedef AccumExponent factor_exponent;
};
template<typename Factor, typename Factoree, typename AccumExponent, typename Quot>
struct extract_factor_impl<Factor, Factoree, AccumExponent, Quot, nat<>>
: extract_factor_impl<Factor, Quot, typename add<AccumExponent, nat<1>>::type> {};

template<typename Factor, typename Factoree>
struct extract_factor_ : extract_factor_impl<Factor, Factoree> {};
template<typename Factor, typename FactoreeNum, typename FactoreeDen>
struct extract_factor_<Factor, rational<FactoreeNum, FactoreeDen>> {
  typedef extract_factor_<Factor, FactoreeNum> extracted_num_;
  typedef extract_factor_<Factor, FactoreeDen> extracted_den_;
  typedef typename make_rational_from_coprime_positive_nats<
    typename extracted_num_::rest_of_factoree,
    typename extracted_den_::rest_of_factoree>::type rest_of_factoree;
  typedef typename subtract<
    typename extracted_num_::factor_exponent,
    typename extracted_den_::factor_exponent>::type factor_exponent;
};
template<typename Factor, typename Factoree>
struct extract_factor_<Factor, negative<Factoree>> {
  typedef extract_factor_<Factor, Factoree> extracted_;
  typedef negative<typename extracted_::rest_of_factoree> rest_of_factoree;
  typedef typename extracted_::factor_exponent factor_exponent;
};


template<uintmax_t NumRemainingDigits, typename RemainingNat,
  typename UInt, typename Base, typename DigitAppender, typename Accum = array<UInt> >
struct convert_impl;

template<typename RemainingNat,
  typename UInt, typename Base, typename DigitAppender, typename Accum>
struct convert_impl<0, RemainingNat, UInt, Base, DigitAppender, Accum> {
  typedef Accum type;
  static const bool overflow = true;
};
template<typename UInt, typename Base, typename DigitAppender, typename Accum>
struct convert_impl<0, nat<>, UInt, Base, DigitAppender, Accum> {
  typedef Accum type;
  static const bool overflow = false;
};
template<uintmax_t NumRemainingDigits, typename RemainingNat,
  typename UInt, typename Base, typename DigitAppender, typename Accum>
struct convert_impl//<NumRemainingDigits, RemainingNat, UInt, Base, DigitAppender, Accum>
: convert_impl<
    (NumRemainingDigits-1),
    typename divide_nat<RemainingNat, Base>::quot,
    UInt, Base, DigitAppender,
    typename DigitAppender::template append<
      to_int_<UInt, typename divide_nat<RemainingNat, Base>::rem>::value(),
      Accum>::type
    > {};

template<uintmax_t NumRemainingDigits,
  typename UInt, typename Base, typename DigitAppender, typename Accum>
struct convert_impl<NumRemainingDigits, nat<>, UInt, Base, DigitAppender, Accum>
: convert_impl<
    (NumRemainingDigits-1),
    nat<>, UInt, Base, DigitAppender,
    typename DigitAppender::template append<UInt(0), Accum>::type
    > {};


template<typename UInt, bool BigEndian> struct convert_impl_append_to_positive;

template<typename UInt> struct convert_impl_append_to_positive<UInt, true> {
  template<UInt Value, typename Array> struct append;
  template<UInt Value, UInt...Values>
  struct append<Value, array<UInt, Values...>> {
    typedef array<UInt, Value, Values...> type;
  };
};
template<typename UInt> struct convert_impl_append_to_positive<UInt, false> {
  template<UInt Value, typename Array> struct append;
  template<UInt Value, UInt...Values>
  struct append<Value, array<UInt, Values...>> {
    typedef array<UInt, Values..., Value> type;
  };
};

template<typename UInt, UInt BaseMinusOne, bool BigEndian> struct convert_impl_append_to_negative;
template<typename UInt, UInt BaseMinusOne> struct convert_impl_append_to_negative<UInt, BaseMinusOne, true> {
  template<UInt Value, typename Array> struct append;
  template<UInt Value, UInt...Values>
  struct append<Value, array<UInt, Values...>> {
    typedef array<UInt, (BaseMinusOne - Value), Values...> type;
  };
};
template<typename UInt, UInt BaseMinusOne> struct convert_impl_append_to_negative<UInt, BaseMinusOne, false> {
  template<UInt Value, typename Array> struct append;
  template<UInt Value, UInt...Values>
  struct append<Value, array<UInt, Values...>> {
    typedef array<UInt, Values..., (BaseMinusOne - Value)> type;
  };
};



// i suspect it should have two interfaces, auto-sized and manually-sized,
// user wise
template<
  typename Integer,
  typename UInt,
  bool BigEndian, // the order of the UInt digits; each UInt is internally the native bit-order
  bool Signed, // represented as radix complement, typically two's complement
  bool AutoSized = true,
  bool ZeroHasADigit = true /*relevant if AutoSized*/,
  bool ModuloIfOverflow = false /*relevant if !AutoSized*/,
  uintmax_t ExactSize = 0 /*used if !AutoSized*/,
  UInt BaseMinusOne = std::numeric_limits<UInt>::max()>
struct convert_to_base {
  static_assert(sizeof(Integer) && false, "convert_to_base not valid for non-integers");
};
template<typename UInt, bool BigEndian, bool Signed,
  bool ModuloIfOverflow, uintmax_t ExactSize, UInt BaseMinusOne>
struct convert_to_base<nat<>, UInt, BigEndian, Signed, true, true,
    ModuloIfOverflow, ExactSize, BaseMinusOne> {
  typedef array<UInt, UInt(0)> type;
};
template<typename UInt, bool BigEndian, bool Signed,
  bool ModuloIfOverflow, uintmax_t ExactSize, UInt BaseMinusOne>
struct convert_to_base<nat<>, UInt, BigEndian, Signed, true, false,
    ModuloIfOverflow, ExactSize, BaseMinusOne> {
  typedef array<UInt> type;
};
template<milliodigit...Milliodigits, typename UInt,
  bool BigEndian, bool Signed, bool AutoSized, bool ZeroHasADigit,
  bool ModuloIfOverflow, uintmax_t ExactSize, UInt BaseMinusOne>
struct convert_to_base<nat<Milliodigits...>, UInt,
    BigEndian, Signed, AutoSized, ZeroHasADigit, ModuloIfOverflow,
    ExactSize, BaseMinusOne> {
  typedef nat<Milliodigits...> Nat;
  typedef typename add<typename uinteger_literal<BaseMinusOne>::type, nat<1>>::type Base;
  typedef typename boost::conditional<
    Signed, typename twice<Nat>::type, Nat>::type space_needed_for_sign_;
  static const uintmax_t num_native_digits_ = to_int_<uintmax_t,
    typename logarithm<Base, typename add<space_needed_for_sign_, nat<1>>::type, rounding_strategy<round_up>>::type
    >::value();
  static_assert(AutoSized || ModuloIfOverflow || num_native_digits_ <= ExactSize,
    "overflow when converting compile-time-integer");
  static const uintmax_t num_digits_ = (AutoSized ? num_native_digits_ : ExactSize);
  typedef typename convert_impl<num_digits_, Nat, UInt, Base,
    convert_impl_append_to_positive<UInt, BigEndian> >::type type;
};
template<milliodigit...Milliodigits, typename UInt,
  bool BigEndian, bool Signed, bool AutoSized, bool ZeroHasADigit,
  bool ModuloIfOverflow, uintmax_t ExactSize, UInt BaseMinusOne>
struct convert_to_base<negative<nat<Milliodigits...>>, UInt,
    BigEndian, Signed, AutoSized, ZeroHasADigit, ModuloIfOverflow,
    ExactSize, BaseMinusOne> {
  static_assert(Signed, "unsigned convert_to_base called on negative integer");
  typedef nat<Milliodigits...> Nat;
  typedef typename add<typename uinteger_literal<BaseMinusOne>::type, nat<1>>::type Base;
  typedef typename subtract<Nat, nat<1>>::type complement_value_;
  typedef typename twice<complement_value_>::type space_needed_for_sign_;
  static const uintmax_t num_native_digits_a_ = to_int_<uintmax_t,
    typename logarithm<Base, typename add<space_needed_for_sign_, nat<1>>::type, rounding_strategy<round_up>>::type
    >::value();
  static const uintmax_t num_native_digits_ =
    ((num_native_digits_a_ == 0) ? 1 : num_native_digits_a_);
  static_assert(AutoSized || ModuloIfOverflow || num_native_digits_ <= ExactSize,
    "overflow when converting compile-time-integer");
  static const uintmax_t num_digits_ = (AutoSized ? num_native_digits_ : ExactSize);
  typedef typename convert_impl<num_digits_, complement_value_, UInt, Base,
    convert_impl_append_to_negative<UInt, BaseMinusOne, BigEndian> >::type type;
};

//if twoscompl,
//when making unsigned pre compl val
//first digit must be 0
//(..probably twoscompl only makes sense for power-of-two bases)
//oh ah compl: (base-1)-val
//and the +1 / -1 overall for negation stays.
//  http://homepage.cs.uiowa.edu/~jones/ternary/numbers.shtml
//how do we decide whether it's negative
//okay, staying with "two's complement" and only power-of-two negative numbers
//handled for now.

template<typename BasePrefixArray, typename Array>
struct show_as_base_impl_nat;
template<char...Prefix, uint8_t...Digit>
struct show_as_base_impl_nat<array<char, Prefix...>, array<uint8_t, Digit...>> {
  typedef array<char, Prefix..., ((Digit<10) ? ('0'+char(Digit)) : ('a'+char(Digit-10)))...> type;
};

template<typename Integer, uintmax_t Base, char...BasePrefix>
struct show_as_base
: show_as_base_impl_nat<
    array<char, BasePrefix...>,
    typename convert_to_base<Integer, uint8_t, true, false, true, true, false, 0, Base-1>::type
> {};
template<typename Nat, uintmax_t Base, char...BasePrefix>
struct show_as_base<negative<Nat>, Base, BasePrefix...>
: concat<array<char, '-'>, typename show_as_base<Nat, Base, BasePrefix...>::type> {};
template<typename Num, typename Den, uintmax_t Base, char...BasePrefix>
struct show_as_base<rational<Num, Den>, Base, BasePrefix...>
: concat<typename show_as_base<Num, Base, BasePrefix...>::type, array<char, '/'>,
         typename show_as_base<Den, Base, BasePrefix...>::type> {};
template<typename Num, typename Den, uintmax_t Base, char...BasePrefix>
struct show_as_base<negative<rational<Num, Den>>, Base, BasePrefix...>
: concat<array<char, '-'>,
         typename show_as_base<Num, Base, BasePrefix...>::type, array<char, '/'>,
         typename show_as_base<Den, Base, BasePrefix...>::type> {};

template<char...Digits> struct parse_nat;

// partially based on open-std document N3531 2013-03-08
// http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2013/n3531.pdf
template <unsigned Base, typename Accumulator, char... Digits> struct parse_digits;
template <char... Digits> struct parse_nat;
template <char... Digits> struct parse_nat<'0','x',Digits...> : parse_digits<16, nat<>, Digits...> {};
template <char... Digits> struct parse_nat<'0','X',Digits...> : parse_digits<16, nat<>, Digits...> {};
template <char... Digits> struct parse_nat<'0',Digits...> : parse_digits<8, nat<>, Digits...> {};
template <char... Digits> struct parse_nat : parse_digits<10, nat<>, Digits...> {};

template <char... Digits> struct parse_integer : parse_nat<Digits...> {};
template <char... Digits> struct parse_integer<'-',Digits...> {
  typedef negative<typename parse_nat<Digits...>::type> type;
};

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

struct make_any_number {
  template<typename T> constexpr inline operator T()const { return T(); }
};
} /* end namespace ct::impl */


// ADL (Argument Dependent Lookup) will find these templates when appropriate.
// operator/ is rational division; div() is integer division.

template<typename A> constexpr inline number<A>
operator+(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline number<typename impl::negate<A>::type>
operator-(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline number<typename impl::abs_<A>::type>
abs(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline number<typename impl::sign_<A>::type>
sign(number<A>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline number<typename impl::add<A, B>::type>
operator+(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline number<typename impl::subtract<A, B>::type>
operator-(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline number<typename impl::multiply<A, B>::type>
operator*(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline number<typename impl::divide_rational<A, B>::type>
operator/(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A> constexpr inline number<typename impl::reciprocal_<A>::type>
reciprocal(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline number<typename impl::numerator_<A>::type>
numerator(number<A>) { return impl::make_any_number(); }
template<typename A> constexpr inline number<typename impl::denominator_<A>::type>
denominator(number<A>) { return impl::make_any_number(); }

template<typename Quot, typename Rem>
struct ctdiv_t {
  Quot quot;
  Rem rem;
};
template<typename A, typename B> constexpr inline
ctdiv_t<
  number<typename impl::divide_integer<A, B>::quot>,
  number<typename impl::divide_integer<A, B>::rem> >
div(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline number<typename impl::power<A, B>::type>
pow(number<A>, number<B>) { return impl::make_any_number(); }

template<typename A, typename B> constexpr inline bool
operator==(number<A>, number<B>) { return impl::compare<A, B>::value == 0; }
template<typename A, typename B> constexpr inline bool
operator!=(number<A>, number<B>) { return impl::compare<A, B>::value != 0; }
template<typename A, typename B> constexpr inline bool
operator<(number<A>, number<B>) { return impl::compare<A, B>::value == -1; }
template<typename A, typename B> constexpr inline bool
operator>(number<A>, number<B>) { return impl::compare<A, B>::value == 1; }
template<typename A, typename B> constexpr inline bool
operator<=(number<A>, number<B>) { return impl::compare<A, B>::value != 1; }
template<typename A, typename B> constexpr inline bool
operator>=(number<A>, number<B>) { return impl::compare<A, B>::value != -1; }

template<typename A> constexpr inline bool
is_integer(number<A>) { return impl::is_integer_<A>::value; }
template<typename A> constexpr inline bool
is_nonnegative_integer(number<A>) { return impl::is_nonnegative_integer_<A>::value; }
template<typename A> constexpr inline bool
is_positive_integer(number<A>) { return impl::is_positive_integer_<A>::value; }
template<typename A> constexpr inline bool
is_nonnegative(number<A>) { return impl::is_nonnegative_<A>::value; }
template<typename A> constexpr inline bool
is_positive(number<A>) { return impl::is_positive_<A>::value; }
template<typename A> constexpr inline bool
is_negative(number<A>) { return impl::is_negative_<A>::value; }
template<typename A> constexpr inline bool
is_zero(number<A>) { return impl::is_zero_<A>::value; }

template<uintmax_t Int>
struct from_uint {
  typedef number<typename impl::literal<uintmax_t, Int>::type> type;
  static constexpr type value = type();
};
template<intmax_t Int>
struct from_int {
  typedef number<typename impl::literal<intmax_t, Int>::type> type;
  static constexpr type value = type();
};

template<typename IntegralType, typename A> constexpr inline IntegralType
to_int(number<A>) {
  static_assert(impl::is_integer_<A>::value, "to_int on non-integer");
  // TODO assert no overflow?
  // TODO choose smallest reasonable int type by default?
  return impl::to_int_<IntegralType, A>::value();
}

template<
  typename UInt,
  bool BigEndian, // the order of the UInt digits; each UInt is internally the native bit-order
  bool Signed, // represented as radix complement, typically two's complement
  bool ZeroHasADigit = true,
  UInt BaseMinusOne = std::numeric_limits<UInt>::max(),
  typename A>
constexpr inline typename
impl::convert_to_base<A, UInt, BigEndian, Signed,
    true, ZeroHasADigit, false, 0, BaseMinusOne>::type
convert_to_base_auto_sized(number<A>) { return impl::make_any_number(); }

template<
  typename UInt,
  bool BigEndian, // the order of the UInt digits; each UInt is internally the native bit-order
  bool Signed, // represented as radix complement, typically two's complement
  uintmax_t ExactSize,
  bool ModuloIfOverflow = false /*non-modulo: error if overflow*/,
  UInt BaseMinusOne = std::numeric_limits<UInt>::max(),
  typename A>
constexpr inline typename
impl::convert_to_base<A, UInt, BigEndian, Signed,
    false, false, ModuloIfOverflow, ExactSize, BaseMinusOne>::type
convert_to_base_manually_sized(number<A>) { return impl::make_any_number(); }
//TODO document these

// TODO allow e.g. 0x prefixes in the interface somehow
// TODO operator<<(std::ostream, number<>)
template<int Base = 10, typename A>
constexpr inline typename
impl::show_as_base<A, Base>::type
to_chars(number<A>) { return impl::make_any_number(); }

template<int Base = 10, typename A>
constexpr inline const char*
to_string(number<A> a) { return to_string(to_chars<Base>(a)); }

template<typename FactorExponent, typename RestOfFactoree>
struct extracted_factor {
  number<FactorExponent> factor_exponent;
  number<RestOfFactoree> rest_of_factoree;
};

template<typename Factor, typename Factoree>
constexpr inline
extracted_factor<
  typename impl::extract_factor_<Factor, Factoree>::factor_exponent,
  typename impl::extract_factor_<Factor, Factoree>::rest_of_factoree
>
extract_factor(number<Factor>, number<Factoree>) { return impl::make_any_number(); }

#if 0
template<
  typename Integer,
  typename UInt,
  bool BigEndian, // the order of the UInt digits; each UInt is internally the native bit-order
  bool Signed, // represented as radix complement, typically two's complement
  bool ZeroHasADigit = true,
  UInt BaseMinusOne = std::numeric_limits<UInt>::max()>
struct convert_to_base_auto_sized {
  typedef typename impl::convert_to_base<Integer, UInt, BigEndian, Signed,
    true, ZeroHasADigit, false, 0, BaseMinusOne>::type type;
  static constexpr type value = type();
};

template<
  typename Integer,
  typename UInt,
  bool BigEndian, // the order of the UInt digits; each UInt is internally the native bit-order
  bool Signed, // represented as radix complement, typically two's complement
  uintmax_t ExactSize,
  bool ModuloIfOverflow = false /*non-modulo: error if overflow*/,
  UInt BaseMinusOne = std::numeric_limits<UInt>::max()>
struct convert_to_base_manually_sized {
  typedef typename impl::convert_to_base<Integer, UInt, BigEndian, Signed,
    false, false, ModuloIfOverflow, ExactSize, BaseMinusOne>::type type;
  static constexpr type value = type();
};
#endif

template<typename A> constexpr inline
number<typename impl::round_<A, rounding_strategy<round_down, negative_continuous_with_positive>>::type>
floor(number<A>) { return impl::make_any_number(); }

template<typename A> constexpr inline
number<typename impl::round_<A, rounding_strategy<round_up, negative_continuous_with_positive>>::type>
ceil(number<A>) { return impl::make_any_number(); }

template<typename A> constexpr inline
number<typename impl::round_<A, rounding_strategy<round_to_nearest_with_ties_rounding_to_even>>::type>
nearbyint(number<A>) { return impl::make_any_number(); }

template<typename A, typename RoundingStrategy> constexpr inline
number<typename impl::round_<A, RoundingStrategy>::type>
round(number<A>, RoundingStrategy) { return impl::make_any_number(); }

template<typename Base, typename Number, typename RoundingStrategy> constexpr inline
number<typename impl::logarithm<Base, Number, RoundingStrategy>::type>
log(number<Base>, number<Number>, RoundingStrategy) { return impl::make_any_number(); }

template<typename Number, typename RoundingStrategy> constexpr inline
number<typename impl::logarithm<nat<2>, Number, RoundingStrategy>::type>
log2(number<Number>, RoundingStrategy) { return impl::make_any_number(); }

template<typename Number, typename RoundingStrategy> constexpr inline
number<typename impl::logarithm<nat<10>, Number, RoundingStrategy>::type>
log10(number<Number>, RoundingStrategy) { return impl::make_any_number(); }

//typedef nat<> zero;

// and operator bool, maybe explicit operator float, etc.
// I suppose the user-visible types should *anyway* be different from the intermediates...

} /* end namespace ct */

#if LASERCAKE_HAVE_USER_DEFINED_LITERALS

template<char...Digits>
constexpr inline ct::number<typename ct::impl::parse_nat<Digits...>::type>
operator "" _integer() {
  return ct::impl::make_any_number();
}

#include <boost/preprocessor/cat.hpp>
#define INTEGER(n) (BOOST_PP_CAT(n,_integer))

#else

// This might not be the fastest possible; our GCC 4.6 support is grudging.
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/stringize.hpp>
#define CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS_AUX(z, n, str) \
  (sizeof(str) > n ? str[n] : '\0')
#define CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS(ident) \
  BOOST_PP_ENUM( \
    BOOST_PP_LIMIT_REPEAT, \
    CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS_AUX, \
    BOOST_PP_STRINGIZE(ident))
#define INTEGER(n) (::ct::impl::_integer<sizeof(BOOST_PP_STRINGIZE(n)), \
                       CONVERT_IDENTIFIER_TO_TEMPLATE_CHAR_ARGUMENTS(n)>())

namespace ct {
namespace impl {
template<size_t NumDigits, char...Digits>
constexpr inline number<typename impl::parse_integer<Digits...>::type>
_integer() {
#if !LASERCAKE_HAVE_USER_DEFINED_LITERALS
  static_assert(NumDigits <= BOOST_PP_LIMIT_REPEAT,
    "ct:: compile-time numeric literal has more digits than supported by this environment.");
#endif
  return impl::make_any_number();
}
}} /* end namespaces ct::impl, ct */

#endif

#define RATIONAL(n, d) (INTEGER(n) / INTEGER(d))
#define INTEGERtype(n) decltype(INTEGER(n))
#define RATIONALtype(n, d) decltype(RATIONAL(n, d))

/*
 * ideas
// NAT could be INT: for operator"" _NAT it is trivial, and for the workaround it requires
// detecting leading '-' in the template which is easy.
// INT(i) and RAT(num, den) ?
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
*/




#endif
