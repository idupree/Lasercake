/*

    Copyright Eli Dupree and Isaac Dupree, 2014

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

#ifndef LASERCAKE_BOUNDED_INT_CALCULUS_HPP__
#define LASERCAKE_BOUNDED_INT_CALCULUS_HPP__


#include <memory>
#include <array>

#include <boost/type_traits.hpp>
#include <boost/range/iterator_range.hpp>

#include "../data_structures/numbers.hpp"

// solely for get_primitive definition
#include "../data_structures/bounds_checked_int.hpp"

// for tests
#include <iostream>

namespace bounded_int_calculus {

// TODO use a more efficient formula
template<typename IntType>
IntType binomial_coefficient(IntType n, IntType k) {
  if ((n == k) || (k == 0)) return 1;
  k = std::min(k, n - k);
  IntType result = 1;
  for (IntType i = 0; i < k; ++i) {
    result = result * (n - i) / (i + 1); // math says this never needs rounding
  }
  return result;
}

// arrayvector<std::vector<float>, 3> is a sequence container
// that stores 3 floats directly in the arrayvector object
// and any floats after the third in a conditionally-present
// std::vector<float> object.  If size() is always <= ArraySize
// then arrayvector never allocates.
//
// value_type must be DefaultConstructible and Assignable (this is
// a conceptual bug).
template<typename VectorType, size_t ArraySize> class arrayvector {
public:
  typedef typename VectorType::value_type value_type;
private:
  std::array<value_type, ArraySize> first_n_values_;
  size_t size_;
  std::unique_ptr<VectorType> real_vector_;
public:
  arrayvector():size_(0),real_vector_(nullptr){}
  size_t size()const { return size_; }
  inline value_type& operator[](size_t idx) {
    if (idx < ArraySize) {
      return first_n_values_[idx];
    }
    else {
      caller_correct_if(real_vector_.get(), "out of bounds arrayvector access");
      return (*real_vector_)[idx - ArraySize];
    }
  }
  inline value_type& back() {
    return (*this)[size_ - 1];
  }
  inline value_type& front() {
    return first_n_values_[0];
  }
  void push_back(value_type value) {
    if (size_ < ArraySize) {
      first_n_values_[size_] = value;
    }
    else {
      if (!real_vector_) {
        real_vector_ = std::unique_ptr<VectorType>(new VectorType());
        //LOG << "arrayvector size exceeded: " << ArraySize << "\n";
      }
      real_vector_->push_back(value);
    }
    ++size_;
  }
  inline void pop_back() {
    if(size_ > ArraySize) {
      real_vector_->pop_back();
    }
    --size_;
  }
  inline bool empty()const { return size_ == 0; }
  class iterator : public boost::iterator_facade<iterator, value_type, boost::bidirectional_traversal_tag> {
    public:
      iterator() : av_(NULL),which_(0) {}
      explicit iterator(arrayvector* av_, size_t which_) : av_(av_),which_(which_) {}
    private:
    friend class boost::iterator_core_access;
      void increment() {++which_;}
      void decrement() {--which_;}
      bool equal(iterator const& other)const { return which_ == other.which_; }
      value_type& dereference()const { return (*av_)[which_]; }
      arrayvector* av_;
      size_t which_;
  };
  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size_); }
};

typedef uint32_t num_terms_t;

template<typename Poly, num_terms_t idx_so_far, typename...Args>
struct polynomial_construct_helper;
template<typename Poly, num_terms_t idx_so_far, typename Head, typename...Tail>
struct polynomial_construct_helper<Poly, idx_so_far, Head, Tail...> {
  static inline void func(Poly& poly, Head head, Tail... tail) {
    poly.terms[idx_so_far] = head;
    polynomial_construct_helper<Poly, idx_so_far+1, Tail...>::func(poly, std::forward<Tail>(tail)...);
  }
};
template<typename Poly, num_terms_t idx_so_far>
struct polynomial_construct_helper<Poly, idx_so_far> {
  static inline void func(Poly& poly) {
    static_assert(idx_so_far <= Poly::max_terms, "Constructing a polynomial with more terms than its max");
    poly.num_terms = idx_so_far;
    poly.prune_leading_zeroes();
    poly.check_coefficients();
  }
};

template<typename DomainType, typename ValueType, num_terms_t MaxTerms>
class polynomial {
public:
  template<typename OtherDomainType, typename OtherValueType, num_terms_t OtherMaxTerms>
  friend class polynomial;
  static_assert(MaxTerms > 0, "We don't support 0-term polynomials");
  typedef ValueType value_type;
  typedef DomainType domain_type;
  static const num_terms_t max_terms = MaxTerms;
  static const value_type lots = std::numeric_limits<value_type>::max();
  static const value_type neglots = std::numeric_limits<value_type>::min();
  static const value_type max_coefficient = lots-1;
  static const value_type min_coefficient = neglots+1;
  typedef decltype(lossless_multiply(value_type(), value_type())) bigger_int_t;
private:
  static inline void check_coefficient(value_type c) {
    caller_correct_if(c >= min_coefficient && c <= max_coefficient, "polynomial coefficient out of bounds");
  }
  inline void check_coefficients() {
    for (num_terms_t i = 0; i < num_terms; ++i) {
      check_coefficient(terms[i]);
    }
  }
  static inline bool add_overflows(value_type c1, value_type c2) {
    return (c1>>1)+(c2>>1)+(c1&c2&  max_coefficient &1)+((c1|c2)&(~max_coefficient)&1) > (max_coefficient>>1);
  }
  static inline bool add_underflows(value_type c1, value_type c2) {
    return (c1>>1)+(c2>>1)+(c1&c2&(~min_coefficient)&1)+((c1|c2)&  min_coefficient &1) < (min_coefficient>>1)+(min_coefficient&1);
  }
  static inline value_type checked_add(value_type c1, value_type c2) {
    caller_error_if(add_overflows (c1, c2),  "overflow in polynomial operation");
    caller_error_if(add_underflows(c1, c2), "underflow in polynomial operation");
    return c1 + c2;
  }
  // Add some numbers without allowing the final result to be out-of-bounds,
  // but tolerating (above_max)+(below_min) if the result is in bounds.
  template<typename ArrayType>
  static value_type lenient_sum(ArrayType& pos_summands, ArrayType& neg_summands, num_terms_t num_pos_summands, num_terms_t num_neg_summands) {
    typedef typename ArrayType::value_type sum_type;
    // TODO reduce duplicate code (id z0IwiMnmm4DpFA)
    sum_type accumulated = 0;
    while (true) {
      if ((accumulated <= 0) && (num_pos_summands > 0)) { accumulated += pos_summands[--num_pos_summands]; }
      else if ((accumulated >= 0) && (num_neg_summands > 0)) { accumulated += neg_summands[--num_neg_summands]; }
      else {
        for (num_terms_t i = 0; i < num_pos_summands; ++i) {
          sum_type const& summand = pos_summands[i];
          caller_error_if(summand > max_coefficient, "overflow in polynomial operation");
          accumulated += summand;
          caller_error_if(accumulated > max_coefficient, "overflow in polynomial operation");
        }
        for (num_terms_t i = 0; i < num_neg_summands; ++i) {
          sum_type const& summand = neg_summands[i];
          caller_error_if(summand < min_coefficient, "underflow in polynomial operation");
          accumulated += summand;
          caller_error_if(accumulated < min_coefficient, "underflow in polynomial operation");
        }
        caller_error_if(accumulated > max_coefficient, "overflow in polynomial operation");
        caller_error_if(accumulated < min_coefficient, "underflow in polynomial operation");
        return value_type(accumulated);
      }
    }
  }
  
public:
  polynomial():num_terms(0){}
  void operator=(value_type c) { terms[0] = c; num_terms = (c != 0); }
  // Note: we start with the constant term, not the "leading" term
  polynomial(value_type c):num_terms(c != 0){ terms[0] = c; }
  
  template<typename Poly, num_terms_t idx_so_far, typename...Args>
  friend struct polynomial_construct_helper;
  template<typename T1, typename T2, typename...T>
  polynomial(T1 data1, T2 data2, T... data){ polynomial_construct_helper<polynomial, 0, T1, T2, T...>::func(*this, std::forward<T1>(data1), std::forward<T2>(data2), std::forward<T>(data)...); }
  
  value_type get_term(uint32_t which)const {
    if (which < num_terms) return terms[which];
    return 0;
  }
  void set_term(uint32_t which, value_type new_value) {
    check_coefficient(new_value);
    if (new_value == 0) {
      if (which+1 < num_terms) {
        terms[which] = 0;
      }
      if (which+1 == num_terms) {
        --num_terms;
        prune_leading_zeroes();
      }
    }
    else {
      if (which < num_terms){
        terms[which] = new_value;
      }
      else {
        resize(which);
        push_back(new_value);
      }
    }
  }
  
  polynomial derivative()const {
    polynomial result;
    if (num_terms > 1) {
      for (num_terms_t i = 1; i < num_terms; ++i) {
        result.push_back(terms[i]*i);
      }
    }
    return result;
  }
  
  template<num_terms_t OtherMaxTerms>
  std::enable_if_t<OtherMaxTerms <= max_terms> operator+=(polynomial<domain_type, value_type, OtherMaxTerms> const& other) {
    if (num_terms < other.num_terms) { resize(other.num_terms); }
    for (num_terms_t i = 0; i < other.num_terms; ++i) {
      terms[i] = checked_add(terms[i], other.terms[i]);
    }
    prune_leading_zeroes();
  }
  template<num_terms_t OtherMaxTerms>
  std::enable_if_t<OtherMaxTerms <= max_terms, polynomial> operator+(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    polynomial result(*this); result += other; return std::move(result);
  }
  template<num_terms_t OtherMaxTerms>
  std::enable_if_t<(OtherMaxTerms > max_terms), polynomial<domain_type, value_type, OtherMaxTerms>> operator+(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    polynomial<domain_type, value_type, OtherMaxTerms> result(other); result += *this; return std::move(result);
  }
  void operator+=(value_type const& other) {
    if (num_terms == 0) { check_coefficient(other); num_terms = 1; terms[0] = other; }
    else {
      terms[0] = checked_add(terms[0], other);
      if (num_terms == 1 && terms[0] == 0) num_terms = 0;
    }
  }
  polynomial operator+(value_type const& other)const {
    polynomial result(*this); result += other; return std::move(result);
  }
  template<num_terms_t OtherMaxTerms>
  std::enable_if_t<OtherMaxTerms <= max_terms> operator-=(polynomial<domain_type, value_type, OtherMaxTerms> const& other) {
    if (num_terms < other.num_terms) { resize(other.num_terms); }
    for (num_terms_t i = 0; i < other.num_terms; ++i) {
      terms[i] = checked_add(terms[i], -other.terms[i]);
    }
    prune_leading_zeroes();
  }
  template<num_terms_t OtherMaxTerms>
  std::enable_if_t<OtherMaxTerms <= max_terms, polynomial> operator-(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    polynomial result(*this); result -= other; return std::move(result);
  }
  template<num_terms_t OtherMaxTerms>
  std::enable_if_t<(OtherMaxTerms > max_terms), polynomial<domain_type, value_type, OtherMaxTerms>> operator-(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    polynomial<domain_type, value_type, OtherMaxTerms> result(other); result -= *this; return std::move(result);
  }
  void operator-=(value_type const& other) {
    if (num_terms == 0) { check_coefficient(-other); num_terms = 1; terms[0] = -other; }
    else {
      terms[0] = checked_add(terms[0], -other);
      if (num_terms == 1 && terms[0] == 0) num_terms = 0;
    }
  }
  polynomial operator-(value_type const& other)const {
    polynomial result(*this); result -= other; return std::move(result);
  }
  template<num_terms_t OtherMaxTerms>
  polynomial<domain_type, value_type, max_terms + OtherMaxTerms - 1> operator*(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    polynomial<domain_type, value_type, max_terms + OtherMaxTerms - 1> result;
    if (num_terms > 0 && other.num_terms > 0) {
      const int result_terms = num_terms+other.num_terms-1;
      result.num_terms = result_terms;
      for (int i = 0; i < result_terms; ++i) {
        // to determine result.terms[i]
        const int first = std::max(0, int(i-(other.num_terms-1)));
        const int  last = std::min(i, int(num_terms-1));
        // TODO reduce duplicate code (id z0IwiMnmm4DpFA)
        std::array<bigger_int_t, (max_terms > OtherMaxTerms ? max_terms : OtherMaxTerms)> pos_summands;
        std::array<bigger_int_t, (max_terms > OtherMaxTerms ? max_terms : OtherMaxTerms)> neg_summands;
        num_terms_t num_pos_summands = 0;
        num_terms_t num_neg_summands = 0;
        for (int j = first; j <= last; ++j) {
          const int k = i-j;
          const bigger_int_t summand = lossless_multiply(terms[j], other.terms[k]);
          ((summand<0) ? neg_summands[num_neg_summands++] : pos_summands[num_pos_summands++]) = summand;
        }
        result.terms[i] = lenient_sum(pos_summands, neg_summands, num_pos_summands, num_neg_summands);
      }
    }
    return std::move(result);
  }
//   void operator*=(polynomial const& other) {
//     *this = *this * other;
//   }
  void operator*=(value_type const& other) {
    if (other == 0) { num_terms = 0; }
    else {
      for (num_terms_t i = 0; i < num_terms; ++i) {
        const auto product = lossless_multiply(terms[i], other);
        caller_error_if(product > max_coefficient, "overflow in polynomial operation");
        caller_error_if(product < min_coefficient, "underflow in polynomial operation");
        terms[i] = value_type(product);
      }
    }
  }
  polynomial operator*(value_type const& other)const {
    polynomial result(*this); result *= other; return std::move(result);
  }
  
  polynomial operator-()const {
    polynomial result(*this);
    for (num_terms_t i = 0; i < result.num_terms; ++i) {
      result.terms[i] = -result.terms[i];
      check_coefficient(result.terms[i]);
    }
    return std::move(result);
  }
  template<num_terms_t OtherMaxTerms>
  bool operator==(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    if (num_terms != other.num_terms) return false;
    for (num_terms_t i = 0; i < num_terms; ++i) {
      if (terms[i] != other.terms[i]) return false;
    }
    return true;
  }
  template<num_terms_t OtherMaxTerms>
  bool operator!=(polynomial<domain_type, value_type, OtherMaxTerms> const& other)const {
    return !((*this) == other);
  }
  
  value_type taylor_coefficient(domain_type const& origin, uint32_t which_term)const {
    // Taylor series entries are (nth derivative at a)/(n!) (x-a)^n
    // (nth derivative at a) = \sum_i=n^{num_terms-1} terms[i]*(i!/(i-n)!) a^(i-n)
    // so, taylor series coefficient n = \sum_i=n^{num_terms-1} terms[i]*(i choose n) a^(i-n)
    if (which_term >= num_terms) return 0;
    domain_type factor = 1;
    // TODO reduce duplicate code (id z0IwiMnmm4DpFA)
    std::array<bigger_int_t, max_terms> pos_summands;
    std::array<bigger_int_t, max_terms> neg_summands;
    num_terms_t num_pos_summands = 0;
    num_terms_t num_neg_summands = 0;
    ((terms[which_term]<0) ? neg_summands[num_neg_summands++] : pos_summands[num_pos_summands++]) = terms[which_term];
    for (uint32_t i = which_term+1; i < num_terms; ++i) {
      factor *= origin;
      const bigger_int_t summand = lossless_multiply(terms[i], binomial_coefficient(i, which_term) * factor);
      ((summand<0) ? neg_summands[num_neg_summands++] : pos_summands[num_pos_summands++]) = summand;
    }
    return lenient_sum(pos_summands, neg_summands, num_pos_summands, num_neg_summands);
  }
  
  void move_origin(domain_type const& new_origin) {
    if (num_terms > 1) {
      for (uint32_t i = 0; i+1 < num_terms; ++i) {
        terms[i] = taylor_coefficient(new_origin, i);
      }
    }
  }
  
  polynomial with_origin(domain_type const& new_origin)const {
    polynomial result(*this); result.move_origin(new_origin); return std::move(result);
  }
  
  // Essentially returns the exact result clamped to the range of the type.
  // Notably, its sign is always correct.
  value_type operator()(domain_type const& input)const {
    // If abs(input) >= 2, then the loop at the bottom of this function
    // would stay out-of-bounds as soon as it was out-of-bounds once.
    // That's why we can return early from the loop.
    // With ones, that doesn't work, so we need a special case.
    if ((input == 1) or (input == -1)) {
      // TODO reduce duplicate code (id z0IwiMnmm4DpFA)
      std::array<value_type, max_terms> pos_summands;
      std::array<value_type, max_terms> neg_summands;
      num_terms_t num_pos_summands = 0;
      num_terms_t num_neg_summands = 0;
      for (num_terms_t i = 0; i < num_terms; ++i) {
        const value_type summand = ((input<0) && (i&1)) ? -terms[i] : terms[i];
        ((summand<0) ? neg_summands[num_neg_summands++] : pos_summands[num_pos_summands++]) = summand;
      }
      value_type result = 0;
      while (true) {
        if ((result <= 0) && (num_pos_summands > 0)) { result += pos_summands[--num_pos_summands]; }
        else if ((result >= 0) && (num_neg_summands > 0)) { result += neg_summands[--num_neg_summands]; }
        else {
          for (num_terms_t i = 0; i < num_pos_summands; ++i) {
            if (add_overflows(result, pos_summands[i])) { return lots; }
            result += pos_summands[i];
          }
          for (num_terms_t i = 0; i < num_neg_summands; ++i) {
            if (add_underflows(result, neg_summands[i])) { return lots; }
            result += neg_summands[i];
          }
          return result;
        }
      }
    }
    value_type result = 0;
    for (int i = num_terms-1; i >= 0; --i) {
      auto i_and_above_part = lossless_multiply(result, input) + terms[i];
      // Suppose i_and_above_part is out of bounds. WLOG assume, i_and_above_part > 0, input > 0 (so input > 1).
      // Then, in the next cycle, 
      //   (result * input) + terms[i]
      //   >= ((max+1) * input) + terms[i]
      //   >= ((max+1) * input) - max
      //   >= ((max+1) * 2) - max
      //   == max + 2
      //   > max
      // so it's still out of bounds in the same direction.
      // The +2 gives leeway for max and -min to differ by 1.
      if (i_and_above_part > max_coefficient) { return ((input<0) && (i&1)) ? neglots : lots; }
      if (i_and_above_part < min_coefficient) { return ((input<0) && (i&1)) ? lots : neglots; }
      result = value_type(i_and_above_part);
    }
    return result;
  }
  
  value_type neginf_sign()const { if (num_terms == 0) return 0; else return sign(back() * ((num_terms & 1) ? 1 : -1)); }
  value_type posinf_sign()const { if (num_terms == 0) return 0; else return sign(back()); }
  
  class sign_interval_boundary_iterator {
  private:
    domain_type input_;
    polynomial const* p_;
    
    std::array<domain_type,6> future;
    uint32_t future_size;
    
    // for order 3+ polys
    polynomial derivative_;
    std::unique_ptr<sign_interval_boundary_iterator> derivative_iterator_;
    bool advanced_to_transition_before_current_derivative_iterator_position_;
    
    void advance_to(domain_type input) {
      input_ = input;
    }
    bool maybe_advance_to(domain_type input) {
      if (input_ < input) {
        advance_to(input);
        return true;
      }
      return false;
    }
  public:
    bool operator==(sign_interval_boundary_iterator const& other)const {
      return (p_ == other.p_) && (!p_ || (input_ == other.input_));
    }
    bool operator!=(sign_interval_boundary_iterator const& other)const {
      return (p_ != other.p_) || ( p_ && (input_ != other.input_));
    }
    domain_type operator*()const {
      caller_correct_if(p_, "can't dereference an empty iterator");
      return input_;
    }
    
    sign_interval_boundary_iterator():p_(nullptr){}
    sign_interval_boundary_iterator(polynomial const* p, domain_type start):p_(p),future_size(0){
      assert(p_);
      advance_to(start);
      
      if (p_->num_terms <= 1) {}
      else if (p_->num_terms == 2) {
        const domain_type at_or_below_the_zero = divide(-p_->terms[0], p_->terms[1], rounding_strategy<round_down, negative_continuous_with_positive>());
        future[0] = at_or_below_the_zero+1;
        future[1] = at_or_below_the_zero;
        future_size = 2;
        if (at_or_below_the_zero * p_->terms[1] == -p_->terms[0]) {
          // We hit the zero exactly, which makes this more complicated.
          future[2] = at_or_below_the_zero-1;
          future_size = 3;
        }
      }
      else if (p_->num_terms == 3) {
        // Quadratics are getting more complicated, but we've got the quadratic formula.
        value_type const& a = p_->terms[2];
        value_type const& b = p_->terms[1];
        value_type const& c = p_->terms[0];
        // "a" probably won't overflow with "4*" because if it did, most normal uses
        // of the quadratic would also overflow.
        const bigger_int_t discriminant_ = lossless_multiply(b,b) - lossless_multiply(4*a,c);
        const value_type discriminant = value_type(discriminant_);
        if (discriminant_ >= 0) {
          if (discriminant_ != discriminant) {
            // Overflow/underflow - fall back to the derivative system
            // TODO clean this stuff up (e.g. value_type(discriminant_) could overflow all the way around to the 'correct' value)
            derivative_ = p_->derivative();
            derivative_iterator_ = std::unique_ptr<sign_interval_boundary_iterator>(new sign_interval_boundary_iterator(&derivative_, start));
          }
          else {
            const value_type sqrt_disc_rounded_down = isqrt(discriminant);
            const domain_type num1 = -b - sqrt_disc_rounded_down*sign(a);
            const domain_type num2 = -b + sqrt_disc_rounded_down*sign(a);
            const domain_type denom = a*2;
            const domain_type at_or_above_the_first_zero =  divide(num1, denom, rounding_strategy<round_up  , negative_continuous_with_positive>());
            const domain_type at_or_below_the_second_zero = divide(num2, denom, rounding_strategy<round_down, negative_continuous_with_positive>());
            // Some math facts allow us to simplify: The sqrt of an integer is always an integer or irrational,
            // and an integer polynomial never has exactly one integer root.
            if ((sqrt_disc_rounded_down*sqrt_disc_rounded_down == discriminant) && (at_or_above_the_first_zero * denom == num1)) {
              future[0] = at_or_below_the_second_zero + 1;
              future[1] = at_or_below_the_second_zero;
              future[2] = at_or_below_the_second_zero - 1;
              future[3] = at_or_above_the_first_zero + 1;
              future[4] = at_or_above_the_first_zero;
              future[5] = at_or_above_the_first_zero - 1;
              future_size = 6;
            }
            else {
              future[0] = at_or_below_the_second_zero + 1;
              future[1] = at_or_below_the_second_zero;
              future[2] = at_or_above_the_first_zero;
              future[3] = at_or_above_the_first_zero - 1;
              future_size = 4;
            }
          }
        }
      }
      else {
        assert(p_->num_terms > 3);
        derivative_ = p_->derivative();
        derivative_iterator_ = std::unique_ptr<sign_interval_boundary_iterator>(new sign_interval_boundary_iterator(&derivative_, start));
      }
    }
    sign_interval_boundary_iterator& operator++() {
      caller_correct_if(p_, "can't increment an empty iterator");
      while (future_size > 0) {
        if (maybe_advance_to(future[--future_size])) { return *this; }
      }
      if (!derivative_iterator_) {
        p_ = nullptr;
      }
      else {
        if (derivative_iterator_->p_ && (input_ == derivative_iterator_->input_)) {
          ++(*derivative_iterator_);
          advanced_to_transition_before_current_derivative_iterator_position_ = false;
        }
        
        if (advanced_to_transition_before_current_derivative_iterator_position_) {
          if (derivative_iterator_->p_) {
            advance_to(derivative_iterator_->input_);
          }
          else {
            p_ = nullptr;
          }
        }
        else {
          domain_type min = input_;
          domain_type max = derivative_iterator_->input_;
          if (!derivative_iterator_->p_) {
            // TODO improve (overflow? efficiency?)
            max = (input_ > 0) ? (input_+1) : 1;
            const value_type ps = p_->posinf_sign();
            while(sign((*p_)(max)) != ps) { max <<= 1; }
          }
          assert(max > min);
          const value_type sign_min = sign((*p_)(min));
          const value_type sign_max = sign((*p_)(max));
          if (sign_min == sign_max) {
            if (derivative_iterator_->p_) {
              advance_to(derivative_iterator_->input_);
            }
            else {
              p_ = nullptr;
            }
          }
          else {
            advanced_to_transition_before_current_derivative_iterator_position_ = true;
                 if (sign_min == 0) { advance_to(min + 1); }
            else if (sign_max == 0) { advance_to(max - 1); }
            else {
              while(min + 1 < max) {
                const domain_type mid = (min>>1)+(max>>1)+(min&max&1);
                const value_type sign_mid = sign((*p_)(mid));
                if (sign_mid == 0) {
                  future[future_size++] = mid + 1;
                  future[future_size++] = mid;
                  future[future_size++] = mid - 1;
                  return ++(*this);
                }
                else if      (sign_mid == sign_min){ min = mid; }
                else { assert(sign_mid == sign_max); max = mid; }
              }
              future[future_size++] = max;
              future[future_size++] = min;
              return ++(*this);
            }
          }
        }
      }
      return *this;
    }
  };
  
  inline sign_interval_boundary_iterator sign_interval_boundaries_upper_bound(domain_type d)const {
    sign_interval_boundary_iterator result(this, d);
    ++result;
    return result;
  }
  inline sign_interval_boundary_iterator sign_interval_boundaries_end()const {
    return sign_interval_boundary_iterator();
  }
  
  friend inline std::ostream& operator<<(std::ostream& os, polynomial const& p) {
    for (num_terms_t i = 0; i < p.num_terms; ++i) {
      os << p.terms[i] << "x^" << i;
      if ((i+1) != p.num_terms) os << " + ";
    }
    return os;
  }
private:
  num_terms_t num_terms;
  std::array<value_type, max_terms> terms;
  void prune_leading_zeroes() {
    while((num_terms > 0) && (back() == 0)) { --num_terms; }
  }
  void resize(num_terms_t new_size) {
    // TODO: can memset instead?
    for (num_terms_t i = num_terms; i < new_size; ++i) {
      terms[i] = 0;
    }
    num_terms = new_size;
  }
  void push_back(value_type v) {
    terms[num_terms++] = v;
  }
  value_type back()const {
    return terms[num_terms-1];
  }
};

template<typename DomainType, typename ValueType, num_terms_t MaxTerms>
class polynomial_with_origin {
public:
  template<typename OtherDomainType, typename OtherValueType, num_terms_t OtherMaxTerms>
  friend class polynomial_with_origin;
  typedef DomainType domain_type;
  typedef ValueType value_type;
  static const num_terms_t max_terms = MaxTerms;
  typedef polynomial<domain_type, value_type, max_terms> without_origin_t;
  polynomial_with_origin(domain_type const& origin, without_origin_t const& p)
    :origin(origin),p(p){}
  // "0" problem: 0 + foo != foo (because it has a different origin). TODO can we do better?
  // (technically they're equal, but the representations being different is weird)
  constexpr polynomial_with_origin(decltype(nullptr)):origin(0),p(){}
  value_type get_term(domain_type const& where, uint32_t which_term)const {
    if (where == origin) return p.get_term(which_term);
    return p.taylor_coefficient(where - origin, which_term);
  }
  void set_term(domain_type const& where, uint32_t which_term, value_type new_value) {
    set_origin(where);
    p.set_term(which_term, new_value);
  }
  value_type operator()(domain_type const& where)const { return p(where - origin); }
  domain_type get_origin()const {
    return origin;
  }
  void set_origin(domain_type const& new_origin) {
    if (new_origin != origin) {
      p.move_origin(new_origin - origin);
      origin = new_origin;
    }
  }
  
  polynomial_with_origin operator-()const { return polynomial_with_origin(origin, -p); }
  
  template<num_terms_t OtherMaxTerms>
  void operator+=(polynomial_with_origin<domain_type, value_type, OtherMaxTerms> const& other) { p += other.p.with_origin(origin - other.origin); }
  template<num_terms_t OtherMaxTerms>
  void operator-=(polynomial_with_origin<domain_type, value_type, OtherMaxTerms> const& other) { p -= other.p.with_origin(origin - other.origin); }
  template<num_terms_t OtherMaxTerms>
  void operator*=(polynomial_with_origin<domain_type, value_type, OtherMaxTerms> const& other) { p *= other.p.with_origin(origin - other.origin); }
  void operator+=(value_type const& other) { p += other; }
  void operator-=(value_type const& other) { p -= other; }
  void operator*=(value_type const& other) { p *= other; }
private:
  
  typedef polynomial_with_origin<domain_type, value_type, max_terms> pwo;
public:
  template<num_terms_t OtherMaxTerms>
  polynomial_with_origin<domain_type, value_type, decltype(without_origin_t() + polynomial<domain_type, value_type, OtherMaxTerms>())::max_terms> operator+(
    polynomial_with_origin<domain_type, value_type, OtherMaxTerms> const& other)const { typedef
    polynomial_with_origin<domain_type, value_type, decltype(p + other.p)::max_terms> ret_t;
    if (origin == other.origin)          { return ret_t(      origin,       p + other.p); }
    else if (max_terms >= OtherMaxTerms) { return ret_t(      origin,       p + other.p.with_origin(origin - other.origin)); }
    else                                 { return ret_t(other.origin, other.p +       p.with_origin(other.origin - origin)); }
  }
  template<num_terms_t OtherMaxTerms>
  polynomial_with_origin<domain_type, value_type, decltype(without_origin_t() - polynomial<domain_type, value_type, OtherMaxTerms>())::max_terms> operator-(
    polynomial_with_origin<domain_type, value_type, OtherMaxTerms> const& other)const { typedef
    polynomial_with_origin<domain_type, value_type, decltype(p - other.p)::max_terms> ret_t;
    if (origin == other.origin)          { return ret_t(      origin,        p - other.p); }
    else if (max_terms >= OtherMaxTerms) { return ret_t(      origin,        p - other.p.with_origin(origin - other.origin)); }
    else                                 { return ret_t(other.origin, -other.p +       p.with_origin(other.origin - origin)); } // TODO can optimize?
  }
  template<num_terms_t OtherMaxTerms>
  polynomial_with_origin<domain_type, value_type, decltype(without_origin_t() * polynomial<domain_type, value_type, OtherMaxTerms>())::max_terms> operator*(
    polynomial_with_origin<domain_type, value_type, OtherMaxTerms> const& other)const { typedef
    polynomial_with_origin<domain_type, value_type, decltype(p * other.p)::max_terms> ret_t;
    if (origin == other.origin)          { return ret_t(      origin,       p * other.p); }
    else if (max_terms >= OtherMaxTerms) { return ret_t(      origin,       p * other.p.with_origin(origin - other.origin)); }
    else                                 { return ret_t(other.origin, other.p *       p.with_origin(other.origin - origin)); }
  }
  pwo operator+(value_type const& other)const { return pwo(origin, p + other); }
  pwo operator-(value_type const& other)const { return pwo(origin, p - other); }
  pwo operator*(value_type const& other)const { return pwo(origin, p * other); }
  
  class sign_interval_boundary_iterator {
  private:
    domain_type origin_;
    typename without_origin_t::sign_interval_boundary_iterator backend_;
  public:
    sign_interval_boundary_iterator(){}
    sign_interval_boundary_iterator(polynomial_with_origin const* p, domain_type start):origin_(p->origin),backend_(&p->p, start - p->origin){}
    sign_interval_boundary_iterator& operator++() { ++backend_; return *this; }
    domain_type operator*()const { return *backend_ + origin_; }
    bool operator==(sign_interval_boundary_iterator const& other) {
      return backend_ == other.backend_;
    }
    bool operator!=(sign_interval_boundary_iterator const& other) {
      return backend_ != other.backend_;
    }
  };
  inline sign_interval_boundary_iterator sign_interval_boundaries_upper_bound(domain_type d)const {
    sign_interval_boundary_iterator result(this, d);
    ++result;
    return result;
  }
  inline sign_interval_boundary_iterator sign_interval_boundaries_end()const {
    return sign_interval_boundary_iterator();
  }
  
  friend inline std::ostream& operator<<(std::ostream& os, polynomial_with_origin const& pwo) {
    return os << "[about " << pwo.origin << ", " << pwo.p << "]";
  }
private:
  domain_type origin;
  without_origin_t p;
};


template<int num_dimensions, typename CoordinateType>
class finite_dimensional_vector {
public:
  finite_dimensional_vector():data(){
    // TODO this is stupid, what is actually guaranteed? how do i do it?
    for(auto &i : data) { i = 0; }
  }
  //finite_dimensional_vector(finite_dimensional_vector const&) = default;
  //finite_dimensional_vector(finite_dimensional_vector&&) = default;

  // Hack: allow to initialize this type as
  // finite_dimensional_vector(dim1value, dim2value, dim3value...):
  template<typename T1>
  explicit finite_dimensional_vector(T1 data):data{ std::move(data) }{}
  template<typename T1, typename T2, typename ...T>
  finite_dimensional_vector(T1&& data1, T2&& data2, T&&...data):data{{ std::forward<T1>(data1), std::forward<T2>(data2), std::forward<T>(data)... }}{}
  // This is chosen in overload resolution in preference to the
  // copy-constructor so we can't use it:
  //template<typename ...T>
  //explicit finite_dimensional_vector(T&&...data):data{ std::forward<T>(data)... }{}
  
  template<typename OtherCoordinateType>
  decltype(std::declval<CoordinateType>()*std::declval<OtherCoordinateType>()) dot(finite_dimensional_vector<num_dimensions, OtherCoordinateType> const& other)const {
    // TODO: 0d vectors? deal with "0" problem with polynomial_with_origin
    decltype(std::declval<CoordinateType>()*std::declval<OtherCoordinateType>()) result(data[0] * other.data[0]);
    for (int i = 1; i < num_dimensions; ++i) result += data[i] * other.data[i];
    return result;
  }
  CoordinateType const& operator()(int which)const { return data[which]; }
  CoordinateType      & operator[](int which)      { return data[which]; }
  
  void operator+=(finite_dimensional_vector const& other) {
    for (int i = 0; i < num_dimensions; ++i) data[i] += other.data[i]; }
  void operator-=(finite_dimensional_vector const& other) {
    for (int i = 0; i < num_dimensions; ++i) data[i] -= other.data[i]; }
  
  // this implements the vector-space operation of scalar multiplication,
  // but is also more general, since the factor type need not be a field.
  // For instance, you may multiply a vector of polynomials by a polynomial,
  // even though polynomials under pointwise addition and multiplication
  // are not a field.
  template<typename Factor>
  void operator*=(Factor const& other) {
    for (int i = 0; i < num_dimensions; ++i) data[i] *= other; }
  
  template <int vector_num_dimensions, typename vector_coordinate_type, typename factor_type>
  friend auto operator*(finite_dimensional_vector <vector_num_dimensions, vector_coordinate_type> const & vector, factor_type const & factor) {
    finite_dimensional_vector <vector_num_dimensions, decltype (std:: declval <vector_coordinate_type> ()*factor) result;
    for (int dimension = 0; dimension <num_dimensions; ++ dimension) {
      result [dimension] = vector (dimension)*factor;
    }
    return result;
  }
  template <int vector_num_dimensions, typename vector_coordinate_type, typename factor_type>
  friend auto operator*(factor_type const & factor, finite_dimensional_vector <vector_num_dimensions, vector_coordinate_type> const & vector) {
    finite_dimensional_vector <vector_num_dimensions, decltype (factor*std:: declval <vector_coordinate_type> ()) result;
    for (int dimension = 0; dimension <num_dimensions; ++ dimension) {
      result [dimension] = factor*vector (dimension);
    }
    return result;
  }
  
  template <int shared_num_dimensions, typename coordinate_type_1, typename coordinate_type_2>
  friend auto operator+ (finite_dimensional_vector <shared_num_dimensions, coordinate_type_1> const & vector_1, finite_dimensional_vector <shared_num_dimensions, coordinate_type_2> const & vector_2) {
    finite_dimensional_vector <shared_num_dimensions, decltype (std::declval <coordinate_type_1> () + std:: declval <coordinate_type_2> ())> result;
    for (int dimension = 0; dimension <num_dimensions; ++ dimension) {
      result [dimension] = vector_1 (dimension) + vector_2 (dimension);
    }
    return result;
  }
  template <int shared_num_dimensions, typename coordinate_type_1, typename coordinate_type_2>
  friend auto operator- (finite_dimensional_vector <shared_num_dimensions, coordinate_type_1> const & vector_1, finite_dimensional_vector <shared_num_dimensions, coordinate_type_2> const & vector_2) {
    finite_dimensional_vector <shared_num_dimensions, decltype (std::declval <coordinate_type_1> () - std:: declval <coordinate_type_2> ())> result;
    for (int dimension = 0; dimension <num_dimensions; ++ dimension) {
      result [dimension] = vector_1 (dimension) - vector_2 (dimension);
    }
    return result;
  }
  
  bool operator==(finite_dimensional_vector const& other)const {
    for (int i = 0; i < num_dimensions; ++i) { if (data[i] != other.data[i]) { return false; }} return true; }
  bool operator!=(finite_dimensional_vector const& other)const {
    for (int i = 0; i < num_dimensions; ++i) { if (data[i] != other.data[i]) { return true; }} return false; }
   
  template<typename ValueType, typename DomainType>
  finite_dimensional_vector<num_dimensions, ValueType> get_term(DomainType const& where, uint32_t which_term)const {
    finite_dimensional_vector<num_dimensions, ValueType> result;
    for (int i = 0; i < num_dimensions; ++i) {
      result[i] = data[i].get_term(where, which_term);
    }
    return result;
  }
  template<typename ValueType, typename DomainType>
  void set_term(DomainType const& where, uint32_t which_term, finite_dimensional_vector<num_dimensions, ValueType> new_value) {
    for (int i = 0; i < num_dimensions; ++i) {
      data[i].set_term(where, which_term, new_value[i]);
    }
  }
  /*
  template<typename DomainType>
  finite_dimensional_vector<num_dimensions, typename CoordinateType::value_type> operator()(DomainType const& where)const {
    finite_dimensional_vector<num_dimensions, typename CoordinateType::value_type> result;
    for (int i = 0; i < num_dimensions; ++i) {
      result[i] = this->data[i]();
    }
    return result;
  }*/
  template<typename DomainType>
  void set_origin(DomainType const& new_origin) {
    for (auto& i : this->data) i.set_origin(new_origin);
  }
  
  friend inline std::ostream& operator<<(std::ostream& os, finite_dimensional_vector const& fd) {
    os << "(";
    for (size_t i = 0; i < fd.data.size(); ++i) {
      os << fd.data[i];
      if ((i+1) != fd.data.size()) os << ", ";
    }
    os << ")";
    return os;
  }
private:
  std::array<CoordinateType, num_dimensions> data;
};

void test_iterate(polynomial<int64_t,int64_t,5> p) {
  std::cerr << "iterating " << p << "\n";
  for (auto i = p.sign_interval_boundaries_upper_bound(std::numeric_limits<int64_t>::min()); i != p.sign_interval_boundaries_end(); ++i) {
    std::cerr << *i << "\n";
  }
}
void test() {
  typedef polynomial<int64_t,int64_t,5> p;
  const int64_t big = 1LL << 20;
  const int64_t huge = isqrt(p::max_coefficient)+2;
#define SHOULD_EXCEPT(expr) \
{ \
  bool excepted = false; \
  try { (expr); } \
  catch(std::logic_error) { excepted = true; } \
  assert(excepted); \
}
  std::cerr << std::hex << huge*huge << "\n";
  std::cerr << std::hex << p::max_coefficient << std::dec << "\n";
  assert(p(0)==p());
  assert(p(0,0)==p());
  assert(p(0,0,0)==p());
  assert(p(1)!=p());
  assert(-p(1,2,3)==p(-1,-2,-3));
  assert(p(p(1)+p(1))==p(2));
  assert(p(1,2,3)+p(1,2)==p(2,4,3));
  assert(p(1,2,3)*p(1,2)==p(1,4,7,6));
  std::cerr << p(big,huge,-big)*p(big,huge,-big) << "\n";
  SHOULD_EXCEPT(p(0,huge,1)*p(0,huge,1))
  SHOULD_EXCEPT(p(p::max_coefficient)+p(1))
  SHOULD_EXCEPT(p(p::min_coefficient)-p(1))
  assert(p()(0)==0);
  assert(p()(1)==0);
  assert(p()(-1)==0);
  assert(p()(2)==0);
  assert(p()(-2)==0);
  assert(p(huge,huge)(huge)==p::lots);
  assert(p(-huge,huge)(-huge)==p::neglots);
  assert(p(0,huge,-1)(huge)==0);
  assert(p().derivative()==p());
  assert(p(-5).derivative()==p());
  assert(p(-5,-5).derivative()==p(-5));
  assert(p(-5,-5,-5).derivative()==p(-5,-10));
  assert(p(-5,-5,-5,-5).derivative()==p(-5,-10,-15));
  assert(p().with_origin(5)==p());
  assert(p(1).with_origin(5)==p(1));
  assert(p(1,1).with_origin(5)==p(1+5,1));
  assert(p(1,1,1).with_origin(5)==p(1+5+5*5,1+5*2,1));
  test_iterate(p(10,0,-1));
  test_iterate(p(10000,0,-100,1));
  test_iterate(p(10000,0,-100,1));
  test_iterate(p(p::max_coefficient,-2));
  test_iterate(p(p::min_coefficient,-2));
  test_iterate(p(p::max_coefficient,-2,-2));
  test_iterate(p(p::min_coefficient,-2,-2));
}

} // namespace bounded_int_calculus
namespace std {
  template<int num_dimensions, typename CoordinateType>
  struct hash<bounded_int_calculus::finite_dimensional_vector<num_dimensions, CoordinateType>> {
    public:
    size_t operator()(bounded_int_calculus::finite_dimensional_vector<num_dimensions, CoordinateType> const& i)const {
      size_t seed = std::hash<CoordinateType>()(i(0));
      for (int j = 1; j < num_dimensions; ++j) {
        seed ^= std::hash<CoordinateType>()(i(j)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };
}

#endif
