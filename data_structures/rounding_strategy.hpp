
#ifndef LASERCAKE_ROUNDING_STRATEGY_HPP__
#define LASERCAKE_ROUNDING_STRATEGY_HPP__

// Zero divided by something is always zero.
// Any other division result is always (before rounding)
// positive or negative.
// We specify a division strategy in terms of positive
// numbers and an indication of how negative rounding is
// related to positive rounding.
enum rounding_strategy_for_positive_numbers {
  // The negative variant makes a difference with these:
  round_down, round_up,
  round_to_nearest_with_ties_rounding_up,
  round_to_nearest_with_ties_rounding_down,
  // The negative variant doesn't make a difference with any
  // after this point.
  round_to_nearest_with_ties_rounding_to_even,
  round_to_nearest_with_ties_rounding_to_odd,
  // signal an error if the answer is not exact
  require_exact_answer,//TODO implement this in divide() and shift_right()
  // Non-integer results become the nearest "X.5".
  // This is only meaningful for types that can hold one-half.
  // It provides information that cannot be directly accessed
  // with the other strategies (namely, whether the result was exact):
  round_inexact_to_halfway,
  // "X.5" results remain "X.5"; the rest round to the nearest integer:
  round_to_nearest_with_ties_rounding_to_halfway
};
enum rounding_strategy_for_negative_numbers {
  // "doesn't make a difference" is true for unsigned arguments
  // and for round-to-even and round-to-odd.  It is a compile
  // error to claim "doesn't make a difference" when it might
  // in fact make a difference.
  negative_variant_doesnt_make_a_difference,
  
  // result invariant under negation; roughly,
  // -divide(-x, y, strat) == divide(x, y, strat)
  negative_mirrors_positive,
  
  // result invariant under addition of a constant; roughly,
  // divide(x + C*y, y, strat) - C == divide(x, y, strat)
  negative_continuous_with_positive,
  
  // these just assert that the numerator and denominator is nonnegative.
  negative_is_forbidden
};
template<
  rounding_strategy_for_positive_numbers PosStrategy,
  rounding_strategy_for_negative_numbers NegStrategy
    = negative_variant_doesnt_make_a_difference>
struct rounding_strategy {
  static const rounding_strategy_for_positive_numbers positive_strategy = PosStrategy;
  static const rounding_strategy_for_negative_numbers negative_strategy = NegStrategy;
};

#endif
