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



#include <array>
#include <memory>
#include <unordered_set>
#include <queue>
#include <boost/optional.hpp>
#include <boost/integer.hpp>
#include <boost/variant.hpp>
#include <boost/compressed_pair.hpp>
#include <boost/range/iterator_range.hpp>

#include "time_steward.hpp"

typedef ptrdiff_t num_bits_type;
typedef ptrdiff_t num_coordinates_type;
using boost::optional;
using boost::none;

template<typename value_type, typename reference = value_type&, typename pointer = value_type*>
class value_as_ptr {
public:
  reference operator*() { return v_; }
  pointer operator->() { return boost::addressof(v_); }
  value_as_ptr(reference v) : v_(v) {}
private:
  value_type v_;
};

template<num_bits_type CoordinateBits, num_coordinates_type NumDimensions>
struct bbox_collision_detector_accessories {
public:
  static const num_bits_type coordinate_bits = CoordinateBits;
  static const num_coordinates_type num_dimensions = NumDimensions;
  typedef typename boost::uint_t<coordinate_bits>::fast coordinate_type;
  typedef std::array<coordinate_type, num_dimensions> coordinate_array;
  
  static coordinate_type safe_left_shift_one(num_bits_type shift) {
    if (shift >= coordinate_bits) return 0;
    return coordinate_type(1) << shift;
  }

  static coordinate_type this_many_low_bits(num_bits_type num_bits) {
    return safe_left_shift_one(num_bits) - 1;
  }
  
  static inline coordinate_type max_in_array_of_unsigned(coordinate_array const& arr) {
    if(num_dimensions == 0) {
      return 0;
    }
    else {
      coordinate_type max_val = arr[0];
      for (size_t i = 1; i < num_dimensions; ++i) {
        if (arr[i] > max_val) max_val = arr[i];
      }
      return max_val;
    }
  }
  
// This bounding_box's coordinates are deliberately modulo.
// This lets us naturally represent either signed or unsigned
// spaces, that themselves wrap or don't (if they don't wrap,
// then the callers just won't make any bounding_boxes that cross
// their self-imposed boundary over which it doesn't wrap).
//
// Example with num_dimensions=1: bbox
//   A = {min=2,size=3} = {min=2, size_minus_one=2}
//   B = {min=5,size=1} = {min=5, size_minus_one=0}
// are adjacent but not intersecting.
//
// All of these bounding-boxes contain some space; they are at least 1x1x...
//
// Their maximum size is the width of the entire coordinate space.  This
// doesn't quite fit in the unsigned integer.  Thus it is only correct
// to use size_minus_one for most arithmetic, rather than the size itself.
//
// These 1..space-size sizes are a design decision.  Degenerate zero-size
// bounding-boxes are not much use to bbox_collision_detector.  Max-size
// bounding-boxes are.  We can't have both because of unsigned integer
// ranges.

class bounding_box {
public:
  // The default-constructor can reasonably be private or public; it creates
  // a bounding_box with undefined values.  We choose public (the default)
  // in case default-constructibility is useful to anyone for their
  // containers or similar.

  // Named "constructors".
  static bounding_box min_and_size_minus_one(coordinate_array min, coordinate_array size_minus_one) {
    bounding_box result;
    result.min_ = min;
    result.size_minus_one_ = size_minus_one;
    return result;
  }
  static bounding_box min_and_max(coordinate_array min, coordinate_array max) {
    bounding_box result;
    result.min_ = min;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      result.size_minus_one_[i] = max[i] - min[i];
    }
    return result;
  }
  static bounding_box size_minus_one_and_max(coordinate_array size_minus_one, coordinate_array max) {
    bounding_box result;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      result.min_[i] = max[i] - size_minus_one[i];
    }
    result.size_minus_one_ = size_minus_one;
    return result;
  }

  // Accessors (with array and indexing variants).
  coordinate_array min()const { return min_; }
  coordinate_array size_minus_one()const { return size_minus_one_; }
  coordinate_array max()const {
    coordinate_array result;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      result[i] = min_[i] + size_minus_one_[i];
    }
    return result;
  }
  coordinate_type min(num_coordinates_type dim)const { return min_[dim]; }
  coordinate_type size_minus_one(num_coordinates_type dim)const { return size_minus_one_[dim]; }
  coordinate_type max(num_coordinates_type dim)const { return min_[dim] + size_minus_one_[dim]; }

  // Utility functions.
  bool overlaps(bounding_box const& other)const {
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      if ( other.size_minus_one(i) <       min(i) - other.min(i)
        &&       size_minus_one(i) < other.min(i) -       min(i)) return false;
    }
    return true;
  }

  bool operator==(bounding_box const& other)const {
    return min_ == other.min_ && size_minus_one_ == other.size_minus_one_;
  }
  bool operator!=(bounding_box const& other)const {
    return !(*this == other);
  }

  friend inline std::ostream& operator<<(std::ostream& os, bounding_box const& bb) {
    os << '[';
    for (num_coordinates_type i = 0; i < bounding_box::num_dimensions; ++i) {
      if(i != 0) os << ", ";
      os << bb.min(i) << '+' << bb.size_minus_one(i);
    }
    os << ']';
    return os;
  }

private:
  coordinate_array min_;
  coordinate_array size_minus_one_;
};

// bbox_collision_detector uses "z-ordering", named such because "z" is a visual for
// the zigzag the ordering creates when drawn.  See
//     https://en.wikipedia.org/wiki/Z-order_curve
//
// A zbox is:
//
// If you interleave the coordinates' bits in the z-ordering way, it is a contiguous
// range of interleaved bits from, for example, binary VVVVV000 through VVVVV111,
// for some VVVVV.  Consider the zbox to consist of VVVVV and the number of low bits
// that vary.
// (This example could be a two-dimensional bbox_collision_detector<foo, 4, 2>,
// if four-bit integers existed -- imagine lots more bits for a real zbox.)
//
// Equivalently, a zbox is a square/cube/etc. on a maximally-bit-aligned grid of
// square/cubes of its size, or a rectangle/etc. that consists of 2/etc. of these
// square/cubes adjacent (in the case that the number of low bits in the first
// explanation is not a multiple of num_dimensions).
class zbox {
private:
  // We ensure that every bit except the ones specifically supposed to be on is off.
  // (Specifically, the "low bits" are zero.)
  coordinate_array coords_;

  //small_num_bits_type can represent at least [0, coordinate_bits*num_dimensions] inclusive.
  //smaller_num_bits_type can represent at least [0, coordinate_bits] inclusive.
  typedef typename boost::uint_t<static_num_bits_in_integer_that_are_not_leading_zeroes<coordinate_bits*num_dimensions>::value>::least small_num_bits_type;
  typedef typename boost::uint_t<static_num_bits_in_integer_that_are_not_leading_zeroes<coordinate_bits>::value>::least smaller_num_bits_type;
  small_num_bits_type num_low_bits_;
  std::array<smaller_num_bits_type, num_dimensions> dim_num_low_bits_;

public:
  zbox():num_low_bits_(coordinate_bits * num_dimensions){}

  // Named constructor idiom
  static zbox smallest_joint_parent(zbox zb1, zbox zb2) {
    zbox new_box;
    num_bits_type dim_low_bits_heuristic[num_dimensions];
    num_bits_type largest_dim_low_bits = 0;
    for (num_coordinates_type i = 0; i != num_dimensions; ++i) {
      const coordinate_type uncommon_bits =
          (zb1.coords_[i] | this_many_low_bits(zb2.num_low_bits_by_dimension(i)))
        ^ (zb2.coords_[i] & ~this_many_low_bits(zb1.num_low_bits_by_dimension(i)));
      const num_bits_type this_dimension_low_bits = num_bits_in_integer_that_are_not_leading_zeroes(uncommon_bits);
      dim_low_bits_heuristic[i] = this_dimension_low_bits;
      if(largest_dim_low_bits < this_dimension_low_bits) largest_dim_low_bits = this_dimension_low_bits;
    }
    const num_bits_type low_bits_minor = largest_dim_low_bits - 1;
    const num_bits_type low_bits_major = largest_dim_low_bits;
    num_bits_type lim_low_bits = low_bits_minor;
    num_bits_type dim_low_bits[num_dimensions];
    num_bits_type total_low_bits = 0;
    for (num_coordinates_type i = num_dimensions - 1; i >= 0; --i) {
      if(dim_low_bits_heuristic[i] ==/*a.k.a.>=*/ low_bits_major) lim_low_bits = low_bits_major;
      dim_low_bits[i] = lim_low_bits;
      new_box.dim_num_low_bits_[i] = lim_low_bits;
      total_low_bits += lim_low_bits;
    }
    new_box.num_low_bits_ = total_low_bits;
    for (num_coordinates_type i = 0; i != num_dimensions; ++i) {
      assert_if_ASSERT_EVERYTHING(
            (zb1.coords_[i] & ~this_many_low_bits(dim_low_bits[i]))
        == (zb2.coords_[i] & ~this_many_low_bits(dim_low_bits[i]))
      );
      new_box.coords_[i] = zb1.coords_[i] & ~this_many_low_bits(dim_low_bits[i]);
    }
    return new_box;
  }

  static zbox box_from_coords(coordinate_array const& coords, num_bits_type num_low_bits) {
    assert(num_low_bits >= 0 && num_low_bits <= coordinate_bits*num_dimensions);
    zbox result;
    result.num_low_bits_ = num_low_bits;
    const num_bits_type base_num_low_bits = num_low_bits / num_dimensions;
    const num_bits_type tweak_num_low_bits = num_low_bits % num_dimensions;
    for (num_coordinates_type i = 0; i != num_dimensions; ++i) {
      result.dim_num_low_bits_[i] = base_num_low_bits + (i < tweak_num_low_bits);
      result.coords_[i] = coords[i] & ~this_many_low_bits(result.num_low_bits_by_dimension(i));
    }
    return result;
  }

  bool subsumes(zbox const& other)const {
    if (other.num_low_bits_ > num_low_bits_) return false;
    for (num_coordinates_type i = 0; i != num_dimensions; ++i) {
      if (coords_[i] != (other.coords_[i] & ~this_many_low_bits(num_low_bits_by_dimension(i)))) return false;
    }
    return true;
  }
  bool overlaps(zbox const& other)const {
    for (num_coordinates_type i = 0; i != num_dimensions; ++i) {
      if ( (coords_[i] & ~this_many_low_bits(other.num_low_bits_by_dimension(i)))
        != (other.coords_[i] & ~this_many_low_bits(num_low_bits_by_dimension(i)))) {
          return false;
      }
    }
    return true;
  }
  bool overlaps(bounding_box const& bbox)const {
    for (num_coordinates_type i = 0; i != num_dimensions; ++i) {
      const coordinate_type this_size_minus_one_i = this_many_low_bits(num_low_bits_by_dimension(i));
      if (bbox.size_minus_one(i) <  coords_[i] - bbox.min(i)
      && this_size_minus_one_i  < bbox.min(i) -  coords_[i]) return false;
    }
    return true;
  }
  bool get_bit(num_bits_type bit)const {
    return coords_[bit % num_dimensions] & safe_left_shift_one(bit / num_dimensions);
  }
  num_bits_type num_low_bits_by_dimension(num_coordinates_type dim)const {
    return dim_num_low_bits_[dim];
  }
  bounding_box get_bbox()const {
    coordinate_array size_minus_one;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      size_minus_one[i] = this_many_low_bits(num_low_bits_by_dimension(i));
    }
    return bounding_box::min_and_size_minus_one(coords_, size_minus_one);
  }
  num_bits_type num_low_bits()const {
    return num_low_bits_;
  }
  bool operator==(zbox const& other)const {
    return num_low_bits_ == other.num_low_bits_ && coords_ == other.coords_;
  }
  bool operator!=(zbox const& other)const {
    return !(*this == other);
  }
  friend inline std::ostream& operator<<(std::ostream& os, zbox const& zb) {
    return os << "0x" << std::hex << zb.get_bbox() << std::dec;
  }
};

}; // struct bbox_collision_detector_accessories

template<typename SpatialEntitySubclass>
class bbox_collision_detector;



template<class TimeSteward, num_bits_type CoordinateBits, num_coordinates_type NumDimensions>
class bbox_collision_detector_spatial_entity : virtual public TimeSteward::entity {
public:
  template<typename T> using entity_id = time_steward_system::entity_id<T>;
  template<typename T> using entity_ref = time_steward_system::entity_ref<T>;
  typedef TimeSteward time_steward;
  typedef typename time_steward::time_type time_type;
  typedef typename time_steward::accessor accessor;
  typedef typename time_steward::event event_type;
  typedef typename time_steward::entity entity;
  typedef typename time_steward::event_entity event_entity;
  static const num_bits_type coordinate_bits = CoordinateBits;
  typedef typename boost::uint_t<coordinate_bits>::fast coordinate_type;
  static const num_coordinates_type num_dimensions = NumDimensions;
  typedef bbox_collision_detector_accessories<coordinate_bits, num_dimensions> accessories;
  typedef typename accessories::bounding_box bounding_box;
  typedef typename accessories::coordinate_array coordinate_array;
  
  //For the subclass to implement:
  //
  //static bounding_box bbox(accessor const* accessor,
  //                         entity_ref<const SpatialEntitySubclass> ent) = 0;
  //
  //static time_type escape_time(accessor const* accessor,
  //                             entity_ref<const SpatialEntitySubclass> ent,
  //                             bounding_box const& bbox) = 0;
  //
  // returns nullptr for "no interaction"
  //static std::unique_ptr<event_type> next_interaction(accessor const* accessor,
  //                       entity_ref<const SpatialEntitySubclass> ent1,
  //                       entity_ref<const SpatialEntitySubclass> ent2) = 0;
private:
  entity_id<event_entity> next_event;
  optional<bounding_box> old_bbox;
  template<typename SpatialEntitySubclass> friend class bbox_collision_detector;
};

template<typename SpatialEntitySubclass>
class bbox_collision_detector : public SpatialEntitySubclass::time_steward::entity {
public:
  bbox_collision_detector* clone()const override { return new bbox_collision_detector(*this); }
private:
  static const num_bits_type coordinate_bits = SpatialEntitySubclass::coordinate_bits;
  static const num_coordinates_type num_dimensions = SpatialEntitySubclass::num_dimensions;
  static_assert(num_dimensions >= 0, "You can't make a space with negative dimensions!");
  typedef ptrdiff_t num_zboxes_type;
  static_assert(num_dimensions < std::numeric_limits<num_zboxes_type>::digits - 1,
    "We don't permit so many dimensions that one bounding_box might need more zboxes than we can count.");
  static_assert(coordinate_bits >= 0, "You can't have an int type with negative bits!");
  typedef bbox_collision_detector_accessories<coordinate_bits, num_dimensions> accessories;
  template<typename T> using entity_id = time_steward_system::entity_id<T>;
  template<typename T> using entity_ref = time_steward_system::entity_ref<T>;
public:
  typedef typename SpatialEntitySubclass::time_steward time_steward;
  typedef typename time_steward::time_type time_type;
  typedef typename accessories::coordinate_type coordinate_type;
  typedef typename accessories::coordinate_array coordinate_array;
  typedef typename accessories::bounding_box bounding_box;
private:
  typedef typename time_steward::accessor accessor;
  typedef typename time_steward::event event_type;
  typedef typename time_steward::entity entity;
  typedef typename time_steward::event_entity event_entity;
  typedef typename accessories::zbox zbox;
  
public:
  // The bbox_collision_detector has a "z-tree".  This
  // is a binary tree.  Keys in the tree are zboxes (see above).
  // Values are entity_ids; each zbox may have any number
  // of entity_ids.  This code is happy for objects to overlap
  // each other, and besides, even non-overlapping objects often
  // have a minimal containing zbox in common.
  //
  // Because of the definition of zboxes, either they are
  // A: the same, and thus the same ztree_node
  // B: one is smaller and fully within the node, and it is a
  //          descendant of the other
  // C: they don't overlap at all, and in the ztree neither is a
  //          descendant of the other.
  // (In particular, zboxes can't partially overlap each other.)
  //
  // An object may need to be in up to (2**num_dimensions) zboxes
  // so that the area covered by its zboxes is only a constant
  // factor larger than the object's regular (non-z-order) bounding
  // box.  Consider an object that's a box of width 2 or 3 with
  // x min-coordinate 10111111 and max 11000001 (binary).  The minimal
  // common prefix there is just a single bit; a single bit means
  // a huge box.  Conceptually, dimensions' bits are interleaved
  // before looking for a common prefix; any dimension has the
  // potential to differ in a high bit and disrupt the common prefix.
  // Splitting the object across two zboxes, for each dimension,
  // is sufficient to avoid this explosion.
  //
  // Specifically, a ztree is a Patricia trie on z-ordered
  // bits (zboxes), where the bits are seen as a sequence with the
  // highest-order bits first and the sequence-length equal to the
  // number of high bits specified by a given key/ztree_node/zbox.
  // The ztree_node happens to contain (unlike typical tries) the
  // entire key that it represents, because the key is small and
  // constant-sized and it's generally easier to do so.
  //
  // What goes into child0 vs. child1?  If trying to insert, say,
  // the zbox B = 10100??? (binary, 5 high bits, 3 low bits) into
  // A = 10??????, B goes at or below A's child1 because B's next bit
  // after A's bits is 1.  (Current: 10, next: 101).
  //
  // If there would be a node with zero entity_ids and
  // only one child, then that child node goes there directly
  // instead of that trivial node.  If the tree, to be correct,
  // needs nodes with two children and zero entity_ids,
  // then it will have them.
  /*

  If there's one zbox in the tree [call it Z]

  tree = ztree_node {
    Z
    nullptr
    nullptr
  }

  Two zboxes that differ at bit B:
  child0 of a node with B ignored bits is the child whose Bth bit is 0.
  tree = ztree_node {
    the common leading bits of Z1 and Z2, with B ignored bits
    ptr to ztree_node {
      Z1
      nullptr
      nullptr
    }
    ptr to ztree_node {
      Z2
      nullptr
      nullptr
    }
  }

  */
  struct ztree_node {
    typedef std::shared_ptr<const ztree_node> ztree_node_ptr;
    
    // Making this a persistent data structure could improve performance
    // in the case where there are many entities at the same node.
    // However, it cannot improve the asymptotic speed, since we always
    // need to review this many entities when one is added or removed.
    // It might improve asymptotic memory use.
    typedef std::unordered_set<entity_id<const SpatialEntitySubclass>> objects_here_type;

    const zbox here;
    ztree_node_ptr child0;
    ztree_node_ptr child1;

    objects_here_type objects_here;

    ztree_node(zbox box):here(box),child0(nullptr),child1(nullptr){}
    ztree_node(ztree_node const& other) :
      here(other.here),
      child0(other.child0),
      child1(other.child1),
      objects_here(other.objects_here)
      {}
    // operator= could exist if we wanted to make zbox non-const.
    // ztree_node& operator=(ztree_node const& other) = delete;
  };
  typedef typename ztree_node::ztree_node_ptr ztree_node_ptr;
  typedef std::shared_ptr<ztree_node> ztree_node_mutable_ptr;
  
private:
  class spatial_entity_escapes_its_zboxes : public event_entity {
  public:
    spatial_entity_escapes_its_zboxes(entity_id<bbox_collision_detector> bbcd_id, entity_id<SpatialEntitySubclass> spatial_entity_id, time_type when) : bbcd_id(bbcd_id),spatial_entity_id(spatial_entity_id),when_(when) {}
    spatial_entity_escapes_its_zboxes* clone()const override { return new spatial_entity_escapes_its_zboxes(*this); }

    void operator()(accessor* accessor)const override {
      auto cd = accessor->get_mut(bbcd_id);
      auto ref = accessor->get_mut(spatial_entity_id);
      erase(accessor, cd, ref);
      insert(accessor, cd, ref);
    }
    time_type when()const override {
      return when_;
    }
    entity_id<bbox_collision_detector> bbcd_id;
    entity_id<SpatialEntitySubclass> spatial_entity_id;
    time_type when_;
  };
  class spatial_entity_interaction : public event_entity {
  public:
    spatial_entity_interaction(entity_id<bbox_collision_detector> bbcd_id, entity_id<SpatialEntitySubclass> id0, entity_id<SpatialEntitySubclass> id1, std::shared_ptr<event_type>&& event) : bbcd_id(bbcd_id),id0(id0),id1(id1),event(event) {}
    spatial_entity_interaction* clone()const override { return new spatial_entity_interaction(*this); }

    void operator()(accessor* accessor)const override {
      auto cd = accessor->get_mut(bbcd_id);
      auto ref0 = accessor->get_mut(id0);
      auto ref1 = accessor->get_mut(id1);
      erase(accessor, cd, ref0);
      erase(accessor, cd, ref1);
      (*event)(accessor);
      insert(accessor, cd, ref0);
      insert(accessor, cd, ref1);
    }
    time_type when()const override {
      return event->when();
    }
    entity_id<bbox_collision_detector> bbcd_id;
    entity_id<SpatialEntitySubclass> id0, id1;
    std::shared_ptr<event_type> event;
  };
  
public:
  bbox_collision_detector():objects_tree_(nullptr){}
  
  static void insert(accessor* accessor, entity_ref<bbox_collision_detector> bbcd, entity_ref<SpatialEntitySubclass> e) {
    if (e) {
      // TODO how to handle double inserts
      e->old_bbox = SpatialEntitySubclass::bbox(accessor, e);
      bbcd->objects_tree_ = insert_object(accessor, bbcd->objects_tree_, e, *(e->old_bbox), bbcd.id());
    }
  }
  static void erase(accessor* accessor, entity_ref<bbox_collision_detector> bbcd, entity_ref<SpatialEntitySubclass> e) {
    if (e->next_event) {
      accessor->delete_entity(e->next_event);
      e->next_event = entity_id<event_entity>();
    }
    bbcd->objects_tree_ = delete_object(bbcd->objects_tree_, e.id(), *(e->old_bbox));
    e->old_bbox = none;
  }
  
private:
  bbox_collision_detector(ztree_node_ptr objects_tree):objects_tree_(objects_tree){}
  ztree_node_ptr objects_tree_;
  
  static ztree_node_ptr insert_object(accessor* accessor, ztree_node_ptr tree, entity_ref<SpatialEntitySubclass> e, bounding_box const& bbox, entity_id<bbox_collision_detector> bbcd_id) {
    const coordinate_type max_width_minus_one = accessories::max_in_array_of_unsigned(bbox.size_minus_one());
    // max_width - 1: power-of-two-sized objects easily squeeze into the next smaller category.
    // i.e., exp = log2_rounding_up(max_width)
    const num_bits_type exp = num_bits_in_integer_that_are_not_leading_zeroes(max_width_minus_one);
    const coordinate_type base_box_size = accessories::safe_left_shift_one(exp);

    // The total number of zboxes we use to cover this bounding_box
    // is a power of two between 1 and 2**num_dimensions.
    // Imagine that we start with a set of one box and that,
    // for each dimension, we start with the previous set of boxes,
    // then zbox-ify this dimension, using either
    //   (A) exactly base_box_size width-in-this-dimension, or
    //   (B) twice base_box_size width-in-this-dimension, or
    // if the bit parity didn't work out so well, it
    //   (C) needs to split each zbox into two.
    num_coordinates_type num_dims_using_one_zbox_of_exactly_base_box_size = 0;
    num_coordinates_type num_dims_using_one_zbox_of_twice_base_box_size = 0;
    num_coordinates_type num_dims_using_two_zboxes_each_of_base_box_size;

    // Given that a coordinate is laid out in bits like XYZXYZXYZ,
    // We're at some exp in there (counted from the right); let's say 3.
    // Given exp 3, Z is a less-significant bit and X is more-significant.
    // ('exp' could also be a non-multiple-of-num_dimensions, in which case
    // the ordering of the dimensions would come out differently.)
    //
    // If the object happens to fit, aligned, in X with width base_box_size,
    // then we can just specify this X bit directly.  If that works for X,
    // we can try Y; if not, we can't try Y because X is already doing
    // something nontrivial (perhaps it could be done; the code would
    // be more complicated).  These are
    // "num_dimensions_that_need_one_zbox_of_exactly_base_box_size".
    // This is the best case.
    //
    // Then, if we can't specify where in all dimensions we are
    // z-box-ly yet in one zbox, we try starting from Z:
    // it's possible that Z (and so forth if Z is) can be included
    // in the zbox's low_bits.  This would make the zbox twice
    // as wide in that dimension (e.g. Z) as it would be in the
    // case of if X fits into a single base_box_size at the scale
    // we're looking at.  But it's better than making two separate
    // zboxes that take up that much space anyway.  If the bounding
    // box didn't happen to be aligned with the right parity,
    // we'll have to make two boxes for it anyway instead of putting
    // it in "num_dimensions_that_need_one_zbox_of_twice_base_box_size".
    //
    // All the dimensions in between will be split into two zboxes,
    // each of width base_box_size, for a total width of
    // twice base_box_size.  This is the worst case,
    // "num_dimensions_that_need_two_zboxes_each_of_base_box_size".
    if(base_box_size == 0) {
      num_dims_using_one_zbox_of_exactly_base_box_size = num_dimensions;
    }
    else {
      for (num_coordinates_type i = num_dimensions - 1; i >= 0; --i) {
        if (bbox.size_minus_one(i) <= (base_box_size - 1) - (bbox.min(i) & (base_box_size - 1))) {
          ++num_dims_using_one_zbox_of_exactly_base_box_size;
        }
        else {
          break;
        }
      }
      for (num_coordinates_type i = 0; i < num_dimensions - num_dims_using_one_zbox_of_exactly_base_box_size; ++i) {
        if (!(bbox.min(i) & base_box_size)) {
          ++num_dims_using_one_zbox_of_twice_base_box_size;
        }
        else {
          break;
        }
      }
    }
    num_dims_using_two_zboxes_each_of_base_box_size = num_dimensions - num_dims_using_one_zbox_of_exactly_base_box_size - num_dims_using_one_zbox_of_twice_base_box_size;

    const num_zboxes_type number_of_zboxes_to_use_if_necessary = num_zboxes_type(1) << num_dims_using_two_zboxes_each_of_base_box_size;

    coordinate_array zboxes_union_min;
    coordinate_array zboxes_union_size_minus_one;
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      zboxes_union_min[i] = bbox.min()[i] & ~(base_box_size-1);
      zboxes_union_size_minus_one[i] = ((i < num_dimensions - num_dims_using_one_zbox_of_exactly_base_box_size) ? (base_box_size<<1) : base_box_size)-1;
    }
    const bounding_box zboxes_union = bounding_box::min_and_size_minus_one(zboxes_union_min, zboxes_union_size_minus_one);
    std::unique_ptr<event_entity> event;
    const time_type escape_time = SpatialEntitySubclass::escape_time(accessor, e, zboxes_union);
    if (escape_time != time_steward::never) {
      assert(escape_time > accessor->now()); // TODO an exception
      event = make_unique<spatial_entity_escapes_its_zboxes>(bbcd_id, e.id(), escape_time);
    }
    for (num_zboxes_type i = 0; i < number_of_zboxes_to_use_if_necessary; ++i) {
      coordinate_array coords = bbox.min();
      for (num_coordinates_type j = num_dims_using_one_zbox_of_twice_base_box_size; j < num_dimensions - num_dims_using_one_zbox_of_exactly_base_box_size; ++j) {
        // By checking this bit of "i" arbitrarily, by the last time
        // we get through the "number_of_zboxes_to_use_if_necessary" loop,
        // we have yielded every combination of possibilities of
        // each dimension that varies varying (between +0 and +base_box_size).
        if (i & (1 << (j - num_dims_using_one_zbox_of_twice_base_box_size))) {
          coords[j] += base_box_size;
        }
      }
      const zbox zb = zbox::box_from_coords(coords, exp * num_dimensions + num_dims_using_one_zbox_of_twice_base_box_size);
      if (zb.overlaps(bbox)) {
        tree = insert_zbox(accessor, event, tree, e, zb, bbcd_id);
      }
    }
    if (event) {
      e->next_event = accessor->create_entity(std::move(event)).id();
    }
    return tree;
  }
  static ztree_node_ptr insert_zbox(accessor* accessor, std::unique_ptr<event_entity>& event, ztree_node_ptr tree, entity_ref<const SpatialEntitySubclass> e, zbox box, entity_id<bbox_collision_detector> bbcd_id) {
    if (!tree) {
      ztree_node_mutable_ptr new_tree(new ztree_node(box));
      new_tree->objects_here.insert(e.id());
      return new_tree;
    }
    else {
      if (tree->here.subsumes(box)) {
        ztree_node_mutable_ptr new_tree(new ztree_node(*tree));
        if (box.num_low_bits() == tree->here.num_low_bits()) {
          gather_events(accessor, event, tree, e, bbcd_id, true);
          new_tree->objects_here.insert(e.id());
        }
        else {
          gather_events(accessor, event, tree, e, bbcd_id, false);
          if (box.get_bit(tree->here.num_low_bits() - 1)) new_tree->child1 = insert_zbox(accessor, event, new_tree->child1, e, box, bbcd_id);
          else                                            new_tree->child0 = insert_zbox(accessor, event, new_tree->child0, e, box, bbcd_id);
        }
        return new_tree;
      }
      else {
        ztree_node_mutable_ptr new_tree(new ztree_node(zbox::smallest_joint_parent(tree->here, box)));

        assert_if_ASSERT_EVERYTHING(new_tree->here.num_low_bits() > tree->here.num_low_bits());
        assert_if_ASSERT_EVERYTHING(new_tree->here.subsumes(tree->here));
        assert_if_ASSERT_EVERYTHING(new_tree->here.subsumes(box));
        assert_if_ASSERT_EVERYTHING(box.subsumes(tree->here) || (tree->here.get_bit(new_tree->here.num_low_bits() - 1) != box.get_bit(new_tree->here.num_low_bits() - 1)));

        if (tree->here.get_bit(new_tree->here.num_low_bits() - 1)) new_tree->child1 = tree;
        else                                                       new_tree->child0 = tree;

        return insert_zbox(accessor, event, new_tree, e, box, bbcd_id);
      }
    }
  }
  static void gather_events(accessor* accessor, std::unique_ptr<event_entity>& event, ztree_node_ptr tree, entity_ref<const SpatialEntitySubclass> e, entity_id<bbox_collision_detector> bbcd_id, bool recurse) {
    if (tree) {
      for (entity_id<const SpatialEntitySubclass> other_id : tree->objects_here) {
        std::unique_ptr<event_type> interaction(SpatialEntitySubclass::next_interaction(accessor, e, accessor->get(other_id)));
        if (interaction && ((!event) || (interaction->when() < event->when()))) { // TODO: fix the fact that the order we collect the events creates a directional bias
          assert(interaction->when() > accessor->now()); // TODO an exception
          // TODO: why did I have to disambiguate with futurestd:: here when I didn't have to elsewhere?
          event = futurestd::make_unique<spatial_entity_interaction>(bbcd_id, 
           time_steward_system::reinterpret_entity_id<SpatialEntitySubclass>(e.id()),
           time_steward_system::reinterpret_entity_id<SpatialEntitySubclass>(other_id),
           std::move(interaction));
        }
      }
      if (recurse) {
        gather_events(accessor, event, tree->child0, e, bbcd_id, true);
        gather_events(accessor, event, tree->child1, e, bbcd_id, true);
      }
    }
  }
  

  static ztree_node_ptr delete_object(ztree_node_ptr tree, entity_id<SpatialEntitySubclass> id, bounding_box const& bbox) {
    if (tree && tree->here.overlaps(bbox)) {
      ztree_node_mutable_ptr new_tree(new ztree_node(*tree));
      new_tree->objects_here.erase(id);
      new_tree->child0 = delete_object(new_tree->child0, id, bbox);
      new_tree->child1 = delete_object(new_tree->child1, id, bbox);

      // collapse nodes with no objects and 0-1 children.
      if (new_tree->objects_here.empty()) {
        if (!new_tree->child0) {
          // (old 'child1' a.k.a. new 'tree' could be nullptr)
          return new_tree->child1;
        }
        else if (!new_tree->child1) {
          return new_tree->child0;
        }
      }
      return new_tree;
    }
    return tree;
  }
  

  
  
  
  
  
  template<typename OrderingFunctor, typename ComparedType>
  struct reverse_first_ordering : private OrderingFunctor {
    typedef bool result_type; typedef ComparedType first_argument_type; typedef ComparedType second_argument_type;
    reverse_first_ordering() : OrderingFunctor() {}
    reverse_first_ordering(OrderingFunctor const& o) : OrderingFunctor(o) {}

    bool operator()(ComparedType const& a, ComparedType const& b) {
      return static_cast<OrderingFunctor&>(*this)(b.first(), a.first());
    }
  };

  template<typename GetCost, typename CostOrdering = std::less<typename GetCost::cost_type> >
  class iterator {
  private:
    typedef typename GetCost::cost_type cost_type;

    typedef boost::variant<ztree_node const*, entity_id<SpatialEntitySubclass>> node_variant_type;
    typedef boost::compressed_pair<cost_type, node_variant_type> queue_value_type_;
    // swap greater/less because pq sorts by greatest and we want least by default (as sorting normally is)
    typedef std::priority_queue<queue_value_type_, std::vector<queue_value_type_>, reverse_first_ordering<CostOrdering, queue_value_type_> > queue_type_;

    struct iteree {
      cost_type cost;
      entity_id<SpatialEntitySubclass> const& id;
      iteree(cost_type const& cost, entity_id<SpatialEntitySubclass> id)
        : cost(cost), id(id) {}
    };
  public:
    typedef iteree value_type;
    typedef value_type reference;
    typedef value_as_ptr<value_type, reference> pointer;
    typedef ptrdiff_t difference_type; //Ha. Ha.
    typedef std::input_iterator_tag iterator_category;

  private:
    struct contents_ : private GetCost {
      contents_(bbox_collision_detector const& bbcd, GetCost const& getcost, CostOrdering const& costordering)
      : GetCost(getcost), queue_(costordering), seen_()
  #ifdef BBOX_COLLISION_DETECTOR_DEBUG
      , bbcd_(&bbcd)
      , revision_count_(bbcd.revision_count_)
  #endif
      {}

      contents_(contents_&&) = default;
      queue_type_ queue_;
      std::unordered_set<entity_id<SpatialEntitySubclass>> seen_;
  #ifdef BBOX_COLLISION_DETECTOR_DEBUG
      bbox_collision_detector const* bbcd_;
      size_t revision_count_;
  #endif

      GetCost& get_get_cost() { return *this; }

      template<typename VariantMember>
      inline void push_cost(cost_type const& cost, VariantMember v) {
        queue_.push(queue_value_type_(cost, v));
      }
      template<typename IndirectCostType, typename VariantMember>
      inline typename boost::disable_if<boost::is_convertible<IndirectCostType, cost_type> >::type
      push_cost(IndirectCostType const& maybe_cost, VariantMember v) {
        if(maybe_cost) {
          queue_.push(queue_value_type_(*maybe_cost, v));
        }
      }

      void add_child(ztree_node const* child) {
        if(child) {
          push_cost(get_get_cost().min_cost(child->here.get_bbox()), child);
        }
      }
      void add_child(entity_id<SpatialEntitySubclass> id) {
        if(!seen_.insert(id)->second) {
          push_cost(get_get_cost().cost(id), id);
        }
      }

      void add_children_of(ztree_node const* node) {
        add_child(node->child0);
        add_child(node->child1);
        for(entity_id<SpatialEntitySubclass> id : node->objects_here) {
          add_child(id);
        }
      }
    };

    //dynamic allocation, pooh.
    boost::shared_ptr<contents_> c_;

    void advance_to_a_returnable_() {
      if(c_) {
#ifdef BBOX_COLLISION_DETECTOR_DEBUG
      caller_correct_if(c_->revision_count_ == c_->detector_->revision_count_,
                        "Error: using a bbox_collision_detector iterator "
                        "after the container has changed!");
#endif
        while(true) {
          if(c_->queue_.empty()) {
            c_.reset();
            break;
          }
          if(ztree_node const*const* node_top_ptr = boost::get<ztree_node const*>(&c_->queue_.top().second())) {
            ztree_node const* node_top = *node_top_ptr;
            c_->queue_.pop();
            c_->add_children_of(node_top);
          }
          else {
            break;
          }
        }
      }
    }

    struct unspecified_bool_{int member; private:unspecified_bool_();};
    typedef int unspecified_bool_::* unspecified_bool_type;

  public:
    iterator() : c_() {}
    explicit iterator(bbox_collision_detector const& bbcd,
                      GetCost const& getcost = GetCost(),
                      CostOrdering const& costordering = CostOrdering())
      : c_(new contents_(bbcd, getcost, costordering)) {}
    template<typename T>
    explicit iterator(bbox_collision_detector const& bbcd,
                      T const& initial,
                      GetCost const& getcost = GetCost(),
                      CostOrdering const& costordering = CostOrdering())
      : c_(new contents_(bbcd, getcost, costordering)) {
      push(initial);
    }

    template<typename T>
    void push(T const& v) { c_->add_child(v); advance_to_a_returnable_(); }
    template<typename InputIterator>
    void push(InputIterator begin, InputIterator end) {
      for( ; begin != end; ++begin) { c_->add_child(*begin); }
      advance_to_a_returnable_();
    }

    reference operator*() const {
      caller_error_if(c_->queue_.empty(), "can't dereference an empty iterator");
      queue_value_type_ const& top = c_->queue_.top();
      return iteree(top.first(), boost::get<entity_id<SpatialEntitySubclass>>(top.second()));
    }
    pointer operator->() const { return pointer(*(*this)); }
    pointer operator++(int) {
      pointer result(*(*this));
      ++*this;
      return result;
    }
    iterator& operator++() {
      c_->queue_.pop(); advance_to_a_returnable_(); return *this;
    }

    operator unspecified_bool_type() const { return c_ ? &unspecified_bool_::member : nullptr; }

    bool operator==(iterator const& other) const { return c_ == other.c_; }
    bool operator!=(iterator const& other) const { return !(*this == other); }

    iterator(iterator const&) = default;
    iterator(iterator&&) = default;
    iterator& operator=(iterator&&) = default;
  };

  template<typename GetCost>
  static inline bool bool_with_error_if_implicit_conversions_to_bool_and_cost_type_are_ambiguous(bool arg) { return arg; }
  template<typename GetCost>
  static inline bool bool_with_error_if_implicit_conversions_to_bool_and_cost_type_are_ambiguous(typename GetCost::cost_type const&) { return true; }

  template<typename GetCostBool>
  static inline void filter_impl(
        ztree_node const* tree,
        std::unordered_set<entity_id<const SpatialEntitySubclass>>& results,
        GetCostBool& getcost) {
    if (tree) {
      if(bool_with_error_if_implicit_conversions_to_bool_and_cost_type_are_ambiguous<GetCostBool>(
            getcost.min_cost(tree->here.get_bbox()))) {
        for(entity_id<const SpatialEntitySubclass> id : tree->objects_here) {
          if(bool_with_error_if_implicit_conversions_to_bool_and_cost_type_are_ambiguous<GetCostBool>(
                getcost.cost(id))) {
            results.insert(id);
          }
        }
        filter_impl(tree->child0.get(), results, getcost);
        filter_impl(tree->child1.get(), results, getcost);
      }
    }
  }

public:
  template<typename GetCost>
  inline boost::iterator_range<iterator<GetCost>> iterate(GetCost const& getcost) const {
    typedef iterator<GetCost> iter;
    boost::iterator_range<iter> result(iter(*this, objects_tree_.get(), getcost), iter());
    return result;
  }
  
  template<typename GetCost>
  inline boost::optional<typename iterator<GetCost>::value_type> find_least(GetCost const& getcost) const {
    typedef iterator<GetCost> iter;
    typedef boost::optional<typename iterator<GetCost>::value_type> result_type;
    iter i(*this, objects_tree_.get(), getcost);
    return i ? result_type(*i) : result_type();
  }

  template<typename GetCostBool>
  inline std::unordered_set<entity_id<const SpatialEntitySubclass>> filter(GetCostBool getcost)const {
    std::unordered_set<entity_id<const SpatialEntitySubclass>> results;
    filter_impl(objects_tree_.get(), results, getcost);
    return results;
  }

public:
  ztree_node const* debug_get_tree()const { return &*objects_tree_; }
};
