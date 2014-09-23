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
template<typename T>
using optional = time_steward_system::optional<T>;
using time_steward_system::none;

template<typename value_type, typename reference = value_type&, typename pointer = value_type*>
class value_as_ptr {
public:
  reference operator*() { return v_; }
  pointer operator->() { return boost::addressof(v_); }
  value_as_ptr(reference v) : v_(v) {}
private:
  value_type v_;
};

template<typename T, typename ...List>
struct type_is_in_list;
template<typename T>
struct type_is_in_list<T> { static const bool value = false; };
template<typename T, typename Head, typename ...Tail>
struct type_is_in_list<T, Head, Tail...> { static const bool value = std::is_same<T, Head>::value || type_is_in_list<T, Tail...>::value; };

template<num_bits_type CoordinateBits, num_coordinates_type NumDimensions>
struct bbox_collision_detector_system {
public:
  static const num_bits_type coordinate_bits = CoordinateBits;
  static const num_coordinates_type num_dimensions = NumDimensions;
private:
  static_assert(num_dimensions >= 0, "You can't make a space with negative dimensions!");
  typedef ptrdiff_t num_zboxes_type;
  static_assert(num_dimensions < std::numeric_limits<num_zboxes_type>::digits - 1,
    "We don't permit so many dimensions that one bounding_box might need more zboxes than we can count.");
  static_assert(coordinate_bits >= 0, "You can't have an int type with negative bits!");
public:
  typedef typename boost::uint_t<coordinate_bits>::fast coordinate_type;
  typedef std::array<coordinate_type, num_dimensions> coordinate_array;
  typedef time_steward_system::entity_id entity_id;
  
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
  bool subsumes(bounding_box const& other)const {
    for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
      // relying on unsignedness
      if (other.min(i) - min(i)                           > size_minus_one(i)) return false;
      if (other.min(i) - min(i) + other.size_minus_one(i) > size_minus_one(i)) return false;
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
      const coordinate_type nonshared_bits =
          (zb1.coords_[i] | this_many_low_bits(zb2.num_low_bits_by_dimension(i)))
        ^ (zb2.coords_[i] & ~this_many_low_bits(zb1.num_low_bits_by_dimension(i)));
      const num_bits_type this_dimension_low_bits = num_bits_in_integer_that_are_not_leading_zeroes(nonshared_bits);
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
// What goes into children[0] vs. children[1]?  If trying to insert, say,
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

root_node = ztree_node {
  Z
  entity_id()
  entity_id()
  entity_id()
}

Two zboxes that differ at bit B:
child0 of a node with B ignored bits is the child whose Bth bit is 0.
root_node = ztree_node at entity_id FOO {
  the common leading bits of Z1 and Z2, with B ignored bits
  entity_id()
  entity_id to ztree_node {
    Z1
    FOO
    entity_id()
    entity_id()
  }
  entity_id to ztree_node {
    Z2
    FOO
    entity_id()
    entity_id()
  }
}

*/
struct ztree_node {
  zbox here;
  entity_id parent;
  std::array<entity_id, 2> children;
  persistent_set<entity_id> objects_here;

  ztree_node(zbox box, entity_id parent):here(box),parent(parent){}
};
struct bbcd_entry_metadata {
  persistent_set<entity_id> nodes;
  bounding_box zboxes_union;
};
struct spatial_entity_metadata {
  persistent_map<entity_id, bbcd_entry_metadata> data;
};
struct bbox_collision_detector_root_node {
  bbox_collision_detector_root_node(entity_id root_node_id_):root_node_id_(root_node_id_){}
  entity_id root_node_id_;
};

typedef time_steward_system::fields_list<ztree_node, spatial_entity_metadata, bbox_collision_detector_root_node> fields;
  
// For the FuncsType to implement:
//
// bounding_box bbox(accessor const* accessor, entity_ref e)
// time_type escape_time(accessor const* accessor, entity_ref e, bounding_box const& bbox)
// void anticipate_interactions(accessor* accessor, entity_ref e0, entity_ref e1)
//
// TODO: if FuncsType::required_field exists,
// FuncsType is constructed with one argument: the given field of the bbox_collision_detector root entity.
// Otherwise, FuncsType is default-constructed.
template<typename TimeSteward, class FuncsType>
class operations {
public:
  typedef TimeSteward time_steward;
  typedef typename time_steward::time_type time_type;
private:
  typedef typename time_steward::accessor accessor;
  typedef typename time_steward::event event;
  typedef typename time_steward::trigger trigger;
  typedef typename time_steward_system::trigger_id trigger_id;
  typedef typename time_steward_system::entity_id entity_id;
  typedef typename accessor::entity_ref entity_ref;
  
private:
  static inline FuncsType get_funcs(accessor*, entity_id) { return FuncsType(); }
  
  class spatial_entity_escapes_its_zboxes : public time_steward::event {
  public:
    spatial_entity_escapes_its_zboxes(entity_id bbcd_id, entity_id spatial_entity_id) : bbcd_id(bbcd_id),spatial_entity_id(spatial_entity_id) {}

    void operator()(accessor* accessor)const override {
      entity_ref e = accessor->get(spatial_entity_id);
      FuncsType funcs = get_funcs(accessor, bbcd_id);
      bbcd_entry_metadata const& metadata = accessor->template get<spatial_entity_metadata>(e)->data.find(bbcd_id)->second;
      update_zboxes(accessor, bbcd_id, e, funcs.bbox(accessor, e), metadata.nodes);
      gather_events(accessor, bbcd_id, funcs, e);
    }
    entity_id bbcd_id;
    entity_id spatial_entity_id;
  };
  class spatial_entity_changes_trigger : public time_steward::event {
  public:
    spatial_entity_changes_trigger(entity_id bbcd_id, entity_id spatial_entity_id) : bbcd_id(bbcd_id),spatial_entity_id(spatial_entity_id) {}

    void operator()(accessor* accessor)const override {
      entity_ref e = accessor->get(spatial_entity_id);
      FuncsType funcs = get_funcs(accessor, bbcd_id);
      bbcd_entry_metadata const& metadata = accessor->template get<spatial_entity_metadata>(e)->data.find(bbcd_id)->second;
      const bounding_box bbox = funcs.bbox(accessor, e);
      if (!metadata.zboxes_union.subsumes(bbox)) {
        update_zboxes(accessor, bbcd_id, e, bbox, metadata.nodes);
      }
      gather_events(accessor, bbcd_id, funcs, e);
    }
    entity_id bbcd_id;
    entity_id spatial_entity_id;
  };
  
public:
  static entity_ref create_bbox_collision_detector(accessor* accessor) {
    entity_ref result = accessor->create_entity();
    accessor->template set<bbox_collision_detector_root_node>(result, bbox_collision_detector_root_node(entity_id()));
    return result;
  }
  static void insert(accessor* accessor, entity_id bbcd_id, entity_ref e, entity_ref hint_object = entity_ref()) {
    // TODO maybe time_steward can treat the different members of the outer_metadata map as different fields,
    // (WRT access/back-in-time-change semantics) since we usually only want to access/modify one of them?
    // Or is that too unlikely to be useful? Some sort of profiling may help later.
    auto& outer_metadata = accessor->template get_mut<spatial_entity_metadata>(e);
    if (!outer_metadata) { outer_metadata = spatial_entity_metadata(); }
    bbcd_entry_metadata& metadata = outer_metadata->data[bbcd_id]; // default-construct
    if (hint_object) {
      // hack: copy over the object's nodes as a hint for where to insert
      metadata.nodes = accessor->template get<spatial_entity_metadata>(e)->data.find(bbcd_id)->second.nodes;
    }
    
    accessor->set_trigger(trigger_id(bbcd_id, e.id()), std::shared_ptr<trigger>(new spatial_entity_changes_trigger(bbcd_id, e.id())));
  }
  static void erase(accessor* accessor, entity_id bbcd_id, entity_ref e) {
    accessor->set_trigger(trigger_id(bbcd_id, e.id()));
    // TODO: what about the outstanding events?
    
    bbcd_entry_metadata& metadata = accessor->template get_mut<spatial_entity_metadata>(e)->data.find(bbcd_id)->second;
    erase_from_nodes(accessor, bbcd_id, e, metadata, metadata->nodes);
    accessor->template get_mut<spatial_entity_metadata>(e)->erase(bbcd_id);
  }
  
private:
  
  static void erase_from_nodes(accessor* accessor, entity_id bbcd_id, entity_ref e,
                                    bbcd_entry_metadata& metadata, persistent_set<entity_id> nodes) {
    for (entity_id node_id : nodes) {
      entity_ref node_ref = accessor->get(node_id);
      auto& node = accessor->template get_mut<ztree_node>(node_ref);
      if (node) {
        node->objects_here.erase(e.id());
        if (node->objects_here.empty()) {
          if (!node->children[0]) {
            squish_node(accessor, bbcd_id, node_ref, 1);
          }
          else if (!node->children[1]) {
            squish_node(accessor, bbcd_id, node_ref, 0);
          }
        }
      }
      metadata.nodes.erase(node_id);
    }
  }
  
  static void squish_node(accessor* accessor, entity_id bbcd_id, entity_ref node_ref, size_t which_child) {
    // (old child a.k.a. new 'node' could be none)
    auto& node = accessor->template get_mut<ztree_node>(node_ref);
    entity_id stolen_node_id = node->children[which_child];
    auto& stolen_node = accessor->template get_mut<ztree_node>(accessor->get(stolen_node_id));
    node = stolen_node;
    stolen_node = none;
    if (node) {
      for (entity_id stolen_child_id : node->children) {
        if (stolen_child_id) {
          accessor->template get_mut<ztree_node>(accessor->get(stolen_child_id))->parent = node_ref.id();
        }
      }
      for (entity_id stolen_object_id : node->objects_here) {
        auto& nodes = accessor->template get_mut<spatial_entity_metadata>(accessor->get(stolen_object_id))->data.find(bbcd_id)->second.nodes;
        nodes.erase(stolen_node_id);
        nodes.insert(node_ref.id());
      }
    }
  }
  
  
  class update_zboxes {
    accessor* accessor_;
    entity_id bbcd_id;
    entity_ref e;
    bbcd_entry_metadata& metadata;
  public:
    update_zboxes(accessor* accessor, entity_id bbcd_id, entity_ref e, bounding_box bbox, persistent_set<entity_id> hint_nodes)
      :
      accessor_(accessor),
      bbcd_id(bbcd_id),
      e(e),
      metadata(accessor->template get_mut<spatial_entity_metadata>(e)->data[bbcd_id])
    {
      optional<bbox_collision_detector_root_node> const* root;
      if (hint_nodes.empty()) {
        root = &accessor->template get<bbox_collision_detector_root_node>(accessor->get(bbcd_id));
      }
      
      const coordinate_type max_width_minus_one = max_in_array_of_unsigned(bbox.size_minus_one());
      // max_width - 1: power-of-two-sized objects easily squeeze into the next smaller category.
      // i.e., exp = log2_rounding_up(max_width)
      const num_bits_type exp = num_bits_in_integer_that_are_not_leading_zeroes(max_width_minus_one);
      const coordinate_type base_box_size = safe_left_shift_one(exp);

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
      metadata.zboxes_union = bounding_box::min_and_size_minus_one(zboxes_union_min, zboxes_union_size_minus_one);
      
      persistent_set<entity_id> old_nodes_to_remove = metadata.nodes;
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
        // while we use the zboxes_union concept, we can't cull this.
        //if (zb.overlaps(bbox)) {
          entity_ref best_hint_ref;
          if (hint_nodes.empty()) {
            if ((*root) && (*root)->root_node_id_) {
              best_hint_ref = accessor->get((*root)->root_node_id_);
            }
            else {
              best_hint_ref = accessor->create_entity();
              accessor->template set<ztree_node>(best_hint_ref, ztree_node(zb, entity_id()));
              accessor->template set<bbox_collision_detector_root_node>(accessor->get(bbcd_id), bbox_collision_detector_root_node(best_hint_ref.id()));
            }
          }
          else {
            // Find the best box to use as hint.
            // TODO: can this be less than O((num zboxes)^2)?
            num_bits_type best_num_low_bits = coordinate_bits + 1;
            for (entity_id node_id : hint_nodes) {
              entity_ref node_ref = accessor->get(node_id);
              auto& node = accessor->template get<ztree_node>(node_ref);
              zbox test = zbox::smallest_joint_parent(zb, node->here);
              if (test.num_low_bits() < best_num_low_bits) {
                best_num_low_bits = test.num_low_bits();
                best_hint_ref = node_ref;
              }
            }
          }
          entity_ref inserted_at_ref = insert_zbox_with_hint(best_hint_ref, zb);
          old_nodes_to_remove.erase(inserted_at_ref.id());
        //}
      }
      erase_from_nodes(accessor, bbcd_id, e, metadata, old_nodes_to_remove);
    }
  private:
    // these functions return where the zbox was inserted
    entity_ref insert_zbox_with_hint(entity_ref hint_node_ref, zbox const& box) {
      while (true) {
        auto& node = accessor_->template get<ztree_node>(hint_node_ref);
        if (node->here.subsumes(box)) {
          break;
        }
        if (node->parent) {
          hint_node_ref = accessor_->get(node->parent);
        }
        else {
          hint_node_ref = add_joint_parent_node(hint_node_ref, box);
        }
      }
      return insert_zbox_downwards(hint_node_ref, box);
    }
    void insert_zbox_at(optional<ztree_node>& node, entity_id node_id) {
      node->objects_here.insert(e.id());
      metadata.nodes.insert(node_id);
    }
    entity_ref insert_zbox_downwards(entity_ref node_ref, zbox const& box) {
      auto& node = accessor_->template get_mut<ztree_node>(node_ref);
      if (node->here.subsumes(box)) {
        if (box.num_low_bits() == node->here.num_low_bits()) {
          insert_zbox_at(node, node_ref.id());
          return node_ref;
        }
        else {
          const size_t which_child = box.get_bit(node->here.num_low_bits() - 1);// ? 1 : 0;
          if (node->children[which_child]) {
            return insert_zbox_downwards(accessor_->get(node->children[which_child]), box);
          }
          else {
            entity_ref new_child = accessor_->create_entity();
            node->children[which_child] = new_child.id();
            insert_zbox_at(accessor_->template set<ztree_node>(new_child, ztree_node(box, node_ref.id())), new_child.id());
            return node_ref;
          }
        }
      }
      else {
        return insert_zbox_downwards(add_joint_parent_node(node_ref, box), box);
      }
    }
    entity_ref add_joint_parent_node(entity_ref node_ref, zbox const& box) {
      auto& node = accessor_->template get_mut<ztree_node>(node_ref);
      entity_ref new_node_ref = accessor_->create_entity();
      auto& new_node = accessor_->template set<ztree_node>(new_node_ref, ztree_node(zbox::smallest_joint_parent(node->here, box), node->parent));
      node->parent = new_node_ref.id();
      if (!new_node->parent) {
        accessor_->template set<bbox_collision_detector_root_node>(accessor_->get(bbcd_id), bbox_collision_detector_root_node(node_ref.id()));
      }
      
      assert_if_ASSERT_EVERYTHING(new_node->here.num_low_bits() > node->here.num_low_bits());
      assert_if_ASSERT_EVERYTHING(new_node->here.subsumes(node->here));
      assert_if_ASSERT_EVERYTHING(new_node->here.subsumes(box));
      assert_if_ASSERT_EVERYTHING(box.subsumes(node->here) || (node->here.get_bit(new_node->here.num_low_bits() - 1) != box.get_bit(new_node->here.num_low_bits() - 1)));

      const size_t which_child = node->here.get_bit(new_node->here.num_low_bits() - 1);
      new_node->children[which_child] = node_ref.id();
      return new_node_ref;
    }
  };
  
  class gather_events {
    // From each of the object's 1-to-2^n nodes, we need to check everything that's
    // an ancestor OR descendant of any of those nodes. The descendants are all different.
    // TODO is there an efficient way to avoid duplicating the search of joint ancestors?
    accessor* accessor_;
    entity_id bbcd_id;
    FuncsType funcs;
    entity_ref e;
    std::unordered_set<entity_id> interaction_possibilities_already_found;
  public:
    gather_events(accessor* accessor, entity_id bbcd_id, FuncsType funcs, entity_ref e)
      :
      accessor_(accessor),
      bbcd_id(bbcd_id),
      funcs(funcs),
      e(e)
    {
      
      bbcd_entry_metadata const& metadata = accessor->template get<spatial_entity_metadata>(e)->data.find(bbcd_id)->second;
      accessor->anticipate_event(funcs.escape_time(accessor, e, metadata.zboxes_union), std::shared_ptr<event>(new spatial_entity_escapes_its_zboxes(bbcd_id, e.id())));
      for (entity_id node_id : metadata.nodes) {
        entity_ref node_ref = accessor->get(node_id);
        auto& node = accessor->template get<ztree_node>(node_ref);
        gather_events_downwards(node);
        if (node->parent) {
          entity_ref node_ref2 = accessor->get(node->parent);
          while (true) {
            auto& node2 = accessor->template get<ztree_node>(node_ref2);
            gather_events_at(node2);
            if (node2->parent) {
              node_ref2 = accessor->get(node2->parent);
            }
            else {
              break;
            }
          }
        }
      }
    }
  private:
    
    void gather_events_at(optional<ztree_node> const& node) {
      for (entity_id other_id : node->objects_here) {
        auto p = interaction_possibilities_already_found.insert(other_id);
        if (p.second) {
          funcs.anticipate_interactions(accessor_, e, accessor_->get(other_id));
        }
      }
    }
    
    void gather_events_downwards(optional<ztree_node> const& node) {
      gather_events_at(node);
      for (entity_id child_id : node->children) {
        if (child_id) {
          gather_events_downwards(accessor_->template get<ztree_node>(accessor_->get(child_id)));
        }
      }
    }
  };
  

  
  
  
  
  
  template<typename OrderingFunctor, typename ComparedType>
  struct reverse_first_ordering : private OrderingFunctor {
    typedef bool result_type; typedef ComparedType first_argument_type; typedef ComparedType second_argument_type;
    reverse_first_ordering() : OrderingFunctor() {}
    reverse_first_ordering(OrderingFunctor const& o) : OrderingFunctor(o) {}

    bool operator()(ComparedType const& a, ComparedType const& b) {
      return static_cast<OrderingFunctor&>(*this)(b.first, a.first);
    }
  };

  template<typename GetCost, typename CostOrdering = std::less<typename GetCost::cost_type> >
  class iterator {
  private:
    typedef typename GetCost::cost_type cost_type;

    typedef std::pair<cost_type, entity_id> queue_value_type_;
    // swap greater/less because pq sorts by greatest and we want least by default (as sorting normally is)
    typedef std::priority_queue<queue_value_type_, std::vector<queue_value_type_>, reverse_first_ordering<CostOrdering, queue_value_type_> > queue_type_;

    typedef queue_value_type_ iteree;
  public:
    typedef iteree value_type;
    typedef value_type reference;
    typedef value_as_ptr<value_type, reference> pointer;
    typedef ptrdiff_t difference_type; //Ha. Ha.
    typedef std::input_iterator_tag iterator_category;

  private:
    struct contents_ : private GetCost {
      contents_(accessor const* accessor, GetCost const& getcost, CostOrdering const& costordering)
        : GetCost(getcost), accessor_(accessor), queue_(costordering), seen_()
        {}

      contents_(contents_&&) = default;
      accessor const* accessor_;
      queue_type_ queue_;
      std::unordered_set<entity_id> seen_;

      GetCost& get_get_cost() { return *this; }

      inline void push_cost(cost_type const& cost, entity_id v) {
        queue_.push(queue_value_type_(cost, v));
      }
      template<typename IndirectCostType>
      inline typename boost::disable_if<boost::is_convertible<IndirectCostType, cost_type> >::type
      push_cost(IndirectCostType const& maybe_cost, entity_id v) {
        if(maybe_cost) {
          queue_.push(queue_value_type_(*maybe_cost, v));
        }
      }

      void add_node(entity_id node_id) {
        if (node_id) {
          push_cost(get_get_cost().min_cost(accessor_->template get<ztree_node>(accessor_->get(node_id))->here.get_bbox()), node_id);
        }
      }
      void add_entity(entity_id id) {
        if(seen_.insert(id).second) {
          push_cost(get_get_cost().cost(id), id);
        }
      }

      void add_children_of(ztree_node const& node) {
        for(entity_id id : node.children) {
          add_node(id);
        }
        for(entity_id id : node.objects_here) {
          add_entity(id);
        }
      }
    };

    //dynamic allocation, pooh.
    boost::shared_ptr<contents_> c_;

    void advance_to_a_returnable_() {
      if(c_) {
        while(true) {
          if(c_->queue_.empty()) {
            c_.reset();
            break;
          }
          auto const& top_as_node = c_->accessor_->template get<ztree_node>(c_->accessor_->get(c_->queue_.top().second));
          if(top_as_node) {
            c_->queue_.pop();
            c_->add_children_of(*top_as_node);
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
    explicit iterator(accessor const* accessor,
                      GetCost const& getcost = GetCost(),
                      CostOrdering const& costordering = CostOrdering())
      : c_(new contents_(accessor, getcost, costordering)) {}
    explicit iterator(accessor const* accessor,
                      entity_id const& initial,
                      GetCost const& getcost = GetCost(),
                      CostOrdering const& costordering = CostOrdering())
      : c_(new contents_(accessor, getcost, costordering)) {
      c_->add_node(initial); advance_to_a_returnable_();
    }

    reference operator*() const {
      caller_error_if(c_->queue_.empty(), "can't dereference an empty iterator");
      return c_->queue_.top();
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
        accessor const* accessor,
        entity_id node_id,
        std::unordered_set<entity_id>& results,
        GetCostBool& getcost) {
    if (node_id) {
      auto const& node = accessor->template get<ztree_node>(accessor->get(node_id));
      if(bool_with_error_if_implicit_conversions_to_bool_and_cost_type_are_ambiguous<GetCostBool>(
            getcost.min_cost(node->here.get_bbox()))) {
        for(entity_id id : node->objects_here) {
          if(bool_with_error_if_implicit_conversions_to_bool_and_cost_type_are_ambiguous<GetCostBool>(
                getcost.cost(id))) {
            results.insert(id);
          }
        }
        for(entity_id id : node->children) {
          filter_impl(accessor, id, results, getcost);
        }
      }
    }
  }

public:
  template<typename GetCost>
  static inline boost::iterator_range<iterator<GetCost>> iterate(accessor const* accessor, entity_id bbcd_id, GetCost const& getcost) {
    typedef iterator<GetCost> iter;
    boost::iterator_range<iter> result(iter(accessor, accessor->template get<bbox_collision_detector_root_node>(accessor->get(bbcd_id))->root_node_id_, getcost), iter());
    return result;
  }
  
  template<typename GetCost>
  static inline boost::optional<typename iterator<GetCost>::value_type> find_least(accessor const* accessor, entity_id bbcd_id, GetCost const& getcost) {
    typedef iterator<GetCost> iter;
    typedef boost::optional<typename iterator<GetCost>::value_type> result_type;
    iter i(accessor, accessor->template get<bbox_collision_detector_root_node>(accessor->get(bbcd_id))->root_node_id_, getcost);
    return i ? result_type(*i) : result_type();
  }

  template<typename GetCostBool>
  static inline std::unordered_set<entity_id> filter(accessor const* accessor, entity_id bbcd_id, GetCostBool getcost) {
    std::unordered_set<entity_id> results;
    filter_impl(accessor, accessor->template get<bbox_collision_detector_root_node>(accessor->get(bbcd_id))->root_node_id_, results, getcost);
    return results;
  }
}; // struct operations
}; // struct bbox_collision_detector_system



