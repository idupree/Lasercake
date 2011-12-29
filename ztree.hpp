
#ifndef ELISAACWATER_ZTREE_HPP__
#define ELISAACWATER_ZTREE_HPP__

#include <boost/integer.hpp>
#include <unordered_map>
#include <unordered_set>

typedef int num_bits_type;
typedef size_t num_coordinates_type;

// ObjectIdentifier needs hash and == and to be freely copiable. So, ints will do, pointers will do...
// coordinate_bits should usually be 32 or 64.
template<typename ObjectIdentifier, num_bits_type coordinate_bits, num_coordinates_type num_dimensions>
class space_with_fast_lookup_of_everything_overlapping_localized_area {
  static_assert(num_dimensions >= 0, "You can't make a space with negative dimensions!");
  static_assert(coordinate_bits >= 0, "You can't have an int type with negative bits!");

public:
  typedef typename boost::uint_t<coordinate_bits>::fast Coordinate;
  struct bounding_box {
    std::array<Coordinate, num_dimensions> min, size;
    bool overlaps(bounding_box const& other)const {
      for (num_coordinates_type i = 0; i < num_dimensions; ++i) {
        if (other.min[i] + other.size[i] <       min[i]) return false;
        if (      min[i] +       size[i] < other.min[i]) return false;
      }
      return true;
    }
  };
  
  space_with_fast_lookup_of_everything_overlapping_localized_area():objects_tree(nullptr){}
  
private:
  static const num_bits_type total_bits = coordinate_bits * num_dimensions;
  
  struct zbox {
    // We ensure that every bit except the ones specifically supposed to be on is off.
    std::array<Coordinate, num_dimensions> coords;
    std::array<Coordinate, num_dimensions> interleaved_bits;
    num_bits_type num_low_bits_ignored;
    
    zbox():num_low_bits_ignored(total_bits){ for (num_coordinates_type i = 0; i < num_dimensions; ++i) interleaved_bits[i] = 0; }
    
    bool subsumes(zbox const& other)const {
      if (other.num_low_bits_ignored > num_low_bits_ignored) return false;
      for (num_coordinates_type i = num_low_bits_ignored / coordinate_bits; i < num_dimensions; ++i) {
        Coordinate mask = ~0;
        if (i == num_low_bits_ignored / coordinate_bits) {
          mask &= ~((num_bits_type(1) << (num_low_bits_ignored % coordinate_bits)) - 1);
        }
        if ((interleaved_bits[i] & mask) != (other.interleaved_bits[i] & mask)) return false;
      }
      return true;
    }
    bool get_bit(num_bits_type bit)const {
      return interleaved_bits[bit / coordinate_bits] & ~(Coordinate(1) << (bit % coordinate_bits));
    }
    num_bits_type num_bits_ignored_by_dimension(num_coordinates_type dim)const {
      return (num_low_bits_ignored + (num_dimensions - 1) - dim) / num_dimensions;
    }
    bounding_box get_bbox()const {
      bounding_box result;
      result.min = coords;
      for (num_coordinates_type i = num_low_bits_ignored / coordinate_bits; i < num_dimensions; ++i) {
        result.size[i] = (1 << num_bits_ignored_by_dimension(i));
      }
    }
  };
  
  num_bits_type idx_of_highest_bit(Coordinate i) {
    int upper_bound = coordinate_bits;
    int lower_bound = -1;
    while(true) {
      int halfway_bit_idx = (upper_bound + lower_bound) / 2;
      if (halfway_bit_idx == lower_bound) return lower_bound;
      
      if (i & ~((Coordinate(1) << halfway_bit_idx) - 1)) lower_bound = halfway_bit_idx;
      else                                               upper_bound = halfway_bit_idx;
    }
  }
  
  struct ztree_node {
    zbox here;
    ztree_node *child0;
    ztree_node *child1;
    unordered_set<ObjectIdentifier> objects_here;
  };
  
  zbox smallest_joint_parent(zbox zb1, zbox zb2) {
    zbox new_box;
    for (num_coordinates_type i = num_dimensions; i >= 0; --i) {
      const num_bits_type highest_bit_idx = idx_of_highest_bit(zb1.interleaved_bits[i] ^ zb2.interleaved_bits[i]);
      assert((zb1.interleaved_bits[i] & ~((1 << (highest_bit_idx + 1)) - 1)) == (zb2.interleaved_bits[i] & ~((1 << (highest_bit_idx + 1)) - 1)));
      new_box.interleaved_bits[i] = zb1.interleaved_bits[i] & ~((1 << (highest_bit_idx + 1)) - 1);
      if (highest_bit_idx > 0) {
        new_box.num_low_bits_ignored = highest_bit_idx + i * coordinate_bits;
        for (num_coordinates_type j = 0; j < num_dimensions; ++j) {
          assert(             zb1.coords[j] & ~((Coordinate(1) << new_box.num_bits_ignored_by_dimension(j)) - 1)
                           == zb2.coords[j] & ~((Coordinate(1) << new_box.num_bits_ignored_by_dimension(j)) - 1));
          new_box.coords[j] = zb1.coords[j] & ~((Coordinate(1) << new_box.num_bits_ignored_by_dimension(j)) - 1);
        }
        return new_box;
      }
    }
    new_box.num_low_bits_ignored = 0;
    assert(zb1.coords == zb2.coords);
    new_box.coords = zb1.coords;
    return new_box;
  }
  
  zbox box_from_coords(std::array<Coordinate, num_dimensions> const& coords, num_bits_type num_low_bits_ignored) {
    zbox result;
    result.coords = coords;
    result.num_low_bits_ignored = num_low_bits_ignored;
    for (num_bits_type bit_within_interleaved_bits = num_low_bits_ignored;
                       bit_within_interleaved_bits < total_bits;
                     ++bit_within_interleaved_bits) {
      const num_bits_type bit_idx_within_coordinates = bit_within_interleaved_bits / num_dimensions;
      const num_coordinates_type which_coordinate    = bit_within_interleaved_bits % num_dimensions;
      const num_bits_type interleaved_bit_array_idx  = bit_within_interleaved_bits / coordinate_bits;
      const num_bits_type interleaved_bit_local_idx  = bit_within_interleaved_bits % coordinate_bits;
      result.interleaved_bits[interleaved_bit_array_idx] |= ((coords[which_coordinate] >> bit_idx_within_coordinates) & 1) << interleaved_bit_local_idx;
    }
    return result;
  }
  
  void insert_box(ztree_node*& tree, ObjectIdentifier obj, zbox box) {
    if (!tree) {
      tree = new ztree_node {
        box,
        nullptr,
        nullptr,
        unordered_set<ObjectIdentifier>()
      };
      tree->objects_here.insert(obj);
    }
    else {
      if (tree->here.subsumes(box)) {
        if (box.num_low_bits_ignored == tree->here.num_low_bits_ignored) {
          tree->objects_here.insert(obj);
        }
        else {
          if (box.get_bit(tree->here.num_low_bits_ignored - 1)) insert_box(tree->child1, obj, box);
          else                                                  insert_box(tree->child0, obj, box);
        }
      }
      else {
        const ztree_node *new_tree = new ztree_node {
          smallest_joint_parent(tree->here, box),
          nullptr,
          nullptr,
          unordered_set<ObjectIdentifier>()
        };
        if (tree->here.get_bit(box.num_low_bits_ignored - 1)) new_tree->child1 = tree;
        else                                                  new_tree->child0 = tree;
      
        tree = new_tree;
        insert_box(tree, obj, box);
      }
    }
  }
  
  void delete_object(ztree_node*& tree, ObjectIdentifier obj, bounding_box const& bbox) {
    if (!tree) return;
    if (tree->here.get_bbox().overlaps(bbox)) {
      tree->objects_here.erase(obj);
      delete_object(tree->child0, obj, bbox);
      delete_object(tree->child1, obj, bbox);
      
      if (tree->objects_here.empty()) {
        if (tree->child0) {
          if (!tree->child1) {
            tree = tree->child0;
          }
        }
        else tree = tree->child1; // which could be null
      }
    }
  }
  
  void zget_objects_overlapping(ztree_node* tree, unordered_set<ObjectIdentifier>& results, bounding_box const& bbox) {
    if (tree && tree->here.get_bbox().overlaps(bbox)) {
      for (const ObjectIdentifier obj : tree->objects_here) {
        auto bbox_iter = bboxes_by_object.find(obj);
        assert(bbox_iter != bboxes_by_object.end());
        if (bbox_iter->overlaps(bbox)) results.insert(obj);
      }
      zget_objects_overlapping(tree->child0, results, bbox);
      zget_objects_overlapping(tree->child1, results, bbox);
    }
  }
  
  unordered_map<ObjectIdentifier, bounding_box> bboxes_by_object;
  ztree_node* objects_tree;
  
public:

  void insert(ObjectIdentifier id, bounding_box const& bbox) {
    bboxes_by_object.insert(make_pair(id, bbox));
    Coordinate max_dim = bbox[0];
    for (num_coordinates_type i = 1; i < num_dimensions; ++i) {
      if (bbox[i] > max_dim) max_dim = bbox.size[i];
    }
    int exp = 0; while ((Coordinate(1) << exp) < max_dim) ++exp;
    for (int i = 0; i < (1 << num_dimensions); ++i) {
      std::array<Coordinate, num_dimensions> coords = bbox.min;
      for (num_coordinates_type j = 0; j < num_dimensions; ++j) {
        if (i & (1<<j)) coords[j] += (Coordinate(1) << exp);
      }
      zbox zb = box_from_coords(coords, exp * num_dimensions);
      if (zb.get_bbox().overlaps(bbox))  // don't add absurdly unnecessary boxes...
        insert_box(objects_tree, id, zb);
    }
  }
  void erase(ObjectIdentifier id) {
    auto bbox_iter = bboxes_by_object.find(id);
    assert(bbox_iter != bboxes_by_object.end());
    delete_object(objects_tree, id, (*bbox_iter));
    bboxes_by_object.erase(bbox_iter);
  }
  
  void get_objects_overlapping(unordered_set<ObjectIdentifier>& results, bounding_box const& bbox)const {
    zget_objects_overlapping(objects_tree, results, bbox);
  }
};

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
  Z1/Z2, with B ignored bits
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

#endif
