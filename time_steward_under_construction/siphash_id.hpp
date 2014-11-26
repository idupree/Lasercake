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

#ifndef LASERCAKE_SIPHASH_ID_HPP__
#define LASERCAKE_SIPHASH_ID_HPP__

#include <boost/integer.hpp>

uint64_t siphash24(const void *src,
                   unsigned long src_sz,
                   const char key[16]);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
// 128 bits of random stuff (presumed unique to this library)
const char siphash_key[16] = { 0xb7, 0xac, 0x3d, 0xf8, 0xc3, 0xa2, 0x8c, 0xd9, 0x3a, 0x10, 0x91, 0x68, 0x09, 0x02, 0x74, 0x0e };
#pragma GCC diagnostic pop

// siphash_id takes anything that's hashable and turns it into a
// fixed-size, statistically almost certainly unique, equality-comparable
// and orderable data structure (i.e. 128 bits of data).  The ==
// is true iff the things are the same; the < is an arbitrary total ordering.
class siphash_id {
  // Implementation note:
  // A reply from the SipHash developers about our idea to get a
  // large enough result from SipHash to be statistically unique
  // by hashing each of N and N+1 and keeping both results:
  //
  // Hi Eli,
  //
  // the development of a 128-bit SipHash is under discussion, see e.g.
  // this Twitter thread:
  //https://twitter.com/Kryptoblog/status/397285808714960897
  //
  // The trick you use seems to work, but it'd perhaps be more elegant to
  // just hash the same info under two distinct keys, or (say) under key K
  // and key (K xor 0xffff..ffff). That's still twice as slow as a single
  // pass though.
  //
  // What you should not do obviously is "hash the hash" to get another 64
  // bits, as any collision would propagate.
  //
  // When/if a SipHash-128 is developed, we'll advertize it on Twitter.
  //
  // Best,
  //
  // JP
  
public:
  typedef std::array<uint64_t,2> data_type;
  
  // Accepts any combination of siphash_ids and types castable to uint64_t.
  // Combines all the entropy from all the arguments.
  template<typename... Args>
  static siphash_id combining(Args... args);
  
  // null siphash
  constexpr siphash_id() : data_({{ 0, 0 }}) {}
  constexpr static siphash_id null() {
    return siphash_id();
  }
  // least and greatest nonnull siphash
  constexpr static siphash_id least() {
    return siphash_id(data_type{{ 1, 0 }});
  }
  constexpr static siphash_id greatest() {
    return siphash_id(data_type{{ uint64_t(0)-1, uint64_t(0)-1 }});
  }

  bool operator<(siphash_id const& o)const {
    if (data_[1] < o.data_[1]) return true;
    if (data_[1] > o.data_[1]) return false;
    return data_[0] < o.data_[0];
  }
  bool operator>(siphash_id const& o)const { return o < *this; }
  bool operator<=(siphash_id const& o)const { return !(o < *this); }
  bool operator>=(siphash_id const& o)const { return !(o > *this); }
  bool operator==(siphash_id const& o)const {
    return (data_[0] == o.data_[0]) && (data_[1] == o.data_[1]);
  }
  bool operator!=(siphash_id const& o)const {
    return (data_[0] != o.data_[0]) || (data_[1] != o.data_[1]);
  }
  explicit operator bool()const {
    return (*this) != null();
  }
  
  friend inline std::ostream& operator<<(std::ostream& os, siphash_id const& s) {
    os << s.data_[0] << "," << s.data_[1];
    return os;
  }
  data_type const& data()const { return data_; }
private:
  constexpr siphash_id(data_type data):data_(data){}
  data_type data_;
  friend struct std::hash<siphash_id>;
};
namespace impl {
  template<typename... Args> struct combining_helper;
  template<> struct combining_helper<> {
    static const size_t next_idx = 0;
    static inline void enter_arg(char[]){}
  };
  template<typename... Tail> struct combining_helper<siphash_id, Tail...> {
    static const size_t next_idx = combining_helper<Tail...>::next_idx + 2*sizeof(uint64_t);
    static inline void enter_arg(char in[], siphash_id head, Tail... tail) {
      combining_helper<Tail...>::enter_arg(in, tail...);
      *(uint64_t*)(&in[combining_helper<Tail...>::next_idx                   ]) = head.data()[0];
      *(uint64_t*)(&in[combining_helper<Tail...>::next_idx + sizeof(uint64_t)]) = head.data()[1];
    }
  };
  template<typename Head, typename... Tail> struct combining_helper<Head, Tail...> {
    static_assert(std::is_integral<Head>::value, "siphash_id::combining only accepts siphash_ids and integral types");
    static const size_t next_idx = combining_helper<Tail...>::next_idx + sizeof(Head);
    static inline void enter_arg(char in[], Head head, Tail... tail) {
      combining_helper<Tail...>::enter_arg(in, tail...);
      *(Head*)(&in[combining_helper<Tail...>::next_idx]) = head;
    }
  };
  template<typename Head, typename... Tail> struct combining_helper<Head&, Tail...> : public combining_helper<Head, Tail...> {};
}
template<typename... Args>
siphash_id siphash_id::combining(Args... args) {
  siphash_id result;
  char in[impl::combining_helper<Args...>::next_idx];
  impl::combining_helper<Args...>::enter_arg(in, args...);
  result.data_[0] = siphash24((void*)in, sizeof(in), siphash_key);
  ++in[0];
  result.data_[1] = siphash24((void*)in, sizeof(in), siphash_key);
  return result;
}

namespace std {
  template<>
  struct hash<siphash_id> {
    public:
    size_t operator()(siphash_id const& i)const {
      // Since the siphash is already a hash, just return some of it.
      return i.data()[0];
    }
  };
}
  

namespace persistent_siphash_id_trie_system {
typedef uint32_t num_bits_type;
template<typename KeyType, typename MappedType, bool IsMap> struct value_typer;
template<typename KeyType, typename MappedType> struct value_typer<KeyType, MappedType, false> {
  typedef KeyType type;
  static KeyType get_key(type const& t) { return t; }
};
template<typename KeyType, typename MappedType> struct value_typer<KeyType, MappedType, true> {
  typedef std::pair<KeyType, MappedType> type;
  static KeyType get_key(type const& t) { return t.first; }
};
template<typename MappedType = void, num_bits_type bits_per_level = 4>
class persistent_siphash_id_trie {
public:
  static_assert((128 % bits_per_level) == 0, "I don't think everything works right if 128 isn't a multiple of bits_per_level");
  typedef siphash_id key_type;
  typedef MappedType mapped_type;
  typedef value_typer<key_type, mapped_type, !std::is_void<mapped_type>::value> value_typer_t;
  typedef typename value_typer_t::type value_type;
private:  
  typedef typename boost::uint_t<(1<<bits_per_level)>::fast one_bit_per_child_t;
  typedef size_t ref_count_t;
  static const ref_count_t is_leaf_bit = 1ULL<<63;
  static const ref_count_t ref_count_mask = is_leaf_bit-1;
  struct child_header {
    child_header(ref_count_t r):ref_count_and_is_leaf(r){}
    ref_count_t ref_count_and_is_leaf;
  };
  struct leaf_t : public child_header {
    template <class... Args>
    leaf_t(Args&&... args):child_header(is_leaf_bit),value(std::forward<Args>(args)...){}
    value_type value;
    key_type id() { return value_typer_t::get_key(value); }
  };
  struct node_header : public child_header {
    node_header():child_header(0),children_exist(0){}
    one_bit_per_child_t children_exist;
  };
  class child_ptr {
  public:
    child_ptr(){}
    child_ptr(decltype(nullptr)):data(nullptr){}
    explicit child_ptr(char* data):data(data){}
    inline bool is_leaf()const { return (((child_header*)(data))->ref_count_and_is_leaf) & is_leaf_bit; }
    leaf_t& leaf() { return (*(leaf_t*)(data)); }
    leaf_t const& leaf()const { return (*(leaf_t const*)(data)); }
    node_header& node() { return (*(node_header*)(data)); }
    node_header const& node()const { return (*(node_header const*)(data)); }
    child_ptr& child_by_idx(num_bits_type i)const { return (*(child_ptr*)(data + sizeof(node_header) + i*sizeof(child_ptr))); }
    child_ptr& child_by_bit(num_bits_type b)const { return child_by_idx(popcount(node().children_exist & (b-1))); }
    num_bits_type num_children()const { return popcount(node().children_exist); }
    explicit operator bool()const { return bool(data); }
    bool operator==(child_ptr o)const { return data==o.data; }
    bool operator!=(child_ptr o)const { return data!=o.data; }
    void inc_ref()const {
      if (data) {
        ++(((child_header*)(data))->ref_count_and_is_leaf);
      }
    }
    void dec_ref()const {
      if (data && ((--(((child_header*)(data))->ref_count_and_is_leaf)) & ref_count_mask) == 0) {
        if (is_leaf()) {
          leaf().~leaf_t();
        }
        else {
          const num_bits_type num_children_ = num_children();
          for (num_bits_type i = 0; i < num_children_; ++i) {
            child_by_idx(i).dec_ref();
            child_by_idx(i).~child_ptr();
          }
          node().~node_header();
        }
        delete[] data;
      }
    }
  private:
    char* data;
  };
  template <class... Args>
  static child_ptr allocate_leaf(Args&&... args) {
    char* ptr = new char[sizeof(leaf_t)];
    new(ptr) leaf_t(std::forward<Args>(args)...);
    return child_ptr(ptr);
  }
  static child_ptr allocate_node(num_bits_type num_children) {
    char* ptr = new char[sizeof(node_header) + num_children*sizeof(child_ptr)];
    new(ptr) node_header();
    for (num_bits_type i = 0; i < num_children; ++i) {
      new(ptr + sizeof(node_header) + i*sizeof(child_ptr)) child_ptr();
    }
    return child_ptr(ptr);
  }
  child_ptr root;
  
  static inline num_bits_type which_child(siphash_id id, num_bits_type which_bits) {
    return (id.data()[which_bits>=64]>>which_bits) & ((1<<bits_per_level)-1);
  }
  
  persistent_siphash_id_trie(child_ptr root):root(root){ root.inc_ref(); }
public:
  persistent_siphash_id_trie():root(nullptr){}
  // Not set equality. Notably, though, trie.erase (something not in the trie) == trie.
  bool operator==(persistent_siphash_id_trie o)const { return root==o.root; }
  bool operator!=(persistent_siphash_id_trie o)const { return root!=o.root; }
  bool empty()const { return !bool(root); }
  
  persistent_siphash_id_trie(persistent_siphash_id_trie const& o):root(o.root){ root.inc_ref(); }
  persistent_siphash_id_trie& operator=(persistent_siphash_id_trie const& o) {
    o.root.inc_ref();
    root.dec_ref();
    root = o.root;
    return *this;
  }
  ~persistent_siphash_id_trie() { root.dec_ref(); }
    
  class iterator : public boost::iterator_facade<iterator, const value_type, boost::bidirectional_traversal_tag> {
    public:
      //iterator() : path_len(0) {}
      explicit operator bool()const { return id && back(); }
    private:
      explicit iterator(child_ptr root):path_len(1) { path[0] = root; }
      friend class boost::iterator_core_access;
      friend class persistent_siphash_id_trie;
      
      // State: We can be
      // 1) an iterator to an element (path is the path ending in the leaf, id redundant (== this->first or *this))
      // 2) an iterator to an empty child (path is the path ending in nullptr, possibly with a leaf before it, id is needed to specify which child)
      // 3) the end iterator (path is the singleton root, id null)
      // --end is the last element, ++end is the first element
      // The empty set has one iterator which is its own next and prior.
      void increment() {
        if (!path[0]) return;
        if (id) {
          num_bits_type child = which_child(id, back_bits());
          
          while (!back() || back().is_leaf() || !(back().node().children_exist & ~((1<<(child+1))-1))) {
            if (back() && !back().is_leaf()) { assert(popcount(back().node().children_exist & ((1U<<child)-1)) + 1U == back().num_children()); }
            bool check = bool(back());
            if (path_len == 1) {
              id = siphash_id::null();
              return;
            }
            --path_len;
            child = which_child(id, back_bits());
            if (check) { assert(back().node().children_exist & (1<<child)); }
          }
          
          num_bits_type idx = popcount(back().node().children_exist & ((1U<<child)-1)) + 1;
          assert (idx < back().num_children());
          assert (idx >= 0);
          path[path_len] = back().child_by_idx(idx);
          ++path_len;
        }
        
        while (!back().is_leaf()) {
          path[path_len] = back().child_by_idx(0);
          ++path_len;
        }
        id = back().leaf().id();
      }
      void decrement() {
        if (!path[0]) return;
        if (id) {
          num_bits_type child = which_child(id, back_bits());
          
          while (!back() || back().is_leaf() || !(back().node().children_exist & ((1<<child)-1))) {
            if (back() && !back().is_leaf()) { assert(popcount(back().node().children_exist & ((1U<<child)-1)) == 0); }
            bool check = bool(back());
            if (path_len == 1) {
              id = siphash_id::null();
              return;
            }
            --path_len;
            child = which_child(id, back_bits());
            if (check) { assert(back().node().children_exist & (1<<child)); }
          }
          
          num_bits_type idx = popcount(back().node().children_exist & ((1U<<child)-1)) + 1;
          assert (idx < back().num_children());
          assert (idx >= 0);
          path[path_len] = back().child_by_idx();
          ++path_len;
        }
        while (!back().is_leaf()) {
          path[path_len] = back().child_by_idx(back().num_children()-1);
          ++path_len;
        }
        id = back().leaf().id();
      }
      bool equal(iterator const& other)const { return id == other.id; }
      value_type const& dereference()const { return back().leaf().value; }
      
      child_ptr& back() { return path[path_len-1]; }
      child_ptr const& back()const { return path[path_len-1]; }
      num_bits_type bits(num_bits_type path_idx)const { return 128 - (path_idx+1)*bits_per_level; }
      num_bits_type back_bits()const { return 128 - path_len*bits_per_level; }
      
      siphash_id id;
      num_bits_type path_len;
      std::array<child_ptr, 128/bits_per_level> path;
  };
  iterator begin()const { return ++iterator(root); }
  iterator end()const { return iterator(root); }
  
  template<typename Hack = void>
  /*std::enable_if_t<std::is_same<key_type, value_type>::value, persistent_siphash_id_trie>*/ persistent_siphash_id_trie insert(key_type k)const {
    return insert(find(k));
  }
  template<typename Hack = void>
  static /*std::enable_if_t<std::is_same<key_type, value_type>::value, persistent_siphash_id_trie>*/ persistent_siphash_id_trie insert(iterator const& i) {
    persistent_siphash_id_trie result = set_leaf(i, allocate_leaf(i.id));
    //assert(result.find(i.id).back().leaf().id() == i.id);
    return result;
  }
  template<class... Args>
  persistent_siphash_id_trie emplace(key_type k, Args&&... args)const {
    return emplace(find(k), std::forward<Args>(args)...);
  }
  template<class... Args>
  static persistent_siphash_id_trie emplace(iterator const& i, Args&&... args) {
    persistent_siphash_id_trie result = set_leaf(i, allocate_leaf(i.id, std::forward<Args>(args)...));
    //assert(result.find(i.id).back().leaf().id() == i.id);
    return result;
  }
  persistent_siphash_id_trie erase(key_type k)const {
    return erase(find(k));
  }
  static persistent_siphash_id_trie erase(iterator const& i) {
//     std::vector<child_ptr> old_leaves;
//     for (iterator j = persistent_siphash_id_trie(i.path[0]).begin(); j != persistent_siphash_id_trie(i.path[0]).end(); ++j) {
//       if (j.id != i.id) old_leaves.push_back(j.back());
//     }
    persistent_siphash_id_trie result = unset_leaf(i);
    //assert(!result.find(i.id));
//     size_t idx = 0;
//     for (iterator j = result.begin(); j != result.end(); ++j) {
//       assert(old_leaves[idx++] == j.back());
//     }
    return result;
  }
  iterator find(key_type k)const {
    caller_correct_if(bool(k), "You can't find() for the null ID");
    iterator i(root);
    i.id = k;
    while (i.back()) {
      if (i.back().is_leaf()) {
        if (i.back().leaf().id() == k) {
          return i;
        }
        else {
          i.path[i.path_len++] = nullptr;
        }
      }
      else {
        const num_bits_type child_bit = 1<<which_child(i.id, i.back_bits());
        if (i.back().node().children_exist & child_bit) {
          i.path[i.path_len] = i.back().child_by_bit(child_bit);
          ++i.path_len;
        }
        else {
          i.path[i.path_len++] = nullptr;
        }
      }
    }
    return i;
  }
  iterator lower_bound(key_type k)const {
    iterator i = find(k);
    if (i) return i;
    return ++i;
  }
  iterator upper_bound(key_type k)const {
    return ++find(k);
  }
private:
  static persistent_siphash_id_trie unset_leaf(iterator const& i) {
    if (i.back()) {
      assert(i.back().is_leaf());
      return walk_change_upwards(i, i.path_len - 1, nullptr);
    }
    else {
      return i.path[0];
    }
  }
  static persistent_siphash_id_trie set_leaf(iterator const& i, child_ptr new_leaf) {
    child_ptr walker;
    num_bits_type replacing_idx = i.path_len - 1;
    if (i.back()) {
      assert(i.back().is_leaf());
      walker = new_leaf;
    }
    else if (i.path_len < 2) {
      return new_leaf;
    }
    else {
      child_ptr parent = i.path[i.path_len - 2];
      if (parent.is_leaf()) {
        --replacing_idx;
        num_bits_type which_bits = i.bits(replacing_idx);
        const siphash_id k1 = parent.leaf().id();
        num_bits_type child0 = which_child(i.id, which_bits);
        num_bits_type child1 = which_child(k1, which_bits);
        child_ptr split_at;
        if (child0 == child1) {
          split_at = walker = allocate_node(1);
          while (true) {
            num_bits_type old_child = child0;
            child_ptr old_split_at = split_at;
            which_bits = which_bits - bits_per_level;
            child0 = which_child(i.id, which_bits);
            child1 = which_child(k1, which_bits);
            old_split_at.node().children_exist = 1<<old_child;
            if (child0 == child1) {
              old_split_at.child_by_idx(0) = split_at = allocate_node(1);
              split_at.inc_ref();
              assert(old_split_at.num_children() == 1);
              assert(old_split_at.child_by_idx(0));
            }
            else {
              old_split_at.child_by_idx(0) = split_at = allocate_node(2);
              split_at.inc_ref();
              assert(old_split_at.num_children() == 1);
              assert(old_split_at.child_by_idx(0));
              break;
            }
          }
        }
        else {
          split_at = walker = allocate_node(2);
        }
        split_at.node().children_exist = (1<<child0) | (1<<child1);
        split_at.child_by_idx(child1 > child0) = parent;
        split_at.child_by_idx(child0 > child1) = new_leaf;
        assert(split_at.num_children() == 2);
        assert(split_at.child_by_idx(0));
        assert(split_at.child_by_idx(1));
        parent.inc_ref();
        new_leaf.inc_ref();
      }
      else {
        walker = new_leaf;
      }
    }
    return walk_change_upwards(i, replacing_idx, walker);
  }
    
  static persistent_siphash_id_trie walk_change_upwards(iterator const& i, num_bits_type replacing_idx, child_ptr walker) {
    child_ptr old_child = i.path[replacing_idx];
    --replacing_idx;
    for (; replacing_idx != num_bits_type(-1); --replacing_idx) {
      child_ptr old_parent = i.path[replacing_idx];
      num_bits_type which_bits = i.bits(replacing_idx);
      assert(!old_parent.is_leaf());
      
      const num_bits_type child_idx = which_child(i.id, which_bits);
      assert (walker != old_child);
      const num_bits_type new_num_children = old_parent.num_children() + (walker && !old_child) - (old_child && !walker);
      if ((new_num_children > 1) || (walker && !walker.is_leaf())) {
        child_ptr walker_parent = allocate_node(new_num_children);
        num_bits_type old_idx = 0;
        num_bits_type new_idx = 0;
        for (num_bits_type i = 0; i < (1<<bits_per_level); ++i) {
          if (i == child_idx) {
            if (old_child) { ++old_idx; }
            if (walker) {
              walker_parent.node().children_exist |= 1 << i;
              walker_parent.child_by_idx(new_idx++) = walker;
              walker.inc_ref();
            }
          }
          else if (old_parent.node().children_exist & (1<<i)) {
            walker_parent.node().children_exist |= 1 << i;
            child_ptr old_child_ptr = old_parent.child_by_idx(old_idx++);
            walker_parent.child_by_idx(new_idx++) = old_child_ptr;
            old_child_ptr.inc_ref();
          }
        }
        assert(new_idx == new_num_children);
        assert(walker_parent.num_children() == new_num_children);
        walker = walker_parent;
      }
      else if ((new_num_children == 1) && (!walker)) {
        assert(old_parent.num_children() == 2);
        num_bits_type other_child_bit = old_parent.node().children_exist & ~(1<<child_idx);
        child_ptr other_child = old_parent.child_by_bit(other_child_bit); // TODO optimize
        if (other_child.is_leaf()) {
          walker = other_child;
        }
        else {
          walker = allocate_node(1);
          walker.node().children_exist = other_child_bit;
          walker.child_by_idx(0) = other_child;
          other_child.inc_ref();
        }
      }
      
      old_child = old_parent;
    }
    
    return walker;
  }
};
} // namespace persistent_siphash_id_trie_system

typedef persistent_siphash_id_trie_system::persistent_siphash_id_trie<> persistent_siphash_id_set;
template<typename T>
using persistent_siphash_id_map = persistent_siphash_id_trie_system::persistent_siphash_id_trie<T>;

#endif
