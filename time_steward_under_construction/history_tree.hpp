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

#ifndef LASERCAKE_HISTORY_TREE_HPP__
#define LASERCAKE_HISTORY_TREE_HPP__

template<class TimeSteward, typename TimeSteward::time_type StandardTimeIncrement>
class history_tree {
public:
  typedef typename TimeSteward::time_type time_type;
  struct fiat_event {
    fiat_event(uint64_t distinguisher, std::shared_ptr<event> e):distinguisher(distinguisher),e(e){}
    uint64_t distinguisher;
    std::shared_ptr<event> e;
  };
  struct node {
    //node* parent;
    node(time_type start_time):start_time(start_time),end_time(start_time){}
    time_type start_time;
    time_type end_time;
    std::multimap<time_type, fiat_event> events;
    std::multimap<time_type, node> children;
  };
  TimeSteward steward;
  node history_root;
  typedef std::vector<node*> history;
  history current_history;
  time_type current_time;
  
  struct spatial_representation_entry {
    spatial_representation_entry(double height, history h):height(height),h(h){}
    double height;
    history h;
  };
  struct spatial_representation_column {
    std::vector<spatial_representation_entry> entries;
    std::vector<double> desired_heights(spatial_representation_column const& prev, time_type time_diff) {
      std::vector<double> results;
      double inc = double(1) / double(current.size() + 1);
      double current_height_asymptote = 0;
      for (auto e : entries) {
        bool in_prev = false;
        double prev_height;
        for (auto f : prev.entries) {
          if (f.h.back() == e.h.back()) {
            in_prev = true;
            prev_height = f.height;
            current_height_asymptote += inc;
            break;
          }
        }
        if (in_prev) {
          // decay towards evenly-spaced
          results.push_back(current_height_asymptote - (current_height_asymptote-prev_height)*std::pow(0.9, double(time_diff)/double(StandardTimeIncrement)));
        }
        else {
          results.push_back(results.back());
        }
      }
      return results;
    }
  };
  std::map<time_type, spatial_representation_column> spatial_representation_columns;
  
  history_tree():
    history_root(-1),
    current_time(-1),
  {
    current_history.push_back(&history_root);
    spatial_representation_column c;
    c.entries.push_back(spatial_representation_entry(0.5, current_history));
    spatial_representation_columns.insert(std::make_pair(time_type(0), c));
    validate();
  }
  
  void update_spatial_representation_heights(time_type start, time_type end) {
    auto i = spatial_representation_columns.upper_bound(start);
    assert (i != spatial_representation_columns.begin());
    auto prev = boost::prior(i);
    for (; i != spatial_representation_columns.end(); ++i) {
      // decay towards desired_heights
      double decay_factor = 0;
      if (i->first > start) { decay_factor = std::pow(0.9, double(end - start)/double(i->first - start)); }
      if (decay_factor > 0.9999) break;
      std::vector<double> desired_heights = i->second.desired_heights(prev->second, i->first - prev->first);
      for (size_t j = 0; j < desired_heights.size(); ++j) {
        i->second.entries[j].height = i->second.entries[j].height*decay_factor - desired_heights[j]*(1-decay_factor);
      }
      prev = i;
    }
  }
  
  double time_coord(time_type time, time_type focus_time) {
    double double_max_time = boost::prior(spatial_representation_columns.end())->first;
    double focus_time_coord = 0.1 + 0.8*(double(focus_time) / double_max_time);
    time_type flat_dist = 2*StandardTimeIncrement;
    double flat_slope = 0.1/double(flat_dist);
    if (focus_time-flat_dist < time && time < focus_time+flat_dist) {
      return focus_time_coord + (time-focus_time)*flat_slope;
    }
    else {
      double s = sign(time-focus_time);
      double far_end_time = (s>0) ? (double_max_time+1.0) : -1.0;
      double far_end_coord = (s>0) ? 1 : 0;
      double middle_end_time = (focus_time + s*flat_dist);
      double middle_end_coord = (focus_time_coord + s*0.1);
      double far_slope = (middle_end_coord-far_end_coord)/(middle_end_time-far_end_time);
      double frac = (double(time)-far_end_time) / double(middle_end_time-far_end_time);
      double from_middle = time - middle_end_time;
      return middle_end_coord + (from_middle*flat_slope)*(frac) + (from_middle*far_slope)*(1-frac);
    }
  }
  double height_coord(history const& hist, time_type time, time_type focus_time) {
    
  }
  
  template<class LineFuncType>
  void draw_tree(LineFuncType& line, time_type focus_time) {
    size_t place_in_current_history = 0;
    auto prev = spatial_representation_columns.begin();
    for (auto i = boost::next(prev); i != spatial_representation_columns.end(); ++i) {
      if ((place_in_current_history+1 < w.current_history.size()) && w.current_history[place_in_current_history+1]->start_time == prev.time) {
        ++place_in_current_history;
      }
      for (auto e : i->second.entries) {
        for (auto f : prev->second.entries) {
          if (f.h.back() == e.h.back()) {
            double ewid = time_coord(i->first, focus_time);
            double fwid = time_coord(prev->first, focus_time);
            bool in_current_history = (w.current_history[place_in_current_history] == e.h.back());
            line(fwid, f.height, ewid, e.height, in_current_history);
          }
        }
      }
      prev = i;
    }
    return metadata;
  }
  
  std::pair<time_type, history_tree::history> hist_from_draw_coords(double time_coord, double height, time_type focus_time) {
    history_tree::history best;
    time_type best_time = never;
    double best_dist = 2.0;
    auto prev = spatial_representation_columns.begin();
    for (auto i = boost::next(prev); i != spatial_representation_columns.end(); ++i) {
      time_type max_time = i->first;
      time_type min_time = prev->first;
      double max_time_coord = time_coord(max_time, focus_time);
      double min_time_coord = time_coord(min_time, focus_time);
      if (min_time_coord <= time_coord && time_coord < max_time_coord) {
        while (max_time-min_time > 1) {
          time_type mid_time = divide(max_time+min_time, 2, rounding_strategy<round_down, negative_continuous_with_positive>());
          double mid_time_coord = time_coord(mid_time, focus_time);
          if (mid_time_coord < time_coord) {
            max_time = mid_time;
            max_time_coord = mid_time_coord;
          }
          else {
            min_time = mid_time;
            min_time_coord = mid_time_coord;
          }
        }
        best_time = min_time;
        double frac = double(min_time-prev->first)/double(i->first-prev->first);
        for (auto e : i->second.entries) {
          for (auto f : prev->second.entries) {
            if (f.h.back() == e.h.back()) {
              const double efh = f.height + (e.height-f.height)*frac;
              const double dist = std::abs(efh - height);
              if (dist < best_dist) {
                best_dist = dist;
                best = e.h;
              }
            }
          }
        }
      }
      prev = i;
    }
    return std::pair<time_type, history_tree::history>(best_time, best);
  }
  
  std::unique_ptr<accessor> accessor_after(time_type time) {
    if (time > current_time) {
      set_time(time);
    }
    return steward.accessor_after(time);
  }
  
  
  void expand_to_time(time_type time) {
    expand_to_time_impl(current_history.back(), time);
  }
  void insert_fiat_event(time_type time, uint64_t distinguisher, std::shared_ptr<event> e) {
    validate();
    if (time <= current_time) {
      set_time(time-1);
    }
    size_t place = where_in_current_history_is(time);
    node* cur_node = current_history[place];
    if ((place+1 >= current_history.size()) && (time >= cur_node->end_time)) {
      cur_node->events.insert(std::make_pair(time, fiat_event(distinguisher, e)));
      if (cur_node->end_time < time) {
        update_spatial_representation_heights(cur_node->end_time, time);
        cur_node->end_time = time;
      }
    }
    else {
      assert (time < cur_node->end_time);
      current_history.resize(place+1);
      expand_to_time_impl(current_history.back(), time);
      auto iter = cur_node->children.insert(std::make_pair(time, node(time)));
      node* new_node = &iter->second;
      current_history.push_back(new_node);
      new_node->events.insert(std::make_pair(time, fiat_event(distinguisher, e)));
      auto p = spatial_representation_columns.insert(std::make_pair(time, spatial_representation_column()));
      if (p.second) {
        assert (p.first != spatial_representation_columns.begin());
        p.first->second.entries = boost::prior(p.first)->second.entries;
      }
      for (size_t i = 0; i < p.first->second.entries.size(); ++i) {
        if (p.first->second.entries[i].back() == current_history[current_history.size()-2]) {
          p.first->second.entries.insert(current.begin()+i+1, current_history);
          break;
        }
      }
    }
    validate();
  }
  
  void expand_to_time_impl(node* node, time_type time) {
    if (time > node->end_time) {
      update_spatial_representation_heights(node->end_time, time);
      auto a = spatial_representation_columns.upper_bound(node->end_time);
      auto b = spatial_representation_columns.upper_bound(time);
      assert (a != spatial_representation_columns.begin());
      auto prev = boost::prior(a);
      node* insert_after_this = nullptr;
      history h;
      for (size_t j = 0; j < prev->second.entries.size(); ++j) {
        if (prev->second.entries[j].h.back() == node) {
          h = prev->second.entries[j].h;
          if (j > 0) { insert_after_this = prev->second.entries[j-1].h.back(); }
          break;
        }
      }
      assert (!h.empty());
      for (auto i = a; i != b; ++i) {
        for (size_t j = 0; j < i->second.entries.size(); ++j) {
          if (i->second.entries[j].h.back() == insert_after_this) {
            i->second.entries.insert(i->second.entries.begin()+j+1, h);
            break;
          }
        }
      }
      node->end_time = time;
      
      while (time > boost::prior(spatial_representation_columns.end())->first + (StandardTimeIncrement>>3)) {
        time_type new_time = boost::prior(spatial_representation_columns.end())->first + (StandardTimeIncrement>>3) + 1;
        auto p = spatial_representation_columns.insert(std::make_pair(new_time, boost::prior(spatial_representation_columns.end())->second));
        assert (p.second);
        for (size_t j = 0; j < p.first->second.entries.size(); ) {
          if (p.first->second.entries[j].h.back()->end_time < new_time) {
            p.first->second.entries.erase(p.first->second.entries.begin()+j);
          }
          else {
            ++j;
          }
        }
      }
    }
  }
  size_t where_in_current_history_is(time_type time) {
    size_t result = 0;
    for (size_t i = 1; i < current_history.size(); ++i) {
      if (time >= current_history[i]->start_time) { result = i; }
      else { break; }
    }
    return result;
  }
  void validate_node(node const& n) {
    for (auto const& p : n.events) {
      assert (p.first >= n.start_time);
    }
    for (auto const& p : n.children) {
      assert (p.first > n.start_time);
      assert (p.second.start_time == p.first);
      assert (p.first <= n.end_time);
      validate_node(p.second);
    }
  }
  void validate() {
    node* n = &history_root;
    for (size_t i = 1; i < current_history.size(); ++i) {
      assert(n->children.find(current_history[i]->start_time) != n->children.end());
      n = current_history[i];
    }
    validate_node(history_root);
  }
  
  void set_time(time_type new_time) {
    validate();
    if (new_time > current_time) {
      size_t place_in_current_history = where_in_current_history_is(current_time);
      while (true) {
        node* cur_node = current_history[place_in_current_history];
        node* next_node = (place_in_current_history+1 < current_history.size()) ? current_history[place_in_current_history+1] : nullptr;
        for (auto i = cur_node->events.upper_bound(current_time);
            i != cur_node->events.end() && (i->first <= new_time) &&
            ((!next_node) || (i->first <= next_node->start_time)); ++i) {
          steward.insert_fiat_event(i->first, i->second.distinguisher, i->second.e);
        }
        if (next_node && new_time >= next_node->start_time) {
          ++place_in_current_history;
        }
        else {
          break;
        }
      }
    }
    if (new_time < current_time) {
      size_t place_in_current_history = where_in_current_history_is(current_time);
      while (true) {
        node* cur_node = current_history[place_in_current_history];
        node* next_node = (place_in_current_history+1 < current_history.size()) ? current_history[place_in_current_history+1] : nullptr;
        for (auto i = cur_node->events.upper_bound(new_time);
            i != cur_node->events.end() && (i->first <= current_time) &&
            ((!next_node) || (i->first <= next_node->start_time)); ++i) {
          steward.erase_fiat_event(i->first, i->second.distinguisher);
        }
        if (place_in_current_history > 0 && new_time < cur_node->start_time) {
          --place_in_current_history;
        }
        else {
          break;
        }
      }
    }
    current_time = new_time;
    validate();
  }
  void set_history(history const& new_history) {
    validate();
    size_t first_difference = 0;
    while (current_history[first_difference] == new_history[first_difference]) {
      ++first_difference;
      if (first_difference >= current_history.size() ||
          first_difference >=     new_history.size()) { break; }
    }
    if (first_difference >= current_history.size() &&
        first_difference >=     new_history.size()) { return; }
    const time_type s0 = (first_difference >= current_history.size()) ? max_time : current_history[first_difference]->start_time;
    const time_type s1 = (first_difference >=     new_history.size()) ? max_time :     new_history[first_difference]->start_time;
    const time_type branch_point = std::min(s0, s1);
    if (branch_point-1 < current_time) {
      set_time(branch_point-1);
    }
    
    current_history = new_history;
    validate();
  }
};

#endif
