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
  typedef typename TimeSteward::event event;
  typedef typename TimeSteward::accessor accessor;
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
  time_type max_reached_time;
  static const time_type max_column_separation = (StandardTimeIncrement>>3) + 1;
  
  struct spatial_representation_entry {
    spatial_representation_entry(double height, history h):height(height),h(h){}
    double height;
    history h;
  };
  struct spatial_representation_column {
    std::vector<spatial_representation_entry> entries;
    std::vector<double> desired_heights(spatial_representation_column const& prev, time_type time_diff) {
      std::vector<double> results;
      double inc = double(1) / double(entries.size() + 1);
      double current_height_asymptote = 0;
      for (auto const& e : entries) {
        bool in_prev = false;
        double prev_height;
        for (auto const& f : prev.entries) {
          if (f.h.back() == e.h.back()) {
            in_prev = true;
            prev_height = f.height;
            current_height_asymptote += inc;
            break;
          }
        }
        if (in_prev) {
          // decay towards evenly-spaced
          double decay_factor = std::pow(0.8, double(time_diff)/double(StandardTimeIncrement));
          assert (decay_factor >= 0.0);
          assert (decay_factor <= 1.0);
          results.push_back(prev_height*decay_factor + current_height_asymptote*(1.0-decay_factor));
          assert (results.back() >= 0.0);
          assert (results.back() <= 1.0);
        }
        else {
          assert (!results.empty());
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
    max_reached_time(0)
  {
    history_root.end_time = 0;
    current_history.push_back(&history_root);
    spatial_representation_column c;
    c.entries.push_back(spatial_representation_entry(0.5, current_history));
    spatial_representation_columns.insert(std::make_pair(time_type(0), c));
    spatial_representation_columns.insert(std::make_pair(time_type(max_column_separation), spatial_representation_column()));
    validate();
  }
  
  void update_spatial_representation_heights(time_type start, time_type end) {
    auto i = spatial_representation_columns.upper_bound(start);
    assert (i != spatial_representation_columns.begin());
    auto prev = boost::prior(i);
    for (; i != spatial_representation_columns.end(); ++i) {
      // decay towards desired_heights
      double decay_factor = 1.0;
      if (i->first > end) { decay_factor = std::pow(0.5, double(end - start)/double(i->first - end)); }
      if (decay_factor < 0.0001) break;
      assert (decay_factor >= 0.0);
      assert (decay_factor <= 1.0);
      std::vector<double> desired_heights = i->second.desired_heights(prev->second, i->first - prev->first);
      for (size_t j = 0; j < desired_heights.size(); ++j) {
        assert (i->second.entries[j].height >= 0.0);
        assert (i->second.entries[j].height <= 1.0);
        i->second.entries[j].height = i->second.entries[j].height*(1.0-decay_factor) + desired_heights[j]*decay_factor;
        assert (i->second.entries[j].height >= 0.0);
        assert (i->second.entries[j].height <= 1.0);
      }
      prev = i;
    }
  }
  
  double time_coord(time_type time, time_type focus_time) {
    time_type flat_dist = 2*StandardTimeIncrement;
    double flat_slope = 0.1/double(flat_dist);
    double double_max_time = max_reached_time;
    if (double_max_time < 0.9/flat_slope) {
      return double(time)*flat_slope;
    }
    
    double focus_time_coord = 0.1 + 0.8*(double(focus_time) / double_max_time);
    if (focus_time-flat_dist < time && time < focus_time+flat_dist) {
      return focus_time_coord + double(time-focus_time)*flat_slope;
    }
    else {
      double s = sign(time-focus_time);
      double far_end_time = (s>0) ? (double_max_time+1.0) : -1.0;
      double far_end_coord = (s>0) ? 1 : 0;
      double middle_end_time = (focus_time + s*flat_dist);
      double middle_end_coord = (focus_time_coord + s*0.1);
      double average_slope = (middle_end_coord-far_end_coord)/(middle_end_time-far_end_time);
      //double frac = (double(time)-far_end_time) / double(middle_end_time-far_end_time);
      double from_middle = double(time) - middle_end_time;
      //double result = middle_end_coord + (from_middle*flat_slope)*(frac) + (from_middle*average_slope)*(1-frac);
      double a = from_middle / double(far_end_time-middle_end_time);
      double k = average_slope / flat_slope;
      assert (k >= 0.0);
      if (k >= 1.0) {
        // just extend the flat area (it won't reach the end of the drawing space,which is OK)
        return focus_time_coord + double(time-focus_time)*flat_slope;
      }
      // wolfram alpha "solve ((1-k)x^2 + kx) = a for x"
      double fancy_func_result = (k - std::sqrt(k*k + 4*a*(1-k)))/(2*(k-1));
      double result = middle_end_coord + fancy_func_result * (far_end_coord-middle_end_coord);
      //std::cerr << time << ", " << focus_time << ", " << result << "\n";
      return result;
    }
  }
  time_type time_from_coord(double drawn_time_coord, time_type min_time, time_type max_time, time_type focus_time) {
    while (max_time-min_time > 1) {
      time_type mid_time = divide(max_time+min_time, 2, rounding_strategy<round_down, negative_continuous_with_positive>());
      double mid_time_coord = time_coord(mid_time, focus_time);
      if (mid_time_coord < drawn_time_coord) {
        max_time = mid_time;
      }
      else {
        min_time = mid_time;
      }
    }
    return min_time;
  }
  double height_coord(history const& hist, time_type time, time_type focus_time) {
    caller_error_if (time > hist.back()->end_time, "getting height beyond end of history");
    auto i = spatial_representation_columns.upper_bound(time);
    assert (i != spatial_representation_columns.begin());
    assert (i != spatial_representation_columns.end());
    auto prev = boost::prior(i);
    double ewid = time_coord(i->first, focus_time);
    double fwid = time_coord(prev->first, focus_time);
    double frac = (time_coord(time, focus_time) - fwid) / (ewid - fwid);
    node* n = hist[where_in_history_is(hist, time)];
    double fallback = -1;
    
    for (auto const& f : prev->second.entries) { if (f.h.back() == n) {
      fallback = f.height;
      for (auto const& e : i->second.entries) { if (e.h.back() == n) {
        return e.height*frac + f.height*(1-frac);
      }}
    }}
    assert (fallback >= 0.0);
    return fallback;
  }
  
  template<class LineFuncType>
  void draw_tree(LineFuncType& line, time_type focus_time) {
    size_t place_in_current_history = 0;
    auto prev = spatial_representation_columns.begin();
    for (auto i = boost::next(prev); i != spatial_representation_columns.end(); ++i) {
      if ((place_in_current_history+1 < current_history.size()) && current_history[place_in_current_history+1]->start_time == prev->first) {
        ++place_in_current_history;
      }
      for (auto const& e : i->second.entries) {
        for (auto const& f : prev->second.entries) {
          if (f.h.back() == e.h.back()) {
            double ewid = time_coord(i->first, focus_time);
            double fwid = time_coord(prev->first, focus_time);
            assert (ewid >= fwid);
            bool in_current_history = (current_history[place_in_current_history] == e.h.back());
            line(fwid, f.height, ewid, e.height, in_current_history);
          }
        }
      }
      prev = i;
    }
  }
  
  std::pair<time_type, history_tree::history> hist_from_draw_coords(double drawn_time_coord, double height, time_type focus_time) {
    history_tree::history best;
    time_type best_time = TimeSteward::never;
    double best_dist = 2.0;
    auto prev = spatial_representation_columns.begin();
    for (auto i = boost::next(prev); i != spatial_representation_columns.end(); ++i) {
      time_type max_time = i->first;
      time_type min_time = prev->first;
      double max_time_coord = time_coord(max_time, focus_time);
      double min_time_coord = time_coord(min_time, focus_time);
      if (min_time_coord <= drawn_time_coord && drawn_time_coord < max_time_coord) {
        best_time = time_from_coord(drawn_time_coord, min_time, max_time, focus_time);
        double frac = double(best_time-prev->first)/double(i->first-prev->first);
        for (auto const& e : i->second.entries) {
          for (auto const& f : prev->second.entries) {
            if (f.h.back() == e.h.back()) {
              const double efh = e.height*frac + f.height*(1-frac);
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
    validate();
    expand_to_time_impl(current_history.back(), time);
    validate();
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
      expand_to_time_impl(cur_node, time);
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
        if (p.first->second.entries[i].h.back() == current_history[current_history.size()-2]) {
          spatial_representation_entry e(
            p.first->second.entries[i].height,
            current_history);
          p.first->second.entries.insert(p.first->second.entries.begin()+i+1, e);
          break;
        }
      }
    }
    validate();
  }
  
  void expand_to_time_impl(node* n, time_type time) {
    if (time > n->end_time) {
      time_type old_end_time = n->end_time;
      auto a = spatial_representation_columns.upper_bound(old_end_time);
      auto b = spatial_representation_columns.upper_bound(time);
      assert (a != spatial_representation_columns.begin());
      auto prev = boost::prior(a);
      history h;
      for (size_t j = 0; j < prev->second.entries.size(); ++j) {
        if (prev->second.entries[j].h.back() == n) {
          h = prev->second.entries[j].h;
          for (auto i = a; i != b; ++i) {
            bool inserted = false;
            for (size_t k = j+1; k < prev->second.entries.size(); ++k) {
              for (size_t l = 0; l < i->second.entries.size(); ++l) {
                if (i->second.entries[l].h.back() == prev->second.entries[k].h.back()) {
                  spatial_representation_entry e(
                    0.5 * (i->second.entries[l].height + ((l > 1) ? i->second.entries[l-1].height : 0.0)),
                    h);
                  i->second.entries.insert(i->second.entries.begin()+l, e);
                  inserted = true;
                  goto doublebreak;
                }
              }
            }
            doublebreak:
            
            if (!inserted) {
              spatial_representation_entry e(
                i->second.entries.empty() ? 0.5 : 0.5 * (i->second.entries.back().height),
                h);
              i->second.entries.push_back(e);
            }
          }
          break;
        }
      }
      
      assert (!h.empty());
      
      n->end_time = time;
      if (max_reached_time < time) {
        max_reached_time = time;
      }
      
      while (time >= boost::prior(spatial_representation_columns.end())->first) {
        time_type new_time = boost::prior(spatial_representation_columns.end())->first + max_column_separation;
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
      update_spatial_representation_heights(old_end_time, time);
    }
  }
  size_t where_in_history_is(history const& h, time_type time) {
    size_t result = 0;
    for (size_t i = 1; i < h.size(); ++i) {
      if (time >= h[i]->start_time) { result = i; }
      else { break; }
    }
    return result;
  }
  size_t where_in_current_history_is(time_type time) {
    return where_in_history_is(current_history, time);
  }
  void validate_node(node const& n) {
    assert (boost::prior(spatial_representation_columns.end())->first > n.end_time);
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
    assert (boost::prior(spatial_representation_columns.end())->second.entries.empty());
    node* n = &history_root;
    for (size_t i = 1; i < current_history.size(); ++i) {
      assert(n->children.find(current_history[i]->start_time) != n->children.end());
      n = current_history[i];
    }
    validate_node(history_root);
    auto i = spatial_representation_columns.begin();
    assert (i != spatial_representation_columns.end());
    for (; i != spatial_representation_columns.end(); ++i) {
      //assert (!i->second.entries.empty());
      for (size_t j = 0; j < i->second.entries.size(); ++j) {
        assert (!i->second.entries[j].h.empty());
      }
    }
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
    const time_type s0 = (first_difference >= current_history.size()) ? TimeSteward::max_time : current_history[first_difference]->start_time;
    const time_type s1 = (first_difference >=     new_history.size()) ? TimeSteward::max_time :     new_history[first_difference]->start_time;
    const time_type branch_point = std::min(s0, s1);
    if (branch_point-1 < current_time) {
      set_time(branch_point-1);
    }
    
    current_history = new_history;
    validate();
  }
};

#endif
