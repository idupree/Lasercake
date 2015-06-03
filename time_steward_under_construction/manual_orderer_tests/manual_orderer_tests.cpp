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

#include <vector>
#include "../manual_orderer.hpp"

struct metadata {
  typedef manual_orderer<metadata> orderer;
  typedef orderer::manually_orderable manually_orderable;
  
  manually_orderable prev;
  manually_orderable next;
};
typedef metadata::orderer orderer;
typedef metadata::manually_orderable manually_orderable;

struct tester {
  
  orderer o;
  std::vector<manually_orderable> all;
  void add() {
    manually_orderable a = manually_orderable::construct();
    if (all.empty()) {
      o.insert_only(a);
    }
    else {
      manually_orderable rel = all[rand()%all.size()];
      bool after = rand()%2;
      if (after) {
        a->prev = rel;
        a->next = rel->next;
        rel->next = a;
        if (a->next) { a->next->prev = a; }
      }
      else {
        a->next = rel;
        a->prev = rel->prev;
        rel->prev = a;
        if (a->prev) { a->prev->next = a; }
      }
      o.insert_adjacent(a, rel, after);
    }
    all.push_back(a);
  }
  void del() {
    size_t d = rand()%all.size();
    manually_orderable a = all[d];
    all[d] = all.back();
    all.pop_back();
    if (a->next) { a->next->prev = a->prev; }
    if (a->prev) { a->prev->next = a->next; }
    o.erase(a);
  }
  void check() {
    for (manually_orderable& a : all) {
      assert (o.contains(a));
    }
    for (manually_orderable& a : all) {
      if (a->next) { assert (a < a->next); }
      if (a->prev) { assert (a->prev < a); }
    }
  }
  
  void test() {
    for (int i = 0; i < 10; ++i) {
      assert (all.empty());
      for (int64_t j = 0; j < (1LL << i); ++j) {
        add();
        check();
      }
      for (int64_t j = 0; j < (1LL << i); ++j) {
        del();
        check();
      }
    }
    for (int i = 0; i < 20; ++i) {
      assert (all.empty());
      for (int64_t j = 0; j < (1LL << i); ++j) {
        add();
      }
      check();
      for (int64_t j = 0; j < (1LL << i); ++j) {
        del();
      }
      assert (all.empty());
      for (int64_t j = 0; j < (1LL << i); ++j) {
        add();
        add();
        add();
        del();
      }
      check();
      for (int64_t j = 0; j < (1LL << i); ++j) {
        add();
        del();
        del();
        del();
      }
    }
    
  }
};

int main(int argc, char **argv)
{
  tester t;
  t.test();
  
  return 0;
}

