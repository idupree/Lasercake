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


#include <iostream>
#include <cstddef>

#include "green_caves.cpp"




// Stuff copied from ../../main.cpp. TODO un-duplicate or something
#ifndef BOOST_CHRONO_HEADER_ONLY
#define BOOST_CHRONO_HEADER_ONLY
#endif
#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <iomanip>
#include <sstream>
#include <locale>
typedef int64_t microseconds_t;

namespace chrono = boost::chrono;
microseconds_t get_this_process_microseconds() {
  return chrono::duration_cast<chrono::microseconds>(chrono::process_real_cpu_clock::now().time_since_epoch()).count();
}

// Usage example:
// LOG << std::setw(6) << (ostream_bundle() << "foo" << 564) << std::endl;
struct ostream_bundle : std::ostream {
  template<typename T> ostream_bundle& operator<<(T const& t) { ss_ << t; return *this; }
  std::string str() { return ss_.str(); }
private:
  std::stringstream ss_;
};
std::ostream& operator<<(ostream_bundle& os, ostream_bundle& b) { return os << b.str(); }
std::ostream& operator<<(std::ostream& os, ostream_bundle& b) { return os << b.str(); }
// show_decimal(1234567, 1000, 10) --> "1234.5"
template<typename Integral, typename Integral2>
std::string show_decimal(Integral us, Integral2 divisor, int places, std::locale const& locale = std::locale()) {
  Integral divisordivisor = 1;
  for(int i = 0; i < places; ++i) { divisordivisor *= 10; }
  
  return (ostream_bundle()
    << (us / divisor)
    << std::use_facet< std::numpunct<char> >(locale).decimal_point()
    << std::setfill('0') << std::setw(places) << std::abs(us / (divisor / divisordivisor) % divisordivisor)
  ).str();
}

std::string show_microseconds(microseconds_t us) {
  return show_decimal(us, 1000, 1);
}
// End copied stuff







struct dont_draw {
  void circle(double, double, double) {}
  void rect(double, double, double, double, bool) {}
  void segment(double, double, double, double, double) {}
};

void do_benchmark() {
  green_caves_ui_backend backend;

  const microseconds_t microseconds_before_benchmark = get_this_process_microseconds();
  
  int64_t ms = 0;
  ms += 1000;
  backend.screen_size = fd_vector(120, 100);
  dont_draw d;
  backend.draw(d);
  for (int i = 0; i < 50000; i += 10) {
    backend.update_to_real_time(i);
    if (i == 1000) { backend.set_key(i, LEFT, true); }
    if (i == 1100) { backend.set_key(i, DOWN, true); }
    if (i == 5000) { backend.mouse_down(i, 10, 10); }
    if (i == 9000) { backend.mouse_up(i, 10, 10); }
    if (i == 9200) { backend.set_key(i, LEFT, false); }
    if (i == 9300) { backend.set_key(i, DOWN, false); }
    if (i == 10000) { backend.mouse_down(i, 110, 8); }
    if (i == 10200) { backend.set_key(i, RIGHT, true); }
    if (i == 10300) { backend.set_key(i, UP, true); }
    if (i == 10400) { backend.mouse_down(i, 64, 60); }
    if (i == 15000) { backend.mouse_up(i, 10, 10); }
    if (i == 15100) { backend.mouse_down(i, 105, 8); }
    if (i == 15150) { backend.mouse_up(i, 10, 10); }
    if (i == 15200) { backend.mouse_down(i, 115, 8); }
    if (i == 15250) { backend.mouse_up(i, 10, 10); }
    if (i == 15300) { backend.mouse_down(i, 105, 8); }
    if (i == 15350) { backend.mouse_up(i, 10, 10); }
    if (i == 15400) { backend.mouse_down(i, 115, 8); }
    if (i == 15450) { backend.mouse_up(i, 10, 10); }
    if (i == 15500) { backend.mouse_down(i, 105, 8); }
    if (i == 15550) { backend.mouse_up(i, 10, 10); }
    if (i == 15600) { backend.mouse_down(i, 115, 8); }
    if (i == 15650) { backend.mouse_up(i, 10, 10); }
  }
  
  const microseconds_t microseconds_after_benchmark = get_this_process_microseconds();
  const microseconds_t microseconds_for_benchmark = microseconds_after_benchmark - microseconds_before_benchmark;
  std::cout << show_microseconds(microseconds_for_benchmark) << " ms\n";
}

int main(int argc, char **argv)
{
  do_benchmark();
}

