/*

    Copyright Eli Dupree and Isaac Dupree, 2011, 2012, 2013

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

#ifndef LASERCAKE_MAIN_HPP__
#define LASERCAKE_MAIN_HPP__

#if LASERCAKE_USE_QT
#include <QtOpenGL/QGLWidget>
#include <QtCore/QThread>
#include <QtCore/QMutex>
#include <QtCore/QWaitCondition>
#include <QtGui/QKeyEvent>
#endif

#include <set>

// https://bugreports.qt.io/browse/QTBUG-22829?focusedCommentId=274496#comment-274496
#ifndef Q_MOC_RUN
#include "gl_data_preparation.hpp"
#include "input_representation.hpp"
#include "world.hpp"
#include "gl_rendering.hpp"
#endif

typedef abstract_gl_data gl_data_t;

//is a pointer to avoid copying around all that data
typedef shared_ptr<gl_data_t> gl_data_ptr_t;

// TODO UNITS use units for this:
typedef int64_t microseconds_t;

struct frame_output_t {
  gl_data_ptr_t gl_data_ptr;
  microseconds_t microseconds_for_drawing;
  microseconds_t microseconds_for_simulation;
  std::string extra_debug_info;
};

using input_representation::input_news_t;

struct config_struct {
  std::string scenario;
  bool crazy_lasers;
  distance view_radius;
  bool have_gui;
  bool run_drawing_code;
  bool initially_drawing_debug_stuff;
  //microseconds_t min_microseconds_per_frame;
  int64_t max_frames_per_seconds;
  bool show_frame_timing;
  int64_t exit_after_frames;
  bool use_opengl_thread;
  bool use_simulation_thread;
};
inline std::ostream& operator<<(std::ostream& os, config_struct const& config) {
  return os << '{'
    << "scenario=" << config.scenario
    << "; crazy_lasers=" << config.crazy_lasers
    << "; view_radius=" << config.view_radius
    << "; have_gui=" << config.have_gui
    << "; run_drawing_code=" << config.run_drawing_code
    << "; initially_drawing_debug_stuff=" << config.initially_drawing_debug_stuff
    << "; exit_after_frames=" << config.exit_after_frames
    << "; use_opengl_thread=" << config.use_opengl_thread
    << "; use_simulation_thread=" << config.use_simulation_thread
    << '}';
}

#if !LASERCAKE_USE_QT
// Hackily define some Qt pieces without needing Qt,
// so that it might be less work to port this code
// to another toolkit such as Javascript/Emscripten.
struct QObject{
  QObject(){}
  QObject(QObject*){}
};
struct QThread : public QObject{
  QThread(){}
  QThread(QObject*){}
};
struct QWidget : public QObject{
  QWidget(){}
  QWidget(QObject*){}
};
struct QGLWidget : public QWidget{
  QGLWidget(){}
  QGLWidget(QObject*){}
};
struct QSize {
  QSize() : width_(-1), height_(-1) {}
  QSize(int w, int h) : width_(w), height_(h) {}
  int width() { return width_; }
  int height() { return height_; }

  int width_;
  int height_;
};
#define Q_OBJECT
#define Q_EMIT
#define Q_SLOTS
#define Q_INVOKABLE
#define Q_SIGNALS protected
template<typename... Args>
class pubsub_topic {
public:
  typedef std::function<void(Args...)> subscriber_type;
  // no multithreadedness

  // too bad std::functions are not equality-comparable so
  // it's more complex to say which one to remove

  // or return std::vector<Ret> but then we need to deal with vector of void
  void publish(Args... args) {
    for(auto& sub : subscribers_) {
      sub(args...);
    }
  }
  void subscribe(subscriber_type&& sub) {
    subscribers_.push_back(std::move(sub));
  }
private:
  std::vector<subscriber_type> subscribers_;
};
#endif

#if LASERCAKE_USE_QT
struct sleeper : private QThread {
public:
  using QThread::sleep;
  using QThread::msleep;
  using QThread::usleep;
};
#endif

typedef shared_ptr<worldgen_type> worldgen_ptr;

// for use in its own QThread
class LasercakeSimulator : public QObject {
  Q_OBJECT

public:
  explicit LasercakeSimulator(QObject* parent = 0);

  Q_INVOKABLE void init(worldgen_ptr worldgen, config_struct config);
  Q_INVOKABLE void new_input_as_of(time_unit moment, input_news_t new_input);
  //actually_prepare_graphics=false still sends frame_output_ready with useful debug info.
  Q_INVOKABLE void prepare_graphics(input_news_t input_since_last_prepare, distance view_radius, bool actually_prepare_graphics);

Q_SIGNALS:
  void sim_frame_done(time_unit moment);
  void frame_output_ready(time_unit moment, frame_output_t output);
#if LASERCAKE_USE_QT /* Qt breaks if I #ifdef its signals */
#else
public:
  pubsub_topic<time_unit /*moment*/> sim_frame_done_topic;
  pubsub_topic<time_unit /*moment*/, frame_output_t /*output*/> frame_output_ready_topic;
#endif

private:
  shared_ptr<world> world_ptr_;
  microseconds_t microseconds_last_sim_frame_took_; //hack
  shared_ptr<view_on_the_world> view_ptr_;
  object_identifier currently_focused_object_;
};

struct gl_thread_data_t {
#if LASERCAKE_USE_QT
  QMutex gl_data_lock;
  QWaitCondition wait_for_instruction;
#endif
  atomic::atomic_bool interrupt;

  // access protected by the mutex:
  bool quit_now;
  uint64_t revision;
  gl_data_ptr_t current_data;
  microseconds_t microseconds_last_gl_render_took;
  QSize viewport_size;
};

class LasercakeGLWidget;
// A Qt-signals-based event loop couldn't do quite what I wanted here
// (insufficient info/manipulation of the sequence of queued signals),
// so use threading abstractions + no event loop.
class LasercakeGLThread : public QThread {
  Q_OBJECT

public:
  shared_ptr<gl_thread_data_t> gl_thread_data_;
  // after thread starts, only thread accesses these:
  LasercakeGLWidget* gl_widget_;
  uint64_t last_revision_;
  microseconds_t microseconds_this_gl_render_took_;
  microseconds_t gl_render(gl_data_ptr_t& gl_data_ptr, LasercakeGLWidget& gl_widget, QSize viewport_size);
  gl_renderer gl_renderer_;
protected:
#if LASERCAKE_USE_QT
  void run() override;
#else
  void run();
#endif
};


class LasercakeGLWidget : public QGLWidget {
  Q_OBJECT

public:
  explicit LasercakeGLWidget(bool use_separate_gl_thread, QWidget* parent = 0);
  void update_gl_data(gl_data_ptr_t data);
  microseconds_t get_last_gl_render_microseconds();
  input_representation::input_news_t get_input_news()const;
  // doesn't clear the keys_currently_pressed, only the history:
  void clear_input_news();
  void toggle_fullscreen();
  void toggle_fullscreen(bool fullscreen);

Q_SIGNALS:
  void key_changed(input_representation::key_change_t);
#if LASERCAKE_USE_QT /* Qt breaks if I #ifdef its signals */
#else
public:
  pubsub_topic<input_representation::key_change_t> key_changed_topic;
#endif

#if LASERCAKE_USE_QT
protected:
  bool event(QEvent*) override;
  //void keyPressEvent(QKeyEvent*) override;
  //void keyReleaseEvent(QKeyEvent*) override;
  void mouseMoveEvent(QMouseEvent* event) override;
  void mousePressEvent(QMouseEvent* event) override;
  void mouseReleaseEvent(QMouseEvent* event) override;
  void focusOutEvent(QFocusEvent*) override;
  void resizeEvent(QResizeEvent*) override;
  void paintEvent(QPaintEvent*) override;
  void closeEvent(QCloseEvent*) override;
#endif

private Q_SLOTS:
  void prepare_to_cleanly_close_();

private:
  // We store QKeyEvent->key() as well as -QMouseEvent->button()
  // in qt_key_type_; a minor hack.
  typedef int qt_key_type_;
  void key_change_(
      qt_key_type_ qkey,
      input_representation::key_type ikey,
      bool pressed);
#if LASERCAKE_USE_QT
  void key_change_(QKeyEvent* event, bool pressed);
#endif
  void invoke_render_(); //precondition: you incremented gl_thread_data_->revision
  void grab_input_();
  void ungrab_input_();

  // e.g. in case there can be multiple shift keys pressed at once
  // (two of them on a regular keyboard... or two USB keyboards plugged in,
  //  for that matter!)
  std::multiset<qt_key_type_> keys_currently_pressed_;
  input_representation::key_activity_t input_rep_key_activity_;
  input_representation::keys_currently_pressed_t input_rep_keys_currently_pressed_;
  input_representation::mouse_displacement_t input_rep_mouse_displacement_;
#if LASERCAKE_USE_QT
  QPoint global_cursor_pos_;
#endif
  bool input_is_grabbed_;
  bool use_separate_gl_thread_;
  LasercakeGLThread thread_;
  shared_ptr<gl_thread_data_t> gl_thread_data_;
  bool has_quit_;
};


class LasercakeController : public QObject {
  Q_OBJECT

public:
  explicit LasercakeController(config_struct config, QObject* parent = 0);
  void killSimulator();

public Q_SLOTS:
  void output_new_frame(time_unit moment, frame_output_t output);
  void key_changed(input_representation::key_change_t);

private Q_SLOTS:
  void invoke_simulation_step_();

private:
  config_struct config_;
  boost::scoped_ptr<LasercakeGLWidget> gl_widget_;
  boost::scoped_ptr<QThread> simulator_thread_;
  boost::scoped_ptr<LasercakeSimulator> simulator_;
  microseconds_t monotonic_microseconds_at_beginning_of_frame_;
  microseconds_t monotonic_microseconds_at_beginning_of_ten_frame_block_;
  microseconds_t monotonic_microseconds_at_beginning_of_hundred_frame_block_;
  // We limit the simulation frame-rate to 60 FPS so that the simulation goes
  // at a somewhat consistent speed (which limits the display to the same
  // frequency in the current implementation, which we would also want to limit
  // but it is already limited by this).
  microseconds_t monotonic_microseconds_at_beginning_of_frame_for_framerate_limiter_;
  int64_t frame_;
  time_unit game_time_;
  bool paused_;
  int64_t steps_queued_to_do_while_paused_;
};

#endif
