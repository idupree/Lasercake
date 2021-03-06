#!/usr/bin/env python
# (Python >= 2.6 including >= 3)
docs = """
Usage: ./dev-build.py [cmake args] [certain make args]

This script wraps `cmake` and `make`; you can use CMake directly
or use this, whichever you prefer.  This script takes `cmake` and
some `make` flags directly on the command-line, puts the build in
a directory under build/ that is determined by which `cmake` flags you
specified (so that you can store separate intermediate compilation results
for various flag/compiler combinations), runs Lasercake's self-tests,
and if they're successful, copies Lasercake to the main directory.
"""

import os, os.path, sys, subprocess, re, shutil, hashlib, datetime

try:
    import multiprocessing
    MAKE_PARALLEL_JOBS_DEFAULT = multiprocessing.cpu_count() + 1
except (ImportError,NotImplementedError):
    MAKE_PARALLEL_JOBS_DEFAULT = 2


ansi_grey = '\033[30m'
ansi_red = '\033[31m'
ansi_green = '\033[32m'
ansi_yellow = '\033[33m'
ansi_blue = '\033[34m'
ansi_magenta = '\033[35m'
ansi_cyan = '\033[36m'
ansi_white = '\033[37m'
ansi_end = '\033[0m'

os.environ['PATH'] = '/usr/lib/colorgcc/bin:'+os.environ['PATH']

def say(string):
    sys.stdout.write(string)
    sys.stdout.flush()

def escaped_control_characters(string):
    return re.sub(r"[\x01-\x1F\x7F]", r"\\", string)

def say_we_are_calling(string):
    say(ansi_cyan+'% '+escaped_control_characters(string)+ansi_end+'\n')

def hash_list(l):
    hashes = (hashlib.sha256(arg.encode('utf8')).digest() for arg in l)
    return hashlib.sha256(b''.join(hashes)).hexdigest()

def replace_app(src, dst):
    """
    Copies src to dst, replacing dst, even if dst is a running executable.
    Complications this deals with:
    On Windows, one can move but not remove or edit a running binary.
    On Linux, one can move and unlink but not edit a running binary.
    On Mac, if using .app bundles, the app is a directory rather than a file.
    """
    old = None
    if os.path.exists(dst):
        unique_str = datetime.datetime.fromtimestamp(
                os.stat(dst).st_mtime).isoformat()
        path, basefile = os.path.split(dst)
        base, ext = os.path.splitext(basefile)
        old = os.path.join(path, base+'.'+unique_str+ext)
        os.rename(dst, old)
    if os.path.isdir(src):
        shutil.copytree(src, dst, symlinks=True)
    else:
        shutil.copy2(src, dst)
    if old != None:
        if os.path.isdir(old):
            shutil.rmtree(old)
        else:
            try: os.unlink(old)
            except OSError:
                say("Couldn't remove old (running?) executable "+os.path.basename(old))

def main():
    os.chdir(os.path.dirname(os.path.join('.', __file__)))
    try: subprocess.call(['cmake', '--version'])
    except OSError:
        say(ansi_red+"Error: 'cmake' not found; please install it."+ansi_end+'\n')
        exit(1)
    try: subprocess.call(['make', '--version'])
    except OSError:
        say(ansi_red+"Error: 'make' not found; please install it."+ansi_end+'\n')
        exit(1)
    try:
        status = subprocess.call(['time', '-f', '', 'true'])
        time_has_f = (status == 0)
        has_time = True
    except OSError:
        has_time = False
    cmake_args = []
    make_args = []
    making_lasercake = True
    running_tests = True
    make_jobs_flags = ['-j'+str(MAKE_PARALLEL_JOBS_DEFAULT)]
    for arg in sys.argv[1:]:
        if arg in ['-?', '-h', '--help', '-help']:
            say(docs)
            sys.exit()
        elif arg == 'no-test':
            running_tests = False
        elif arg == 'VERBOSE=1':
            make_args.append(arg)
        elif re.search(r'^-j[0-9]*$', arg):
            make_jobs_flags = [arg]
        else:
            cmake_args.append(arg)
        if arg == '-DBUILD_SELF_TESTS=OFF':
            running_tests = False
    # Make it more likely the directories are named non-conflictingly
    # even though we delete some special characters from their paths.
    hash_for_this_config = hash_list(cmake_args)
    name_for_this_config = hash_for_this_config[:10] + '_' + (','.join([re.sub(r'[^-+.a-zA-Z0-9]', r'_', arg) for arg in cmake_args]) or 'default')
    build_dir = 'build/'+name_for_this_config
    try: os.remove('CMakeCache.txt')
    except OSError: pass
    try:
        os.makedirs(build_dir)
    except OSError:
        say(ansi_green+'This build dir already exists!'+ansi_end+'\n')
        cmake_args = []
    to_call_cmake = ['cmake', '../../'] + cmake_args
    say_we_are_calling('cd '+build_dir+'; '+'   '.join(to_call_cmake))
    say(ansi_cyan+'''  (^^ not escaped properly in these info messages - doin' it right in python)'''+ansi_end+'\n')
    os.chdir(build_dir)
    subprocess.check_call(to_call_cmake)
    if not has_time:
        time_args = []
    elif time_has_f:
        time_args = ['time', '-f', (ansi_magenta+'`make` took %E'+ansi_end)]
    else:
        time_args = ['time']
    to_call_make = time_args + ['make'] + make_jobs_flags + make_args
    say_we_are_calling('   '.join(to_call_make))
    make_status = subprocess.call(to_call_make)
    if has_time and not time_has_f:
        say(ansi_magenta+'  (time to run `make`).'+ansi_end+'\n')
    # How to print the compiler & flags here? Can we get it from cmake
    # somehow? Add to the cmakelists to write those?
    # but with the -l's and all... hmm.
    if make_status != 0:
        # skip make_jobs_flags (-jN) to get just one file's error
        # messages not mixed up with the others
        to_call_make_again = ['make'] + make_args
        say_we_are_calling('   '.join(to_call_make_again))
        subprocess.call(to_call_make_again)
        say(ansi_red+'build failed'+ansi_end+'\n')
        exit(1)
    is_windows_exe = is_apple_app = is_js = False
    app_name = exe_name = None
    if os.path.exists('Lasercake'):
        exe_name = 'Lasercake'
        app_name = [exe_name]
    if os.path.exists('Lasercake.exe'):
        exe_name = 'Lasercake.exe'
        app_name = [exe_name]
        is_windows_exe = True
    if os.path.exists('Lasercake.js'):
        exe_name = 'Lasercake.js'
        app_name = [exe_name, exe_name+'.mem']
        is_js = True
    if os.path.exists('Lasercake.app'):
        exe_name = 'Lasercake.app/Contents/MacOS/Lasercake'
        app_name = ['Lasercake.app']
        is_apple_app = True
        # The bundle isn't complete without calling `make install`
        subprocess.check_call(['make', 'install'])
    if exe_name == None:
        say(ansi_red+"couldn't find Lasercake binary?!"+ansi_end+'\n')
        exit(1)
    # poor estimate that assumes people only cross-compile
    # if it involves Windows (or Javascript, but we can run
    # tests even for Javascript):
    is_cross_compiling = (is_windows_exe != (os.name == 'nt'))
    if making_lasercake:
        say(ansi_green+'and you got:\n./'+build_dir+'/'+exe_name+'\n(etc.)'+ansi_end+'\n')
    if is_cross_compiling:
        say(ansi_yellow+'Cross-compiling, so not running tests or copying the binary.'+ansi_end+'\n')
    else:
        if running_tests:
            say(ansi_cyan+'Testing...\n')
            say_we_are_calling('(cd '+build_dir+'; '+
                ('node ' if is_js else '')+'./'+exe_name+' --run-self-tests)')
            test_status = subprocess.call((['node'] if is_js else []) +
                ['./'+exe_name, '--run-self-tests'])
            if test_status != 0:
                say(ansi_red+'Tests failed.\n')
                exit(test_status)
            say(ansi_green+'success')
        else:
            say(ansi_yellow+'NOT RUNNING TESTS')
        if making_lasercake:
            say('; ')
            for f in app_name:
                say('copying '+build_dir+'/'+f+' to ./'+f+'\n')
        say(ansi_end+'\n')
        if making_lasercake:
            for f in app_name:
                replace_app(f, '../../'+f)
    exit()

if __name__ == '__main__':
    main()

