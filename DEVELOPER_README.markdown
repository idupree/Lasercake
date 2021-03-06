

Localization
------------

Currently we use Qt Linguist.
(TODO: mark strings with tr(), and check that it works.)
We could switch to gettext if you prefer.

Fonts
-----

To add (or remove) a font, you need to make the following changes:
* Add a .ttf font file to resources/
* Remember that Lasercake is GPLed so you shall also:
  Include your most-preferred format for editing the font next to the .ttf.
* List this file in fonts.qrc
* Add another QFontDatabase::addApplicationFont line in main.cpp
* Find out the "font name" of the font so you can refer to it in your code
* Rebuild Lasercake (these fonts are embedded in the executable)


Cross-compiling
---------------

Notes:

### emscripten ###

This might not entirely work yet, but for starters:

- Extra CMake arguments (Arch Linux, Nov 2015):
    -DCMAKE_TOOLCHAIN_FILE="$EMSCRIPTEN/cmake/Modules/Platform/Emscripten.cmake" -DCMAKE_CXX_FLAGS='-s DEMANGLE_SUPPORT=1 -DBOOST_NO_EXCEPTIONS=1' -DUSE_QT=NO

The path to the toolchain file might be different for you.
Example values of $EMSCRIPTEN are /usr/lib/emscripten or
.../emsdk_portable/emscripten/master.

See https://github.com/kripken/emscripten/blob/master/src/settings.js for
more "-s" Emscripten flags you might be interested in. If you want
to turn exceptions back on (at a performance cost), remove
`-DBOOST_NO_EXCEPTIONS=1` and add `-s DISABLE_EXCEPTION_CATCHING=0` in its
place.

You could run it from the command line by
`(cd [the directory where Lasercake.js is]; node Lasercake.js --help)`
although that's kind of silly because the hope is it will have
a browser UI.  The main use for running it from the command line
is that --run-self-tests works.  You may have to use `nodejs` instead
of `node` on some versions of Debian.  (You may have to install
Node.js for this.  Node probably isn't necessary for the build process,
though, just the tests.)

### from Linux to Windows via MinGW: ###

Fedora: install mingw{32,64}-{gcc,gcc-c++,qt}.
  Use mingw32-cmake or mingw64-cmake instead of cmake.

Arch Linux: After installing several mingw32 packages,

- Extra CMake arguments:
    -DCMAKE_TOOLCHAIN_FILE=cmake/Toolchain-ArchLinux-mingw32.cmake

    (use/make a different toolchain file if your distro puts mingw files in
    a different place, or if you want to use a different toolchain)

- CMake might say

        CMake Warning: Manually-specified variables were not used by the project:
            CMAKE_TOOLCHAIN_FILE

  (it does for me); ignore that warning, it's incorrect.

On any distro, then copy the necessary DLLs next to Lasercake.exe
(the path for those dlls will depend on your distro)
(Anyone: do you know whether you're supposed to do it this way,
 or whether CMake or Windows or such has a better way to deal with
 pulling in these DLLs?)
/usr/i486-mingw32/{bin/{zlib1.dll,libpng15-15.dll},lib/{libgcc_s_sjlj-1.dll,libstdc++-6.dll,QtCore4.dll,QtGui4.dll,QtOpenGL4.dll}}
or
/usr/i686-w64-mingw32/{bin/{QtCore4.dll,QtGui4.dll,QtOpenGL4.dll,zlib1.dll,libpng15-15.dll},lib/{libgcc_s_sjlj-1.dll,libstdc++-6.dll}}
or
/usr/x86_64-w64-mingw32/{bin/{QtCore4.dll,QtGui4.dll,QtOpenGL4.dll,zlib1.dll,libpng15-15.dll},lib/{libgcc_s_sjlj-1.dll,libstdc++-6.dll}}
or
/usr/i686-w64-mingw32/sys-root/mingw/bin/{QtGui4.dll,QtCore4.dll,QtOpenGL4.dll,libstdc++-6.dll,zlib1.dll,libpng16-16.dll,libgcc_s_sjlj-1.dll}
or
/usr/x86_64-w64-mingw32/sys-root/mingw/bin/{QtGui4.dll,QtCore4.dll,QtOpenGL4.dll,libstdc++-6.dll,zlib1.dll,libpng16-16.dll,libgcc_s_seh-1.dll}

Then run
`wine lasercake.exe`
and see if it runs!


Releasing
---------

### process ###

0: Update CHANGELOG.markdown
1: Set Lasercake_VERSION_* to the correct numbers in CMakeLists.txt; commit.
2: Git tag a release candidate.
git tag -u 17062391 Lasercake-[version]-rcN -m'Lasercake-[version]-rcN'
git push --tags
3: Build binaries:
On Fedora Linux,
./release-build.py --update-host-fedora --source --linux --mingw Lasercake-[version]-rcN
On OS X 10.6,
./release-build.py --osx-bare Lasercake-[version]-rcN
4: Upload the release candidate:
Add the release candidate to website-source/downloadable/
for f in downloadable/Lasercake-[version]-rcN*; do gpg --detach-sign "$f"; done
./build.py
./upload.py
5: Get the release candidate tested on several platforms.
6: If there are any problems, repeat starting at step 2
7: Otherwise:
8: git tag and rebuild with new non-rc version number [steps 2-4]
9: update website-source/src/{index.html,downloads.html} to point to
the new version number + (downloads page) have the right megabytes
[re-finish step 4]
10: Update Lasercake_VERSION_* in git CMakeLists.txt to the next odd number
to indicate dev version; even numbers are releases.

### notes and rationale ###

#### LTO (link-time optimization)

Using link-time optimization (tested with GCC 4.7, Feb 2013) cuts the
Lasercake binary size in half, though it doesn't make the program run
significantly faster (probably because we already put the performance
critical code in templates in the headers).  Linking with -flto takes
about a minute on my modern CPU so you probably don't want to do this
while developing.
-DLTO=ON (which is shorthand for
    -DCMAKE_CXX_FLAGS=-flto -DCMAKE_C_FLAGS=-flto -DCMAKE_EXE_LINKER_FLAGS=-fwhole-program
)

#### PGO (profile-guided optimization)

Profile-guided optimization did not currently appear to profitable enough
to bother with, but it has been before, so perhaps it will in the future. How to:
% cmake path/to/lasercake -DCMAKE_CXX_FLAGS=-fprofile-generate -DCMAKE_C_FLAGS=-fprofile-generate -DCMAKE_EXE_LINKER_FLAGS=-fprofile-generate
% make -j3
run the generated lasercake in various ways until you've exercised its various
performance-critical code paths.  This generates .gcda profiles next to the .o
binaries in CMakeFiles/lasercake.dir; each run adds to the existing .gcda files.
You probably want to keep this binary in case you didn't like the profile and want
to run it some more, so rename it.
% cmake path/to/lasercake -DCMAKE_CXX_FLAGS='-fprofile-use -fprofile-correction' -DCMAKE_C_FLAGS='-fprofile-use -fprofile-correction' -DCMAKE_EXE_LINKER_FLAGS='-fprofile-use -fprofile-correction'
% make -j3
enjoy your new binary.
(-fprofile-correction is necessary because Lasercake is multi-threaded.
Alternatively, you can run the profile-generating Lasercake with --no-threads.)

#### Mac OS X

Uses CPack DragNDrop to create a DMG containing the .app (in a version
that should start on other users' systems), a ReadMe, and an
alias to /Applications for the user.

Compilers: The paths in release-build.py are right if using Macports Clang.
If you have a new enough XCode, it should have a new enough Clang already.
(I don't have the latest OS X to check.  Apple keep including snapshot
Clang versions, so I don't know which one is the first to fix all the bugs
that make Clang crash when compiling Lasercake.  Upstream Clang 3.1 is
too old and upstream Clang 3.2 is new enough.)  Recent GCC from Macports
probably works too, but Apple are moving towards Clang so it's probably
better to prefer Clang on Mac (all else equal).

We target 10.6 because a good fraction of OS X users are still on it.
Qt has abandoned any official support 10.5.  Targeting 10.5 while
building on 10.6 worked fine for everything except running on 10.5.
Currently, we build on 10.6 anyway, so the attempt at explicitly targeting
10.6 is moot.  If you have too new an OS X / XCode setup, you might not have
the 10.6 SDK.

USE_BOOST_CXX11_LIBS: OS X before 10.7 uses GCC 4.2's libstdc++, which doesn't
support any C++11 library features.  This flag makes us use equivalent Boost
libraries instead.  Alternatives: OS X 10.7 has a copy of libc++ (C++ library
created by LLVM/Clang folks).  I hear 10.7's copy supports (some? all?) C++11
library features, but I have not tested their version (I don't have 10.7, and
in any case using it for releases would prevent the binary from running on
10.6).  It might also work to build your own recent lib(std)c++ and link it
into the binary, although it is problematic (in theory, and in practice if you
are using system Boost) to have more than one C++ runtime in the same
executable.

#### Linux

What is the best file format for Linux ReadMe:s?
Does the ReadMe contain any instruction to install Qt?

#### Windows (cross-compiled from Linux)

We use dynamic linking because I so far failed at getting Lasercake to
cross-compile with static linking (using Fedora Linux).

Currently we ship a .zip, but CMake can produce NSIS installers
that might be worth using instead.

#### Source

Full version: with .git/ and bundled_libs/ (which are similar sizes both
much larger than the rest of the code)

Minimal version: without .git/ or bundled_libs/.  Suitable for Linux-distro
packagers if the distro has a compatible Boost version.

TODO: figure something out regarding Windows line endings for the .zip.
Currently I do nothing special for the line endings and release using Linux,
so they are presumably Unix-style LF line endings in all source releases.
