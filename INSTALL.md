iBWA uses CMake which is a cross-platform build tool. Basically it will
generate a Makefile so you can use `make`. The requirements are gcc, gmake,
cmake 2.8+.

Here's the steps to build:

    $ git clone https://github.com/genome/ibwa.git
    ...
    Resolving deltas: 100% (38/38), done.

    $ cd ibwa
    $ mkdir build
    $ cd build
    $ cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local
    ...
    -- Build files have been written to: .../ibwa/build

    $ make
    ...
    Linking CXX executable ibwa-0.5 
    [100%] Built target ibwa-0.5

    $ sudo make install
