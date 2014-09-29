Feel++ tutorial
===============


# Compile

To compile the code on irma-atlas, type
```
cmake .
make
```

# Adding a new code

Let `myapp.cpp` a new c++ file using Feel++.
Edit `CMakeLists.txt` and add a line with
```cmake
feelpp_add_application(myapp SRCS myapp.cpp)
```
