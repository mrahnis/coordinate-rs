==========
Coordinate
==========
Coordinate converts between LLA and ECEF coordinates. It is something I made in order to try the Rust programming language. It is not the best choice (or even a good one) for coordinate conversion. Other libraries offer implementations with more options, better numerical precision, and better handling of corner cases. You might look to:

* [GDAL]
* [Proj4]
* [GeographicLib]

=========
Compiling
=========
Install Rust compiler from:
https://www.rust-lang.org/

In your terminal or command window cd into the source directory and type:

```shell
cargo build
```

To run tests type:

```shell
cargo test
```

The build process produces a binary executable in ./target/debug/ or in ./target/release/ if building with the --release option.

[GDAL]: http://www.gdal.org/
[Proj4]: https://github.com/OSGeo/proj.4
[GeographicLib]: http://geographiclib.sourceforge.net/