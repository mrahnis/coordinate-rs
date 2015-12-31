#[macro_use]
extern crate clap;
extern crate csv;
extern crate rustc_serialize;

pub mod algorithm;
pub mod convert;
pub mod ellipsoid;

use std::path::Path;
use rustc_serialize::{Encodable};

/// Constructs a Point.
///
#[derive(RustcEncodable)]
pub struct Point {
    /// x coordinate
    x: f64,
    /// y coordinate
    y: f64,
    /// z coordinate
    z: f64,
}

/// A commandline utility to convert coordinate triples between LLA and ECEF.
///
fn main() {

    // cli
    let matches = clap_app!(ecef2lla =>
        (version: "0.1.0")
        (author: "Mike R <mike.rahnis@gmail.com>")
        (about: "Does awesome things")
        (@arg INPUT: -i --input +required "Input filename")
        (@arg ELLIPSOID: -e --ellipsoid +takes_value "Ellipsoid to use")
        (@arg OUTPUT: -o --output +takes_value "Output filename")
        (@subcommand lla2ecef =>
            (about: "Converts geodetic lla to geocentric xyz")
            (version: "0.0.1")
            (author: "Mike <mike@other.com>")
        )
        (@subcommand ecef2lla =>
            (about: "Converts geocentric xyz to geodetic lla")
            (version: "0.0.1")
            (author: "Mike <mike@other.com>")
            (@arg METHOD: -m --method +takes_value "Specify which algorithm to use. Choices are bowring, olson, ublox.")
        )
    ).get_matches();

    // reading from csv converting and printing to console
    let mut rdr = csv::Reader::from_file(matches.value_of("INPUT").unwrap()).unwrap();
    let path = Path::new("output.csv");
    //let mut wtr = csv::Writer::from_file(path);
    let mut wtr = csv::Writer::from_memory();

    let ellipsoid = matches.value_of("ELLIPSOID").unwrap_or("WGS84");

    let ellps:ellipsoid::Ellipsoid =
        match ellipsoid {
            "WGS84" => ellipsoid::Ellipsoid::new(6378137.000, 298.257223563),
            _ => ellipsoid::Ellipsoid::new(6378137.000, 298.257223563)
        };

    if let Some(matches) = matches.subcommand_matches("ecef2lla") {
        if matches.is_present("METHOD") {
            let method = matches.value_of("METHOD").unwrap();
            println!("Using method...{}", method);
            for record in rdr.decode() {
                let (x, y, z): (f64, f64, f64) = record.unwrap();
                let (lon,lat,hae) = convert::ecef2lla(x,y,z,ellps,method);
                println!("({}, {}, {})", lon, lat, hae);
                let pt = Point { x: lat, y: lon, z: hae };
                wtr.encode(pt).ok().expect("CSV writer error");
            }

        } else {
            println!("Using default method...");
        }
    }
    if let Some(matches) = matches.subcommand_matches("lla2ecef") {
        for record in rdr.decode() {
            let (lat, lon, hae): (f64, f64, f64) = record.unwrap();
            let lambda:f64 = lon.to_radians();
            let phi:f64 = lat.to_radians();
            let (x,y,z) = convert::lla2ecef(lambda,phi,hae,ellps);
            println!("({}, {}, {})", x, y, z);
            let pt = Point { x: x, y: y, z: z };
            wtr.encode(pt).ok().expect("CSV writer error");
        }
    }
}
