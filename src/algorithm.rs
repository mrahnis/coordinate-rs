use std::process;

use convert;
use ellipsoid;

/// NOT YET IMPLEMENTED: Convert ECEF coordinates to LLA coordinates by direct solution of Vermielle.
///
/// Adapted from GeographicLib source code available from:
///
/// http://geographiclib.sourceforge.net/
///
/// References:
/// -----------
/// Vermeille, H. (2002). Direct transformation from geocentric coordinates
/// to geodetic coordinates. J Geod 76(9):451–454.
///
pub fn vermielle(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {
    println!("vermielle not yet implemented");
    process::exit(0);
}

/// NOT YET IMPLEMENTED: Convert ECEF coordinates to LLA coordinates by direct solution of Borokowski.
///
/// Adapted from XyzWin source code (file xyz2plh.c) available from NOAA NGS:
///
/// http://www.ngs.noaa.gov/PC_PROD/XYZWIN/
///
/// References:
/// -----------
/// Borkowski, K. M. (1989).  "Accurate algorithms to transform geocentric
/// to geodetic coordinates", *Bulletin Geodesique*, v. 63, pp. 50-56.
///
/// Borkowski, K. M. (1987).  "Transformation of geocentric to geodetic
/// coordinates without approximations", *Astrophysics and Space Science*,
/// v. 139, n. 1, pp. 1-4.  Correction in (1988), v. 146, n. 1, p. 201.
///
pub fn borokowski(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {
    println!("borokowski not yet implemented");
    process::exit(0);
}

/// Convert ECEF coordinates to LLA coordinates by direct solution of uBlox.
///
/// References:
/// -----------
/// https://microem.ru/files/2012/08/GPS.G1-X-00006.pdf
///
/// http://www.u-blox.com/customersupport/docs/GPS.G1-X-00006.pdf
///
pub fn ublox(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {

    let p:f64 = (x*x + y*y).sqrt();
    let theta:f64 = ((z*ellps.a)/(p*ellps.b)).atan();

    if p < 1.0e-10 {
        let (lambda, phi, h) = convert::handle_polar(x,y,z,ellps);
        return (lambda.to_degrees(),phi.to_degrees(),h)     
    }

    let lambda:f64 = convert::longitude(x, y);

    let phi:f64 = ((z + ellps.e2.powi(2) * ellps.b * theta.sin().powi(3)) / (p - ellps.e1.powi(2) * ellps.a * theta.cos().powi(3))).atan();
    let (r, rn, rm) = ellps.radii(phi);
    let h:f64 = ellps.hae(p, rn, phi);

    (lambda.to_degrees(),phi.to_degrees(),h)
}

/// Convert ECEF coordinates to LLA coordinates by the iterative method of Bowring.
///
/// References:
/// -----------
///
pub fn bowring(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {

    let p:f64 = (x*x + y*y).sqrt(); 

    if p < 1.0e-10 {
        let (lambda, phi, h) = convert::handle_polar(x,y,z,ellps);
        return (lambda.to_degrees(),phi.to_degrees(),h)        
    }

    let lambda:f64 = convert::longitude(x, y);

    // initial estimated values
    let mut phi:f64 = (z/(p*(1f64-ellps.e1.powi(2)))).atan();
    let mut phi_current:f64 = phi;

    let (mut r, mut rn, mut rm) = ellps.radii(phi);
    let mut h:f64 = 0f64;

    // or should it iterate while phi - phi_current?
    while (h - ellps.hae(p, rn, phi)).abs() > 0.000001 {

        let (r, rni, rm) = ellps.radii(phi);
        rn = rni;  // don't like this
        h = ellps.hae(p, rn, phi);
        phi_current = phi;
        phi = (z/(p*(1f64-ellps.e1.powi(2)*rn/(rn+h)))).atan();
    }
    (lambda.to_degrees(),phi_current.to_degrees(),h)
}

/// Convert ECEF coordinates to LLA coordinates by direct solution of Olson.
///
/// Adapted from source code available from Eric Farmer:
///
/// https://possiblywrong.wordpress.com/2014/02/14/when-approximate-is-better-than-exact/
///
/// References:
/// -----------
/// Olson, D. K., Converting Earth-Centered, Earth-Fixed Coordinates to
/// Geodetic Coordinates, IEEE Transactions on Aerospace and Electronic
/// Systems, 32 (1996) 473-476.
///
pub fn olson(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {
    // let e2:f64 = ellps.f * (2.0 - ellps.f);

    // derived parameters
    let e2:f64 = ellps.e1.powi(2);

    let a1:f64 = ellps.a * e2;
    let a2:f64 = a1 * a1;
    let a3:f64 = a1 * e2 / 2.0;
    let a4:f64 = 2.5 * a2;
    let a5:f64 = a1 + a3;
    let a6:f64 = 1.0 - e2;
 
    let mut lat:f64;
    let s:f64;
    let ss:f64;
    let c:f64;

    let w:f64 = (x * x + y * y).sqrt();
    let zp:f64 = z.abs();
    let w2:f64 = w * w;
    let r2:f64 = z * z + w2;
    let r:f64  = r2.sqrt();
    let s2:f64 = z * z / r2;
    let c2:f64 = w2 / r2;
    let u:f64 = a2 / r;
    let v:f64 = a3 - a4 / r;
    if c2 > 0.3 {
        s = (zp / r) * (1.0 + c2 * (a1 + u + s2 * v) / r);
        lat = s.asin();
        ss = s * s;
        c = (1.0-ss).sqrt();
    } else {
        c = (w / r) * (1.0 - s2 * (a5 - u - c2 * v) / r);
        lat = c.acos();
        ss = 1.0 - c * c;
        s = ss.sqrt();
    }
    let g:f64 = 1.0 - e2 * ss;
    let rg:f64 = ellps.a / g.sqrt();
    let rf:f64 = a6 * rg;
    let u:f64 = w - rg * c;
    let v:f64 = zp - rf * s;
    let f:f64 = c * u + s * v;
    let m:f64 = c * v - s * u;
    let p:f64 = m / (rf / g + f);
    lat = lat + p;
    if z < 0.0 {
        lat = -lat;
    }
    let lon = convert::longitude(x, y);
    let h = f + m * p /2.0;
    (lon.to_degrees(), lat.to_degrees(), h)
}

 
#[test]
fn test_bowring() {
    let expected = (-76.26108657, 39.99277705, 230.920);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 1162172.3971573876;
    let y:f64 = -4753390.239612344;
    let z:f64 = 4077519.584501206;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = bowring(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_bowring_pole() {
    let expected = (0.0, 90.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 0.0;
    let y:f64 = 0.0;
    let z:f64 = 6356852.3141;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = bowring(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_bowring_meridian() {
    let expected = (0.0, 40.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 4892784.2046;
    let y:f64 = 0.0;
    let z:f64 = 4078049.8509;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = bowring(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_bowring_equator() {
    let expected = (40.0, 0.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 2181485.5329;
    let y:f64 = -5993582.2425;
    let z:f64 = 0.0000;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = bowring(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}
#[test]
fn test_ublox() {
    let expected = (-76.26108657, 39.99277705, 230.920);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 1162172.3971573876;
    let y:f64 = -4753390.239612344;
    let z:f64 = 4077519.584501206;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = ublox(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);

    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_ublox_pole() {
    let expected = (0.0, 90.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 0.0;
    let y:f64 = 0.0;
    let z:f64 = 6356852.3141;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = ublox(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_ublox_meridian() {
    let expected = (0.0, 40.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 4892784.2046;
    let y:f64 = 0.0;
    let z:f64 = 4078049.8509;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = ublox(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_ublox_equator() {
    let expected = (40.0, 0.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 2181485.5329;
    let y:f64 = -5993582.2425;
    let z:f64 = 0.0000;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = ublox(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_olson() {
    let expected = (-76.26108657, 39.99277705, 230.920);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 1162172.3971573876;
    let y:f64 = -4753390.239612344;
    let z:f64 = 4077519.584501206;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = olson(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);

    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_olson_pole() {
    let expected = (0.0, 90.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 0.0;
    let y:f64 = 0.0;
    let z:f64 = 6356852.3141;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = olson(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_olson_meridian() {
    let expected = (0.0, 40.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 4892784.2046;
    let y:f64 = 0.0;
    let z:f64 = 4078049.8509;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = olson(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}

#[test]
fn test_olson_equator() {
    let expected = (40.0, 0.0, 100.000);
    let ll_tolerance:f64 = 0.00000001;
    let h_tolerance:f64 = 0.000001;

    let x:f64 = 2181485.5329;
    let y:f64 = -5993582.2425;
    let z:f64 = 0.0000;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let (lon,lat,h) = olson(x, y, z, wgs84);
    println!("lon {0} lat {1} h {2}", lon,lat,h);
    println!("lon-diff: {}, lat-diff: {}", (lon-expected.0).abs(), (lat-expected.1).abs());
    println!("h-diff: {}", (h-expected.2).abs());
    if ((lon-expected.0).abs(), (lat-expected.1).abs()) > (ll_tolerance, ll_tolerance) {
        assert!(false);
    }
    if (h-expected.2).abs() > h_tolerance {
        assert!(false);
    }
}