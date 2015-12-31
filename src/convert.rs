use ellipsoid;
use algorithm;

/// Calculate geodetic longitude.
///
pub fn longitude(x:f64, y:f64) -> f64 {
    // if x are near zero then we are at the pole; use zero longitude
    let lambda:f64 = if (x.abs() + y.abs()) < 1.0e-10 {
        0.0
    } else {
        (y/x).atan()
    };
    return lambda;
}

/// Handle special cases for ecef2lla algorithms needing help at poles.
///
pub fn handle_polar(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {
    let lambda:f64 = longitude(x,y);
    let phi:f64 = if z < 0.0 {
        -90f64.to_radians()
    } else {
        90f64.to_radians()
    };
    let (r, rn, rm) = ellps.radii(phi);
    let h:f64 = ellps.hae2(x, y, z, r);

    return (lambda, phi, h);
}

/// Convert LLA coordinates to ECEF coordinates.
///
/// # Examples
///
pub fn lla2ecef(lambda:f64,phi:f64,h:f64,ellps:ellipsoid::Ellipsoid) -> (f64, f64, f64) {
    // LLA to ECEF
    // φ : geodetic latitude, phi, in radians
    // ƛ : geodetic longitude, lambda, in radians
    // h : height above the ellipsoid in meters
    // N : radius of curvature of the ellipsoid in meters
    let (r, rn, rm) = ellps.radii(phi);
    let x:f64 = (rn + h)*phi.cos()*lambda.cos();
    let y:f64 = (rn + h)*phi.cos()*lambda.sin();
    let z:f64 = ((ellps.b.powi(2)/ellps.a.powi(2))*rn + h)*phi.sin();

    (x,y,z)
}

/// Convert ECEF coordinates to LLA coordinates.
///
/// # Examples
///
pub fn ecef2lla(x:f64,y:f64,z:f64,ellps:ellipsoid::Ellipsoid,algo:&str) -> (f64, f64, f64) {
    let f:fn(x:f64,y:f64,z:f64,ellipsoid::Ellipsoid) -> (f64, f64, f64) =
        match algo {
            "olson" => algorithm::olson,
            "bowring" => algorithm::bowring,
            "ublox" => algorithm::ublox,
            "vermielle" => algorithm::vermielle,
            "borokowski" => algorithm::borokowski,
            _ => algorithm::olson
        };
    f(x,y,z,ellps)
}

#[test]
fn test_lla2ecef() {
    let expected = (1162172.3971573876, -4753390.239612344, 4077519.584501206);

    let lat:f64 = 39.99277705;
    let lon:f64 = -76.26108657;
    let h:f64 = 230.920;
    let wgs84 = ellipsoid::Ellipsoid::new(6378137.000, 298.257223563);

    let lambda:f64 = lon.to_radians();
    let phi:f64 = lat.to_radians();

    let (x,y,z) = lla2ecef(lambda, phi, h, wgs84);
    if (x,y,z) != expected {
        assert!(false);
    }
}
