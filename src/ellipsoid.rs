/// Ellipsoid parameters.
///
#[derive(Copy, Clone)]
pub struct Ellipsoid {
    /// equatorial axis (m)
    pub a : f64,
    /// flattening 
    pub f : f64, 
    /// polar sxis (m)
    pub b : f64, 
    /// first eccentricity
    pub e1 : f64, 
    /// second eccentricity
    pub e2 : f64, 

    pub asq : f64,
    pub bsq : f64,
    pub eccsq : f64,
    pub ecc : f64
}

/// Construct an Ellipsoid from parameters a (equatorial axis) and f (inverse flattening).
///
impl Ellipsoid {

    pub fn new(a:f64, invf:f64) -> Ellipsoid {
        let f:f64 = invf.recip();
        let b:f64 = a*(1f64-f);
        let e1:f64 = ((a*a - b*b)/(a*a)).sqrt();
        let e2:f64 = ((a*a - b*b)/(b*b)).sqrt();

        let asq = a.powi(2);
        let bsq = b.powi(2);
        let eccsq = 1.0 - bsq/asq;
        let ecc = eccsq.sqrt();

        Ellipsoid {a: a, f: f, b: b, e1: e1, e2: e2, asq: asq, bsq: bsq, eccsq: eccsq, ecc: ecc}
    }

    /// Calculate ellipsoid radii (geocentric, normal to ellipsoid, and meridional) at a latitude, phi, given in radians.
    ///
    pub fn radii(self, phi:f64) -> (f64, f64, f64) {
        // ellps.a / (1f64 - ellps.e1.powi(2) * phi.sin().powi(2)).sqrt()
        let clat  =  phi.cos();
        let slat  =  phi.sin();

        let dsq   =  1.0 - self.eccsq * slat * slat;
        let d     =  dsq.sqrt();

        let rn    =  self.a/d;
        let rm    =  rn * (1.0 - self.eccsq ) / dsq;

        let rho   =  rn * clat;
        let z     =  (1.0 - self.eccsq ) * rn * slat;
        let rsq   =  rho*rho + z*z;
        let r     =  rsq.sqrt();

        (r, rn, rm)
    }
    /// Calculate height above ellipsoid.
    ///
    pub fn hae(self, p:f64, rn:f64, phi:f64) -> f64 {
        p/phi.cos() - rn
    }

    /// Calculate height above ellipsoid.
    ///
    pub fn hae2(self, x:f64, y:f64, z:f64, r:f64) -> f64 {
    (x*x + y*y + z*z).sqrt() - r
    }    
}
