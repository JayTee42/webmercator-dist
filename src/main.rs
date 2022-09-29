use std::f64::consts::{E, PI};

const SEMI_MAJOR: f64 = 6378137.0;
const SEMI_MINOR: f64 = 6356752.314245;

fn spherical_mercator(lat: f64) -> f64 {
    SEMI_MAJOR * ((PI / 4.0) + (lat / 2.0)).tan().ln()
}

fn ellipsoidal_mercator(lat: f64) -> f64 {
    let e = (1.0 - ((SEMI_MINOR * SEMI_MINOR) / (SEMI_MAJOR * SEMI_MAJOR))).sqrt();
    let esin = e * lat.sin();
    let f = ((1.0 - esin) / (1.0 + esin)).powf(e / 2.0);

    SEMI_MAJOR * (f * ((PI / 4.0) + (lat / 2.0)).tan()).ln()
}

fn inverse_ellipsoidal_mercator(y: f64) -> f64 {
    let e = (1.0 - ((SEMI_MINOR * SEMI_MINOR) / (SEMI_MAJOR * SEMI_MAJOR))).sqrt();
    let e2 = e * e;
    let e4 = e2 * e2;
    let e6 = e2 * e4;
    let e8 = e4 * e4;

    let chi = (PI / 2.0) - (2.0 * E.powf(-y / SEMI_MAJOR).atan());

    let l2 = (e2 / 2.0) + ((5.0 * e4) / 24.0) + (e6 / 12.0) + ((13.0 * e8) / 360.0);
    let r2 = (2.0 * chi).sin();
    let l4 = ((7.0 * e4) / 48.0) + ((29.0 * e6) / 240.0) + ((811.0 * e8) / 11520.0);
    let r4 = (4.0 * chi).sin();
    let l6 = ((7.0 * e6) / 120.0) + ((81.0 * e8) / 1120.0);
    let r6 = (6.0 * chi).sin();
    let l8 = (4279.0 * e8) / 161280.0;
    let r8 = (8.0 * chi).sin();

    chi + (l2 * r2) + (l4 * r4) + (l6 * r6) + (l8 * r8)
}

fn meridian_distance(lat: f64) -> f64 {
    let n = (SEMI_MAJOR - SEMI_MINOR) / (SEMI_MAJOR + SEMI_MINOR);
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n2 * n2;

    let l0 = 1.0 + (n2 / 4.0) + (n4 / 64.0);
    let l2 = (-n / 2.0) + (n3 / 16.0);
    let r2 = (2.0 * lat).sin();
    let l4 = (-n2 / 16.0) + (n4 / 64.0);
    let r4 = (4.0 * lat).sin();
    let l6 = -n3 / 48.0;
    let r6 = (6.0 * lat).sin();
    let l8 = (-5.0 * n4) / 512.0;
    let r8 = (8.0 * lat).sin();
    let pre = (SEMI_MAJOR + SEMI_MINOR) / 2.0;
    let post = (2.0 * n * r2) / (1.0 + (2.0 * n * (2.0 * lat).cos()) + (n * n)).sqrt();

    pre * ((l0 * lat) - (l2 * r2) + (l4 * r4) - (l6 * r6) + (l8 * r8) - post)
}

fn main() {
    for lat_deg in 0..90 {
        let lat = (lat_deg as f64).to_radians();
        let y_spherical = spherical_mercator(lat);
        let y_ellipsoidal = ellipsoidal_mercator(lat);
        let lat_spherical = inverse_ellipsoidal_mercator(y_spherical);
        let dist = meridian_distance(lat_spherical) - meridian_distance(lat);

        println!(
            "At {lat_deg}°: map: {:.3} km, ground: {:.3} km",
            (y_spherical - y_ellipsoidal) / 1000.0,
            dist / 1000.0
        );
    }
}
