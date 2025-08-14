use std::f64::consts::PI;

fn main() {
    // Calculate chord angle for 30 degrees
    let angle_30_rad = 30.0 * PI / 180.0;
    let chord_length_30 = 2.0 * (angle_30_rad / 2.0).sin();
    let chord_length2_30 = chord_length_30 * chord_length_30;
    
    // Calculate chord angle for 90 degrees (separation between centers)
    let angle_90_rad = 90.0 * PI / 180.0;
    let chord_length_90 = 2.0 * (angle_90_rad / 2.0).sin();
    let chord_length2_90 = chord_length_90 * chord_length_90;
    
    println\!("30° angle:");
    println\!("  radians: {}", angle_30_rad);
    println\!("  chord_length: {}", chord_length_30);
    println\!("  chord_length²: {}", chord_length2_30);
    
    println\!("\n90° angle:");
    println\!("  radians: {}", angle_90_rad);
    println\!("  chord_length: {}", chord_length_90);
    println\!("  chord_length²: {}", chord_length2_90);
    
    println\!("\nIntersection test:");
    println\!("  radius1 + radius2: {} + {} = {}", chord_length2_30, chord_length2_30, chord_length2_30 + chord_length2_30);
    println\!("  distance_between_centers: {}", chord_length2_90);
    println\!("  Should intersect? {} >= {} = {}", 
             chord_length2_30 + chord_length2_30, 
             chord_length2_90, 
             (chord_length2_30 + chord_length2_30) >= chord_length2_90);
}
