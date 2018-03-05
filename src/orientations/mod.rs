use super::*;

pub mod bunge;
pub mod rmat;
pub mod ang_axis;
pub mod rod_vec;
pub mod quat;

pub use self::bunge::*;
pub use self::rmat::*;
pub use self::ang_axis::*;
pub use self::rod_vec::*;
pub use self::quat::*;


#[derive(Clone, Debug)]
pub struct Homochoric{
    ori: Array2<f64>,
}

