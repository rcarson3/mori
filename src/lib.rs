#[macro_use]
extern crate ndarray;

use ndarray::prelude::*;

pub mod orientations{

    use super::*;

    pub struct Bunge{
        ori: Array2<f64>,
    }

    pub struct RMat{
        ori: Array3<f64>,
    } 

    pub struct AngAxis{
        ori: Array2<f64>,
    }

    pub struct AngAxisComp{
        ori: Array2<f64>
    }

    pub struct RodVec{
        ori: Array2<f64>,
    }

    pub struct RodVecComp{
        ori: Array2<f64>,
    }

    pub struct Quat{
        ori: Array2<f64>,
    }

    impl Bunge{

        //Creates an array of zeros for the initial Bunge angles when data is not fed into it
        pub fn new(size: usize) -> Bunge{
            assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

            let ori = Array2::<f64>::zeros((3, size).f());

            Bunge{
                ori,
            }
        }//End of new

        //Creates a Bunge type with the supplied data as long as the supplied data is in the following format
        //shape (3, nelems), memory order = fortran/column major.
        //If it doesn't fit those standards it will fail.
        pub fn new_init(ori: Array2<f64>) -> Bunge{

            let nrow = ori.rows();

            assert!(nrow == 3, "Number of rows of array was: {}, which is not equal to 3", nrow);
            //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
            {
                let strides = ori.strides();

                assert!(strides[0] == 1, "The memory stride is not column major (f order)");
            }

            Bunge{
                ori,
            }
        }//End of new_init

        //Just returns self if called
        pub fn to_bunge(self) -> Bunge{
            self
        }//End of to_bunge

        //Converts the Bunge angles over to a rotation matrix which has the following properties
        //shape (3, 3, nelems), memory order = fortran/column major.
        pub fn to_rmat(&self) -> RMat{

            let nelems = self.ori.len_of(Axis(1));
            
            let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

            //Individual cosines for each Bunge angle
            //Defining the types here so that the collect down below knows exactly what type they should be
            let mut cs:Array1<f64>;
            let mut ss:Array1<f64>;

            for i in 0..nelems{
                //Collecting the cosines and sines for a single Bunge angle
                cs = self.ori.slice(s![.., i]).iter().map(|&x| f64::cos(x)).collect();
                ss = self.ori.slice(s![.., i]).iter().map(|&x| f64::sin(x)).collect();

                ori[(0, 0, i)] = cs[0] * cs[2] - ss[0] * ss[2] * cs[1];
                ori[(0, 1, i)] = -cs[0] * ss[2] - ss[0] * cs[1] * cs[2];
                ori[(0, 2, i)] = ss[0] * ss[1];

                ori[(1, 0, i)] = ss[0] * cs[2] + cs[0] * cs[1] * ss[2];
                ori[(1, 1, i)] = -ss[0] * ss[2] + cs[0] * cs[1] * cs[2];
                ori[(1, 2, i)] = -cs[0] * ss[1];

                ori[(2, 0, i)] = ss[1] * ss[2];
                ori[(2, 1, i)] = ss[1] * cs[2];
                ori[(2, 2, i)] = cs[1];
            }

            RMat{
                ori,
            }
        }//End of to_rmat

        //Converts the Bunge angles over to an angle-axis representation which has the following properties
        //shape (4, nelems), memory order = fortran/column major.
        pub fn to_ang_axis(&self) -> AngAxis{

            let nelems = self.ori.len_of(Axis(1));
            
            let mut ori = Array2::<f64>::zeros((4, nelems).f());

            //Going ahead and defining the types for these loop variables
            //This should help and ensure that they are what types we expect them to be
            //down below even though the compiler should be able to infer them correctly.
            let mut t:f64;
            let mut sigma:f64;
            let mut delta:f64;
            let mut tau:f64;
            let mut alpha:f64;
            let mut itau:f64;
            let mut p:f64;
            //This is a constant factor used in the loop.
            let inv2 = 1.0_f64/2.0_f64;


            for i in 0..nelems{

                t       = f64::tan(self.ori[(1, i)]/2.0_f64);
                sigma   = inv2 * (self.ori[(0, i)] + self.ori[(2, i)]);
                delta   = inv2 * (self.ori[(0, i)] - self.ori[(2, i)]);
                tau     = f64::sqrt(t * t + f64::sin(sigma) * f64::sin(sigma));
                alpha   = 2.0_f64 * f64::atan(tau / (f64::cos(sigma)));
                itau = 1.0_f64/tau;

                //If alpha is greater than pi we need to set p equal to 1 and
                //set alpha equal to 2*pi - alpha. If it isn't then p is set
                //to -1. Afterwards everything else is the same.
                if alpha > std::f64::consts::PI{
                    p = 1.0_f64;
                    alpha = 2.0_f64 * std::f64::consts::PI - alpha;
                }else{
                    p = -1.0_f64;
                }

                ori[(0, i)] = p * itau * t * f64::cos(delta);
                ori[(1, i)] = p * itau * t * f64::sin(delta);
                ori[(2, i)] = p * itau * t * f64::cos(sigma);
                ori[(3, i)] = alpha;
            }

            AngAxis{
                ori,
            }
        }//End of to_ang_axis

        //Converts the Bunge angles over to a compact angle-axis representation which has the following properties
        //shape (3, nelems), memory order = fortran/column major.
        pub fn to_ang_axis_comp(&self) -> AngAxisComp{
            //We first convert to a angle axis representation. Then we scale our normal vector by our the rotation
            //angle which is the fourth component of our angle axis vector.
            let ang_axis = self.to_ang_axis();

            let nelems = ang_axis.ori.len_of(Axis(1));

            let mut ori = Array2::<f64>::zeros((3, nelems).f());

            for i in 0 .. nelems{
                ori[(0, i)] = ang_axis.ori[(0, i)] * ang_axis.ori[(3, i)];
                ori[(1, i)] = ang_axis.ori[(1, i)] * ang_axis.ori[(3, i)];
                ori[(2, i)] = ang_axis.ori[(2, i)] * ang_axis.ori[(3, i)];
            }

            AngAxisComp{
                ori,
            } 
        }//End of to_ang_axis_comp

        //Converts the Bunge angles over to a rodrigues vector representation which has the following properties
        //shape (4, nelems), memory order = fortran/column major.
        pub fn to_rod_vec(&self) -> RodVec{
            //We first convert to a angle axis representation. Then we just need to change the last component
            //of our angle axis representation to be tan(phi/2) instead of phi
            let ang_axis = self.to_ang_axis();

            let nelems = ang_axis.ori.len_of(Axis(1));

            let mut ori = Array2::<f64>::zeros((4, nelems).f());

            let inv2 = 1.0_f64/2.0_f64;


            for i in 0 .. nelems{
                ori[(0, i)] = ang_axis.ori[(0, i)];
                ori[(1, i)] = ang_axis.ori[(1, i)];
                ori[(2, i)] = ang_axis.ori[(2, i)];
                ori[(3, i)] = f64::tan(inv2 * ang_axis.ori[(3, i)]);
            }

            RodVec{
                ori,
            }
        }//End of to_rod_vec

        //Converts the Bunge angles over to a compact rodrigues vector representation which has the following properties
        //shape (3, nelems), memory order = fortran/column major.
        pub fn to_rod_vec_comp(&self) -> RodVecComp{
            //We first convert to a rodrigues vector representation. Then we scale our normal vector by our the rotation
            //angle which is the fourth component of our angle axis vector.
            //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
            //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
            //code in a helper function that isn't seen.
            let rod_vec = self.to_rod_vec();

            let nelems = rod_vec.ori.len_of(Axis(1));

            let mut ori = Array2::<f64>::zeros((3, nelems).f());

            for i in 0 .. nelems{
                ori[(0, i)] = rod_vec.ori[(0, i)] * rod_vec.ori[(3, i)];
                ori[(1, i)] = rod_vec.ori[(1, i)] * rod_vec.ori[(3, i)];
                ori[(2, i)] = rod_vec.ori[(2, i)] * rod_vec.ori[(3, i)];
            }

            RodVecComp{
                ori,
            }
        }//End of to_rod_vec_comp

        //Converts the Bunge angles over to a unit quaternion representation which has the following properties
        //shape (4, nelems), memory order = fortran/column major.
        pub fn to_quat(&self) -> Quat{

            let nelems = self.ori.len_of(Axis(1));
            
            let mut ori = Array2::<f64>::zeros((4, nelems).f());

            //Going ahead and defining the types for these loop variables
            //This should help and ensure that they are what types we expect them to be
            //down below even though the compiler should be able to infer them correctly.
            let mut sigma:f64;
            let mut delta:f64;
            let mut c:f64;
            let mut s:f64;
            let mut q0:f64;
            let mut p:f64;
            //This is a constant factor used in the loop.
            let inv2 = 1.0_f64/2.0_f64;

            for i in 0..nelems{
                sigma   = inv2 * (self.ori[(0, i)] + self.ori[(2, i)]);
                delta   = inv2 * (self.ori[(0, i)] - self.ori[(2, i)]);
                c       = f64::cos(inv2 * self.ori[(1, i)]);
                s       = f64::sin(inv2 * self.ori[(1, i)]);
                q0      = c * f64::cos(sigma);

                if q0 < 0.0_f64{
                    p = -1.0_f64;
                }else{
                    p = 1.0_f64;
                }

                ori[(0, i)] = p * q0;
                ori[(1, i)] = -p * s * f64::cos(delta);
                ori[(2, i)] = -p * s * f64::sin(delta);
                ori[(3, i)] = -p * c * f64::cos(sigma);
            }

            Quat{
                ori,
            }            
        }//End of to_quat
    }//End of Impl Bunge

    impl RMat{

        //Creates a series of identity matrices for the initial rotation matrix when the data is not fed into it
        pub fn new(size: usize) -> RMat{
            assert!(size > 0, "Size inputted: {}, was not greater than 0", size);
            let mut ori = Array3::<f64>::zeros((3, 3, size).f());
            let eye = Array2::<f64>::eye(3);

            for i in 0..size{
                ori.slice_mut(s![.., .., i]).assign(&eye);
            }

            RMat{
                ori,
            }
        }//End of new

        //Creates a rotation matrix type with the supplied data as long as the supplied data is in the following format
        //shape (3, 3, nelems), memory order = fortran/column major.
        //If it doesn't fit those standards it will fail.
        //The data is also being assummed to be orthogonal, but it's currently not checked to see if that is the case.
        //I might add a check later on to make sure that's the case.
        pub fn new_init(ori: Array3<f64>) -> RMat{

            let nrow = ori.len_of(Axis(0));
            let ncol = ori.len_of(Axis(1));

            assert!(nrow == 3 && ncol == 3, "Number of rows and cols of array was: {:?}, which is not equal to 3", [nrow, ncol]);
            //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
            {
                let strides = ori.strides();

                assert!(strides[0] == 1, "The memory stride is not column major (f order)");
            }

            RMat{
                ori,
            }
        }//End of new_init

        pub fn to_bunge(&self) -> Bunge{
            
            let nelems = self.ori.len_of(Axis(2));

            let mut ori = Array2::<f64>::zeros((3, nelems));

            //Temporary variable that's most likely going to be used often
            let mut eta:f64;
            //We need to check the R_33 component to see if it's near 1.0 
            let tol = std::f64::EPSILON;

            //Eventually, I might split this up into two arrays where the conditional is avoided.
            for i in 0..nelems{

                if f64::abs(self.ori[(2, 2, i)] - 1.0_f64) < tol{

                    ori[(0, i)] = f64::atan2(self.ori[(0, 1, i)], self.ori[(0, 0, i)]);
                    ori[(1, i)] = std::f64::consts::FRAC_PI_2 * (1.0_f64 - self.ori[(2, 2, 0)]);
                    ori[(2, i)] = 0.0_f64;

                }else{

                    eta         = 1.0_f64 / f64::sqrt(1.0_f64 - self.ori[(2, 2, i)] * self.ori[(2, 2, i)]);
                    ori[(0, i)] = f64::atan2(self.ori[(2, 0, i)] * eta, -self.ori[(2, 1, i)] * eta);
                    ori[(1, i)] = f64::acos(self.ori[(2, 2, i)]);
                    ori[(2, i)] = f64::atan2(self.ori[(0, 2, i)] * eta, self.ori[(1, 2, i)] * eta);

                }
            }

            Bunge{
                ori,
            }
        }//End of to_bunge

        //If to_rmat is called it just returns self
        pub fn to_rmat(self) -> RMat{
            self
        }//End of to_rmat

        //Converts the rotation matrix over to an angle-axis representation which has the following properties
        //shape (4, nelems), memory order = fortran/column major.
        pub fn to_ang_axis(&self) -> AngAxis{

            let nelems = self.ori.len_of(Axis(2));

            let mut ori = Array2::<f64>::zeros((4, nelems));

            //This is the angle of rotation about our normal axis.
            let mut phi:f64;
            //This is going to be the tr(RMat)
            let mut trR:f64;
            //This is going to be the inverse of sin phi
            let mut inv_sin:f64;
            let inv2 = 1.0_f64/2.0_f64;

            for i in 0..nelems{

                trR = self.ori[(0, 0, i)] + self.ori[(1, 1, i)] + self.ori[(2, 2, i)];
                phi = f64::acos(inv2 * (trR - 1.0_f64));

                if f64::abs(phi) < std::f64::EPSILON{

                    ori[(0, i)] = 0.0_f64;
                    ori[(1, i)] = 0.0_f64;
                    ori[(2, i)] = 1.0_f64;
                    ori[(3, i)] = 0.0_f64;

                }else{
                    inv_sin = 1.0_f64 / f64::sin(phi);
                    //The first three terms are the axial vector of RMat times (1/(2*sin(phi)))
                    ori[(0, i)] = inv_sin * inv2 * (self.ori[(2, 1, i)] - self.ori[(1, 2, i)]);
                    ori[(1, i)] = inv_sin * inv2 * (self.ori[(0, 2, i)] - self.ori[(2, 0, i)]);
                    ori[(2, i)] = inv_sin * inv2 * (self.ori[(1, 0, i)] - self.ori[(0, 1, i)]);
                    ori[(3, i)] = phi;

                }
            }

            AngAxis{
                ori,
            }
        }//End of to_ang_axis

        //Converts the rotation matrix over to a compact angle-axis representation which has the following properties
        //shape (3, nelems), memory order = fortran/column major.
        pub fn to_ang_axis_comp(&self) -> AngAxisComp{
            //We first convert to a angle axis representation. Then we scale our normal vector by our the rotation
            //angle which is the fourth component of our angle axis vector.
            let ang_axis = self.to_ang_axis();

            let nelems = ang_axis.ori.len_of(Axis(1));

            let mut ori = Array2::<f64>::zeros((3, nelems).f());

            for i in 0 .. nelems{
                ori[(0, i)] = ang_axis.ori[(0, i)] * ang_axis.ori[(3, i)];
                ori[(1, i)] = ang_axis.ori[(1, i)] * ang_axis.ori[(3, i)];
                ori[(2, i)] = ang_axis.ori[(2, i)] * ang_axis.ori[(3, i)];
            }

            AngAxisComp{
                ori,
            }
        }//End of to_ang_axis_comp

        //Converts the rotation matrix over to a rodrigues vector representation which has the following properties
        //shape (4, nelems), memory order = fortran/column major.
        pub fn to_rod_vec(&self) -> RodVec{
            //We first convert to a angle axis representation. Then we just need to change the last component
            //of our angle axis representation to be tan(phi/2) instead of phi
            let ang_axis = self.to_ang_axis();

            let nelems = ang_axis.ori.len_of(Axis(1));

            let mut ori = Array2::<f64>::zeros((4, nelems).f());

            let inv2 = 1.0_f64/2.0_f64;


            for i in 0 .. nelems{
                ori[(0, i)] = ang_axis.ori[(0, i)];
                ori[(1, i)] = ang_axis.ori[(1, i)];
                ori[(2, i)] = ang_axis.ori[(2, i)];
                ori[(3, i)] = f64::tan(inv2 * ang_axis.ori[(3, i)]);
            }

            RodVec{
                ori,
            }
        }//End of to_rod_vec

        //Converts the rotation matrix over to a compact rodrigues vector representation which has the following properties
        //shape (3, nelems), memory order = fortran/column major.
        pub fn to_rod_vec_comp(&self) -> RodVecComp{
            //We first convert to a rodrigues vector representation. Then we scale our normal vector by our the rotation
            //angle which is the fourth component of our angle axis vector.
            //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
            //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
            //code in a helper function that isn't seen.
            let rod_vec = self.to_rod_vec();

            let nelems = rod_vec.ori.len_of(Axis(1));

            let mut ori = Array2::<f64>::zeros((3, nelems).f());

            for i in 0 .. nelems{
                ori[(0, i)] = rod_vec.ori[(0, i)] * rod_vec.ori[(3, i)];
                ori[(1, i)] = rod_vec.ori[(1, i)] * rod_vec.ori[(3, i)];
                ori[(2, i)] = rod_vec.ori[(2, i)] * rod_vec.ori[(3, i)];
            }

            RodVecComp{
                ori,
            }
        }//End of to_rod_vec_comp

        //Converts the rotation matrix over to a unit quaternion representation which has the following properties
        //shape (4, nelems), memory order = fortran/column major.
        pub fn to_quat(&self) -> Quat{

            let nelems = self.ori.len_of(Axis(2));

            let mut ori = Array2::<f64>::zeros((4, nelems));

            let inv2 = 1.0_f64/2.0_f64;

            //We're creating individual variables for each component of the unit quaternion
            let mut q0:f64;
            let mut q1:f64;
            let mut q2:f64;
            let mut q3:f64;
            let mut inv_norm:f64;

            for i in 0..nelems{

                q0 = inv2 * f64::sqrt(1.0_f64 + self.ori[(0, 0, i)] + self.ori[(1, 1, i)] + self.ori[(2, 2, i)]);
                q1 = inv2 * f64::sqrt(1.0_f64 + self.ori[(0, 0, i)] - self.ori[(1, 1, i)] - self.ori[(2, 2, i)]);
                q2 = inv2 * f64::sqrt(1.0_f64 - self.ori[(0, 0, i)] + self.ori[(1, 1, i)] - self.ori[(2, 2, i)]);
                q3 = inv2 * f64::sqrt(1.0_f64 - self.ori[(0, 0, i)] - self.ori[(1, 1, i)] + self.ori[(2, 2, i)]);

                if self.ori[(1, 2, i)] > self.ori[(2, 1, i)]{
                    q1 = -q1;
                }
                if self.ori[(2, 0, i)] > self.ori[(0, 2, i)]{
                    q2 = -q2;
                }
                if self.ori[(0, 1, i)] > self.ori[(1, 0, i)]{
                    q3 = -q3;
                }

                //We need to normalize our quaternion when we store it.
                inv_norm = 1.0_f64 / f64::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);

                ori[(0, i)] = q0 * inv_norm;
                ori[(1, i)] = q1 * inv_norm;
                ori[(2, i)] = q2 * inv_norm;
                ori[(3, i)] = q3 * inv_norm;

            }

            Quat{
                ori,
            }
        }//End of to_quat
    }//End of Impl RMat


}
