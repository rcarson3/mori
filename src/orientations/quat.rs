use super::*;

#[derive(Clone, Debug)]
pub struct Quat{
    pub ori: Array2<f64>,
}

impl Quat{
    //Creates an array of zeros for the initial unit quaternion when data is not fed into it
    pub fn new(size: usize) -> Quat{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        ori.slice_mut(s![0, ..]).mapv_inplace(|_x| 1.0_f64);

        Quat{
            ori,
        }
    }//End of new

    //Creates a unit quaternion type with the supplied data as long as the supplied data is in the following format
    //shape (4, nelems), memory order = fortran/column major.
    //If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> Quat{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        Quat{
            ori,
        }
    }//End of new_init

    //Converts the unit quaternion representation over to Bunge angles which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_bunge(&self) -> Bunge{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let mut xi:f64;
        let mut q03:f64;
        let mut q12:f64;
        let mut inv_xi:f64;

        let mut t1:f64;
        let mut t2:f64;

        let tol = std::f64::EPSILON;

        for i in 0..nelems{

            q03 = self.ori[(0, i)] * self.ori[(0, i)] + self.ori[(3, i)] * self.ori[(3, i)];
            q12 = self.ori[(1, i)] * self.ori[(1, i)] + self.ori[(2, i)] * self.ori[(2, i)];
            xi = f64::sqrt(q03 * q12);

            if f64::abs(xi) < tol && f64::abs(q12) < tol {

                ori[(0, i)] = f64::atan2(-2.0_f64 * self.ori[(0, i)] * self.ori[(3, i)], self.ori[(0, i)] * self.ori[(0, i)] - self.ori[(3, i)] * self.ori[(3, i)]);
                //All of the other values are zero

            }else if f64::abs(xi) < tol && f64::abs(q03) < tol{

                ori[(0, i)] = f64::atan2(2.0_f64 * self.ori[(1, i)] * self.ori[(2, i)], self.ori[(1, i)] * self.ori[(1, i)] - self.ori[(2, i)] * self.ori[(2, i)]);
                ori[(1, i)] = std::f64::consts::PI;
                //The other value is zero

            }else{

                inv_xi      = 1.0_f64 / xi;
                t1          = inv_xi * (self.ori[(1, i)] * self.ori[(3, i)] - self.ori[(0, i)] * self.ori[(2, i)] );
                t2          = inv_xi * (-self.ori[(0, i)] * self.ori[(1, i)] - self.ori[(2, i)] * self.ori[(3, i)] );

                ori[(0, i)] = f64::atan2(t1, t2);
                ori[(1, i)] = f64::atan2(2.0_f64 * xi, q03 - q12);

                t1          = inv_xi * (self.ori[(0, i)] * self.ori[(2, i)] - self.ori[(1, i)] * self.ori[(3, i)] );
                t2          = inv_xi * (self.ori[(2, i)] * self.ori[(3, i)] - self.ori[(0, i)] * self.ori[(1, i)] );

                ori[(2, i)] = f64::atan2(t1, t2);

            }
        }

        Bunge{
            ori,
        }
    }//End of to_Bunge

    //Converts the unit quaternion representation over to rotation matrix which has the following properties
    //shape (3, 3, nelems), memory order = fortran/column major.
    pub fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        let mut qbar:f64;

        for i in 0..nelems{

            qbar =  self.ori[(0, i)] * self.ori[(0, i)] - (self.ori[(1, i)] * self.ori[(1, i)] + self.ori[(2, i)] * self.ori[(2, i)] + self.ori[(3, i)] * self.ori[(3, i)]);

            ori[(0, 0, i)] = qbar + self.ori[(1, i)] * self.ori[(1, i)];
            ori[(1, 0, i)] = 2.0_f64 * (self.ori[(1, i)] * self.ori[(2, i)] + self.ori[(0, i)] * self.ori[(3, i)]);
            ori[(2, 0, i)] = 2.0_f64 * (self.ori[(1, i)] * self.ori[(3, i)] - self.ori[(0, i)] * self.ori[(2, i)]);

            ori[(1, 1, i)] = 2.0_f64 * (self.ori[(1, i)] * self.ori[(2, i)] - self.ori[(0, i)] * self.ori[(3, i)]);
            ori[(1, 1, i)] = qbar + self.ori[(2, i)] * self.ori[(2, i)];
            ori[(2, 1, i)] = 2.0_f64 * (self.ori[(2, i)] * self.ori[(3, i)] + self.ori[(0, i)] * self.ori[(1, i)]);

            ori[(0, 2, i)] = 2.0_f64 * (self.ori[(1, i)] * self.ori[(3, i)] + self.ori[(0, i)] * self.ori[(2, i)]);
            ori[(1, 2, i)] = 2.0_f64 * (self.ori[(2, i)] * self.ori[(3, i)] - self.ori[(0, i)] * self.ori[(1, i)]);
            ori[(2, 2, i)] = qbar + self.ori[(3, i)] * self.ori[(3, i)];

        }

        RMat{
            ori,
        }
    }//End of to_rmat

    //Converts the unit quaternion representation over to angle-axis representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems));

        let tol = std::f64::EPSILON;
        let mut inv_t:f64;
        let mut s:f64;
        let mut phi:f64;

        for i in 0..nelems{

            if f64::abs(self.ori[(0, i)]) < tol{

                ori[(0, i)] = self.ori[(1, i)];
                ori[(1, i)] = self.ori[(2, i)];
                ori[(2, i)] = self.ori[(3, i)];
                ori[(3, i)] = std::f64::consts::PI;

            }else{

                phi = 2.0_f64 * f64::acos(self.ori[(0, i)]);

                inv_t = f64::sqrt(self.ori[(1, i)] * self.ori[(1, i)] + self.ori[(2, i)] * self.ori[(2, i)] + self.ori[(3, i)] * self.ori[(3, i)]);
                s   = self.ori[(0, i)].signum() / inv_t;

                ori[(0, i)] = s * self.ori[(1, i)];
                ori[(1, i)] = s * self.ori[(2, i)];
                ori[(2, i)] = s * self.ori[(3, i)];
                ori[(3, i)] = phi;

            }
        }

        AngAxis{
            ori,
        }
    }

    //Converts the unit quaternion over to a compact angle-axis representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_ang_axis_comp(&self) -> AngAxisComp{
        //We first convert to a angle axis representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    //Converts the unit quaternion over to a rodrigues vector representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_rod_vec(&self) -> RodVec{
        //We first convert to a angle axis representation. Then we just need to change the last component
        //of our angle axis representation to be tan(phi/2) instead of phi
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rod_vec()
    }//End of to_rod_vec

    //Converts the unit quaternion over to a compact rodrigues vector representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
        //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
        //code in a helper function that isn't seen.
        let rod_vec = self.to_rod_vec();
        rod_vec.to_rod_vec_comp()
    }//End of to_rod_vec_comp

    //This returns a clone of itself
    pub fn to_quat(&self) -> Quat{
        self.clone()
    }

    pub fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }


}//End of impl of unit Quaternion