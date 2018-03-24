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

        azip!(mut quat (ori.axis_iter_mut(Axis(1))) in {quat[0] = 1.0_f64});

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

    //Return a view of ori
    pub fn ori_view(&self) -> ArrayView2<f64>{
        self.ori.view()
    }

    //Return a mutable view of ori
    pub fn ori_view_mut(&mut self) -> ArrayViewMut2<f64>{
        self.ori.view_mut()
    }
}//End of Impl of Quat

impl OriConv for Quat{
    //Converts the unit quaternion representation over to Bunge angles which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut bunge (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let q03 = quat[0] * quat[0] + quat[3] * quat[3];
            let q12 = quat[1] * quat[1] + quat[2] * quat[2];
            let xi = f64::sqrt(q03 * q12);
            //We get to now go through all of the different cases that this might break down into
            if xi.abs() < tol && q12.abs() < tol {
                bunge[0] = f64::atan2(-2.0_f64 * quat[0] * quat[3], quat[0] * quat[0] - quat[3] * quat[3]);
                //All of the other values are zero
            }else if xi.abs() < tol && q03.abs() < tol{
                bunge[0] = f64::atan2(2.0_f64 * quat[1] * quat[2], quat[1] * quat[1] - quat[2] * quat[2]);
                bunge[1] = std::f64::consts::PI;
                //The other value is zero
            }else{
                let inv_xi = 1.0_f64 / xi;
                //The atan2 terms are pretty long so we're breaking it down into a couple of temp terms
                let t1 = inv_xi * (quat[1] * quat[3] - quat[0] * quat[2]);
                let t2 = inv_xi * (-quat[0] * quat[1] - quat[2] * quat[3]);
                //We can now assign the first two bunge angles
                bunge[0] = t1.atan2(t2);
                bunge[1] = f64::atan2(2.0_f64 * xi, q03 - q12);
                //Once again these terms going into the atan2 term are pretty long
                let t1 = inv_xi * (quat[0] * quat[2] - quat[1] * quat[3]);
                let t2 = inv_xi * (quat[2] * quat[3] - quat[0] * quat[1]);
                //We can finally find the final bunge angle
                bunge[2] = t1.atan2(t2);
            }
        });

        Bunge::new_init(ori)
    }//End of to_bunge

    //Converts the unit quaternion representation over to rotation matrix which has the following properties
    //shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref quat (self.ori.axis_iter(Axis(1))) in {
            let qbar =  quat[0] * quat[0] - (quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

            rmat[[0, 0]] = qbar + quat[1] * quat[1];
            rmat[[1, 0]] = 2.0_f64 * (quat[1] * quat[2] + quat[0] * quat[3]);
            rmat[[2, 0]] = 2.0_f64 * (quat[1] * quat[3] - quat[0] * quat[2]);

            rmat[[1, 1]] = 2.0_f64 * (quat[1] * quat[2] - quat[0] * quat[3]);
            rmat[[1, 1]] = qbar + quat[2] * quat[2];
            rmat[[2, 1]] = 2.0_f64 * (quat[2] * quat[3] + quat[0] * quat[1]);

            rmat[[0, 2]] = 2.0_f64 * (quat[1] * quat[3] + quat[0] * quat[2]);
            rmat[[1, 2]] = 2.0_f64 * (quat[2] * quat[3] - quat[0] * quat[1]);
            rmat[[2, 2]] = qbar + quat[3] * quat[3];
        });

        RMat::new_init(ori)
    }//End of to_rmat

    //Converts the unit quaternion representation over to angle-axis representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems));

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref quat (self.ori.axis_iter(Axis(1))) in {
            if quat[0].abs() < tol{
                angaxis[0] = quat[1];
                angaxis[1] = quat[2];
                angaxis[2] = quat[3];
                angaxis[3] = std::f64::consts::PI;
            }else{
                //This is our angle of rotation
                let phi = 2.0_f64 * quat[0].acos();
                let s   = quat[0].signum() / f64::sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);

                angaxis[0] = s * quat[1];
                angaxis[1] = s * quat[2];
                angaxis[2] = s * quat[3];
                angaxis[3] = phi;
            }
        });

        AngAxis::new_init(ori)
    }//End of to_ang_axis

    //Converts the unit quaternion over to a compact angle-axis representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        //We first convert to a angle axis representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    //Converts the unit quaternion over to a rodrigues vector representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{
        //We first convert to a angle axis representation. Then we just need to change the last component
        //of our angle axis representation to be tan(phi/2) instead of phi
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rod_vec()
    }//End of to_rod_vec

    //Converts the unit quaternion over to a compact rodrigues vector representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        let rod_vec = self.to_rod_vec();
        rod_vec.to_rod_vec_comp()
    }//End of to_rod_vec_comp

    //This returns a clone of the original Quaternion structure
    fn to_quat(&self) -> Quat{
        self.clone()
    }//End of to_quat

    //Converts the quaternion representation over to a homochoric representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric
}//End of impl of unit Quaternion