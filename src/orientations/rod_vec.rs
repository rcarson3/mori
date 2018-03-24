// This file is a part of the mori - Material Orientation Library in Rust
// Copyright 2018 Robert Carson <>

//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at

//        http://www.apache.org/licenses/LICENSE-2.0

//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
// limitations under the License.

use super::*;

#[derive(Clone, Debug)]
pub struct RodVec{
    pub ori: Array2<f64>,
}

#[derive(Clone, Debug)]
pub struct RodVecComp{
    ori: Array2<f64>,
}

impl RodVec{

    //Creates an array of zeros for the initial rodrigues vector parameterization when data is not fed into it
    pub fn new(size: usize) -> RodVec{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))) in {rodvec[3] = 1.0_f64});

        RodVec{
            ori,
        }
    }//End of new

    //Creates a rodrigues vector parameterization  type with the supplied data as long as the supplied data is in the following format
    //shape (4, nelems), memory order = fortran/column major.
    //If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> RodVec{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        RodVec{
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
}//End of Impl of RodVec

impl OriConv for RodVec{
    //Converts the rodrigues vector representation over to Bunge angles which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        //When a pure conversion doesn't exist we just use the already existing ones in other orientation
        //representations 
        rmat.to_bunge()   
    }//End of to_bunge


    //Converts the rodrigues vector representation over to a rotation matrix which has the following properties
    //shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{
        let ang_axis = self.to_ang_axis();
        //We could convert this to a pure converesion if we wanted to save on memory usage later on
        ang_axis.to_rmat()
    }//End of to_rmat

    //Converts the rodrigues vector representation over to axis-angle representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref rodvec (self.ori.axis_iter(Axis(1))) in {
            angaxis[0] = rodvec[0];
            angaxis[1] = rodvec[1];
            angaxis[2] = rodvec[2];
            angaxis[3] = 2.0_f64 * rodvec[3].atan();
        });

        AngAxis::new_init(ori)
    }//End of to_ang_axis

    //Converts the rodrigues vector representation over to a compact axial vector representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    //It returns a clone of itself.
    fn to_rod_vec(&self) -> RodVec{
        self.clone()
    }//End of to_rod_vec

    //Converts the rodrigues vector representation over to a compact rodrigues which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        //We first convert to a rodrigues vector representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        //If we want to be more efficient about this in the future with out as many copies used we can reuse a lot of the code
        //used in the to_ang_axis code. However, we will end up with a lot of similar/repeated code then. We could put that
        //code in a helper function that isn't seen.

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        azip!(mut rodvec_comp (ori.axis_iter_mut(Axis(1))), ref rodvec (self.ori.axis_iter(Axis(1))) in {
            rodvec_comp[0] = rodvec[0] * rodvec[3];
            rodvec_comp[1] = rodvec[1] * rodvec[3];
            rodvec_comp[2] = rodvec[2] * rodvec[3];
        });

        RodVecComp::new_init(ori)
    }//End of to_rod_vec_comp

    //Converts the rodrigues vector representation over to a unit quaternion which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{
        //Will replace this with a more direct conversion later on
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat()
    }//End of to_quat

    //Converts the rodrigues vector representation over to a homochoric representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric
}//End of impl Ori_Conv of Rodrigues Vector

impl RodVecComp{

    //Creates an array of zeros for the initial rodrigues vector parameterization when data is not fed into it
    pub fn new(size: usize) -> RodVecComp{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))) in {rodvec[3] = 1.0_f64});

        RodVecComp{
            ori,
        }
    }//End of new

    //Creates a rodrigues vector parameterization  type with the supplied data as long as the supplied data is in the following format
    //shape (4, nelems), memory order = fortran/column major.
    //If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> RodVecComp{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        RodVecComp{
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
}//End of Impl of RodVecComp
