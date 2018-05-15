// This file is a part of the mori - Material Orientation Library in Rust
// Copyright 2018 Robert Carson 
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

use super::*;

///A structure that holds an array of compact Rodrigues vectors
#[derive(Clone, Debug)]
pub struct RodVecComp{
    ori: Array2<f64>,
}

impl RodVecComp{

    ///Creates an array of zeros for the initial compact Rodrigues vector parameterization when data is not fed into it
    pub fn new(size: usize) -> RodVecComp{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let ori = Array2::<f64>::zeros((3, size).f());

        RodVecComp{
            ori,
        }
    }//End of new

    ///Creates a compact Rodrigues vector parameterization type with the supplied data as long as the supplied data is in the following format
    ///shape (4, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> RodVecComp{

        let nrow = ori.rows();

        assert!(nrow == 3, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        RodVecComp{
            ori,
        }
    }//End of new_init

    ///Return a ndarray view of the orientation data
    pub fn ori_view(&self) -> ArrayView2<f64>{
        self.ori.view()
    }

    ///Return a ndarray mutable view of the orientation data
    pub fn ori_view_mut(&mut self) -> ArrayViewMut2<f64>{
        self.ori.view_mut()
    }

    ///Returns a new RodVecComp that is equal to the equivalent of transposing a rotation matrix.
    ///It turns out this is simply the negative of the normal vector due to the vector being formed
    ///from an axial vector of the rotation matrix --> Rmat\^T = -Rx where Rx is the axial vector.
    pub fn transpose(&self) -> RodVecComp{
        let nelems = self.ori.len_of(Axis(1));
        let mut ori = Array2::<f64>::zeros((3, nelems).f());
        ori.assign(&(-1.0 * &self.ori));

        RodVecComp::new_init(ori)
    }

    ///Performs the equivalent of transposing a rotation matrix on the internal orientations.
    ///It turns out this is simply the negative of the normal vector due to the vector being formed
    ///from an axial vector of the rotation matrix --> Rmat\^T = -Rx where Rx is the axial vector.
    pub fn transpose_inplace(&mut self){
        self.ori.mapv_inplace(|x| {-1.0_f64 * x});
    }
}//End of Impl of RodVecComp

///The orientation conversions of a series of compact Rodrigues vectors to a number of varying different orientation
///representations commonly used in material orientation processing. It should be noted that the compact
///Rodrigues vector is much more likely to be numerically unstable compared to the other orientation representations.
///This is due to the fact that the tan(phi\2) term is included into the compacted vector. 
impl OriConv for RodVecComp{
    ///Converts the compact Rodrigues vector representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        //When a pure conversion doesn't exist we just use the already existing ones in other orientation
        //representations 
        rmat.to_bunge()   
    }//End of to_bunge

    ///Converts the compact Rodrigues vector representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{
        let ang_axis = self.to_ang_axis();
        //We could convert this to a pure converesion if we wanted to save on memory usage later on
        ang_axis.to_rmat()
    }//End of to_rmat

    ///Converts the compact Rodrigues vector representation over to axis-angle representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{
        let rod_vec = self.to_rod_vec();
        rod_vec.to_ang_axis()
    }//End of to_ang_axis

    ///Converts the compact Rodrigues vector representation over to a compact axial vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp()
    }//End of to_ang_axis_comp

    ///Converts the compact Rodrigues vector representation over to a Rodrigues which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref rodvec_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_rodvec = f64::sqrt({
                rodvec_comp[0] * rodvec_comp[0] 
                + rodvec_comp[1] * rodvec_comp[1] 
                + rodvec_comp[2] * rodvec_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_rodvec.abs() < tol {
                rodvec[2] = 1.0_f64;
            }else if norm_rodvec == std::f64::INFINITY {
                rodvec[3] = norm_rodvec;
            }else {
                let inv_norm_rodvec = 1.0_f64 / norm_rodvec;
                rodvec[0] = rodvec_comp[0] * inv_norm_rodvec;
                rodvec[1] = rodvec_comp[1] * inv_norm_rodvec;
                rodvec[2] = rodvec_comp[2] * inv_norm_rodvec;
                rodvec[3] = norm_rodvec;
            }
        });

        RodVec::new_init(ori)
    }//End of to_rod_vec

    ///Returns a clone of the compact Rodrigues vector structure
    fn to_rod_vec_comp(&self) -> RodVecComp{
        self.clone()
    }//End of to_rod_vec_comp

    ///Converts the compact Rodrigues vector representation over to a unit quaternion which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{
        //Will replace this with a more direct conversion later on
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat()
    }//End of to_quat

    ///Converts the compact Rodrigues vector representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) -> Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric

    ///Converts the compact Rodrigues vector representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let rmat = self.to_rmat();
        rmat.to_bunge_inplace(bunge);
    }

    ///Converts the compact Rodrigues vector representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rmat_inplace(&self, rmat: &mut RMat){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_rmat_inplace(rmat);
    }

    ///Converts the compact Rodrigues vector representation over to axis-angle representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){
        let rod_vec = self.to_rod_vec();
        rod_vec.to_ang_axis_inplace(ang_axis);
    }

    ///Converts the compact Rodrigues vector representation over to compact axis-angle representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_comp_inplace(&self, ang_axis_comp: &mut AngAxisComp){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_ang_axis_comp_inplace(ang_axis_comp);
    }

    ///Converts the compact Rodrigues vector representation over to a Rodrigues which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_inplace(&self, rod_vec: &mut RodVec){
        let mut ori = rod_vec.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);


        let tol = std::f64::EPSILON;

        azip!(mut rodvec (ori.axis_iter_mut(Axis(1))), ref rodvec_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_rodvec = f64::sqrt({
                rodvec_comp[0] * rodvec_comp[0] 
                + rodvec_comp[1] * rodvec_comp[1] 
                + rodvec_comp[2] * rodvec_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_rodvec.abs() < tol {
                rodvec[2] = 1.0_f64;
            }else if norm_rodvec == std::f64::INFINITY {
                rodvec[3] = norm_rodvec;
            }else {
                let inv_norm_rodvec = 1.0_f64 / norm_rodvec;
                rodvec[0] = rodvec_comp[0] * inv_norm_rodvec;
                rodvec[1] = rodvec_comp[1] * inv_norm_rodvec;
                rodvec[2] = rodvec_comp[2] * inv_norm_rodvec;
                rodvec[3] = norm_rodvec;
            }
        });
    }

    ///Returns a clone of the compact Rodrigues vector structure
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_comp_inplace(&self, rod_vec_comp: &mut RodVecComp){
        let mut ori = rod_vec_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);
    }

    ///Converts the compact Rodrigues vector representation over to a unit quaternion which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_quat_inplace(&self, quat: &mut Quat){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_quat_inplace(quat);
    }

    ///Converts the compact Rodrigues vector representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric_inplace(homochoric);
    }


}//End of impl Ori_Conv of Compact Rodrigues Vector 