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

///A structure that holds an array of compact axis-angle representation of a rotation
#[derive(Clone, Debug)]
pub struct AngAxisComp{
    ori: Array2<f64>
}

impl AngAxisComp{

    ///Creates an array of zeros for the initial compact axis-angle parameterization when data is not fed into it
    pub fn new(size: usize) -> AngAxisComp{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((3, size).f());

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))) in {angaxis[2] = 1.0_f64});

        AngAxisComp{
            ori,
        }
    }//End of new

    ///Creates a compact axis-angle parameterization  type with the supplied data as long as the supplied data is in the following format
    ///shape (3, nelems), memory order = fortran/column major.
    ///If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> AngAxisComp{

        let nrow = ori.rows();

        assert!(nrow == 3, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        AngAxisComp{
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
}//End of AngAxisComp impl

///The orientation conversions of a series of compact axis-angle representation to a number of varying different orientation
///representations commonly used in material orientation processing. 
impl OriConv for AngAxisComp{
    ///Converts the compact axis-angle representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    fn to_rmat(&self) -> RMat{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let mut angaxis = Array1::<f64>::zeros((4).f());
            
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                angaxis[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
                angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
                angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
                angaxis[3] = norm_angaxis;
            }


            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);
        });

        RMat::new_init(ori) 
    }//End of to_rmat

    ///Converts the compact axis-angle representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_bunge(&self) -> Bunge{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut bunge (ori.axis_iter_mut(Axis(2))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let mut angaxis = Array1::<f64>::zeros((4).f());
            
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                angaxis[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
                angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
                angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
                angaxis[3] = norm_angaxis;
            }


            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            let mut rmat = Array2::<f64>::zeros((3, 3).f());

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);

            if f64::abs(rmat[[2, 2]]) > (1.0_f64 - tol){
                bunge[0] = f64::atan2(rmat[[0, 1]], rmat[[0, 0]]);
                bunge[1] = std::f64::consts::FRAC_PI_2 * (1.0_f64 - rmat[[2, 2]]);
                bunge[2] = 0.0_f64;
            }else{
                let eta  = 1.0_f64 / f64::sqrt(1.0_f64 - rmat[[2, 2]] * rmat[[2, 2]]);
                bunge[0] = f64::atan2(rmat[[2, 0]] * eta, -rmat[[2, 1]] * eta);
                bunge[1] = rmat[[2, 2]].acos();
                bunge[2] = f64::atan2(rmat[[0, 2]] * eta, rmat[[1, 2]] * eta);
            }
        });

        Bunge::new_init(ori)
    }//End of to_bunge

    ///Converts the compact axis-angle representation over to an angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_ang_axis(&self) -> AngAxis{
        //We first convert to a axis-angle representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our axis-angle vector.
        // let ang_axis = self.to_ang_axis();

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                angaxis[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
                angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
                angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
                angaxis[3] = norm_angaxis;
            }
        });

        AngAxis::new_init(ori)
    }

    ///Returns a clone of the compact axis-angle structure
    fn to_ang_axis_comp(&self) -> AngAxisComp{
        self.clone()
    }//End of to_ang_axis_comp

    ///Converts the compact axis-angle representation over to a Rodrigues vector representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_rod_vec(&self) -> RodVec{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let tol = std::f64::EPSILON;

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                rod_vec[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                let tan2 = f64::tan(inv2 * norm_angaxis);

                rod_vec[0] = angaxis_comp[0] * inv_norm_angaxis;
                rod_vec[1] = angaxis_comp[1] * inv_norm_angaxis;
                rod_vec[2] = angaxis_comp[2] * inv_norm_angaxis;
                rod_vec[3] = tan2;
            }
        });

        RodVec::new_init(ori)
    }//End of to_rod_vec

    ///Converts the compact axis-angle representation over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    fn to_rod_vec_comp(&self) -> RodVecComp{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;
        let tol = std::f64::EPSILON;

        azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });

            if norm_angaxis.abs() > tol{ 

                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                let tan2 = f64::tan(inv2 * norm_angaxis);

                rod_vec[0] = angaxis_comp[0] * inv_norm_angaxis * tan2;
                rod_vec[1] = angaxis_comp[1] * inv_norm_angaxis * tan2;
                rod_vec[2] = angaxis_comp[2] * inv_norm_angaxis * tan2;
            }
        });

        RodVecComp::new_init(ori)
    }//End of to_rod_vec_comp

    ///Converts the compact axis-angle representation over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_quat(&self) -> Quat{
        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64 / 2.0_f64;
        let tol = std::f64::EPSILON;

        azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });

            let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

            let s = f64::sin(inv2 * norm_angaxis); 

            if norm_angaxis.abs() > tol{
                quat[0] = f64::cos(inv2 * norm_angaxis);
                quat[1] = s * angaxis_comp[0] * inv_norm_angaxis;
                quat[2] = s * angaxis_comp[1] * inv_norm_angaxis;
                quat[3] = s * angaxis_comp[2] * inv_norm_angaxis;
            }else{
                quat[0] = 1.0_f64;
            }
        });

        Quat::new_init(ori)
    }//End of to_quat

    ///Converts the compact axis-angle representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    fn to_homochoric(&self) ->Homochoric{
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric()
    }//End of to_homochoric

    ///Converts the compact axis-angle representation over to Bunge angles which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_bunge_inplace(&self, bunge: &mut Bunge){
        let mut ori = bunge.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut bunge (ori.axis_iter_mut(Axis(2))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let mut angaxis = Array1::<f64>::zeros((4).f());
            
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                angaxis[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
                angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
                angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
                angaxis[3] = norm_angaxis;
            }


            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            let mut rmat = Array2::<f64>::zeros((3, 3).f());

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);

            if f64::abs(rmat[[2, 2]]) > (1.0_f64 - tol){
                bunge[0] = f64::atan2(rmat[[0, 1]], rmat[[0, 0]]);
                bunge[1] = std::f64::consts::FRAC_PI_2 * (1.0_f64 - rmat[[2, 2]]);
                bunge[2] = 0.0_f64;
            }else{
                let eta  = 1.0_f64 / f64::sqrt(1.0_f64 - rmat[[2, 2]] * rmat[[2, 2]]);
                bunge[0] = f64::atan2(rmat[[2, 0]] * eta, -rmat[[2, 1]] * eta);
                bunge[1] = rmat[[2, 2]].acos();
                bunge[2] = f64::atan2(rmat[[0, 2]] * eta, rmat[[1, 2]] * eta);
            }
        });
    }

    ///Converts the compact axis-angle representation over to a rotation matrix which has the following properties
    ///shape (3, 3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rmat_inplace(&self, rmat: &mut RMat){
        let mut ori = rmat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(2));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut rmat (ori.axis_iter_mut(Axis(2))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let mut angaxis = Array1::<f64>::zeros((4).f());
            
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                angaxis[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
                angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
                angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
                angaxis[3] = norm_angaxis;
            }


            let c = angaxis[3].cos();
            let s = angaxis[3].sin();

            rmat[[0, 0]] = c + (1.0_f64 - c) * (angaxis[0] * angaxis[0]);
            rmat[[1, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) + s * angaxis[2];
            rmat[[2, 0]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) - s * angaxis[1];

            rmat[[0, 1]] = (1.0_f64 - c) * (angaxis[0] * angaxis[1]) - s * angaxis[2];
            rmat[[1, 1]] = c + (1.0_f64 - c) * (angaxis[1] * angaxis[1]);
            rmat[[2, 1]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) + s * angaxis[0];

            rmat[[0, 2]] = (1.0_f64 - c) * (angaxis[0] * angaxis[2]) + s * angaxis[1];
            rmat[[1, 2]] = (1.0_f64 - c) * (angaxis[1] * angaxis[2]) - s * angaxis[0];
            rmat[[2, 2]] = c + (1.0_f64 - c) * (angaxis[2] * angaxis[2]);
        }); 
    }

    ///Converts the compact axis-angle representation over to an angle-axis representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major. 
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_inplace(&self, ang_axis: &mut AngAxis){
        let mut ori = ang_axis.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let tol = std::f64::EPSILON;

        azip!(mut angaxis (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                angaxis[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                angaxis[0] = angaxis_comp[0] * inv_norm_angaxis;
                angaxis[1] = angaxis_comp[1] * inv_norm_angaxis;
                angaxis[2] = angaxis_comp[2] * inv_norm_angaxis;
                angaxis[3] = norm_angaxis;
            }
        });
    }

    ///Returns a clone of the compact axis-angle structure
    ///This operation is done inplace and does not create a new structure
    fn to_ang_axis_comp_inplace(&self, ang_axis_comp: &mut AngAxisComp){
        let mut ori = ang_axis_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        ori.assign(&self.ori);

    }

    ///Converts the compact axis-angle representation over to a Rodrigues vector representation which has the following properties
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

        let inv2 = 1.0_f64/2.0_f64;

        azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });
            //If we follow the same convention that we use with quaternions for cases with no rotation
            //then we set it equal to the following vector with the no rotation ([0, 0, 1], 0)
            if norm_angaxis.abs() < tol{
                rod_vec[2] = 1.0_f64; 
            }else{
                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                let tan2 = f64::tan(inv2 * norm_angaxis);

                rod_vec[0] = angaxis_comp[0] * inv_norm_angaxis;
                rod_vec[1] = angaxis_comp[1] * inv_norm_angaxis;
                rod_vec[2] = angaxis_comp[2] * inv_norm_angaxis;
                rod_vec[3] = tan2;
            }
        });
    }

    ///Converts the compact axis-angle representation over to a compact Rodrigues vector representation which has the following properties
    ///shape (3, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_rod_vec_comp_inplace(&self, rod_vec_comp: &mut RodVecComp){
        let mut ori = rod_vec_comp.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64/2.0_f64;
        let tol = std::f64::EPSILON;

        azip!(mut rod_vec (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });

            if norm_angaxis.abs() > tol{ 

                let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

                let tan2 = f64::tan(inv2 * norm_angaxis);

                rod_vec[0] = angaxis_comp[0] * inv_norm_angaxis * tan2;
                rod_vec[1] = angaxis_comp[1] * inv_norm_angaxis * tan2;
                rod_vec[2] = angaxis_comp[2] * inv_norm_angaxis * tan2;
            }
        });
    }
    
    ///Converts the compact axis-angle representation over to a unit quaternion representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_quat_inplace(&self, quat: &mut Quat){
        let mut ori = quat.ori_view_mut();

        let new_nelem = ori.len_of(Axis(1));
        let nelem = self.ori.len_of(Axis(1));

        assert!(new_nelem == nelem, 
        "The number of elements in the original ori field do no match up with the new field.
        The old field had {} elements, and the new field has {} elements",
        nelem, new_nelem);

        let inv2 = 1.0_f64 / 2.0_f64;
        let tol = std::f64::EPSILON;

        azip!(mut quat (ori.axis_iter_mut(Axis(1))), ref angaxis_comp (self.ori.axis_iter(Axis(1))) in {
            let norm_angaxis = f64::sqrt({
                angaxis_comp[0] * angaxis_comp[0] 
                + angaxis_comp[1] * angaxis_comp[1] 
                + angaxis_comp[2] * angaxis_comp[2]
                });

            let inv_norm_angaxis = 1.0_f64 / norm_angaxis;

            let s = f64::sin(inv2 * norm_angaxis); 

            quat[0] = f64::cos(inv2 * norm_angaxis);

            if norm_angaxis.abs() > tol{
                quat[1] = s * angaxis_comp[0] * inv_norm_angaxis;
                quat[2] = s * angaxis_comp[1] * inv_norm_angaxis;
                quat[3] = s * angaxis_comp[2] * inv_norm_angaxis;
            }
        });
    }

    ///Converts the compact axis-angle representation over to a homochoric representation which has the following properties
    ///shape (4, nelems), memory order = fortran/column major.
    ///This operation is done inplace and does not create a new structure
    fn to_homochoric_inplace(&self, homochoric: &mut Homochoric){
        let ang_axis = self.to_ang_axis();
        ang_axis.to_homochoric_inplace(homochoric);
    }

}//End of Impl OriConv for AngAxisComp