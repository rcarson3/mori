// This file is a part of the mori - Material Orientation Library in Rust
// Copyright 2018 Robert Carson

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
pub struct AngAxis{
    pub ori: Array2<f64>,
}

#[derive(Clone, Debug)]
pub struct AngAxisComp{
    pub ori: Array2<f64>
}

impl AngAxis{

    //Creates an array of zeros for the initial angle axis parameterization when data is not fed into it
    pub fn new(size: usize) -> AngAxis{
        assert!(size > 0, "Size inputted: {}, was not greater than 0", size);

        let mut ori = Array2::<f64>::zeros((4, size).f());

        ori.slice_mut(s![3, ..]).mapv_inplace(|_x| 1.0_f64);

        AngAxis{
            ori,
        }
    }//End of new

    //Creates an angle axis parameterization  type with the supplied data as long as the supplied data is in the following format
    //shape (3, nelems), memory order = fortran/column major.
    //If it doesn't fit those standards it will fail.
    pub fn new_init(ori: Array2<f64>) -> AngAxis{

        let nrow = ori.rows();

        assert!(nrow == 4, "Number of rows of array was: {}, which is not equal to 4", nrow);
        //We need to deal with a borrowing of ori here, so we need to have strides dropped at one point.
        {
            let strides = ori.strides();

            assert!(strides[0] == 1, "The memory stride is not column major (f order)");
        }

        AngAxis{
            ori,
        }
    }//End of new_init

    //Converts the angle axis representation over to a rotation matrix which has the following properties
    //shape (3, 3, nelems), memory order = fortran/column major.
    pub fn to_rmat(&self) -> RMat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array3::<f64>::zeros((3, 3, nelems).f());

        //The cosines and sines of the function
        let mut c:f64;
        let mut s:f64;

        for i in 0..nelems{

            c = f64::cos(self.ori[(3, i)]);
            s = f64::sin(self.ori[(3, i)]);

            ori[(0, 0, i)] = c + (1.0_f64 - c) * (self.ori[(0, i)] * self.ori[(0, i)]);
            ori[(1, 0, i)] = (1.0_f64 - c) * (self.ori[(0, i)] * self.ori[(1, i)]) + s * self.ori[(2, i)];
            ori[(2, 0, i)] = (1.0_f64 - c) * (self.ori[(0, i)] * self.ori[(2, i)]) - s * self.ori[(1, i)];

            ori[(0, 1, i)] = (1.0_f64 - c) * (self.ori[(0, i)] * self.ori[(1, i)]) - s * self.ori[(2, i)];
            ori[(1, 1, i)] = c + (1.0_f64 - c) * (self.ori[(1, i)] * self.ori[(1, i)]);
            ori[(2, 1, i)] = (1.0_f64 - c) * (self.ori[(1, i)] * self.ori[(2, i)]) + s * self.ori[(0, i)];

            ori[(0, 2, i)] = (1.0_f64 - c) * (self.ori[(0, i)] * self.ori[(2, i)]) + s * self.ori[(1, i)];
            ori[(1, 2, i)] = (1.0_f64 - c) * (self.ori[(1, i)] * self.ori[(2, i)]) - s * self.ori[(0, i)];
            ori[(2, 2, i)] = c + (1.0_f64 - c) * (self.ori[(2, i)] * self.ori[(2, i)]);

        }

        RMat{
            ori,
        }
    }//End of to_RMat

    //Converts the angle axis representation over to Bunge angles which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_bunge(&self) -> Bunge{
        let rmat = self.to_rmat();
        //When a pure conversion doesn't exist we just use the already existing ones in other orientation
        //representations 
        rmat.to_bunge()
    }

    //This should be fixed to return a clone of self
    pub fn to_ang_axis(&self) -> AngAxis{
        self.clone()
    }

    //Converts the angle axis representation over to a compact angle-axis representation which has the following properties
    //shape (3, nelems), memory order = fortran/column major.
    pub fn to_ang_axis_comp(&self) -> AngAxisComp{
        //We first convert to a angle axis representation. Then we scale our normal vector by our the rotation
        //angle which is the fourth component of our angle axis vector.
        // let ang_axis = self.to_ang_axis();

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        for i in 0 .. nelems{
            ori[(0, i)] = self.ori[(0, i)] * self.ori[(3, i)];
            ori[(1, i)] = self.ori[(1, i)] * self.ori[(3, i)];
            ori[(2, i)] = self.ori[(2, i)] * self.ori[(3, i)];
        }

        AngAxisComp{
            ori,
        }
    }//End of to_ang_axis_comp

    //Converts the angle axis representation over to a rodrigues vector representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_rod_vec(&self) -> RodVec{
        //We first convert to a angle axis representation. Then we just need to change the last component
        //of our angle axis representation to be tan(phi/2) instead of phi
        // let ang_axis = self.to_ang_axis();

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv2 = 1.0_f64/2.0_f64;


        for i in 0 .. nelems{
            ori[(0, i)] = self.ori[(0, i)];
            ori[(1, i)] = self.ori[(1, i)];
            ori[(2, i)] = self.ori[(2, i)];
            ori[(3, i)] = f64::tan(inv2 * self.ori[(3, i)]);
        }

        RodVec{
            ori,
        }
    }//End of to_rod_vec

    //Converts the angle axis representation over to a compact rodrigues vector representation which has the following properties
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

    //Converts the angle axis representation over to a unit quaternion representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_quat(&self) -> Quat{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((3, nelems).f());

        let inv2 = 1.0_f64 / 2.0_f64;
        let mut s:f64;

        for i in 0..nelems{

            s = f64::sin(inv2 * self.ori[(3, i)]); 

            ori[(0, i)] = f64::cos(inv2 * self.ori[(3, i)]);
            ori[(1, i)] = s * self.ori[(0, i)];
            ori[(2, i)] = s * self.ori[(1, i)];
            ori[(3, i)] = s * self.ori[(2, i)];

        }

        Quat{
            ori,
        }
    }//End of to_quat

    //Converts the angle axis representation over to a homochoric representation which has the following properties
    //shape (4, nelems), memory order = fortran/column major.
    pub fn to_homochoric(&self) ->Homochoric{

        let nelems = self.ori.len_of(Axis(1));

        let mut ori = Array2::<f64>::zeros((4, nelems).f());

        let inv3  = 1.0_f64 / 3.0_f64;
        let inv34 = 3.0_f64 / 4.0_f64;
        //Inner calculation used in the power term 
        let mut pow_term:f64;

        for i in 0 .. nelems{

            pow_term    = inv34 * (self.ori[(3, i)] - f64::sin(self.ori[(3, i)])); 

            ori[(0, i)] = self.ori[(0, i)];
            ori[(1, i)] = self.ori[(1, i)];
            ori[(2, i)] = self.ori[(2, i)];
            ori[(3, i)] = pow_term.powf(inv3);
        }

        Homochoric{
            ori,
        }
    }
}//End of Impl AngAxis

