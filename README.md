# ToTem NRSfM: Object-wise Non-Rigid Structure-from-Motion with a Topological Template

We present a Non-Rigid Structure-from-Motion (NRS*f*M) method to reconstruct an object whose topology is known.
   We represent the topology by a 3D shape that weakly resembles the object, which we call a **To**pological **Tem**plate (ToTem).

   ![teaser](https://github.com/agnivsen/ToTem-NRSfM/assets/5153445/dadb67c5-d1e6-4eec-a830-a6e87ef3599d)

   The ToTem has two main differences with the template used in Shape-from-Template (S*f*T).
   First, the shape in the ToTem is not necessarily feasible for the object, whereas it must be in the S*f*T's template.
   Second, the ToTem only models shape, excluding the classical texture map representing color in the S*f*T's template. 
   These two differences greatly alleviate the practical difficulty of constructing a template.
   However, they make the reconstruction problem challenging, as they preclude the use of strong deformation constraints between the template shape and the reconstruction and the possibility of directly establishing correspondences between the template and the images. 
   
   
   Our method uses an isometric deformation prior and proceeds in four steps:
   1. It reconstructs point clouds from the images
   2. It aligns the ToTem to the point clouds
   3. It creates a coherent surface parameterization
   4. It performs a global refinement, posed as a Non-Rigid Bundle Adjustment (NRBA)
     
   We show experimentally that our method outperforms the existing methods for its isolated steps and NRS*f*M methods overall, in terms of 3D accuracy, ability to reconstruct the object's visible surface, and ability to approximate the object's invisible surface. 




   ---

   
