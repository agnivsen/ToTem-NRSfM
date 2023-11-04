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


**Pre-print:** [shorturl.at/zDUV0](https://encov.ip.uca.fr/publications/pubfiles/2023_Sengupta_etal_IJCV_topological.pdf) (also provided in the 'Docs' folder)


   ---

## Dependencies

The following are the dependencies for the code (included in the 'Data' folder):

- [CrustOpen](mathworks.com/matlabcentral/fileexchange/63731-surface-reconstruction-from-scattered-points-cloud-open-surfaces): based on [Bernardini et al, 1999]
- Generalized Procrustes Analysis (GPA): based on [[Bai et al., 2022]](https://encov.ip.uca.fr/publications/pubfiles/2021_Bai_etal_IJCV_defgpa.pdf)
- [Gloptipoly](https://homepages.laas.fr/henrion/software/gloptipoly3/) on [SeDuMi](https://github.com/sqlp/sedumi)
  
 ---

 ## Scripts

 The following scripts are included in this repository:
 
- _IsometricNRSfM.m_: replicates section 4 of our article
- _Parameterisation.m_: replicates section 6 of our article
- _SurfaceBasedNonRigidBundleAdjustment.m_: replicates section 7 of our article

The code provided in this repository has been tested on the following platforms:
- Matlab 2022a on Windows 11
- Matlab 2022b on Ubuntu 20.04 LTS (Focal Fossa)


## Citation

Our article has been recently accepted by the International Journal of Computer Vision (IJCV) and is still in press. If you find this code useful, you may cite the pre-print available at:

```
@article{sengupta2023totem,
  title={ToTem NRSfM: Object-wise Non-Rigid Structure-from-Motion with a Topological Template},
  author={Sengupta, Agniva and Bartoli, Adrien},
  journal={International Journal of Computer Vision},
  URL = {https://encov.ip.uca.fr/publications/pubfiles/2023_Sengupta_etal_IJCV_topological.pdf}
  year={2023}
}
```
 ---

 **References**

[Bernardini et al, 1999]: Bernardini, F., Mittleman, J., Rushmeier, H., Silva, C., & Taubin, G. (1999). The ball-pivoting algorithm for surface reconstruction. IEEE transactions on visualization and computer graphics, 5(4), 349-359.

[Bai et al., 2022]: Bai, Fang, and Adrien Bartoli. "Procrustes analysis with deformations: A closed-form solution by eigenvalue decomposition." International Journal of Computer Vision (2022): 1-27.
