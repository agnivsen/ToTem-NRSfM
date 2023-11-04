# Deformable Procrustes Analysis

The Matlab package for the generalized Procrustes analysis (GPA) problem with deformation models:

- [linear basis warps (LBWs) in IJCV 2022](https://link.springer.com/article/10.1007/s11263-021-01571-8), e.g., the affine transformation and the thin-plate spline (TPS) warp

- [kernel based transformation models in RSS 2022](http://www.roboticsproceedings.org/rss18/p002.html)

Maintainer: Fang Bai (fang.bai@yahoo.com)

---
## Content

[**DefGPA**](./DefGPA/). This folder contains the implementation of Generalzied Procrustes Analysis (GPA) with Linear Basis Warps (LBWs). The implementation specializes two commonly used LBWs: 1) the affine transformation, and 2) the thin-plate spline (TPS) warp.

- [DefGPA.m](./DefGPA/DefGPA.m) The class implements the ***DefGPA*** method. 

- [ThinPlateSpline.m](./DefGPA/ThinPlateSpline.m) The implementation of static class ***ThinPlateSpline*** which includes all supporting functions for the TPS warp

- [ImageWarp.m](./DefGPA/ImageWarp.m) A function *ImageWarp* used in experimental visualization. It warps an image using a given transformation function *warpFunc*.

[**@KernelGPA**](@KernelGPA/KernelGPA.m) This class implements the ***KernelGPA*** method.

[**@RigidGPA**](@RigidGPA/RigidGPA.m) This class implements the ***RigidGPA*** method.


---
## How to use the code

An example on how to use the code is provided in [run_examples.m](./run_examples.m)

**STEP1. Create the problem handle**

We provide four GPA methods.

>- GPA = DefGPA('AFFINE')

>- GPA = DefGPA('TPS')

>- GPA = KernelGPA

>- GPA = RigidGPA

**STEP2. Add data**

An example of the desired data format is given in [readDataset.m](./readDataset.m).

>- GPA.addDataByArray(DataCell, VisVecCell)

The i-th cell contains the information for the i-th point-cloud.

**i_th_point_cloud = DataCell{i}** is a 3byN (or 2byN) matrix, with its columns arranged according to correspondences.

**i_th_vis = VisVecCell{i}** is a 1byN vector, indicating the visibility of each points. In particular, if **i_th_vis(k) == 0**, then the k-th point does not occur in the **i_th_point-cloud**, which means the values in the k-th column of **i_th_point_cloud** are meanless.


**STEP3. Set tuning parameters**

- **The tuning parameters of the TPS warp**. For TPS warp, we need to assign the *number of control points*, and the *smooth paramter* indicating the regularization strength by the bending energy. Typically, we suggest use 5^d (d=2 or 3) control points.

>- GPA.dimFeatureSpace = 5^3

>- GPA.smoothParam = 0.01; % for 3D data

>- GPA.smoothParam = 10; % for 2D data

- **The tuning paramters of the kernel model**. We implement the Gaussian kernel. So we need to set the kernel bandwidth decided by the *p-quantile of the pairwise Euclidean distances*, and the regularization strength (i.e., the *smooth paramter*).

>- GPA.smoothParam = 0.05;

>- GPA.quantilePercentParam = 0.2;

**STEP4. Run the data **

>- GPA.run()
>- GPA.run_N_Fold_CrossValidation(N)

**GPA.run()** solves the data once, using the specified transformation model, i.e., rigid, affine, TPS, or Kernel.

**GPA.run_N_Fold_CrossValidation(N)** gives the result of N-fold cross-validation. We suggest using N = 10 or N = 20 to search for the best tuning parameters.

*Note*: If N = inf, then GPA.run_N_Fold_CrossValidation(inf) gives the result of leave-one cross-validation.

**STEP5. The output and results **

>- GPA.mShape
>- transformed_point_cloud = GPA.transformPoints(input_point_cloud, point_cloud_index)
>- GPA.rsdError()

**GPA.mShape**. A 3byN (or 2byN) matrix with each column giving the coordinates of the estimated reference point-cloud.

**transformed_point_cloud = GPA.transformPoints(input_point_cloud, point_cloud_index)**. The function **transformPoints()** is the estimated transformation model. Its first argument is the input point-cloud, i.e., the point-cloud you want to transform. Its second argument expects a scalar indicating which transformation you want to use, e.g., 2 means you want to use the transformation of the second point-cloud.

**GPA.rsdError()**. The optimal cost of the problem. For DefGPA, we have two costs, **reference space cost** and **datum space cost**. For KernelGPA, we only have **reference space cost**.

*Note*. The above three are the most important outputs for GPA. However, it might be good to just type GPA **without ;** in the command window to see all the computed properties of your GPA call.


---
## Dataset

[readDataset.m](./readDataset.m) is the interface file to read datasets. The original datasets are provided at the [EnCoV website](http://igt.ip.uca.fr/~ab/code_and_datasets/index.php):

Below is the list of datasets used for evaluation.

- **HandBag(2D)**. Deformable Surfaces and Objects. http://igt.ip.uca.fr/~ab/code_and_datasets/datasets/handbag_v1p0.zip

- **PillowCover(2D)**. Deformable Surfaces and Objects. http://igt.ip.uca.fr/~ab/code_and_datasets/datasets/pillow_cover_v1p0.zip

- **StereoToyRug(3D)**. Deformable Surfaces and Objects. http://igt.ip.uca.fr/~ab/code_and_datasets/datasets/stereo_toy_rug_v1p0.rar

- **Semisynthetic Liver (3D)**. A simulated 3D Liver. Laparoscopy. http://igt.ip.uca.fr/~ab/code_and_datasets/datasets/semisynthetic_liver_v1p0.zip

- **TOPACS (3D)**. Real CT data with around 1320 landmarks. Initial correspondences are found by matching feature descriptors and then refined by an ICP algorithm (two shape registeration). The global correspondences across multiple shapes are found by a graph matching algorithm, and the ambiguous ones are removed based on distances. http://igt.ip.uca.fr/~ab/code_and_datasets/datasets/topacs6_v1p0.zip

- **LiTS (3D)**. Real CT data with 54 fiducial landmarks. http://igt.ip.uca.fr/~ab/code_and_datasets/datasets/lits_landmarks_v1p0.zip

The support for downloadable links is not stable. ***Please copy and paste the link in a new brower tab, then it should start to download.***



---
## Experiment script

[**IJCV_experiment_scripts**](https://bitbucket.org/clermontferrand/deformableprocrustes/src/master/IJCV_experiment_scripts/)

This contains the experimental scripts to reproduce of the results of the IJCV paper:

[*Procrustes Analysis with Deformations: A Closed-Form Solution by Eigenvalue Decomposition*](https://link.springer.com/article/10.1007/s11263-021-01571-8)

Run [run.m](./IJCV_experiment_scripts/run.m) to check all the results.

*Note:* The leave-one cross-validation process is very slow, in particular for the synthetic Liver dataset. Please comment this if you're not interested in this stattistics.

[**RSS_experiment_scripts**](https://bitbucket.org/clermontferrand/deformableprocrustes/src/master/RSS_experiment_scripts/)

This contains the experimental scripts to reproduce of the results of the RSS paper:

[*KernelGPA: A Deformable SLAM Back-end*](http://www.roboticsproceedings.org/rss18/p002.html)

Run [run.m](./RSS_experiment_scripts/run.m) to check all the results.

*Note:* As the leave-one cross-validation is very slow, we suggest to use N-fold (N=20) cross-validation to find the best tuning parameter, and use leave-one cross-validation to visualize the fitness only.



---
## Reference
**DefGPA**

```
@article{bai2022defgpa,
	title={{Procrustes Analysis with Deformations: A Closed-Form Solution by Eigenvalue Decomposition}},
	author={Bai, Fang and Bartoli, Adrien},
	journal={International Journal of Computer Vision},
	year={2022},
	volume={130},
	number={2},
	pages={567-593},
	doi={10.1007/s11263-021-01571-8}
}
```

| [Springer Official](https://link.springer.com/article/10.1007/s11263-021-01571-8) | [EnCov (one column pdf)](http://igt.ip.uca.fr/encov/publications/pubfiles/2021_Bai_etal_IJCV_defgpa.pdf) | [arXiv (double column pdf)](https://arxiv.org/pdf/2206.14528.pdf) |

**KernelGPA**
```
@INPROCEEDINGS{Bai-RSS-22, 
    AUTHOR    = {Fang Bai AND Adrien Bartoli}, 
    TITLE     = {{KernelGPA: A Deformable SLAM Back-end}}, 
    BOOKTITLE = {Proceedings of Robotics: Science and Systems}, 
    YEAR      = {2022}, 
    ADDRESS   = {New York City, NY, USA}, 
    MONTH     = {June}, 
    DOI       = {10.15607/RSS.2022.XVIII.002} 
}
```

| [RSS Official (open access pdf)](http://www.roboticsproceedings.org/rss18/p002.pdf) |

---
## License

The code is relased under the [Apache License](./LICENSE).
