# VUMAT-Implementation-in-ABAQUS-Explicit
A finite deformation fortran code is shared in the repository along with files associated with various test cases. Zip files are shared due to some uploading issues to ascertain that all files are available to the user. User requires a Intel fortran compiler & VS code linked with ABAQUS. User can supply the shared fortran code while doing simulations. All the required material properties are defined in the Fortran code.
Various Standard test cases were considered for validating the code with ABAQUS inbuilt models.
# Uniaxial Tension
A single element is subjected to uni-axial load resulting in deformation as shown in the figure below.

![Uniaxial_tension](https://user-images.githubusercontent.com/88960574/129486840-44a366f3-1f65-43e0-860e-1e6fb6bc8087.png)

**Comparison between Inbuilt model with VUMAT result**

![stress_vs_strain](https://user-images.githubusercontent.com/88960574/129487649-5e90d17c-53c8-451c-a171-965d56b5f329.png)

# Simple Shear
A single element is subjected to shear load resulting in deformation as shown in the figure below.

![Simple_shear](https://user-images.githubusercontent.com/88960574/129486844-4525c41c-e9cd-4248-a177-83565bcea323.png)

**Comparison between Inbuilt model with VUMAT result**

![stress_vs_strain](https://user-images.githubusercontent.com/88960574/129487708-187a7b20-ead8-40b4-9f0f-192d7d686e95.png)

# Quarter plate with Hole
A quarter-plate is subjected to load in x-direction and fixed along y-axis direction. The resulting deformation can be seen in the image below. Sufficient numbers of elements are considered to qualify mesh convergence study.

![Quater_plate with hole](https://user-images.githubusercontent.com/88960574/129486847-69213198-1433-4ce3-8629-af3d1b4d9648.png)

**Comparison between Inbuilt model with VUMAT result**

![stress_vs_strain](https://user-images.githubusercontent.com/88960574/129487678-daddafb5-3a0d-4fee-a297-44de7c157a2a.png)

