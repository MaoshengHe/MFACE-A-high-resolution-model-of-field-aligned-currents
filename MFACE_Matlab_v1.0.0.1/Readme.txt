
------------------
Description
------
MFACE is a high-resolution model for auroral field-aligned currents derived from ten years of CHAMP data using a novel technique based on EOF decomposition. Predictor variables include MLat, MLT, DoY, IMF conditions, Vsw, and AE index. Compared to existing models, MFACE offers significantly better spatial resolution, accurately reproduces observed FAC thickness and intensity, improves the magnetic local time distribution, and provides the seasonal dependence of FAC latitudes and the NBZ current signature. Our EOF decomposition reveals that EOF1 primarily represents the Bz-controlled large-scale R1/R2 current pattern, while EOF2 represents the By-controlled cusp current signature.

------------------
v1.0.0.1 Update
------
The license has been changed from v1.0.0.0 to the Apache License, Version 2.0.
The Matlab code has been translated to Python 3 and released.

For details, refer to:
- https://doi.org/10.1029/2012GL053168
- https://doi.org/10.1002/2014JA019776

------------------
Examples of the output are available at
------
-FAC evolution at active levels of geomagnetic activity
https://www.youtube.com/watch?v=XV5xiRcsd4Q&ab_channel=M.He
-FAC evolution at moderate levels of geomagnetic activity
https://www.youtube.com/watch?v=VKMcKda0Khs&ab_channel=M.He
-FAC evolution at quiet levels of geomagnetic activity
https://www.youtube.com/watch?v=wxuOZqJwjMs&ab_channel=M.He
To request similar movies or figures for any specific periods, please contact Maosheng He.
------------------
Contents
------
- MFACE_v1.m: Code of MFACE v1.0
- MFACE_30min_Lag_N_v1.0: MFACE v1.0 coefficients for the Northern Hemisphere
- MFACE_30min_Lag_S_v1.0: MFACE v1.0 coefficients for the Southern Hemisphere
- AEmodel.mat: Model coefficients for AE index

--
- MFACE_Example.m: Example code
- exampleData.mat: Data for the example
- MFACE_Example.png: Output of the example

-----------------
Communication to
-------
Maosheng He  
Jacobs University Bremen,  
Campus Ring 1,  
28725 Bremen, Germany  
+49-421-200-3209  
m.he@jacobs-university.de / hmq512@gmail.com 
Sept. 20, 2013
------------------------
Copyright (c) 2012, Jacobs University Bremen
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

-----------------
Cite As 
-------
He, M., J. Vogt, H. LŸhr, E. Sorbalo, A. Blagau, G. Le, and G. Lu (2012), A high-resolution model of field-aligned currents through empirical orthogonal functions analysis (MFACE), Geophys. Res. Lett., 39, L18105, doi:10.1029/2012GL053168.

He, M., J. Vogt, H. LŸhr, and E. Sorbalo (2014), Local time resolved dynamics of field-aligned currents and their response to solar wind variability, J. Geophys. Res. Space Physics, 119, 5305Ð5315, doi:10.1002/2014JA019776.
