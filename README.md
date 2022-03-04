## Shoulder check

This repository contains all the matlab code used for data analysis and visualisation of the research article: *Turning the head while biking makes older people lose cycling direction and balance*. This article is currently under peer review and is already availlable as a pre-rint (https://doi.org/10.1101/2022.03.01.481993 ).



## Data analysis

1) Clone or download this repository 

2) Download the IMU sensor data from (). and unpack the zip file in a subfolder Data of this repository. *(Note: if you save the data at a different folder you will have to adjust the path to the datafiles in the first lines of the two matlab scripts)*.

3) Run the script *GetRotationAxisSteer_Subjects.m* to compute the axis of rotation between the frame and the handlebars in each subject for each type of bike. This axis of rotation is saved in a .mat file in the datafolder.

4) Run the script *GetDataMatrix_ShoulderCheck.m* to create a large matrix with all the outcomes and information on subject number, type of bike and so on. This datamatrix will be saved in a subfolder called *Outcomes*.

5) Run the scripts in the subfolder *PlotFiguresPaper* to do the statistics an visualise the results.










