This file contains 6 sequences of images and RF signals acquired on a phantom consisting of a background medium 
traversed by a cylindrical inclusion. Each sequence consists of a serie of US images and RF signals acquired with 
an indentation ranging from 0 to 7 mm in steps of 0.2 mm. The data (images and RF signals) are contained in the 6 
folders Acqui1 → Acqui6 corresponding to 6 acquired image sequences at different positions of the phantom. 
Each Acqui_i folder contains two subfolders: the RF folder containing RF signals in .mat file format (readable in Matlab), 
and the folder US_Image containing the images of the sequence in PNG format. Both images and RF signals are numbered, 
with each index corresponding to an indentation level and a force measured by the force sensor as outlined in the Excel file (Tabulation 1). 
Each RFi.mat file comprises 3152 rows representing the signal along the temporal axis, and 256 columns corresponding to the number of A lines 
in the image.

The Excel file has 5 tabs: The first tab contains, for each of the 6 acquired image sequences, 
the frame number in the sequence, the corresponding indentation of the probe (in mm), 
the recorded voltage value on the force sensor (V), and the corresponding calculated force value (N). 
Thus, each image in the 6 sequences is identified by a frame number, an indentation, and a force value. 
The second tab provide the acquisition parameters of the ultrasound images (frequency, depth, gain, etc.) 
performed using the SonixTablet ultrasound system from Ultrasonix (now bk medical). 
The third tab contains stress-strain curves, the mean, and the standard deviation of the Young's modulus 
for both the inclusion and the background of the phantom. 
The Young's modulus is obtained through compression tests conducted by an electromechanical testing machine (Bose, Electroforce 3200)
 on 10 small cylindrical samples taken from the background and the inclusion. The fourth tab contains the geometry and dimensions of the phantom.
The fifth tab contains the recipe used to make the gelatin phantom.


