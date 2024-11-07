%%%% Enhancing Photoacoustic Imaging Quality Using Virtual Sensor Points

% Overview

Photoacoustic imaging (PAI) is an advanced technique that combines the high contrast of optical imaging 
with the spatial resolution of ultrasound imaging. Despite its potential, 
challenges such as limited sensor points and measurement noise can adversely affect image quality. 
This project introduces a novel algorithm that generates data 
for virtual sensors positioned between real physical ones, 
utilizing the Kalman filter. 
By incorporating these virtual sensors,
we aim to reduce costs while significantly enhancing image quality in PAI, 
particularly in scenarios with a limited number of physical sensors.

%% Requirements

To run this code, you will need:
- MATLAB: Ensure you have MATLAB installed on your machine.
- k-Wave Toolbox: This toolbox is essential for simulating photoacoustic imaging. 
You can download it from [k-Wave Toolbox](http://www.k-wave.org/).

%%% Code Overview

The provided MATLAB code allows you to simulate photoacoustic imaging 
with different sensor configurations. You can change the type of sensor (linear or circular) 
and adjust the number of sensors based on your requirements.
