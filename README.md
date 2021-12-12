# QRS-detection-project

## _Quadratic filter algorithm for ECG signal processing and detecting the signal’s key features, in this case the QRS points._

The QRS complex representing ventricular depolarization often, but not always, consists of three separate waves: the Q-wave, the R-wave, and the S-wave. The most important step for detecting the QRS complex, the
segmentation of the heartbeat and the interval between two beats is precisely the QRS detection within the ECG (electrocardiogram) signal. The problems that complicate the mentioned detection are different types of noise and different types of abnormal morphologies. 

The algorithm that improves the ratio of QRS signals and noise, which is described here, is based on a square filter, developed by Pornchai Phukpattaranont. Improvements in the QRS signal-to-noise ratio allow us to use a single fixed threshold without any additional post-processing techniques in the beat detection step. This algorithm is implemented within MATLAB (Matrix laboratory). A set of data available on the PhysioBank ATM page, within the MIT-BIH Arrhythmia Database (mitdb), was used to check the operation of the algorithm itself

Tools:
- MATLAB

## _Result_

![alt text](https://github.com/smuminovic/QRS-detection-project/blob/main/result.png?raw=true)

