# Mobility_RIS
These codes correspond to the following paper accepted in the ICC 2023 Conference. 

IEEE
* M. Haghshenas, P. Ramezani and E. Björnson, "Efficient LOS Channel Estimation for RIS-Aided Communications Under Non-Stationary Mobility," ICC 2023 - IEEE International Conference on Communications, Rome, Italy, 2023, pp. 2007-2012, doi: 10.1109/ICC45041.2023.10279469. keywords: {Training;Maximum likelihood estimation;Base stations;Tracking;Wireless networks;Refining;Channel estimation;Reconfigurable intelligent surface;parametric channel estimation;maximum likelihood estimator},


arXiv (Open Access)
* Haghshenas, Mehdi, Parisa Ramezani, and Emil Björnson. "Efficient LOS Channel Estimation for RIS-Aided Communications Under Non-Stationary Mobility." arXiv preprint arXiv:2303.16544 (2023).
# List of main files #
* **SmartInitRun**: It generates Fig. 2, where it loads specific data (fast13.mat), including the large-scale parameters of the channel.
* **RandomWalkRun**: It generates and extracts the large-scale parameters of the channel. 
* **PlotTrajectory**: It results in Fig. 1. It also has two configurations to demonstrate the user's random walk in the room.
* **MobilityMLTracking**: It generates Fig. 3. One can configure how often the channel should be re-estimated and RIS be reconfigured.


# List of secondary files:
* **OldPaper.m**: it includes the code for the following paper: 
  * Björnson, Emil, and Parisa Ramezani. "Maximum Likelihood Channel Estimation for RIS-Aided Communications With LOS Channels." arXiv preprint arXiv:2209.10867 (2022).
* **ML_ULA_direct_test.m**: To test ULA structured RIS + direct channel
* **search2D.m**: To test the search over azimuth and elevation angles
* **search2D_direct.m** To test the algorithm when RIS is equipped with UPA and the direct channel exists. 
