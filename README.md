# Mobility_RIS
These codes correspond to the following paper accepted in the ICC 2023 Conference. 
* Haghshenas, Mehdi, Parisa Ramezani, and Emil Björnson. "Efficient LOS Channel Estimation for RIS-Aided Communications Under Non-Stationary Mobility." arXiv preprint arXiv:2303.16544 (2023).
# List of main files #
* **SmartInitRun**: It generates Fig. 2 where it loads specific data including the large-scale parameters of the channel.
* **RandomWalkRun**: It generates and extracts the large-scale parameters of the channel. 
* **PlotTrajectory**: It results in Fig. 1. It also has two configurations to demonstrate the user's random walk in the room.
* **MobilityMLTracking**: It generates Fig. 3. One can configure how often the channel should be re-estimated and RIS be reconfigured.


# List of secondary files:
* **OldPaper**: it includes the code for the following paper: 
  * Björnson, Emil, and Parisa Ramezani. "Maximum Likelihood Channel Estimation for RIS-Aided Communications With LOS Channels." arXiv preprint arXiv:2209.10867 (2022).
