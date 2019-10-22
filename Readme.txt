To run the code, run "ems_ekf_ins_gnss_main_TC.m"

The format of the INS-GNSS Tightly-Coupled dataset, "ems_data_TC_simulated.mat", is as follows:

time: time measurement syncing IMU (@100hz) and GNSS (@1Hz) data
	  the first IMU sample and the first GNSS sample are synced

error_free_accel: pseudo error-free accelerometer data
error_free_gyro: pseudo error-free gyroscope data
sim_accel: simulated noisy accelerometer data
sim_gyro: simulated noisy gyroscope data

roll: ground truth roll orientation data
pitch: ground pitch orientation data
heading ground truth heading orientation data

lat: ground truth latitude data
lon: ground truth longitude data
h: ground truth height data
lat_noisy: noisy latitude data
lon_noisy: noisy longitude data
h_noisy: noisy height data

north: ground truth north position data
east: ground truth east position data
north_noisy: noisy north position data
east_noisy: noisy east position data

vel: ground truth velocity data in NED frame
vel_N: ground truth velocity data in ENU frame
vel_noisy: noisy velocity data in NED frame
vel_N_noisy: noisy velocity data in ENU frame

no_of_sat: number of satellites observed at each measurement
SVN: corresponding IDs of observed satellites
sat_Px: satellite position data in ECEF frame along x axis in the same order as SVN
sat_Py: satellite position data in ECEF frame along y axis in the same order as SVN
sat_Pz: satellite position data in ECEF frame along z axis in the same order as SVN
sat_Vx: satellite velocity data in ECEF frame along x axis in the same order as SVN
sat_Vy: satellite velocity data in ECEF frame along y axis in the same order as SVN
sat_Vz: satellite velocity data in ECEF frame along z axis in the same order as SVN
sim_sat_range: simulated satellite range measurement data in the same order as SVN

noiseInfo: noise parameters including: 
		   Accelerometer (XYZ) bias 
		   Accelerometer (XYZ) bias random noise std
           Gyrsocope (XYZ) bias
		   Gyrsocope (XYZ) random noise std
		   Position (NED) random noise std
		   Velocity (NED) random noise std
		   Receiver clock bias
		   Range measurement random noise std
		   