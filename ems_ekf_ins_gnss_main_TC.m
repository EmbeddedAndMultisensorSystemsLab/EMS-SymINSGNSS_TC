function ems_ekf_ins_gnss_main_TC()
%   main module
%
%   This is a basic INS/GNSS Tighly-Coupled system using EKF. This is only 
%   a DEMO with basic updates (range measurements) are applied. More advanced
%   updates such as range rate measurements, nonholonomic constraints, 
%   zero-speed, and adaptive EKF are NOT implemented in this DEMO. The purpose
%   of the DEMO is to demonstrate the basics of EKF in a basic INS/GNSS 
%   Tightly-Coupled fusion scenario. 
% 
%   For commercial use or embdded C/C++ versions, please contact mohamed.atia@carleton.ca.
% 
%   Copyright (C) 2019, Mohamed Atia, all rights reserved.
%   The software is given under GNU Lesser General Public License
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this program. If not, see
%   <http://www.gnu.org/licenses/>.


set(0,'DefaultFigureWindowStyle','docked');
opengl('save', 'software');
clear;
close all;
fclose all;
clc;


%--> Load EKF parameters
load('EKF_parameters_simulated');

%--> Call symbolic mechanization equation definition script
ems_symbolic_engine_TC;

%--> Enable/Disable animation of results
enable_animation = 1;

%--> Converting symbolic functions into real-valued functions
C_LB_from_Euler_fun         = matlabFunction(C_LB_from_Euler);
w_L_IL_fun                  = matlabFunction(w_L_IL);
a_dot_fun                   = matlabFunction(a_dot);
b_dot_fun                   = matlabFunction(b_dot);
c_dot_fun                   = matlabFunction(c_dot);
d_dot_fun                   = matlabFunction(d_dot);
C_LB_from_quat_fun          = matlabFunction(C_LB_from_quat);
C_EN_fun                    = matlabFunction(C_EN);
C_EL_fun                    = matlabFunction(C_EL);
w_N_EN_fun                  = matlabFunction(w_N_EN);
w_N_IE_fun                  = matlabFunction(w_N_IE);
V_N_dot_fun                 = matlabFunction(V_N_dot);
pos_quat_dot_fun_EN         = matlabFunction(pos_quat_dot_EN);
pos_quat_dot_fun_EL         = matlabFunction(pos_quat_dot_EL);
C_EN_from_pos_quat_fun      = matlabFunction(C_EN_from_pos_quat);
C_EL_from_pos_quat_fun      = matlabFunction(C_EL_from_pos_quat);
pos_quat_from_lat_lon_fun   = matlabFunction(pos_quat_from_lat_lon);
w_L_IE_fun                  = matlabFunction(w_L_IE);
gravity_fun                 = matlabFunction(g);
F_fun                       = matlabFunction(F);
G_fun                       = matlabFunction(G);
H_row_range_fun             = matlabFunction(H_row_range);

%--> Earth ellipsoid shape parameters
earth_a           = 6378137;
earth_f           = 1/298.257223563;
earth_b           = earth_a*(1-earth_f);
earth_e2          = 1-(earth_b^2)/(earth_a^2);
we_value          = 2*pi/(24*60*60);

%--> Load the dataset
load('ems_data_TC_simulated');

%--> Plot IMU and reference navigation data
figure;
plot(ems_data.time,ems_data.sim_gyro(:,1)*R2D,'r', 'LineWidth',1); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,ems_data.sim_gyro(:,2)*R2D,'g', 'LineWidth',1);
plot(ems_data.time,ems_data.sim_gyro(:,3)*R2D,'b', 'LineWidth',1); legend({'gyro x','gyro y','gyro z'}, 'FontSize', 12); title('Simulated IMU - Gyroscope'); xlabel('Time (s)');ylabel('Rotation Rate (deg/s)');

figure;
plot(ems_data.time,ems_data.sim_accel(:,1),'r', 'LineWidth',1); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,ems_data.sim_accel(:,2),'g', 'LineWidth',1);
plot(ems_data.time,ems_data.sim_accel(:,3),'b', 'LineWidth',1); legend({'accel x','accel y','accel z'}, 'FontSize', 12); title('Simulated IMU - Accelerometer'); xlabel('Time (s)');ylabel('Acceleration (m/s^2)');

figure;
plot(ems_data.time,ems_data.roll*R2D,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,ems_data.pitch*R2D,'g', 'LineWidth',2);
plot(ems_data.time,ems_data.heading*R2D,'b', 'LineWidth',2); legend({'Roll','Pitch','Heading'}, 'FontSize', 12);  title('Reference 3D Orientation'); xlabel('Time (s)');ylabel('Orientation (deg)');

figure;
plot(ems_data.time,ems_data.vel(:,1),'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,ems_data.vel(:,2),'g', 'LineWidth',2);
plot(ems_data.time,ems_data.vel(:,3),'b', 'LineWidth',2); legend({'V_n','V_e','V_d'}, 'FontSize', 12);  title('Reference 3D Velocity'); xlabel('Time (s)');ylabel('Velocity (m/s)');

figure;
plot(ems_data.time,ems_data.north,'r', 'LineWidth',2);hold on; grid on;set(gca,'FontSize',12);
plot(ems_data.time,ems_data.east,'g', 'LineWidth',2); 
plot(ems_data.time,-ems_data.h,'b', 'LineWidth',2); legend({'North','East','Down'}, 'FontSize', 12);  title('Reference 3D Position'); xlabel('Time (s)'); ylabel('Position (m)');

figure;
plot(ems_data.east, ems_data.north, 'LineWidth',2); title('Reference 2D Position'); xlabel('East (m)'); ylabel('North (m)'); grid on; set(gca,'FontSize',12);

%--> Initialization
num_of_samples = length(ems_data.time);
sampling_period_sec = mean(diff(ems_data.time));
sampling_freq = round(1/sampling_period_sec);

%--> Initialize vectors sizes
time_vector_sec                         = zeros(num_of_samples,1);
accel_x                               = zeros(num_of_samples,1);
accel_y                               = zeros(num_of_samples,1);
accel_z                               = zeros(num_of_samples,1);
gyro_x                              = zeros(num_of_samples,1);
gyro_y                              = zeros(num_of_samples,1);
gyro_z                              = zeros(num_of_samples,1);
a_value                                 = zeros(num_of_samples,1);
b_value                                 = zeros(num_of_samples,1);
c_value                                 = zeros(num_of_samples,1);
d_value                                 = zeros(num_of_samples,1);
Euler_pitch_value                       = zeros(num_of_samples,1);
Euler_roll_value                        = zeros(num_of_samples,1);
Euler_heading_value                     = zeros(num_of_samples,1);
ve_value                                = zeros(num_of_samples,1);
vn_value                                = zeros(num_of_samples,1);
vu_value                                = zeros(num_of_samples,1);
a_pos_value                             = zeros(num_of_samples,1);
b_pos_value                             = zeros(num_of_samples,1);
c_pos_value                             = zeros(num_of_samples,1);
d_pos_value                             = zeros(num_of_samples,1);
lat_value                               = zeros(num_of_samples,1);
lon_value                               = zeros(num_of_samples,1);
alt_value                               = zeros(num_of_samples,1);
Rn_value                                = zeros(num_of_samples,1);
Rm_value                                = zeros(num_of_samples,1);
EP_value                                = zeros(num_of_samples,1);
NP_value                                = zeros(num_of_samples,1);

gyro_bias_x_value                       = zeros(num_of_samples,1);
gyro_bias_y_value                       = zeros(num_of_samples,1);
gyro_bias_z_value                       = zeros(num_of_samples,1);
acc_bias_x_value                        = zeros(num_of_samples,1);
acc_bias_y_value                        = zeros(num_of_samples,1);
acc_bias_z_value                        = zeros(num_of_samples,1);

gyro_sf_x_value                         = zeros(num_of_samples,1);
gyro_sf_y_value                         = zeros(num_of_samples,1);
gyro_sf_z_value                         = zeros(num_of_samples,1);
acc_sf_x_value                          = zeros(num_of_samples,1);
acc_sf_y_value                          = zeros(num_of_samples,1);
acc_sf_z_value                          = zeros(num_of_samples,1);

g_value                                 = zeros(num_of_samples,1);

gyro_rw_stdv_x_value                    = zeros(num_of_samples,1);
gyro_rw_stdv_y_value                    = zeros(num_of_samples,1);
gyro_rw_stdv_z_value                    = zeros(num_of_samples,1);
gyro_bias_gauss_markov_stdv_x_value     = zeros(num_of_samples,1);
gyro_bias_gauss_markov_stdv_y_value     = zeros(num_of_samples,1);
gyro_bias_gauss_markov_stdv_z_value     = zeros(num_of_samples,1);
gyro_bias_x_time_cnst_value             = zeros(num_of_samples,1);
gyro_bias_y_time_cnst_value             = zeros(num_of_samples,1);
gyro_bias_z_time_cnst_value             = zeros(num_of_samples,1);
gyro_sf_gauss_markov_stdv_x_value       = zeros(num_of_samples,1);
gyro_sf_gauss_markov_stdv_y_value       = zeros(num_of_samples,1);
gyro_sf_gauss_markov_stdv_z_value       = zeros(num_of_samples,1);
gyro_sf_x_time_cnst_value               = zeros(num_of_samples,1);
gyro_sf_y_time_cnst_value               = zeros(num_of_samples,1);
gyro_sf_z_time_cnst_value               = zeros(num_of_samples,1);

acc_rw_stdv_x_value                     = zeros(num_of_samples,1);
acc_rw_stdv_y_value                     = zeros(num_of_samples,1);
acc_rw_stdv_z_value                     = zeros(num_of_samples,1);
acc_bias_gauss_markov_stdv_x_value      = zeros(num_of_samples,1);
acc_bias_gauss_markov_stdv_y_value      = zeros(num_of_samples,1);
acc_bias_gauss_markov_stdv_z_value      = zeros(num_of_samples,1);
acc_bias_x_time_cnst_value              = zeros(num_of_samples,1);
acc_bias_y_time_cnst_value              = zeros(num_of_samples,1);
acc_bias_z_time_cnst_value              = zeros(num_of_samples,1);
acc_sf_gauss_markov_stdv_x_value        = zeros(num_of_samples,1);
acc_sf_gauss_markov_stdv_y_value        = zeros(num_of_samples,1);
acc_sf_gauss_markov_stdv_z_value        = zeros(num_of_samples,1);
acc_sf_x_time_cnst_value                = zeros(num_of_samples,1);
acc_sf_y_time_cnst_value                = zeros(num_of_samples,1);
acc_sf_z_time_cnst_value                = zeros(num_of_samples,1);

a_pos_rw_stdv_value                     = zeros(num_of_samples,1);
b_pos_rw_stdv_value                     = zeros(num_of_samples,1);
c_pos_rw_stdv_value                     = zeros(num_of_samples,1);
d_pos_rw_stdv_value                     = zeros(num_of_samples,1);
alt_rw_stdv_value                       = zeros(num_of_samples,1);

receiver_clk_bias_value                 = zeros(num_of_samples,1);
receiver_clk_bias_stdv_value            = zeros(num_of_samples,1);
receiver_clk_drift_value                = zeros(num_of_samples,1);
receiver_clk_drift_stdv_value           = zeros(num_of_samples,1);

wg_noise_value                          = zeros(num_of_samples,1);

%--> Initialize the state
Euler_roll_value(1)     = ems_data.roll(1);
Euler_pitch_value(1)    = ems_data.pitch(1);
Euler_heading_value(1)  = ems_data.heading(1);
attitude_quat_value = angle2quat(Euler_heading_value(1),Euler_pitch_value(1),Euler_roll_value(1),'ZYX');
C_BL_value = angle2dcm(Euler_heading_value(1),Euler_pitch_value(1),Euler_roll_value(1),'ZYX');
C_LB_value = C_BL_value';
a_value(1) = attitude_quat_value(1);
b_value(1) = attitude_quat_value(2);
c_value(1) = attitude_quat_value(3);
d_value(1) = attitude_quat_value(4);
lat_value(1)  = ems_data.lat(1)*R2D;
lon_value(1)  = ems_data.lon(1)*R2D;
alt_value(1)  = ems_data.h(1);
Rn_value(1)   = earth_a/sqrt(1-earth_e2*sin(lat_value(1)*D2R)*sin(lat_value(1)*D2R));
Rm_value(1)   = earth_a*(1-earth_e2)/((1-earth_e2*sin(lat_value(1)*D2R)*sin(lat_value(1)*D2R))^(1.5));
EP_value(1) = ems_data.east(1);
NP_value(1) = ems_data.north(1);
g_value(1) = gravity_fun(Rm_value(1),Rn_value(1),alt_value(1),lat_value(1));
C_EN_value = C_EN_fun(lat_value(1)*D2R , lon_value(1)*D2R);
pos_quat_vector = dcm2quat (C_EN_value');
a_pos_value(1) = pos_quat_vector(1)/sqrt(sum(pos_quat_vector.^2));
b_pos_value(1) = pos_quat_vector(2)/sqrt(sum(pos_quat_vector.^2));
c_pos_value(1) = pos_quat_vector(3)/sqrt(sum(pos_quat_vector.^2));
d_pos_value(1) = pos_quat_vector(4)/sqrt(sum(pos_quat_vector.^2));
ve_value(1) = ems_data.vel_N(1,1);
vn_value(1) = ems_data.vel_N(1,2);
vu_value(1) = ems_data.vel_N(1,3);

%--> Initialize biases
gyro_bias_x_value(1) = 0.0;
gyro_bias_y_value(1) = 0.0;
gyro_bias_z_value(1) = 0.0;
acc_bias_x_value(1)  = 0.0;
acc_bias_y_value(1)  = 0.0;
acc_bias_z_value(1)  = 0.0;

%--> Initialize scale factors
gyro_sf_x_value(1) = 0.0;
gyro_sf_y_value(1) = 0.0;
gyro_sf_z_value(1) = 0.0;
acc_sf_x_value(1)  = 0.0;
acc_sf_y_value(1)  = 0.0;
acc_sf_z_value(1)  = 0.0;

%--> Noise is to set to zero in INS prediction
wg_noise_value(1) = 0.0;

%--> Position system noise in Q matrix
a_pos_rw_stdv_value(1)  = 0.0;
b_pos_rw_stdv_value(1)  = 0.0;
c_pos_rw_stdv_value(1)  = 0.0;
d_pos_rw_stdv_value(1)  = 0.0;
alt_rw_stdv_value(1)    = 0.0;

%--> Accelerometer System noise and Gauss-Markov model parameters
acc_rw_stdv_x_value(1) = parameters(1);
acc_rw_stdv_y_value(1) = parameters(2);
acc_rw_stdv_z_value(1) = parameters(3);
acc_bias_gauss_markov_stdv_x_value (1) = parameters(4);
acc_bias_gauss_markov_stdv_y_value (1) = parameters(5);
acc_bias_gauss_markov_stdv_z_value (1) = parameters(6);
acc_bias_x_time_cnst_value(1) = parameters(7);
acc_bias_y_time_cnst_value(1) = parameters(8);
acc_bias_z_time_cnst_value(1) = parameters(9);
acc_sf_gauss_markov_stdv_x_value (1) = parameters(10);
acc_sf_gauss_markov_stdv_y_value (1) = parameters(11);
acc_sf_gauss_markov_stdv_z_value (1) = parameters(12);
acc_sf_x_time_cnst_value(1) = parameters(13);
acc_sf_y_time_cnst_value(1) = parameters(14);
acc_sf_z_time_cnst_value(1) = parameters(15);

%--> Gyroscope System noise and Gauss-Markov model parameters
gyro_rw_stdv_x_value(1) = parameters(16);
gyro_rw_stdv_y_value(1) = parameters(17);
gyro_rw_stdv_z_value(1) = parameters(18);
gyro_bias_gauss_markov_stdv_x_value (1) = parameters(19);
gyro_bias_gauss_markov_stdv_y_value (1) = parameters(20);
gyro_bias_gauss_markov_stdv_z_value (1) = parameters(21);
gyro_bias_x_time_cnst_value(1) = parameters(22);
gyro_bias_y_time_cnst_value(1) = parameters(23);
gyro_bias_z_time_cnst_value(1) = parameters(24);
gyro_sf_gauss_markov_stdv_x_value (1) = parameters(25);
gyro_sf_gauss_markov_stdv_y_value (1) = parameters(26);
gyro_sf_gauss_markov_stdv_z_value (1) = parameters(27);
gyro_sf_x_time_cnst_value(1) = parameters(28);
gyro_sf_y_time_cnst_value(1) = parameters(29);
gyro_sf_z_time_cnst_value(1) = parameters(30);

%--> Receiver clock errors parameters
receiver_clk_bias_stdv_value(1) = parameters(31);
receiver_clk_drift_stdv_value(1) = parameters(32);

%--> Initialize state error covariances
east_pos_error_covariance(1) = parameters(33);
north_pos_error_covariance(1) = parameters(34);
alt_error_covariance(1) = parameters(35);
ve_error_covariance(1) = parameters(36);
vn_error_covariance(1) = parameters(37);
vu_error_covariance(1) = parameters(38);
b_error_covariance(1) = parameters(39);
c_error_covariance(1) = parameters(40);
d_error_covariance(1) = parameters(41);
gyro_bias_x_error_covariance(1) = parameters(42);
gyro_bias_y_error_covariance(1) = parameters(43);
gyro_bias_z_error_covariance(1) = parameters(44);
acc_bias_x_error_covariance(1) = parameters(45);
acc_bias_y_error_covariance(1) = parameters(46);
acc_bias_z_error_covariance(1) = parameters(47);
gyro_sf_x_error_covariance(1) = parameters(48);
gyro_sf_y_error_covariance(1) = parameters(49);
gyro_sf_z_error_covariance(1) =parameters(50);
acc_sf_x_error_covariance(1) = parameters(51);
acc_sf_y_error_covariance(1) = parameters(52);
acc_sf_z_error_covariance(1) = parameters(53);
receiver_clk_bias_error_covariance(1) =parameters(54);
receiver_clk_drift_error_covariance(1) =parameters(55);

P = diag([
    east_pos_error_covariance(1);
    north_pos_error_covariance(1);
    alt_error_covariance(1);
    ve_error_covariance(1);
    vn_error_covariance(1);
    vu_error_covariance(1);
    b_error_covariance(1);
    c_error_covariance(1);
    d_error_covariance(1);
    gyro_bias_x_error_covariance(1);
    gyro_bias_y_error_covariance(1);
    gyro_bias_z_error_covariance(1);
    acc_bias_x_error_covariance(1);
    acc_bias_y_error_covariance(1);
    acc_bias_z_error_covariance(1);
    gyro_sf_x_error_covariance(1);
    gyro_sf_y_error_covariance(1);
    gyro_sf_z_error_covariance(1);
    acc_sf_x_error_covariance(1);
    acc_sf_y_error_covariance(1);
    acc_sf_z_error_covariance(1);
    receiver_clk_bias_error_covariance(1);
    receiver_clk_drift_error_covariance(1)
    ]);

%--> Draw a 3D Aircraft for animation
X = 60;Y = 15;Z = 5;
origin = [X/2 Y/2 Z/2];
initial_vert = ...
    [X Y 0;             %(1)
    0 Y 0;              %(2)
    0 Y Z;              %(3)
    X Y Z;              %(4)
    0 0 Z;              %(5)
    X 0 Z;              %(6)
    X 0 0;              %(7)
    0 0 0;              %(8)
    1.3*X Y/2 0;        %(9)
    1.3*X Y/2 0.5*Z;    %(10)
    ];

for p = 1:length(initial_vert)
    initial_vert(p,:) = initial_vert(p,:) - origin;
end

CubePoints = [initial_vert(:,1),initial_vert(:,2),initial_vert(:,3)];
[ CubeXData , CubeYData , CubeZData ] = get_cube_axis_data(initial_vert);
faces = [1 2 3 4; 4 3 5 6; 6 7 8 5; 1 2 8 7; 6 7 1 4; 2 3 5 8;1 9 10 4;4 10 6 6;6 10 9 7;1 9 7 7];

%--> Main processing loop
for index = 1:length(ems_data.time)-1
    
    time_vector_sec(index) = ems_data.time(index);
    
    %--> Read IMU data
    accel_x(index) = ems_data.sim_accel(index,1);
    accel_y(index) = ems_data.sim_accel(index,2);
    accel_z(index) = ems_data.sim_accel(index,3);

    gyro_x(index) = ems_data.sim_gyro(index,1);
    gyro_y(index) = ems_data.sim_gyro(index,2);
    gyro_z(index) = ems_data.sim_gyro(index,3);      

    %--> Calculate Earth Rotation Rate, Transport Rate, and Corriollis Effect
    w_N_EN_value = w_N_EN_fun(Rm_value(index),Rn_value(index),alt_value(index),lat_value(index)*D2R,ve_value(index),vn_value(index));
    w_N_IE_value = w_N_IE_fun(lat_value(index)*D2R,we_value);
    w_L_IL_value = C_LN*(w_N_EN_value + w_N_IE_value);
    w_B_IL_value = C_LB_value'*w_L_IL_value;
    corriollis_effect_vector_in_B = C_LB_value'*C_LN*cross((w_N_EN_value + 2*w_N_IE_value),[ve_value(index);vn_value(index);vu_value(index)]);
    
    %--> Perform attitude quaternion mechanization
    a_dot_value = a_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),gyro_x(index),gyro_y(index),gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    b_dot_value = b_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),gyro_x(index),gyro_y(index),gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    c_dot_value = c_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),gyro_x(index),gyro_y(index),gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    d_dot_value = d_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),gyro_x(index),gyro_y(index),gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    
    if (~isreal([a_dot_value b_dot_value c_dot_value d_dot_value]))
        disp('imaginary quaternion rate');return;
    end
    
    %--> Advance quaternion
    a_value(index+1) = a_value(index) + a_dot_value*sampling_period_sec;
    b_value(index+1) = b_value(index) + b_dot_value*sampling_period_sec;
    c_value(index+1) = c_value(index) + c_dot_value*sampling_period_sec;
    d_value(index+1) = d_value(index) + d_dot_value*sampling_period_sec;
    
    %--> Normalize the quaternion
    normalized_quat = quatnormalize([a_value(index+1) b_value(index+1) c_value(index+1) d_value(index+1)]);
    a_value(index+1) = normalized_quat(1);
    b_value(index+1) = normalized_quat(2);
    c_value(index+1) = normalized_quat(3);
    d_value(index+1) = normalized_quat(4);
    
    %--> Calculate Euler angles from Quaternion
    [A, p, r] = quat2angle([a_value(index+1) b_value(index+1) c_value(index+1) d_value(index+1)],'ZYX');
    Euler_pitch_value(index+1)     = p;
    Euler_roll_value(index+1)      = r;
    Euler_heading_value(index+1)   = A;
    
    if (~isreal([Euler_roll_value(index+1) Euler_pitch_value(index+1) Euler_heading_value(index+1)]))
        disp('imaginary Euler angles');return;
    end
    
    %-->Reset quaternion in case -180/180 boundaries crossed
    quat_vector = angle2quat(Euler_heading_value(index+1),Euler_pitch_value(index+1),Euler_roll_value(index+1),'ZYX');
    a_value(index+1) = quat_vector(1);
    b_value(index+1) = quat_vector(2);
    c_value(index+1) = quat_vector(3);
    d_value(index+1) = quat_vector(4);
    
    normalized_quat = quatnormalize([a_value(index+1) b_value(index+1) c_value(index+1) d_value(index+1)]);
    a_value(index+1) = normalized_quat(1);
    b_value(index+1) = normalized_quat(2);
    c_value(index+1) = normalized_quat(3);
    d_value(index+1) = normalized_quat(4);
    
    C_BL_value = angle2dcm(Euler_heading_value(index+1),Euler_pitch_value(index+1),Euler_roll_value(index+1),'ZYX');
    C_LB_value = C_BL_value';
    
    %--> Advance velocity
    V_N_dot_value = V_N_dot_fun(Rm_value(index),Rn_value(index),accel_x(index),accel_y(index),accel_z(index),...
                                acc_sf_x_value(index),acc_sf_y_value(index),acc_sf_z_value(index),...
                                acc_bias_x_value(index),acc_bias_y_value(index),acc_bias_z_value(index),0.0,0.0,0.0,...
                                alt_value(index),b_value(index),c_value(index),d_value(index),lat_value(index)*D2R,ve_value(index),vn_value(index),vu_value(index),we_value,wg_noise_value(index));
    
    if (~isreal(V_N_dot_value))
        disp('imaginary velocity rate');return;
    end
    
    ve_value(index+1) = ve_value(index) + V_N_dot_value(1)*sampling_period_sec;
    vn_value(index+1) = vn_value(index) + V_N_dot_value(2)*sampling_period_sec;
    vu_value(index+1) = vu_value(index) + V_N_dot_value(3)*sampling_period_sec;
    
    %--> Advance position
    pos_quat_dot_value = pos_quat_dot_fun_EN(Rm_value(index),Rn_value(index),a_pos_rw_stdv_value(index), alt_value(index),...
                                             b_pos_value(index),b_pos_rw_stdv_value(index), c_pos_value(index),c_pos_rw_stdv_value(index), d_pos_value(index),d_pos_rw_stdv_value(index),...
                                             lat_value(index)*D2R,ve_value(index),vn_value(index), wg_noise_value(index));
    
    if (~isreal(pos_quat_dot_value))
        disp('imaginary position quaternion rate');return;
    end
    
    a_pos_value(index+1) = a_pos_value(index) + pos_quat_dot_value(1)*sampling_period_sec;
    b_pos_value(index+1) = b_pos_value(index) + pos_quat_dot_value(2)*sampling_period_sec;
    c_pos_value(index+1) = c_pos_value(index) + pos_quat_dot_value(3)*sampling_period_sec;
    d_pos_value(index+1) = d_pos_value(index) + pos_quat_dot_value(4)*sampling_period_sec;
    
    normalization_factor = sqrt(a_pos_value(index+1)^2 + b_pos_value(index+1)^2 + c_pos_value(index+1)^2 + d_pos_value(index+1)^2);
    a_pos_value(index+1) = a_pos_value(index+1)/normalization_factor;
    b_pos_value(index+1) = b_pos_value(index+1)/normalization_factor;
    c_pos_value(index+1) = c_pos_value(index+1)/normalization_factor;
    d_pos_value(index+1) = d_pos_value(index+1)/normalization_factor;
    
    %--> Calculate position DCM from position quaternion
    C_EN_value = C_EN_from_pos_quat_fun(b_pos_value(index+1),c_pos_value(index+1),d_pos_value(index+1));
    
    if (~isreal(C_EN_value))
        disp('imaginary position DCM');return;
    end
    
    %--> Calculate lat and lon from position DCM
    lon_from_C_EN_value = atan2(C_EN_value(1,3),C_EN_value(3,3));
    lat_from_C_EN_value = atan2(C_EN_value(2,3),sqrt(C_EN_value(2,1)^2+C_EN_value(2,2)^2));
    lat_value(index+1) = lat_from_C_EN_value*R2D;
    lon_value(index+1) = lon_from_C_EN_value*R2D;
    
    %--> Advance altitude
    alt_value(index+1) = alt_value(index) + vu_value(index)*sampling_period_sec;
    
    %--> Calculate new earth parameters and east, north position components
    Rn_value(index+1)  = earth_a/sqrt(1-earth_e2*sin(lat_value(index+1)*D2R)*sin(lat_value(index+1)*D2R));
    Rm_value(index+1)  = earth_a*(1-earth_e2)/((1-earth_e2*sin(lat_value(index+1)*D2R)*sin(lat_value(index+1)*D2R))^(1.5));
    EP_value(index+1) = (lon_value(index+1)-lon_value(1))*D2R*(Rn_value(index+1)+alt_value(index+1))*cos(lat_value(index+1)*D2R);
    NP_value(index+1) = (lat_value(index+1)-lat_value(1))*D2R*(Rm_value(index+1)+alt_value(index+1));
    
    %--> Advance biases
    gyro_bias_x_value(index+1) = gyro_bias_x_value(index);
    gyro_bias_y_value(index+1) = gyro_bias_y_value(index);
    gyro_bias_z_value(index+1) = gyro_bias_z_value(index);
    acc_bias_x_value(index+1) = acc_bias_x_value(index);
    acc_bias_y_value(index+1) = acc_bias_y_value(index);
    acc_bias_z_value(index+1) = acc_bias_z_value(index);
    
    %--> Advance scale factors
    gyro_sf_x_value(index+1) = gyro_sf_x_value(index);
    gyro_sf_y_value(index+1) = gyro_sf_y_value(index);
    gyro_sf_z_value(index+1) = gyro_sf_z_value(index);
    acc_sf_x_value(index+1) = acc_sf_x_value(index);
    acc_sf_y_value(index+1) = acc_sf_y_value(index);
    acc_sf_z_value(index+1) = acc_sf_z_value(index);
    
    %--> Advance receiver clock bias and drift
    receiver_clk_bias_value(index+1) = receiver_clk_bias_value(index);
    receiver_clk_drift_value(index+1) = receiver_clk_drift_value(index);
    
    %--> Kalman Filter model matrices calculation
    F_value = F_fun(Rm_value(index),Rn_value(index),accel_x(index),accel_y(index),accel_z(index),...
                    acc_sf_x_value(index),acc_sf_y_value(index),acc_sf_z_value(index),...
                    acc_bias_x_value(index),acc_bias_y_value(index),acc_bias_z_value(index),...
                    acc_rw_stdv_x_value(index),acc_rw_stdv_y_value(index),acc_rw_stdv_z_value(index),...
                    acc_sf_x_time_cnst_value(index),acc_sf_y_time_cnst_value(index),acc_sf_z_time_cnst_value(index),...
                    acc_bias_x_time_cnst_value(index),acc_bias_y_time_cnst_value(index),acc_bias_z_time_cnst_value(index),...
                    alt_value(index),b_value(index),c_value(index),d_value(index),...
                    gyro_x(index),gyro_y(index),gyro_z(index),...
                    gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),...
                    gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),...
                    gyro_rw_stdv_x_value(index),gyro_rw_stdv_y_value(index),gyro_rw_stdv_z_value(index),...
                    gyro_sf_x_time_cnst_value(index),gyro_sf_y_time_cnst_value(index),gyro_sf_z_time_cnst_value(index),...
                    gyro_bias_x_time_cnst_value(index),gyro_bias_y_time_cnst_value(index),gyro_bias_z_time_cnst_value(index),...
                    lat_value(index)*R2D,ve_value(index),vn_value(index),vu_value(index),we_value,wg_noise_value(index));
    
    if (~isreal(F_value))
        disp('imaginary transition matrix');return;
    end
    
    Phi = eye(size(F_value)) + F_value*sampling_period_sec;
    
    G_value = G_fun(acc_rw_stdv_x_value(index),acc_rw_stdv_y_value(index),acc_rw_stdv_z_value(index),...
                    acc_sf_x_time_cnst_value(index),acc_sf_y_time_cnst_value(index),acc_sf_z_time_cnst_value(index),...
                    acc_bias_x_time_cnst_value(index),acc_bias_z_time_cnst_value(index),acc_bias_y_time_cnst_value(index),...
                    acc_sf_gauss_markov_stdv_x_value(index),acc_sf_gauss_markov_stdv_y_value(index),acc_sf_gauss_markov_stdv_z_value(index),...
                    acc_bias_gauss_markov_stdv_x_value(index),acc_bias_gauss_markov_stdv_y_value(index),acc_bias_gauss_markov_stdv_z_value(index),....
                    alt_rw_stdv_value(index),...
                    b_value(index), c_value(index),d_value(index),...
                    gyro_rw_stdv_x_value(index),gyro_rw_stdv_y_value(index),gyro_rw_stdv_z_value(index),...
                    gyro_sf_x_time_cnst_value(index),gyro_sf_y_time_cnst_value(index),gyro_sf_z_time_cnst_value(index),...
                    gyro_bias_x_time_cnst_value(index),gyro_bias_y_time_cnst_value(index),gyro_bias_z_time_cnst_value(index),...
                    gyro_sf_gauss_markov_stdv_x_value(index),gyro_sf_gauss_markov_stdv_y_value(index),gyro_sf_gauss_markov_stdv_z_value(index),...
                    gyro_bias_gauss_markov_stdv_x_value(index),gyro_bias_gauss_markov_stdv_y_value(index),gyro_bias_gauss_markov_stdv_z_value(index),...
                    receiver_clk_bias_stdv_value(index), receiver_clk_drift_stdv_value(index));

    if (~isreal(G_value))
        disp('imaginary noise shaping matrix');return;
    end
    
    Qd = sampling_period_sec^2*(G_value*G_value');
    
    %--> Advance state error covariance matrix
    P = Phi*P*Phi' + Qd;
    
    if (~isreal(P))
        disp('imaginary noise covariance matrix');return;
    end
    
    %--> Apply updates from obervations
    if (mod(index,sampling_freq) == 0)
        %--> Initialize error vector and matrices
        clear z H R;
        
        no_of_sats = ems_data.no_of_sat{((index)/sampling_freq)+1};        
        
        %--> Convert reciever position from Geodetic frame to ECEF frame
        x = (Rn_value(index+1) + alt_value(index+1))*cos(lat_from_C_EN_value)*cos(lon_from_C_EN_value);
        y = (Rn_value(index+1) + alt_value(index+1))*cos(lat_from_C_EN_value)*sin(lon_from_C_EN_value);
        z = ((Rn_value(index+1)*(1-earth_e2))+alt_value(index+1))*sin(lat_from_C_EN_value);        
        
        %--> Set the innovation sequence
        deltaPx = (x.*ones(no_of_sats,1)) - (ems_data.sat_Px{((index)/sampling_freq)+1});
        deltaPy = (y.*ones(no_of_sats,1)) - (ems_data.sat_Py{((index)/sampling_freq)+1});
        deltaPz = (z.*ones(no_of_sats,1)) - (ems_data.sat_Pz{((index)/sampling_freq)+1});

        range_from_INS = sqrt(deltaPx.^2+deltaPy.^2+deltaPz.^2);      
        range_from_GNSS = cell2mat(ems_data.sim_sat_range(((index)/sampling_freq)+1)) + (ems_data.noiseInfo.clk_bias - receiver_clk_bias_value(index+1)).*ones(no_of_sats,1);
        
        z = range_from_GNSS - range_from_INS;
        
        %--> Set the measurment matrix jacobian
        H = zeros(no_of_sats, 23);
        
        for m=1:no_of_sats
            H(m,:) = H_row_range_fun(Rm_value(index+1),Rn_value(index+1),alt_value(index+1),...
                                    (lon_from_C_EN_value - ems_data.lon(1))*(Rn_value(index+1)+alt_value(index+1))*cos(lat_from_C_EN_value)+ems_data.east(1),...
                                    ems_data.east(1),ems_data.lat(1),ems_data.lon(1),(lat_from_C_EN_value - ems_data.lat(1))*(Rm_value(index+1)+alt_value(index+1)),...
                                    ems_data.north(1),ems_data.sat_Px{(index/sampling_freq)+1}(m),ems_data.sat_Py{(index/sampling_freq)+1}(m),ems_data.sat_Pz{(index/sampling_freq)+1}(m));
        end
        
        %--> Adjust measurment noise matrix
        R = diag(ones(1*no_of_sats,1));   
        
        range_error_std = parameters(56);
       
        R(1:no_of_sats,1:no_of_sats) = R(1:no_of_sats,1:no_of_sats).*(range_error_std.^2);
        
        
        if (~isreal(P))
            disp('imaginary updated error covariance matrix ');return;
        end
        
        %--> EKF update
        K = (P*H')/(H*P*H'+R);
        state_correction_vector = K*z;
        P = P - K*H*P;
        
        %--> Normalize P
        P = (P+P')/2;

        if (~isreal(P))
            disp('imaginary error covariance matrix ');return;
        end
        
        %--> Correct states
        correct_states_TC;
    end
    

    gyro_rw_stdv_x_value(index+1) = gyro_rw_stdv_x_value(index);
    gyro_rw_stdv_y_value(index+1) = gyro_rw_stdv_y_value(index);
    gyro_rw_stdv_z_value(index+1) = gyro_rw_stdv_z_value(index);
    gyro_bias_gauss_markov_stdv_x_value (index+1) = gyro_bias_gauss_markov_stdv_x_value (index);
    gyro_bias_gauss_markov_stdv_y_value (index+1) = gyro_bias_gauss_markov_stdv_y_value (index);
    gyro_bias_gauss_markov_stdv_z_value (index+1) = gyro_bias_gauss_markov_stdv_z_value (index);
    gyro_bias_x_time_cnst_value(index+1) = gyro_bias_x_time_cnst_value(index);
    gyro_bias_y_time_cnst_value(index+1) = gyro_bias_y_time_cnst_value(index);
    gyro_bias_z_time_cnst_value(index+1) = gyro_bias_z_time_cnst_value(index);
    gyro_sf_gauss_markov_stdv_x_value (index+1) = gyro_sf_gauss_markov_stdv_x_value (index);
    gyro_sf_gauss_markov_stdv_y_value (index+1) = gyro_sf_gauss_markov_stdv_y_value (index);
    gyro_sf_gauss_markov_stdv_z_value (index+1) = gyro_sf_gauss_markov_stdv_z_value (index);
    gyro_sf_x_time_cnst_value(index+1) = gyro_sf_x_time_cnst_value(index);
    gyro_sf_y_time_cnst_value(index+1) = gyro_sf_y_time_cnst_value(index);
    gyro_sf_z_time_cnst_value(index+1) = gyro_sf_z_time_cnst_value(index);
    
    acc_rw_stdv_x_value(index+1) = acc_rw_stdv_x_value(index);
    acc_rw_stdv_y_value(index+1) = acc_rw_stdv_y_value(index);
    acc_rw_stdv_z_value(index+1) = acc_rw_stdv_z_value(index);
    acc_bias_gauss_markov_stdv_x_value (index+1) = acc_bias_gauss_markov_stdv_x_value (index);
    acc_bias_gauss_markov_stdv_y_value (index+1) = acc_bias_gauss_markov_stdv_y_value (index);
    acc_bias_gauss_markov_stdv_z_value (index+1) = acc_bias_gauss_markov_stdv_z_value (index);
    acc_bias_x_time_cnst_value(index+1) = acc_bias_x_time_cnst_value(index);
    acc_bias_y_time_cnst_value(index+1) = acc_bias_y_time_cnst_value(index);
    acc_bias_z_time_cnst_value(index+1) = acc_bias_z_time_cnst_value(index);
    acc_sf_gauss_markov_stdv_x_value (index+1) = acc_sf_gauss_markov_stdv_x_value (index);
    acc_sf_gauss_markov_stdv_y_value (index+1) = acc_sf_gauss_markov_stdv_y_value (index);
    acc_sf_gauss_markov_stdv_z_value (index+1) = acc_sf_gauss_markov_stdv_z_value (index);
    acc_sf_x_time_cnst_value(index+1) = acc_sf_x_time_cnst_value(index);
    acc_sf_y_time_cnst_value(index+1) = acc_sf_y_time_cnst_value(index);
    acc_sf_z_time_cnst_value(index+1) = acc_sf_z_time_cnst_value(index);
    
    a_pos_rw_stdv_value(index + 1) = a_pos_rw_stdv_value(index);
    b_pos_rw_stdv_value(index + 1) = b_pos_rw_stdv_value(index);
    c_pos_rw_stdv_value(index + 1) = c_pos_rw_stdv_value(index);
    d_pos_rw_stdv_value(index + 1) = d_pos_rw_stdv_value(index);
    alt_rw_stdv_value(index + 1)   = alt_rw_stdv_value(index);
    
    wg_noise_value(index+1) = 0.0;
    
    g_value(index + 1) = gravity_fun(Rm_value(index),Rn_value(index),alt_value(index),lat_value(index));
    
    %--> Disply results in 3D
    C_LB_value_plus_90 = C_LB_from_Euler_fun(Euler_roll_value(index),...
        -Euler_pitch_value(index),...
        pi/2 - Euler_heading_value(index));
    
    updated_vert = C_LB_value_plus_90*initial_vert';
    
    updated_vert = updated_vert';
    
    for p = 1:length(updated_vert)
        updated_vert(p,:) = updated_vert(p,:) + [EP_value(index), NP_value(index), alt_value(index)];
    end
    
    [ CubeXData , CubeYData , CubeZData ] = get_cube_axis_data(updated_vert);
    CubePoints = [updated_vert(:,1),updated_vert(:,2),updated_vert(:,3)];
    
    %--> Create first figure after the 5th IMU sample
    if(index == 5)
        figure;
        hold on;
        plot3(ems_data.east,ems_data.north,ems_data.h, 'LineWidth', 1, 'color', [0.93 .69 .13]);grid on;xlabel('east');ylabel('north');zlabel('alt');
        view(-1,50);title('3D Trajectory of Ground Truth, Noisy Updates, and EKF Solution');
        h_noisy_updates = plot3(ems_data.east_noisy(1:20:index),ems_data.north_noisy(1:20:index),ems_data.h_noisy(1:20:index), '-*','MarkerSize', 5,'LineWidth', 0.5, 'color','k');grid on;xlabel('east');ylabel('north');zlabel('alt');
        h_ekf_position  = plot3(EP_value(1:index), NP_value(1:index), alt_value(1:index), 'color', 'r','LineWidth', 2);hold on;grid on;
        set(h_ekf_position,'YDataSource','NP_value(1:index)');
        set(h_ekf_position,'XDataSource','EP_value(1:index)');
        set(h_ekf_position,'ZDataSource','alt_value(1:index)');
        set(h_noisy_updates,'YDataSource','ems_data.north_noisy(1:20:index)');
        set(h_noisy_updates,'XDataSource','ems_data.east_noisy(1:20:index)');
        set(h_noisy_updates,'ZDataSource','ems_data.h_noisy(1:20:index)');
        h_vehicle = patch('Faces',faces,'Vertices',CubePoints,'FaceVertexCData',hsv(10),'FaceColor','flat');
        set(h_vehicle,'XData',CubeXData,'YData',CubeYData,'ZData',CubeZData);
        set(gca,'FontSize',12);
        xlabel('East (m)');ylabel('North (m)');zlabel('Altitude(m)');
        axis([min(ems_data.east) max(ems_data.east) min(ems_data.north) max(ems_data.north) -100 300]);
        legend({'Ground Truth','Noisy GPS','EKF Solution'}, 'FontSize', 12); 
    else
        if enable_animation == 1
            if (mod(index,sampling_freq) == 0)
                refreshdata(h_ekf_position, 'caller');
                refreshdata(h_noisy_updates, 'caller');
                set(h_vehicle,'XData',CubeXData,'YData',CubeYData,'ZData',CubeZData);
                pause(0.05);
            end
        end
    end
    fprintf('processing epoch %d/%d\n',index,num_of_samples);
end

%--> Refresh figures
refreshdata(h_ekf_position, 'caller');
refreshdata(h_noisy_updates, 'caller');
set(h_vehicle,'XData',CubeXData,'YData',CubeYData,'ZData',CubeZData);
pause(0.05);

%--> plot reference and EKF result
figure;
subplot (3,1,1); 
plot(ems_data.time,ems_data.roll.*R2D,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,Euler_roll_value.*R2D,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('Roll (deg)');title('Roll');
subplot (3,1,2); 
plot(ems_data.time,ems_data.pitch.*R2D,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,Euler_pitch_value.*R2D,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('Pitch (deg)');title('Pitch');
subplot (3,1,3); 
plot(ems_data.time,ems_data.heading.*R2D,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,Euler_heading_value.*R2D,'b', 'LineWidth',2);legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('Heading (deg)');title('Heading');

figure;
subplot (3,1,1); 
plot(ems_data.time,ems_data.vel(:,1),'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,vn_value,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('V_n (m/s)');title('V_n');
subplot (3,1,2); 
plot(ems_data.time,ems_data.vel(:,2),'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,ve_value,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('V_e (m/s)');title('V_e');
subplot (3,1,3); 
plot(ems_data.time,ems_data.vel(:,3),'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,-vu_value,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('V_d (m/s)');title('V_d');

figure;
subplot (3,1,1); 
plot(ems_data.time,ems_data.north,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,NP_value,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('North Position (m)');title('North Position');
subplot (3,1,2); 
plot(ems_data.time,ems_data.east,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,EP_value,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('East Position (m)');title('East Position');
subplot (3,1,3); 
plot(ems_data.time,-ems_data.h,'r', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot(ems_data.time,-alt_value,'b', 'LineWidth',2); legend({'Reference','EKF Solution'}, 'FontSize', 12);  xlabel('Time (s)');ylabel('Down Position (m)');title('Down Position');

figure;
plot(ems_data.time,receiver_clk_bias_value, 'b', 'LineWidth',2); hold on; grid on; set(gca,'FontSize',12);
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.clk_bias*ones(1,2), '--b', 'LineWidth', 2);
legend({'Estimated clock bias', sprintf('True clock bias(= %1.2f)',  ems_data.noiseInfo.clk_bias)}, 'FontSize', 12);xlabel('Time (s)');ylabel('Bias (m)');title('Receiver Clock Bias');

figure;
plot(ems_data.time, gyro_bias_x_value*R2D, 'r', 'LineWidth',2);grid on;hold on;set(gca,'FontSize',12);
plot(ems_data.time, gyro_bias_y_value*R2D, 'g', 'LineWidth',2);hold on;
plot(ems_data.time, gyro_bias_z_value*R2D, 'b', 'LineWidth',2);hold on;
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.gyro_bias(1)*ones(1,2), '--r', 'LineWidth', 2);hold on;
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.gyro_bias(2)*ones(1,2), '--g', 'LineWidth', 2);hold on;
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.gyro_bias(3)*ones(1,2), '--b', 'LineWidth', 2);hold on;
legend({'Estimated bias X','Estimated bias Y', 'Estimated bias Z', ...
sprintf('True bias X (= %1.2f)', ems_data.noiseInfo.gyro_bias(1)),...
sprintf('True bias Y (= %1.2f)',  ems_data.noiseInfo.gyro_bias(2)),...
sprintf('True bias Z (= %1.2f)',  ems_data.noiseInfo.gyro_bias(3))}, 'FontSize', 12); xlabel('Time (s)');ylabel('Bias (deg/s)');title('Gyroscope Biases');

figure;
plot(ems_data.time, acc_bias_x_value, 'r', 'LineWidth',2);grid on;hold on; set(gca,'FontSize',12);
plot(ems_data.time, acc_bias_y_value, 'g', 'LineWidth',2);hold on;
plot(ems_data.time, acc_bias_z_value, 'b', 'LineWidth',2);hold on;
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.accel_bias(1)*ones(1,2), '--r', 'LineWidth', 2);hold on;
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.accel_bias(2)*ones(1,2), '--g', 'LineWidth', 2);hold on;
plot([ems_data.time(1), ems_data.time(end)], ems_data.noiseInfo.accel_bias(3)*ones(1,2), '--b', 'LineWidth', 2);hold on;
legend({'Estimated bias X','Estimated bias Y', 'Estimated bias Z', ...
sprintf('True bias X (= %1.2f)', ems_data.noiseInfo.accel_bias(1)),...
sprintf('True bias Y (= %1.2f)',  ems_data.noiseInfo.accel_bias(2)),...
sprintf('True bias Z (= %1.2f)',  ems_data.noiseInfo.accel_bias(3))}, 'FontSize', 12); xlabel('Time (s)');ylabel('Bias (m/s^2)');title('Accelerometer Biases');


%--> Plot errors
v_east_ref_vector       = ems_data.vel_N(:,1)';
v_north_ref_vector      = ems_data.vel_N(:,2)';
v_up_ref_vector         = ems_data.vel_N(:,3)';
p_east_ref_vector       = ems_data.east';
p_north_ref_vector      = ems_data.north';
alt_ref_vector          = ems_data.h';
roll_ref_vector         = ems_data.roll';
pitch_ref_vector        = ems_data.pitch';
heading_ref_vector      = ems_data.heading';

figure;title('Orientation Errors');
subplot(3,1,1);
plot(ems_data.time, (roll_ref_vector' - Euler_roll_value)*R2D,'b', 'LineWidth',2);title('Roll Error'); grid on; xlabel('Time (s)');ylabel('Error (deg)'); set(gca,'FontSize',12);
subplot(3,1,2);
plot(ems_data.time, (pitch_ref_vector' - Euler_pitch_value)*R2D,'b', 'LineWidth',2);title('Pitch Error'); grid on; xlabel('Time (s)');ylabel('Error (deg)'); set(gca,'FontSize',12);
subplot(3,1,3);
plot(ems_data.time, (heading_ref_vector' - Euler_heading_value)*R2D,'b', 'LineWidth',2);title('Heading Error'); grid on; xlabel('Time (s)');ylabel('Error (deg)'); set(gca,'FontSize',12);

figure;title('Velocity Errors');
subplot(3,1,2);
plot(ems_data.time, v_east_ref_vector' - ve_value,'b', 'LineWidth',2);title('V_e Error'); grid on; xlabel('Time (s)');ylabel('Error (m/s)'); set(gca,'FontSize',12);
subplot(3,1,1);
plot(ems_data.time, v_north_ref_vector' - vn_value,'b', 'LineWidth',2);title('V_n Error'); grid on; xlabel('Time (s)');ylabel('Error (m/s)'); set(gca,'FontSize',12);
subplot(3,1,3);
plot(ems_data.time, -v_up_ref_vector' + vu_value,'b', 'LineWidth',2);title('V_d Error'); grid on; xlabel('Time (s)');ylabel('Error (m/s)'); set(gca,'FontSize',12);

figure;title('Position Errors');
subplot(3,1,2);
plot(ems_data.time, p_east_ref_vector' - EP_value,'b', 'LineWidth',2);title('East Position Error'); grid on; xlabel('Time (s)');ylabel('Error (m)'); set(gca,'FontSize',12);
subplot(3,1,1);
plot(ems_data.time, p_north_ref_vector' - NP_value,'b', 'LineWidth',2);title('North Position Error'); grid on; xlabel('Time (s)');ylabel('Error (m)'); set(gca,'FontSize',12);
subplot(3,1,3);
plot(ems_data.time, -alt_ref_vector' + alt_value,'b', 'LineWidth',2);title('Down Position Error'); grid on; xlabel('Time (s)');ylabel('Error (m)'); set(gca,'FontSize',12);

%--> Print orientation, posiiton, and velocity errors
pct = round(length(v_north_ref_vector)*.1);
fprintf('\n\nError mean \n');
fprintf('V north error(m/s) = %.10f\n', sqrt(mean((v_north_ref_vector(pct:end)'-vn_value(pct:end)).^2)));
fprintf('V east error(m/s) = %.10f\n', sqrt(mean((v_east_ref_vector(pct:end)'-ve_value(pct:end)).^2)));
fprintf('V up error(m/s) = %.10f\n', sqrt(mean((v_up_ref_vector(pct:end)'-vu_value(pct:end)).^2)));
fprintf('North error(m) = %.10f\n', sqrt(mean((p_north_ref_vector(pct:end)'-NP_value(pct:end)).^2)));
fprintf('East error(m) = %.10f\n', sqrt(mean((p_east_ref_vector(pct:end)'-EP_value(pct:end)).^2)));
fprintf('Alt error(m) = %.10f\n', sqrt(mean((alt_ref_vector(pct:end)'-alt_value(pct:end)).^2)));
fprintf('Roll error(deg) = %.10f\n', sqrt(mean((roll_ref_vector(pct:end)'-Euler_roll_value(pct:end)).^2))*R2D);
fprintf('Pitch error(deg) = %.10f\n', sqrt(mean((pitch_ref_vector(pct:end)'-Euler_pitch_value(pct:end)).^2))*R2D);
fprintf('Heading error(deg) = %.10f\n', sqrt(mean((heading_ref_vector(pct:end)'-Euler_heading_value(pct:end)).^2))*R2D);

fprintf('\n\nError std\n');
fprintf('V north error std(m/s) = %.10f\n', std(v_north_ref_vector(pct:end)'-vn_value(pct:end)));
fprintf('V east error std(m/s) = %.10f\n', std(v_east_ref_vector(pct:end)'-ve_value(pct:end)));
fprintf('V up error std(m/s) = %.10f\n', std(v_up_ref_vector(pct:end)'-vu_value(pct:end)));
fprintf('North error std(m) = %.10f\n', std(p_north_ref_vector(pct:end)'-NP_value(pct:end)));
fprintf('East error std(m) = %.10f\n', std(p_east_ref_vector(pct:end)'-EP_value(pct:end)));
fprintf('Alt error std(m) = %.10f\n', std(alt_ref_vector(pct:end)'-alt_value(pct:end)));
fprintf('Roll error std(deg) = %.10f\n', std(roll_ref_vector(pct:end)'-Euler_roll_value(pct:end))*R2D);
fprintf('Pitch error std(deg) = %.10f\n', std(pitch_ref_vector(pct:end)'-Euler_pitch_value(pct:end))*R2D);
fprintf('Heading error std(deg) = %.10f\n', std(heading_ref_vector(pct:end)'-Euler_heading_value(pct:end))*R2D);

end
